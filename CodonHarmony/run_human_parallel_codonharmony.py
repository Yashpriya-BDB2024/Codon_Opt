"""
run_human_parallel_codonharmony.py - full human Codon Harmony grid, in parallel, with a
progress bar and resume support. 

Imports codonharmony_grid.py. Reads the 8 human sequences from the
trials workbook, runs every grid cell (default 2 x 3 x 3 x 3 = 54) for each sequence and
writes ONE Excel with the same sheets as the sequential script
(Summary | Results | Grid | Distinct per protein).

Why parallel: one Codon Harmony run is slow (~30 min per sequence), and every
(sequence, cell) run is independent. Seeds are base_seed + cell index (set inside each
run), so the results are IDENTICAL to a sequential run, whatever the worker count or
finishing order.

Usage:
    python run_human_parallel_codonharmony.py                  # full run
    python run_human_parallel_codonharmony.py --test           # tiny test: 2 seqs x 2 cells
    python run_human_parallel_codonharmony.py --resume         # continue an interrupted run
    python run_human_parallel_codonharmony.py --homopolymer 4,5,6   # grid overrides work too
"""

import argparse
import json
import os
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from tqdm import tqdm

import codonharmony_grid as chg          # the grid script, imported as-is

HERE = Path(__file__).resolve().parent
TRIALS = next(HERE.glob("Codon*ptimization*trials*.xlsx"))      # spaces or underscores
TABLE = HERE / "human_codon.csv"
# TABLE = HERE.parent / "host_files" / "Human" / "human_codon.csv"        # same CSV as the DNA Chisel run
SPECIES = "human"
CHUNK = 1                                # cells per parallel task (runs are long -> 1 = best load balance)
# Default workers: 40 (leaves a few cores free on a 44-core box). Windows cannot
# run more than 61 worker processes, so the value is always capped at 61.
DEFAULT_WORKERS = min(40, max(1, (os.cpu_count() or 2) - 1))

_STATE = {}                              # per-worker-process objects (set by _init_worker)


def load_human_sequences():
    """Trials sheet: headers on row 3, data from column B. Keep Homo sapiens rows."""
    d = pd.read_excel(TRIALS, header=2).iloc[:, 1:5]
    d.columns = ["Name", "Protein", "Host", "DNA"]
    d = d.dropna(subset=["Name", "Protein"])
    return d[d["Host"].astype(str).str.lower().str.contains("homo")].reset_index(drop=True)


def _init_worker(profile_json, host_threshold, nested_table, motifs, splice, start,
                 species, seed, domain):
    """Runs once per worker process: loads Codon Harmony's backend (slow, not picklable)."""
    _STATE.update(
        backend=chg.load_backend(profile_json, chg.DEFAULT_RESTRICTION_ENZYMES, host_threshold),
        nested_table=nested_table, motifs=motifs, splice=splice, start=start,
        species=species, seed=seed, domain=domain)


def run_task(name, protein, kind, items):
    """One parallel task: a slice of the grid cells [(cell_index, cell), ...] for ONE sequence."""
    s = _STATE
    rows = []
    for ci, cell in items:
        row, _ = chg.run_cell_row(name, protein, kind, ci, cell, s["backend"], s["nested_table"],
                                  s["motifs"], s["splice"], s["start"], s["species"],
                                  s["seed"], s["domain"])
        rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", action="store_true", help="quick test (2 seqs x 2 cells)")
    ap.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    ap.add_argument("--output", default=str(HERE / "outputs" / "human_codonharmony_grid.xlsx"))
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--host-threshold", type=float, default=chg.HOST_THRESHOLD)
    ap.add_argument("--resume", action="store_true",
                    help="reuse finished runs from <output>.partial.jsonl")
    ap.add_argument("--inner-cycles", default=None, help="Comma-separated override")
    ap.add_argument("--cycles", default=None)
    ap.add_argument("--max-relax", default=None)
    ap.add_argument("--homopolymer", default=None)
    args = ap.parse_args()

    # ---- grid ---------------------------------------------------------------
    grid = {k: list(v) for k, v in chg.GRID.items()}
    if args.test:
        grid = {"inner_cycles": [5], "cycles": [5], "max_relax": [0.2, 0.3],
                "local_homopolymer_threshold": [5]}
    if args.inner_cycles:
        grid["inner_cycles"] = [int(x) for x in args.inner_cycles.split(",")]
    if args.cycles:
        grid["cycles"] = [int(x) for x in args.cycles.split(",")]
    if args.max_relax:
        grid["max_relax"] = [float(x) for x in args.max_relax.split(",")]
    if args.homopolymer:
        grid["local_homopolymer_threshold"] = [int(x) for x in args.homopolymer.split(",")]
    cells = chg.build_cells(grid)

    # ---- sequences ----------------------------------------------------------
    seqs = load_human_sequences()
    if args.test:
        seqs = seqs.head(2)
    records = []
    for si, r in seqs.iterrows():
        try:
            prot, kind = chg.protein_from_cell(r["Protein"], "protein")
            records.append((si, r["Name"], prot, kind))
        except Exception as e:
            print(f"{r['Name']}: INPUT ERROR: {e}", flush=True)

    print(f"Trials file: {TRIALS.name}\nHuman sequences used ({len(records)}): "
          + ", ".join(f"{n} ({len(p)} aa)" for _, n, p, _ in records), flush=True)
    print(f"Codon table: {TABLE}", flush=True)

    # ---- codon table -> Codon Harmony profile (written once, read by each worker)
    nested_table = chg.load_dnachisel_table(str(TABLE))
    missing = sorted(chg.STANDARD_AAS - set(nested_table))
    if missing:
        print(f"WARNING: table missing amino-acid families {missing}")
    profile_json = chg.to_codon_harmony_profile(
        nested_table, os.path.join(tempfile.gettempdir(), f"_ch_profile_{SPECIES}.json"))
    domain = chg.DOMAIN[SPECIES]
    motifs = chg.PROKARYOTE_MOTIFS if domain == "prokaryote" else chg.EUKARYOTE_MOTIFS
    splice, start = domain == "eukaryote", domain == "prokaryote"

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    partial = Path(str(out) + ".partial.jsonl")

    # ---- resume -------------------------------------------------------------
    results = {}                                     # (si, ci) -> row
    if args.resume and partial.exists():
        for line in partial.read_text().splitlines():
            d = json.loads(line)
            results[(d["si"], d["ci"])] = d["row"]
        print(f"Resuming: {len(results)} finished runs loaded from {partial.name}", flush=True)
    elif partial.exists():
        partial.unlink()                             # fresh run: discard stale checkpoint

    todo = []
    for si, name, prot, kind in records:
        pending = [(ci, c) for ci, c in enumerate(cells) if (si, ci) not in results]
        for i in range(0, len(pending), CHUNK):
            todo.append((si, name, prot, kind, pending[i:i + CHUNK]))

    workers = max(1, min(args.workers, 61, max(1, len(todo))))
    total = len(records) * len(cells)
    print(f"{len(records)} sequences x {len(cells)} cells = {total} runs "
          f"({len(results)} already done), {workers} workers", flush=True)

    initargs = (profile_json, args.host_threshold, nested_table, motifs, splice, start,
                SPECIES, args.seed, domain)
    if todo:
        with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker,
                                 initargs=initargs) as pool, \
                tqdm(total=total, initial=len(results), unit="run", desc="Harmonizing") as bar, \
                open(partial, "a") as ck:
            futures = {pool.submit(run_task, name, prot, kind, items): (si, items)
                       for si, name, prot, kind, items in todo}
            for f in as_completed(futures):
                si, items = futures[f]
                for (ci, _), row in zip(items, f.result()):
                    results[(si, ci)] = row
                    ck.write(json.dumps({"si": si, "ci": ci, "row": row}, default=str) + "\n")
                ck.flush()
                bar.update(len(items))

    rows = [results[k] for k in sorted(results)]     # original order: sequence, then cell
    per_protein = chg.write_workbook(out, rows, cells, SPECIES)
    n_ok = sum(r["Status"] == "OK" for r in rows)
    print(f"Saved {len(rows)} rows ({n_ok} OK) to {out}")
    for p in per_protein:
        print(f"  {p[chg.NAME_COL]}: {p['Distinct sequences']} distinct of {p['Grid cells']} cells")
    if partial.exists():
        partial.unlink()


if __name__ == "__main__":
    main()