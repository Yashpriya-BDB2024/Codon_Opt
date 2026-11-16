"""
run_human_parallel_codontransformer.py - human CodonTransformer grid for the 8 human genes,
in parallel (one task per protein x (temperature, top_p) cell), with a progress bar.

Imports codontransformer_grid.py AS-IS (the lead's file is not edited) and calls its own
run_grid() for every cell. Generation, diversity and per-sequence metrics are entirely the
lead's code: nothing from DNA Chisel / Codon Harmony is used, and no sequence is altered.

The only things added here are the sheet/column layout that score_codontransformer.py reads
(like the other two workbooks): sheets Summary | Results | Grid, the DNA column named
"Optimized DNA Sequence", and Variant / Host / Domain labels. All metrics for the comparison
are computed later by the scoring script from the DNA.

Reproducibility: codontransformer_grid seeds torch before EVERY call, so a cell's draws do not
depend on worker count or finishing order.

Usage:
    python run_human_parallel_codontransformer.py --test          # 2 proteins x 2 cells
    python run_human_parallel_codontransformer.py                 # full: 8 x 16 cells x 5 draws
    python run_human_parallel_codontransformer.py --workers 40
"""
import argparse
import contextlib
import io
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from tqdm import tqdm

import codontransformer_grid as ctg          # the grid script, imported as-is

HERE = Path(__file__).resolve().parent
TRIALS = next(HERE.glob("Codon*ptimization*trials*.xlsx"))      # spaces or underscores
TABLE = HERE / "human_codon.csv"                                # same CSV as the other tools (CAI column of the grid only)
# TABLE = HERE.parent / "host_files" / "Human" / "human_codon.csv"
ORGANISM = "Homo sapiens"
# Default workers: 40 (leaves a few cores free on a 44-core box). Windows cannot run more
# than 61 worker processes, so the value is always capped at 61.
DEFAULT_WORKERS = min(40, max(1, (os.cpu_count() or 2) - 1))

_STATE = {}                                  # per-worker-process objects
STOPS = ("TAA", "TAG", "TGA")


def load_human_sequences():
    """Trials sheet: headers on row 3, data from column B. Keep Homo sapiens rows."""
    d = pd.read_excel(TRIALS, header=2).iloc[:, 1:5]
    d.columns = ["Name", "Protein", "Host", "DNA"]
    d = d.dropna(subset=["Name", "Protein"])
    return d[d["Host"].astype(str).str.lower().str.contains("homo")].reset_index(drop=True)


def _init_worker(model, tokenizer, organism, device, attention_type, num_sequences, seed,
                 table_path, threads):
    """Runs once per worker process: loads the model once and the codon table."""
    try:
        import torch
        if threads:
            torch.set_num_threads(threads)   # avoid 40 procs x all cores oversubscription
    except Exception:
        pass
    predict_fn, _ = ctg.build_predict_fn(model, tokenizer, organism, device,
                                         attention_type, num_sequences, seed)
    _STATE.update(predict_fn=predict_fn, organism=organism,
                  table=ctg.load_codon_table(table_path))


def run_cell(si, ci, name, protein, kind, combo):
    """One task: ONE protein x ONE (temperature, top_p) cell, via the lead's run_grid()."""
    s = _STATE
    with contextlib.redirect_stdout(io.StringIO()):          # keep the progress bar clean
        seq_rows, combo_rows = ctg.run_grid([(name, protein, kind)], s["predict_fn"], [combo],
                                            s["table"], ctg.EUKARYOTE_MOTIFS, s["organism"])
    return si, ci, seq_rows, combo_rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", action="store_true", help="quick test (2 proteins x 2 cells)")
    ap.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    ap.add_argument("--threads-per-worker", type=int, default=1,
                    help="torch threads per worker (1 is right for many CPU workers)")
    ap.add_argument("--output", default=str(HERE / "outputs" / "human_codontransformer_grid.xlsx"))
    ap.add_argument("--organism", default=ORGANISM)
    ap.add_argument("--model", default="adibvafa/CodonTransformer")
    ap.add_argument("--tokenizer", default="adibvafa/CodonTransformer")
    ap.add_argument("--temperatures", default="0.2,0.4,0.6,0.8")
    ap.add_argument("--top-ps", default="0.8,0.9,0.95,1.0")
    ap.add_argument("--num-sequences", type=int, default=ctg.DEFAULT_NUM_SEQUENCES)
    ap.add_argument("--attention-type", default="original_full",
                    choices=["original_full", "block_sparse"])
    ap.add_argument("--device", default=None, help="cuda / cpu (default: auto)")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    temps, tops = ctg.parse_floats(args.temperatures), ctg.parse_floats(args.top_ps)
    seqs = load_human_sequences()
    if args.test:
        temps, tops, seqs = temps[:1], tops[:2], seqs.head(2)
    combos = [(t, p) for t in temps for p in tops]

    records = []
    for si, r in seqs.iterrows():
        try:
            prot, kind = ctg.protein_from_cell(r["Protein"], "protein")
            records.append((si, r["Name"], prot, kind))
        except Exception as e:
            print(f"{r['Name']}: INPUT ERROR: {e}", flush=True)
    print(f"Trials file: {TRIALS.name}\nHuman sequences ({len(records)}): "
          + ", ".join(f"{n} ({len(p)} aa)" for _, n, p, _ in records), flush=True)
    print(f"Codon table (CAI column of the grid only): {TABLE}", flush=True)

    tasks = [(si, ci, n, p, k, c) for si, n, p, k in records for ci, c in enumerate(combos)]
    workers = max(1, min(args.workers, 61, len(tasks)))
    print(f"{len(records)} proteins x {len(combos)} cells = {len(tasks)} tasks "
          f"({args.num_sequences} draws each), {workers} workers, organism={args.organism}",
          flush=True)
    if workers > 4:
        try:
            import torch
            if (args.device or ("cuda" if torch.cuda.is_available() else "cpu")).startswith("cuda"):
                print(f"NOTE: {workers} workers will each load the model onto the GPU; "
                      f"on one GPU a few workers (1-4) are usually as fast.", flush=True)
        except Exception:
            pass

    initargs = (args.model, args.tokenizer, args.organism, args.device, args.attention_type,
                args.num_sequences, args.seed, str(TABLE), args.threads_per_worker)
    results = {}
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker,
                             initargs=initargs) as pool, \
            tqdm(total=len(tasks), unit="cell", desc="Generating") as bar:
        futures = [pool.submit(run_cell, *t) for t in tasks]
        for f in as_completed(futures):
            si, ci, seq_rows, combo_rows = f.result()
            results[(si, ci)] = (seq_rows, combo_rows)
            bar.update(1)

    seq_rows = [r for k in sorted(results) for r in results[k][0]]       # protein, then cell
    combo_rows = [r for k in sorted(results) for r in results[k][1]]
    for r in seq_rows:                                   # labels only; Name + Variant is unique
        r["Variant"] = f"T{r['Temperature']:g}_p{r['top_p']:g}_d{r['Draw']}"
        r["Host"], r["Domain"] = "human_codon", "eukaryote"

    res = pd.DataFrame(seq_rows).rename(columns={"Generated DNA Sequence": "Optimized DNA Sequence"})
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(out) as w:
        pd.DataFrame(combo_rows).to_excel(w, sheet_name="Summary", index=False)
        res.to_excel(w, sheet_name="Results", index=False)
        pd.DataFrame([{"Temperature": t, "top_p": p} for t, p in combos]).to_excel(
            w, sheet_name="Grid", index=False)
    dna = res["Optimized DNA Sequence"].dropna().astype(str)
    n_ok = int((res["Status"] == "OK").sum())
    n_stop = int(sum(len(s) >= 3 and s[-3:] in STOPS for s in dna))
    print(f"Saved {len(res)} rows ({n_ok} OK) to {out}")
    print(f"Sequences ending in a stop codon (left untouched; the scoring script's "
          f"--stop-policy handles it): {n_stop} of {len(dna)}")


if __name__ == "__main__":
    main()
