"""
Run_CHO.py  -  full CHO grid run, in parallel, with a progress bar.

Imports codon_grid_search.py (not modified). Reads the 80 CHO sequences from the
trials workbook, runs all 2,592 cartesian variants for each, and writes ONE Excel.

Usage (from the DNAChisel folder):
    python run_cho_parallel.py                        # full run
    python run_cho_parallel.py --test                 # tiny test: 2 seqs x 17 variants
"""

import argparse
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

import codon_grid_search as cgs          # the grid script, imported as-is

HERE = Path(__file__).resolve().parent
TRIALS = next(HERE.parent.glob("Codon*ptimization*trials*.xlsx"))   # spaces or underscores
TABLE = HERE.parent / "host_files" / "CHO" / "cho_codon.csv"
CHUNK = 25                               # variants per parallel task
# Default workers: 40 (leaves a few cores free on a 44-core box). Windows cannot
# run more than 61 worker processes, so the value is always capped at 61.
DEFAULT_WORKERS = min(40, max(1, (os.cpu_count() or 2) - 1))


def load_cho_sequences():
    """Trials sheet: headers on row 3, data from column B. Keep CHO rows."""
    d = pd.read_excel(TRIALS, header=2).iloc[:, 1:5]
    d.columns = ["Name", "Protein", "Host", "DNA"]
    d = d.dropna(subset=["Name", "Protein"])
    return d[d["Host"].astype(str).str.lower().str.contains("cho")].reset_index(drop=True)


def run_chunk(name, protein, chunk, seeds):
    """One parallel task: a slice of the variants for ONE sequence."""
    table = cgs.load_codon_table(str(TABLE))
    # KTX24h ends in "--" in the trials sheet = two stop codons (TGA TGA in its DNA),
    # not amino acids. Drop them for this one sequence (same as a trailing "*").
    if name == "KTX24h" and protein.endswith("--"):
        protein = protein[:-2]
    start, kind = cgs.to_start_dna(protein, "protein")
    rows = []
    for vname, params in chunk:
        t0 = time.time()
        seq, status, score, seed = cgs.run_variant(
            start, params, table, cgs.EUKARYOTE_MOTIFS, cgs.DEFAULT_SITES, seeds=seeds)
        row = {cgs.NAME_COL: name, "Variant": vname, "Host": "cho_codon",
               "Domain": "eukaryote", "Input kind": kind, "Status": status,
               "Seed used": seed, "Objective score": score,
               "Runtime (s)": round(time.time() - t0, 2)}
        row.update(cgs.flatten_params(params))
        if status.startswith("FAILED"):
            row["Optimized DNA Sequence"] = ""
        else:
            row.update(cgs.metrics_row(start, seq, table, cgs.EUKARYOTE_MOTIFS))
            row["Optimized DNA Sequence"] = seq
        rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", action="store_true", help="quick test (2 seqs, OFAT grid)")
    ap.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    ap.add_argument("--output", default=str(HERE / "outputs" / "cho_grid.xlsx"))
    args = ap.parse_args()

    seqs = load_cho_sequences()
    config = cgs.DEFAULT_CONFIG
    if args.test:
        seqs = seqs.head(2)
        combos = cgs.generate_ofat(config)
    else:
        combos = cgs.generate_cartesian(config, 0, 42)        # 0 = full grid
    variants = [(cgs.name_variant(i + 1, p, config["baseline"]), p)
                for i, p in enumerate(combos)]
    workers = max(1, min(args.workers, 61))
    # same seed logic as codon_grid_search.main(): seed 42 + 3 fallback seeds
    rng = np.random.RandomState(42)
    seeds = [42] + [int(x) for x in rng.randint(1, 10 ** 6, size=3)]

    chunks = [variants[i:i + CHUNK] for i in range(0, len(variants), CHUNK)]
    total = len(seqs) * len(variants)
    print(f"{len(seqs)} sequences x {len(variants)} variants = {total} runs, "
          f"{workers} workers", flush=True)

    results = {}
    with ProcessPoolExecutor(max_workers=workers) as pool, \
            tqdm(total=total, unit="variant", desc="Optimizing") as bar:
        futures = {}
        for si, r in seqs.iterrows():
            for ci, chunk in enumerate(chunks):
                f = pool.submit(run_chunk, r["Name"], r["Protein"], chunk, seeds)
                futures[f] = (si, ci, len(chunk))
        for f in as_completed(futures):
            si, ci, n = futures[f]
            results[(si, ci)] = f.result()
            bar.update(n)

    rows = [row for key in sorted(results) for row in results[key]]   # original order
    res = pd.DataFrame(rows)
    ok = res["Status"] == "OK"
    # Summary: identical to codon_grid_search.main() (same metrics, same order, same rounding)
    metric_cols = [c for c in (
        "Nucleotides changed", "% identity", "GC% global",
        f"GC% min ({cgs.GC_WINDOW_AUDIT}nt win)", f"GC% max ({cgs.GC_WINDOW_AUDIT}nt win)",
        "Longest homopolymer", f"Repeated {cgs.AUDIT_KMER}-mers", "CAI",
        "Codon-match score", f"5' Tm ({cgs.TM_AUDIT_WINDOW}nt)",
        f"3' Tm ({cgs.TM_AUDIT_WINDOW}nt)", "Forbidden motifs left", "Runtime (s)")
        if c in res.columns]
    summary = (res.assign(_ok=ok.astype(int)).groupby("Variant")
               .agg(n_sequences=("Status", "size"), n_passed=("_ok", "sum"),
                    **{c: (c, "mean") for c in metric_cols})
               .reset_index())
    summary[metric_cols] = summary[metric_cols].round(2)
    grid = pd.DataFrame([{"Variant": v, **cgs.flatten_params(p)} for v, p in variants])

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(out) as w:
        summary.to_excel(w, sheet_name="Summary", index=False)
        res.to_excel(w, sheet_name="Results", index=False)
        grid.to_excel(w, sheet_name="Grid", index=False)
    print(f"Saved {len(res)} rows ({ok.sum()} OK) to {out}")


if __name__ == "__main__":
    main()
