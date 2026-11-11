# -*- coding: utf-8 -*-
"""
score_codonharmony.py  -  score table: every Status-OK Codon Harmony sequence + the GeneArt sequences
====================================================================================================

Rows    : every grid row with Status == "OK" (Results sheet of human_codonharmony_grid.xlsx)
          + the GeneArt sequence of each gene (trials workbook) = the benchmark baseline.
          No DNA Chisel rows: the optimisers are not mixed.
Columns : Name, Source ("CodonHarmony" or "GeneArt"), Variant (the Codon Harmony cell, e.g.
          C017_ic5_c10_rel0.2_hp5, or "GeneArt") + the 22 metrics of codon_compare.COMPARE_METRICS
          (without Kozak) - the layout of score_all_sequences.py.

SCORING IS UNCHANGED. codon_features.py, codon_compare.py and score_all_sequences.py are imported, not
modified; every metric comes from score_all_sequences.score_one -> codon_features.compute_features().

Terminal stop codon (input handling only, not scoring)
------------------------------------------------------
compute_features strips a trailing stop codon for the codon-level metrics (CAI, tAI, ENC ...) but computes
the nucleotide-level metrics (GC, CpG, structure, hairpins, repeats, restriction/splice counts, Tm
fragments) on the DNA exactly as given. Codon Harmony writes no stop codon; a GeneArt sequence may end in
one. To make every sequence cover the SAME nucleotides, --stop-policy decides what is handed to the scorer:
    strip (default)  remove one terminal stop codon (if present) from every sequence of both sources.
                     A no-op when no sequence has one. After this, every sequence must be exactly
                     3 x len(protein) nt, which is checked and stops the run if violated.
    asis             hand the DNA to the scorer exactly as it is (what score_all_sequences.py does);
                     a warning is printed if the two sources differ in stop-codon handling.
The number of stops found/removed per source is always printed.

Run from Sequence_scoring/ (pip install numpy pandas openpyxl seqfold tqdm):
    python score_codonharmony.py --grid ../CodonHarmony/outputs/human_codonharmony_grid.xlsx
    python score_codonharmony.py --workers 40        # default: min(40, cores-1), max 61

Outputs (outputs/scores/): scores_codonharmony.csv and scores_codonharmony.xlsx (sheets CodonHarmony, GeneArt).
The DNA Chisel files (scores_all.*) are not touched.
"""

from __future__ import annotations

import argparse
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

import codon_features as cf              # given, unchanged
import score_all_sequences as sas        # given, unchanged: score_one, METRICS, translate, load_geneart, paths

HERE = Path(__file__).resolve().parent
DNA = sas.DNA                            # "Optimized DNA Sequence" - same column name in the CH workbook
METRICS = sas.METRICS
SOURCE = "CodonHarmony"
STOPS = ("TAA", "TAG", "TGA")


def has_stop(s):
    return len(s) >= 3 and s[-3:] in STOPS


def main():
    ap = argparse.ArgumentParser(description="Score table: Status-OK Codon Harmony sequences + GeneArt.")
    ap.add_argument("--grid", default=str(sas.ROOT / "CodonHarmony" / "outputs" / "human_codonharmony_grid.xlsx"))
    ap.add_argument("--trials", default=None, help="the ORIGINAL trials workbook (GeneArt sequences + input proteins)")
    ap.add_argument("--out-dir", default=str(HERE / "outputs" / "scores"))
    ap.add_argument("--workers", type=int, default=min(40, max(1, (os.cpu_count() or 2) - 1)))
    ap.add_argument("--stop-policy", choices=["strip", "asis"], default="strip",
                    help="terminal stop codon handling before scoring (see module docstring)")
    args = ap.parse_args()

    try:
        import seqfold  # noqa: F401
    except Exception:
        warnings.warn("seqfold is not installed: TIS_accessibility_5p, MFE_min_window and "
                      "stable_structure_count will be empty. Run: pip install seqfold")

    res = pd.read_excel(args.grid, sheet_name="Results")
    n_all = len(res)
    res = res[res["Status"] == "OK"].copy()
    res[DNA] = res[DNA].fillna("").astype(str).map(cf.clean_seq)
    res = res[res[DNA] != ""].reset_index(drop=True)
    print("Grid rows: {} | Status OK: {}".format(n_all, len(res)))

    ga = sas.load_geneart(args.trials)            # exits with a clear message if the workbook is not found
    prot = dict(zip(ga["Name"], ga["Protein"]))
    missing = sorted(set(res["Name"]) - set(prot))
    if missing:
        sys.exit("Genes in the grid but not in the trials workbook: {}".format(missing))

    # ---- terminal stop codon: report, then apply the policy (input handling only) -------------------
    n_ch, n_ga = sum(map(has_stop, res[DNA])), sum(map(has_stop, ga["GeneArt"]))
    print("Terminal stop codon found: CodonHarmony {}/{} | GeneArt {}/{}".format(n_ch, len(res), n_ga, len(ga)))
    if args.stop_policy == "strip":
        res["_dna"] = res[DNA].map(cf.strip_trailing_stop)
        ga["_dna"] = ga["GeneArt"].map(cf.strip_trailing_stop)
        print("stop-policy strip: removed {} (CodonHarmony) + {} (GeneArt) terminal stop codon(s); "
              "scorer input = coding sequence only, for every sequence.".format(n_ch, n_ga))
    else:
        res["_dna"], ga["_dna"] = res[DNA], ga["GeneArt"]
        if (n_ch == 0) != (n_ga == 0) or (0 < n_ch < len(res)) or (0 < n_ga < len(ga)):
            warnings.warn("stop-policy asis: stop codons are not handled the same way in the two sources; "
                          "nucleotide-level metrics cover 3 nt more on sequences that carry one.")

    # ---- safeguards: translation, and (strip policy) identical nucleotide span ---------------------
    bad = [n for n, s in zip(res["Name"], res["_dna"]) if sas.translate(cf.strip_trailing_stop(s)) != prot[n]]
    bad += [n for n, s in zip(ga["Name"], ga["_dna"]) if sas.translate(cf.strip_trailing_stop(s)) != prot[n]]
    if bad:
        sys.exit("{} sequence(s) do not translate to the input protein, e.g. {}".format(len(bad), bad[:5]))
    print("Translation check: all {} sequences encode their input protein.".format(len(res) + len(ga)))
    if args.stop_policy == "strip":
        off = [(n, len(s), 3 * len(prot[n])) for n, s in list(zip(res["Name"], res["_dna"])) +
               list(zip(ga["Name"], ga["_dna"])) if len(s) != 3 * len(prot[n])]
        if off:
            sys.exit("{} sequence(s) are not exactly 3 x protein length after stop handling "
                     "(name, nt, expected), e.g. {}".format(len(off), off[:5]))
        print("Span check: every sequence is exactly 3 x protein length - all scored over the same nucleotides.")

    seqs = list(dict.fromkeys(list(res["_dna"]) + list(ga["_dna"])))
    print("Status OK rows: {} | GeneArt: {} | distinct DNA to score: {}".format(len(res), len(ga), len(seqs)))

    workers = max(1, min(args.workers, 61))
    try:
        from tqdm import tqdm
    except Exception:
        tqdm = lambda x, **k: x                                    # noqa: E731
    scores = {}
    if workers == 1:
        for s in tqdm(seqs, desc="Scoring", unit="seq"):
            scores[s] = sas.score_one(s)[1]
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for s, sc in tqdm(ex.map(sas.score_one, seqs, chunksize=2), total=len(seqs), desc="Scoring", unit="seq"):
                scores[s] = sc

    grid = pd.concat([res[["Name", "Variant"]].assign(Source=SOURCE),
                      pd.DataFrame([scores[s] for s in res["_dna"]], columns=METRICS)], axis=1)
    geneart = pd.concat([ga[["Name"]].assign(Variant="GeneArt", Source="GeneArt"),
                         pd.DataFrame([scores[s] for s in ga["_dna"]], columns=METRICS)], axis=1)
    cols = ["Name", "Source", "Variant"] + METRICS
    grid, geneart = grid[cols], geneart[cols]

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    pd.concat([geneart, grid], ignore_index=True).to_csv(out / "scores_codonharmony.csv", index=False)
    with pd.ExcelWriter(out / "scores_codonharmony.xlsx", engine="openpyxl") as w:
        grid.to_excel(w, sheet_name="CodonHarmony", index=False)
        geneart.to_excel(w, sheet_name="GeneArt", index=False)
        for ws in w.book.worksheets:
            ws.freeze_panes = "D2"
            for col in ws.columns:
                ws.column_dimensions[col[0].column_letter].width = max(9, len(str(col[0].value)) + 2)
    print("Rows per gene (CodonHarmony, Status OK):")
    print(grid.groupby("Name").size().to_string())
    print("wrote {} and scores_codonharmony.xlsx".format(out / "scores_codonharmony.csv"))


if __name__ == "__main__":
    main()
