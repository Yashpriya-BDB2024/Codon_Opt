# -*- coding: utf-8 -*-
"""
score_all_sequences.py  -  score table: every Status-OK DNA Chisel sequence + the GeneArt sequences
==================================================================================================

Rows    : every grid row with Status == "OK" (Results sheet of the grid output)
          + the GeneArt sequence of each gene ("DNA optimized sequence" in the trials workbook)
Columns : Name, Source, Variant + the 22 metrics of codon_compare.COMPARE_METRICS
          (the comparison panel, without Kozak). See README_score_table.md for what each one means.

codon_features.py and codon_compare.py are imported, NOT modified. Every value is taken from
codon_features.compute_features()["scores"], with the same human_codon.csv / human_tgcn.csv for every
sequence. Identical DNA gives identical scores, so each distinct DNA is scored once and the result is
copied to every row that has that DNA (same numbers, less run time).

Run from Sequence_scoring/ (pip install numpy pandas openpyxl seqfold tqdm):
    python score_all_sequences.py --grid ../DNAChisel/outputs/human_grid_1.xlsx
    python score_all_sequences.py --workers 40        # default: min(40, cores-1), max 61

Outputs (outputs/scores/): scores_all.csv  (all rows)   scores_all.xlsx  (sheets DNAChisel, GeneArt)
"""

from __future__ import annotations

import argparse
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

import codon_features as cf          # given, unchanged
import codon_compare as cc           # given, unchanged (only for its list of compared metrics)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent if (HERE.parent / "host_files").exists() else HERE
CODON_CSV = ROOT / "host_files" / "Human" / "human_codon.csv"
TGCN_CSV = ROOT / "host_files" / "Human" / "human_tgcn.csv"
DNA = "Optimized DNA Sequence"
METRICS = [m for m, _ in cc.COMPARE_METRICS if m != "kozak_score"]      # 22 metrics, Kozak excluded


def translate(seq):
    return "".join(cf.GENETIC_CODE.get(seq[i:i + 3], "X") for i in range(0, len(seq) - len(seq) % 3, 3))


def load_geneart(trials=None):
    """Trials workbook: headers on row 3, data from column B = Name, protein, host, GeneArt DNA."""
    if trials is None:
        found = sorted(ROOT.glob("Codon*ptimization*trials*.xlsx"))
        if len(found) != 1:
            sys.exit("Found {} trials workbooks in {}; pass the original one with --trials".format(len(found), ROOT))
        trials = found[0]
    d = pd.read_excel(trials, header=2).iloc[:, 1:5]
    d.columns = ["Name", "Protein", "Host", "GeneArt"]
    d = d.dropna(subset=["Name", "Protein"])
    d = d[d["Host"].astype(str).str.lower().str.contains("homo")].reset_index(drop=True)
    d["Protein"] = d["Protein"].astype(str).str.strip().str.upper()
    d["GeneArt"] = d["GeneArt"].map(cf.clean_seq)
    return d


def score_one(seq):
    r = cf.compute_features(seq, host="human", host_label="human",
                            codon_table=str(CODON_CSV), tgcn=str(TGCN_CSV))
    return seq, {m: r["scores"].get(m) for m in METRICS}


def main():
    ap = argparse.ArgumentParser(description="Score table: Status-OK DNA Chisel sequences + GeneArt.")
    ap.add_argument("--grid", default=str(ROOT / "DNAChisel" / "outputs" / "human_grid_1.xlsx"))
    ap.add_argument("--trials", default=None, help="the ORIGINAL trials workbook (GeneArt sequences)")
    ap.add_argument("--out-dir", default=str(HERE / "outputs" / "scores"))
    ap.add_argument("--workers", type=int, default=min(40, max(1, (os.cpu_count() or 2) - 1)))
    args = ap.parse_args()

    try:
        import seqfold  # noqa: F401
    except Exception:
        warnings.warn("seqfold is not installed: TIS_accessibility_5p, MFE_min_window and "
                      "stable_structure_count will be empty. Run: pip install seqfold")

    res = pd.read_excel(args.grid, sheet_name="Results")
    res = res[res["Status"] == "OK"].copy()
    res[DNA] = res[DNA].fillna("").astype(str).map(cf.clean_seq)
    res = res[res[DNA] != ""].reset_index(drop=True)
    ga = load_geneart(args.trials)
    prot = dict(zip(ga["Name"], ga["Protein"]))
    missing = sorted(set(res["Name"]) - set(prot))
    if missing:
        sys.exit("Genes in the grid but not in the trials workbook: {}".format(missing))

    # safeguard: every sequence must encode its input protein
    bad = [n for n, s in zip(res["Name"], res[DNA]) if translate(cf.strip_trailing_stop(s)) != prot[n]]
    bad += [n for n, s in zip(ga["Name"], ga["GeneArt"]) if translate(cf.strip_trailing_stop(s)) != prot[n]]
    if bad:
        sys.exit("{} sequence(s) do not translate to the input protein, e.g. {}".format(len(bad), bad[:5]))

    seqs = list(dict.fromkeys(list(res[DNA]) + list(ga["GeneArt"])))
    print("Status OK rows: {} | GeneArt: {} | distinct DNA to score: {}".format(len(res), len(ga), len(seqs)))

    workers = max(1, min(args.workers, 61))
    try:
        from tqdm import tqdm
    except Exception:
        tqdm = lambda x, **k: x                                    # noqa: E731
    scores = {}
    if workers == 1:
        for s in tqdm(seqs, desc="Scoring", unit="seq"):
            scores[s] = score_one(s)[1]
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for s, sc in tqdm(ex.map(score_one, seqs, chunksize=4), total=len(seqs), desc="Scoring", unit="seq"):
                scores[s] = sc

    grid = pd.concat([res[["Name", "Variant"]].assign(Source="DNAChisel"),
                      pd.DataFrame([scores[s] for s in res[DNA]], columns=METRICS)], axis=1)
    geneart = pd.concat([ga[["Name"]].assign(Variant="GeneArt", Source="GeneArt"),
                         pd.DataFrame([scores[s] for s in ga["GeneArt"]], columns=METRICS)], axis=1)
    cols = ["Name", "Source", "Variant"] + METRICS
    grid, geneart = grid[cols], geneart[cols]

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    pd.concat([geneart, grid], ignore_index=True).to_csv(out / "scores_all.csv", index=False)
    with pd.ExcelWriter(out / "scores_all.xlsx", engine="openpyxl") as w:
        grid.to_excel(w, sheet_name="DNAChisel", index=False)
        geneart.to_excel(w, sheet_name="GeneArt", index=False)
        for ws in w.book.worksheets:
            ws.freeze_panes = "D2"
            for col in ws.columns:
                ws.column_dimensions[col[0].column_letter].width = max(9, len(str(col[0].value)) + 2)
    print("Rows per gene (DNAChisel, Status OK):")
    print(grid.groupby("Name").size().to_string())
    print("wrote {} and scores_all.xlsx".format(out / "scores_all.csv"))


if __name__ == "__main__":
    main()
