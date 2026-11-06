# -*- coding: utf-8 -*-
"""
score_ncbi_sequences.py  -  scores for the NCBI sequences only (same table layout as scores_all.xlsx)
=====================================================================================================

Imports score_all_sequences.py and re-uses its scorer unchanged (score_one(), METRICS, the same
human_codon.csv / human_tgcn.csv), so the NCBI values are directly comparable with scores_all.xlsx.

Input : the workbook written by remove_ncbi_stop.py (NCBI sequences without the stop codon)
Output: outputs/scores/scores_ncbi.xlsx (sheet NCBI) and scores_ncbi.csv
        columns: Name, Source (= NCBI), Variant (= NCBI) + the 22 metrics

Each NCBI sequence is checked first: whole codons, no stop codon, translates exactly to the protein of
its row. A sequence that fails is NOT scored; it is listed and skipped.

Run (from Sequence_scoring/, after remove_ncbi_stop.py; seqfold installed):
    python score_ncbi_sequences.py --trials trials_NCBI_stop_removed.xlsx
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import pandas as pd

import score_all_sequences as sas      # given scorer, unchanged

cf = sas.cf


def main():
    ap = argparse.ArgumentParser(description="Score the NCBI sequences only.")
    ap.add_argument("--trials", default=str(sas.ROOT / "trials_NCBI_stop_removed.xlsx"),
                    help="workbook written by remove_ncbi_stop.py")
    ap.add_argument("--out-dir", default=str(sas.HERE / "outputs" / "scores"))
    args = ap.parse_args()

    try:
        import seqfold  # noqa: F401
    except Exception:
        warnings.warn("seqfold is not installed: TIS_accessibility_5p, MFE_min_window and "
                      "stable_structure_count will be empty. Run: pip install seqfold")

    d = pd.read_excel(args.trials, header=2).iloc[:, 1:6]
    d.columns = ["Name", "Protein", "Host", "GeneArt", "NCBI"]
    d = d.dropna(subset=["Name", "Protein", "NCBI"])
    d = d[d["Host"].astype(str).str.lower().str.contains("homo")].reset_index(drop=True)

    rows, skipped = [], []
    for _, r in d.iterrows():
        seq = cf.clean_seq(r["NCBI"])
        prot = str(r["Protein"]).strip().upper()
        if len(seq) % 3:
            skipped.append((r["Name"], "length {} is not a multiple of 3".format(len(seq))))
            continue
        t = sas.translate(seq)
        if "*" in t:
            skipped.append((r["Name"], "contains a stop codon: run remove_ncbi_stop.py first"))
        elif t != prot:
            skipped.append((r["Name"], "does not translate to the protein of its row"))
        else:
            rows.append({"Name": r["Name"], "Source": "NCBI", "Variant": "NCBI", **sas.score_one(seq)[1]})

    for name, why in skipped:
        print("SKIPPED {}: {}".format(name, why))
    if not rows:
        sys.exit("No NCBI sequence passed the checks.")
    out = pd.DataFrame(rows, columns=["Name", "Source", "Variant"] + sas.METRICS)

    o = Path(args.out_dir)
    o.mkdir(parents=True, exist_ok=True)
    out.to_csv(o / "scores_ncbi.csv", index=False)
    with pd.ExcelWriter(o / "scores_ncbi.xlsx", engine="openpyxl") as w:
        out.to_excel(w, sheet_name="NCBI", index=False)
        ws = w.book["NCBI"]
        ws.freeze_panes = "D2"
        for col in ws.columns:
            ws.column_dimensions[col[0].column_letter].width = max(9, len(str(col[0].value)) + 2)
    print("scored {} NCBI sequence(s): {}".format(len(out), ", ".join(out["Name"])))
    print("wrote {}".format(o / "scores_ncbi.xlsx"))


if __name__ == "__main__":
    main()
