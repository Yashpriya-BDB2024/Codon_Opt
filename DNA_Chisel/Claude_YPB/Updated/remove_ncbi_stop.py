# -*- coding: utf-8 -*-
"""
remove_ncbi_stop.py  -  remove the terminal stop codon from the NCBI sequences and verify the protein
=====================================================================================================

Reads the trials workbook, and for every row that has an NCBI sequence:
  1. the sequence must be a whole number of codons;
  2. if its LAST codon is TAA / TAG / TGA, that one codon is removed (nothing else is touched);
  3. VERIFY: the remaining sequence is a whole number of codons, has no stop codon, and translates to
     exactly the amino-acid sequence in the same row.
Only rows that pass are written. A row that fails is left exactly as it was and is reported as FAILED.

The GeneArt column ("DNA optimized sequence") and every other cell are never changed: the script
re-opens the input and the output and checks that the ONLY cells that differ are the NCBI cells of the
rows reported as "stop removed". The input file is not modified; a new workbook is written.

Run (from Sequence_scoring/):
    python remove_ncbi_stop.py --input "../Codon_optimization_trials_27102026_1.xlsx"
    -> trials_NCBI_stop_removed.xlsx (same folder as the input) + a printed report
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import openpyxl

import codon_features as cf          # given, unchanged (genetic code only)

STOPS = {"TAA", "TAG", "TGA"}


def translate(seq):
    return "".join(cf.GENETIC_CODE[seq[i:i + 3]] for i in range(0, len(seq), 3))


def find_columns(ws):
    """Locate the header row and the Name / protein / GeneArt / NCBI columns by their header text."""
    want = {"name": "Name", "amino acid sequence": "Protein", "dna optimized sequence": "GeneArt",
            "ncbi sequence": "NCBI"}
    for row in ws.iter_rows(min_row=1, max_row=15):
        found = {want[str(c.value).strip().lower()]: c.column for c in row
                 if c.value and str(c.value).strip().lower() in want}
        if "NCBI" in found:
            missing = [v for v in want.values() if v not in found]
            if missing:
                sys.exit("Header row {} is missing columns: {}".format(row[0].row, missing))
            return row[0].row, found
    sys.exit("No header row with 'NCBI sequence' found in the first 15 rows.")


def check_and_trim(raw, protein):
    """Returns (status, new_sequence, note). status: 'stop removed' | 'no stop, unchanged' | 'FAILED'."""
    seq = "".join(str(raw).split()).upper()
    protein = "".join(str(protein).split()).upper()
    if set(seq) - set("ACGT"):
        return "FAILED", None, "non-ACGT characters: {}".format(sorted(set(seq) - set("ACGT")))
    if len(seq) % 3:
        msg = "length {} is not a multiple of 3 (expected {} nt for the protein, {} with stop)".format(
            len(seq), 3 * len(protein), 3 * len(protein) + 3)
        part = translate(seq[:len(seq) - len(seq) % 3])
        if protein.startswith(part):
            msg += "; its {} codons match the start of the protein exactly, so the sequence is cut short at the 3' end ({} nt missing)".format(
                len(part), 3 * len(protein) - len(seq))
        return "FAILED", None, msg
    removed = seq[-3:] in STOPS
    new = seq[:-3] if removed else seq
    if len(new) % 3 or not new:
        return "FAILED", None, "empty or partial codon after trimming"
    prot = translate(new)
    if "*" in prot:
        return "FAILED", None, "stop codon remains inside the sequence (position {})".format(prot.index("*") + 1)
    if prot != protein:
        n = next((i for i, (a, b) in enumerate(zip(prot, protein)) if a != b), min(len(prot), len(protein)))
        return "FAILED", None, ("translation differs from the protein: {} vs {} residues, first difference at residue {}"
                                .format(len(prot), len(protein), n + 1))
    return ("stop removed" if removed else "no stop, unchanged"), new, "{} nt -> {} nt".format(len(seq), len(new))


def main():
    ap = argparse.ArgumentParser(description="Remove the terminal stop codon from the NCBI sequences and verify the protein.")
    ap.add_argument("--input", required=True, help="trials workbook with the NCBI sequence column")
    ap.add_argument("--output", default=None, help="default: trials_NCBI_stop_removed.xlsx next to the input")
    args = ap.parse_args()
    src = Path(args.input)
    out = Path(args.output) if args.output else src.with_name("trials_NCBI_stop_removed.xlsx")

    wb = openpyxl.load_workbook(src)
    ws = wb.active
    hdr, col = find_columns(ws)
    results, changed = [], set()
    for r in range(hdr + 1, ws.max_row + 1):
        name, prot, ncbi = (ws.cell(r, col[k]).value for k in ("Name", "Protein", "NCBI"))
        if not ncbi or not str(ncbi).strip():
            continue
        status, new, note = check_and_trim(ncbi, prot)
        if status == "stop removed":
            ws.cell(r, col["NCBI"]).value = new
            changed.add(r)
        results.append((name, status, note))
    wb.save(out)

    # integrity: only the NCBI cells of the 'stop removed' rows may differ from the input
    a, b = openpyxl.load_workbook(src).active, openpyxl.load_workbook(out).active
    diffs = [(c.row, c.column) for row in a.iter_rows() for c in row if c.value != b.cell(c.row, c.column).value]
    allowed = {(r, col["NCBI"]) for r in changed}
    if set(diffs) != allowed or a.max_row != b.max_row or a.max_column != b.max_column:
        sys.exit("INTEGRITY CHECK FAILED: cells other than the trimmed NCBI cells differ. Output not trusted.")
    ga_same = all(a.cell(r, col["GeneArt"]).value == b.cell(r, col["GeneArt"]).value for r in range(1, a.max_row + 1))

    print("{:<10} {:<20} {}".format("Gene", "Status", "Detail"))
    for name, status, note in results:
        print("{:<10} {:<20} {}".format(str(name), status, note))
    n_fail = sum(s == "FAILED" for _, s, _ in results)
    print("\nstop removed: {} | no stop, unchanged: {} | FAILED: {}".format(
        sum(s == "stop removed" for _, s, _ in results), sum(s == "no stop, unchanged" for _, s, _ in results), n_fail))
    print("GeneArt column identical to the input: {} | cells changed in total: {} (all NCBI)".format(ga_same, len(diffs)))
    print("wrote {}".format(out))
    if n_fail:
        print("FAILED rows were left unchanged; fix the NCBI sequence in the workbook and run again.")
        sys.exit(1)


if __name__ == "__main__":
    main()
