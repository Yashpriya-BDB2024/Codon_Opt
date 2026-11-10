"""
check_restriction_sites.py - diagnostic: WHY do Codon Harmony outputs still contain
restriction sites? Compares, for every OK sequence in the result workbook:
  (a) plain-text count of each enzyme's recognition site (both strands), and
  (b) Bio.Restriction's own detection (the library Codon Harmony relies on).
Reading the result:
  - (b) finds ~nothing while (a) finds many  -> Codon Harmony cannot SEE those sites
    (detection problem: Type IIS handling / Biopython version). Not an optimisation failure.
  - (b) finds about the same as (a)          -> it sees them but fails to remove them
    (the tool simply does not guarantee removal).
Usage:  python check_restriction_sites.py --input outputs/human_codonharmony_grid.xlsx
"""
import argparse
import re

import pandas as pd
import Bio.Restriction as BR
from Bio.Restriction import RestrictionBatch, Analysis
from Bio.Seq import Seq

try:
    import codonharmony_grid as chg
except ImportError:
    import codonharmony_grid_new as chg

ap = argparse.ArgumentParser()
ap.add_argument("--input", required=True)
ap.add_argument("--max-seqs", type=int, default=100, help="sample size (speed)")
a = ap.parse_args()

res = pd.read_excel(a.input, sheet_name="Results")
seqs = res[res["Status"] == "OK"]["Optimized DNA Sequence"].dropna().tolist()[:a.max_seqs]
names = chg.DEFAULT_RESTRICTION_ENZYMES
rb = RestrictionBatch(names)
rcomp = lambda s: s[::-1].translate(str.maketrans("ACGT", "TGCA"))

rows = {n: {"enzyme": n, "site": getattr(BR, n).site, "regex_sites": 0, "bio_detects": 0} for n in names}
for dna in seqs:
    found = Analysis(rb, Seq(dna), linear=True).full()
    for enz, cuts in found.items():
        rows[str(enz)]["bio_detects"] += len(cuts)
    for n in names:
        site = rows[n]["site"]
        rows[n]["regex_sites"] += len(re.findall(f"(?={site})", dna)) + \
            (len(re.findall(f"(?={rcomp(site)})", dna)) if rcomp(site) != site else 0)
out = pd.DataFrame(rows.values()).sort_values("regex_sites", ascending=False)
print(f"{len(seqs)} OK sequences checked\n")
print(out[out.regex_sites > 0].to_string(index=False))
print(f"\nTOTAL  text-search sites: {out.regex_sites.sum()}   Bio.Restriction detects: {out.bio_detects.sum()}")
