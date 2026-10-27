"""
dnachisel_variants.py

Runs 10 DNA Chisel variants on the "Generated DNA Sequence" column of
human_codon_sequences.xlsx (made by main.py / reverse_translation.py).
GeneArt is NOT used here; this script only produces the DNA Chisel sequences.

Usage:
    python dnachisel_variants.py --input human_codon_sequences.xlsx --output dnachisel_variants.xlsx
    python dnachisel_variants.py --input human_codon_sequences.xlsx --output test.xlsx --limit 1   (quick test)
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from dnachisel import (
    AvoidPattern,
    AvoidRareCodons,
    CodonOptimize,
    DnaOptimizationProblem,
    EnforceGCContent,
    EnforceTranslation,
    UniquifyAllKmers,
)

NAME_COL = "Name"
INPUT_COL = "Generated DNA Sequence"
SPECIES = "h_sapiens"

# ----------------------------------------------------------------------------
# Cutoffs. "geneart" level = values read from the part-2 report.
# "lenient" and "strict" are looser / tighter versions to see how cutoffs matter.
# ----------------------------------------------------------------------------
MOTIFS = ["AWTAAA", "ATTTA", "TATAAA"]         # polyA signal (W = A/T), AU-rich, TATA
SITES = ["BsaI", "BsmBI", "BbsI", "SapI", "EcoRI", "XhoI", "NotI",
         "XbaI", "SpeI", "KpnI", "SalI", "AscI", "PacI"]   # absent in GeneArt
GC_TARGET = (50, 0.57)        # (window, target) - overall GeneArt GC is ~57 %

LEVELS = {
    # looser than GeneArt
    "lenient": {"gc": (50, 0.40, 0.72), "runs": ["7xA", "6xC", "6xG", "6xT"],
                "kmer": 15, "rare": 0.05, "motifs": True, "sites": True},
    # exactly what GeneArt stayed within (GC 44-68 %, max runs A5/C4/G4/T4,
    # no repeated 10-mer, lowest codon frequency used = 0.11)
    "geneart": {"gc": (50, 0.44, 0.68), "runs": ["6xA", "5xC", "5xG", "5xT"],
                "kmer": 10, "rare": 0.11, "motifs": True, "sites": True},
    # tighter than GeneArt (k=9 not tried: GeneArt itself repeats 9-mers)
    "strict": {"gc": (50, 0.46, 0.66), "runs": ["5xA", "4xC", "4xG", "4xT"],
               "kmer": 10, "rare": 0.15, "motifs": True, "sites": True},
}

# ----------------------------------------------------------------------------
# The 10 variants. EnforceTranslation is always on (protein must not change).
# "codon" = CodonOptimize method (None = no CodonOptimize objective).
# ----------------------------------------------------------------------------
VARIANTS = {
    # 1-2: codon strategy only (baselines, no other constraints)
    "V01_best_only": {"codon": "use_best_codon"},
    "V02_match_only": {"codon": "match_codon_usage"},

    # 3-8: both codon strategies x three cutoff levels
    "V03_best_lenient": {"codon": "use_best_codon", **LEVELS["lenient"]},
    "V04_match_lenient": {"codon": "match_codon_usage", **LEVELS["lenient"]},
    "V05_best_geneart": {"codon": "use_best_codon", **LEVELS["geneart"]},
    "V06_match_geneart": {"codon": "match_codon_usage", **LEVELS["geneart"]},
    "V07_best_strict": {"codon": "use_best_codon", **LEVELS["strict"]},
    "V08_match_strict": {"codon": "match_codon_usage", **LEVELS["strict"]},

    # 9: "light" optimization - only rare codons replaced, no CodonOptimize
    "V09_light_geneart": {"codon": None, **LEVELS["geneart"]},

    # 10: GC as a soft objective (target) instead of a hard range
    "V10_match_geneart_softGC": {"codon": "match_codon_usage",
                                 **{**LEVELS["geneart"], "gc": None},
                                 "gc_target": GC_TARGET},
}


def build_problem(sequence, cfg):
    """Turn one variant definition into a DnaOptimizationProblem."""
    constraints = [EnforceTranslation()]                       # protein stays identical
    objectives = []
    if cfg.get("codon"):
        objectives.append(CodonOptimize(species=SPECIES, method=cfg["codon"]))

    if cfg.get("gc"):                                          # hard GC range per window
        window, low, high = cfg["gc"]
        constraints.append(EnforceGCContent(mini=low, maxi=high, window=window))

    constraints += [AvoidPattern(p) for p in cfg.get("runs", [])]   # homopolymer limits

    if cfg.get("kmer"):                                        # no repeated k-mers
        constraints.append(UniquifyAllKmers(cfg["kmer"]))

    if cfg.get("rare"):                                        # no rare human codons
        constraints.append(AvoidRareCodons(min_frequency=cfg["rare"], species=SPECIES))

    if cfg.get("motifs"):                                      # unwanted motifs
        constraints += [AvoidPattern(m) for m in MOTIFS]

    if cfg.get("sites"):                                       # restriction sites
        constraints += [AvoidPattern(f"{enzyme}_site") for enzyme in SITES]

    if cfg.get("gc_target"):                                   # soft GC objective
        window, target = cfg["gc_target"]
        objectives.append(EnforceGCContent(target=target, window=window))

    return DnaOptimizationProblem(
        sequence=sequence, constraints=constraints, objectives=objectives, logger=None
    )


def run_variant(sequence, cfg):
    """Resolve constraints, optimize, and return (sequence, status)."""
    try:
        problem = build_problem(sequence, cfg)
        problem.resolve_constraints()
        problem.optimize()
        status = "OK" if problem.all_constraints_pass() else "Constraints not met"
        return problem.sequence, status
    except Exception as error:   # e.g. no solution for the given constraints
        return sequence, f"FAILED: {str(error)[:100]}"


def count_changes(before, after):
    """Number of nucleotides that differ between two sequences."""
    return sum(a != b for a, b in zip(before, after))


def main():
    parser = argparse.ArgumentParser(description="Run 10 DNA Chisel variants.")
    parser.add_argument("--input", required=True, help="human_codon_sequences.xlsx")
    parser.add_argument("--output", required=True, help="Output Excel file")
    parser.add_argument("--limit", type=int, default=None, help="Only first N genes (test)")
    args = parser.parse_args()

    data = pd.read_excel(args.input)
    for col in (NAME_COL, INPUT_COL):
        if col not in data.columns:
            raise ValueError(f"Missing column: {col}")
    if args.limit:
        data = data.head(args.limit)

    rows = []
    wide = data.copy()
    for variant in VARIANTS:
        wide[variant] = ""

    for idx, record in data.iterrows():
        start_seq = "".join(str(record[INPUT_COL]).split()).upper()
        for variant, cfg in VARIANTS.items():
            np.random.seed(42)                      # same random start every time
            new_seq, status = run_variant(start_seq, cfg)
            rows.append({
                NAME_COL: record[NAME_COL],
                "Variant": variant,
                "Status": status,
                "Nucleotides changed": count_changes(start_seq, new_seq),
                "Optimized DNA Sequence": new_seq,
            })
            wide.at[idx, variant] = new_seq
            print(f"{record[NAME_COL]:<10} {variant:<30} {status}")

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output_path) as writer:
        pd.DataFrame(rows).to_excel(writer, sheet_name="Results", index=False)
        wide.to_excel(writer, sheet_name="Sequences", index=False)
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
