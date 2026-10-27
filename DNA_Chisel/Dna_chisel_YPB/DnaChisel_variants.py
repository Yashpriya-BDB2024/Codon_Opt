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
# Cutoffs taken from the part-2 report (what GeneArt sequences stayed within)
# ----------------------------------------------------------------------------
GC_50 = (50, 0.44, 0.68)      # (window, mini, maxi)  50-bp windows: 44-68 %
GC_100 = (100, 0.49, 0.66)    # 100-bp windows: 49-66 %
GC_TARGET = (50, 0.57)        # (window, target) - overall GeneArt GC is ~57 %
HOMOPOLYMERS = ["6xA", "5xC", "5xG", "5xT"]    # runs longer than GeneArt's max
MOTIFS = ["AWTAAA", "ATTTA", "TATAAA"]         # polyA signal (W = A/T), AU-rich, TATA
SITES = ["BsaI", "BsmBI", "BbsI", "SapI", "EcoRI", "XhoI", "NotI",
         "XbaI", "SpeI", "KpnI", "SalI", "AscI", "PacI"]   # absent in GeneArt
KMER = 10                     # no repeated 10-mers in GeneArt
RARE = 0.11                   # lowest codon frequency GeneArt ever used
RARE_STRICT = 0.15            # a stricter rare-codon cutoff to try

# ----------------------------------------------------------------------------
# The 10 variants. Every variant always has EnforceTranslation (protein must
# not change) and one CodonOptimize objective. Other keys switch constraints on.
# ----------------------------------------------------------------------------
VARIANTS = {
    # 1-2: codon strategy only (baselines)
    "V01_best_codon": {"codon": "use_best_codon"},
    "V02_match_usage": {"codon": "match_codon_usage"},

    # 3-4: + local GC window + homopolymer limits
    "V03_best_GC_runs": {"codon": "use_best_codon", "gc": GC_50, "runs": True},
    "V04_match_GC_runs": {"codon": "match_codon_usage", "gc": GC_50, "runs": True},

    # 5-6: + unique 10-mers + rare codon cutoff
    "V05_best_GC_runs_kmer_rare": {"codon": "use_best_codon", "gc": GC_50, "runs": True,
                                   "kmer": KMER, "rare": RARE},
    "V06_match_GC_runs_kmer_rare": {"codon": "match_codon_usage", "gc": GC_50, "runs": True,
                                    "kmer": KMER, "rare": RARE},

    # 7-8: everything (also unwanted motifs + restriction sites)
    "V07_best_full": {"codon": "use_best_codon", "gc": GC_50, "runs": True, "kmer": KMER,
                      "rare": RARE, "motifs": True, "sites": True},
    "V08_match_full": {"codon": "match_codon_usage", "gc": GC_50, "runs": True, "kmer": KMER,
                       "rare": RARE, "motifs": True, "sites": True},

    # 9: stricter cutoffs (100-bp GC window, higher rare-codon cutoff)
    "V09_best_full_strict": {"codon": "use_best_codon", "gc": GC_100, "runs": True,
                             "kmer": KMER, "rare": RARE_STRICT, "motifs": True, "sites": True},

    # 10: GC as a soft objective (target) instead of a hard range
    "V10_match_soft_GC": {"codon": "match_codon_usage", "runs": True, "motifs": True,
                          "sites": True, "gc_target": GC_TARGET},
}


def build_problem(sequence, cfg):
    """Turn one variant definition into a DnaOptimizationProblem."""
    constraints = [EnforceTranslation()]                       # protein stays identical
    objectives = [CodonOptimize(species=SPECIES, method=cfg["codon"])]

    if "gc" in cfg:                                            # hard GC range per window
        window, low, high = cfg["gc"]
        constraints.append(EnforceGCContent(mini=low, maxi=high, window=window))

    if cfg.get("runs"):                                        # homopolymer limits
        constraints += [AvoidPattern(p) for p in HOMOPOLYMERS]

    if "kmer" in cfg:                                          # no repeated k-mers
        constraints.append(UniquifyAllKmers(cfg["kmer"]))

    if "rare" in cfg:                                          # no rare human codons
        constraints.append(AvoidRareCodons(min_frequency=cfg["rare"], species=SPECIES))

    if cfg.get("motifs"):                                      # unwanted motifs
        constraints += [AvoidPattern(m) for m in MOTIFS]

    if cfg.get("sites"):                                       # restriction sites
        constraints += [AvoidPattern(f"{enzyme}_site") for enzyme in SITES]

    if "gc_target" in cfg:                                     # soft GC objective
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
