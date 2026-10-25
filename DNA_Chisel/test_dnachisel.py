"""
============================================================
DNA Chisel — Starter & Test Script
License: MIT | Commercial use: YES
Algorithm: Sliding-window, deterministic, multi-parameter

What this script does:
  1. Back-translates a protein sequence to DNA using
     human/CHO codon tables
  2. Optimizes simultaneously for:
     - CAI (codon usage for human/CHO)
     - GC content (40-65% per 60bp window)
     - Restriction site removal (BsaI, BsmBI)
     - Hairpin avoidance
  3. Reports before/after comparison

Run with:
    source ~/codon_optimization/dnachisel/venv/bin/activate
    python3 test_dnachisel.py
============================================================
"""

from dnachisel import (
    DnaOptimizationProblem,
    CodonOptimize,
    EnforceGCContent,
    AvoidPattern,
    AvoidHairpins,
    EnforceTranslation,
)
from python_codon_tables import get_codons_table

# ── Replace this with your actual protein sequence ──────────
PROTEIN = (
    "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLT"
    "LKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKS"
    "AMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDF"
)

# ── Helper: GC content calculation ──────────────────────────
def gc_content(seq):
    seq = seq.upper()
    return round(sum(1 for b in seq if b in "GC") / len(seq) * 100, 1)

def gc_range(seq, window=60):
    windows = [gc_content(seq[i:i+window])
               for i in range(0, len(seq) - window, window)]
    return min(windows), max(windows)

# ── Step 1: Back-translate using human codon table ──────────
print("=" * 60)
print("  DNA Chisel — Codon Optimization for CHO")
print("=" * 60)
print(f"\nProtein length: {len(PROTEIN)} amino acids")

codon_table = get_codons_table("c_griseus_10029")

def backtranslate(protein, codon_table):
    dna = ""
    for aa in protein:
        if aa not in codon_table:
            continue
        codons = codon_table[aa]
        best = max(codons, key=lambda c: codons[c])
        dna += best
    return dna + "TAA"

start_seq = backtranslate(PROTEIN, codon_table)
print(f"\nStarting sequence ({len(start_seq)} bp):")
print(f"  {start_seq[:60]}...")
print(f"  GC content:  {gc_content(start_seq)}%")
gc_lo, gc_hi = gc_range(start_seq)
print(f"  GC range:    {gc_lo}% – {gc_hi}%  (lower range = smoother = better)")

# ── Step 2: Set up optimization problem ─────────────────────
print("\nOptimization parameters:")
print("  → Maximize CAI for Cricetulus griseus (CHO cells)")
print("  → Enforce GC content 40–65% per 60bp window")
print("  → Remove BsaI restriction sites")
print("  → Remove BsmBI restriction sites")
print("  → Avoid hairpin loops (stem≥6, window=200bp)")

problem = DnaOptimizationProblem(
    sequence=start_seq,
    constraints=[
        EnforceTranslation(),
        EnforceGCContent(mini=0.40, maxi=0.65, window=60),
        AvoidPattern("BsaI_site"),
        AvoidPattern("BsmBI_site"),
        AvoidHairpins(stem_size=6, hairpin_window=200),
    ],
    objectives=[
        CodonOptimize(
            species="c_griseus_10029",
            method="use_best_codon"
        ),
    ]
)

# ── Step 3: Run optimization ─────────────────────────────────
print("\nRunning optimization...")
problem.resolve_constraints()
problem.optimize()

# ── Step 4: Results ──────────────────────────────────────────
opt_seq = problem.sequence
print("\n" + "=" * 60)
print("  Results")
print("=" * 60)
print(f"\nOptimized sequence ({len(opt_seq)} bp):")
print(f"  {opt_seq[:60]}...")

print(f"\n  {'Metric':<25} {'Before':>10} {'After':>10}")
print(f"  {'-'*45}")
print(f"  {'GC content (overall)':<25} {gc_content(start_seq):>9}% {gc_content(opt_seq):>9}%")
gc_lo_s, gc_hi_s = gc_range(start_seq)
gc_lo_o, gc_hi_o = gc_range(opt_seq)
print(f"  {'GC min (per window)':<25} {gc_lo_s:>9}% {gc_lo_o:>9}%")
print(f"  {'GC max (per window)':<25} {gc_hi_s:>9}% {gc_hi_o:>9}%")
print(f"  {'Constraints passed':<25} {'—':>10} {str(problem.all_constraints_pass()):>10}")

print(f"\nObjectives summary:")
print(f"  {problem.objectives_text_summary()}")

print("\n" + "=" * 60)
print("  Next steps:")
print("  1. Replace PROTEIN with your actual sequence")
print("  2. Add more constraints as needed:")
print("     AvoidPattern('EcoRI_site')")
print("     EnforceGCContent(mini=0.45, maxi=0.60, window=30)")
print("  3. Change species for non-CHO hosts:")
print("     species='e_coli_316407'  # E. coli")
print("     species='s_cerevisiae_4932'  # Yeast")
print("=" * 60)
