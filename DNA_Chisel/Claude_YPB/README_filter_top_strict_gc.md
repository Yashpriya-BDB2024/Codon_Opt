# Selecting the top 10 optimized sequences per gene (strict GC window)

**Script:** `DNAChisel/filter_top_strict_gc.py`
**Input:** the `Results` sheet of the grid-search output (`human_grid.xlsx`: 8 genes × 2,592 variants)
**Output:** `human_strictGC_top10.xlsx`, a `Top` sheet with 10 sequences per gene (80 rows)

The grid script and the run scripts are not modified or imported. The filter reads their output only, and uses **only columns that the grid output already contains**. It does not compute any scoring metric.

---

## 1. Summary

| Question | Decision | Basis |
|---|---|---|
| Which sequences are eligible? | Valid, made with the rare-codon constraint on, distinct, clean (no motifs, repeats or long homopolymers), local GC in range, overall GC 45–55% | Project lead's rules (option 1) |
| Is there a CAI cut-off? | No | Project lead |
| Is there a codon-method restriction? | No (`use_best_codon` is not forced) | Project lead |
| How are the survivors ranked? | By CAI, highest first, **as the last step** | Project lead |
| How many per gene? | 10 | Project lead |

---

## 2. How to run

From the `DNAChisel` folder, after the full grid run:

```bash
python filter_top_strict_gc.py                    # reads outputs/human_grid.xlsx
python filter_top_strict_gc.py --top-n 10 --input outputs/human_grid.xlsx
```

It needs `pandas` and `openpyxl`. It prints the funnel (Section 5) and a per-gene summary, and gives a `WARNING` for any gene with fewer than N survivors. The output is `outputs/<name>_strictGC_top<N>.xlsx`.

The `Top` sheet has the rank, how many variants gave the same DNA, and then every column of the grid output, with the DNA last.

---

## 3. The order of the steps and the reason for each

The order matters in two places only: step 2 must come before step 3, and CAI ranking must come last. Steps 4 to 8 are independent pass/fail tests on the DNA, so their order does not change the final set. It only changes the funnel counts.

| # | Step | Rule | Why |
|---|---|---|---|
| 1 | Valid result | `Status == OK` | `Constraints not met` and `FAILED` rows have no usable sequence. Across the run, 19,764 of 20,736 were OK, 207 were `Constraints not met` and 765 were `FAILED` |
| 2 | Rare-codon constraint on | `rare_codon_min_freq` is set (0.11) | See Section 4 |
| 3 | Remove duplicates | One row per distinct DNA per gene | Many settings converge on identical DNA, so keeping copies would fill the top 10 with repeats. The first variant is kept and the others are listed in the output |
| 4 | No forbidden motifs | `Forbidden motifs left = 0` | Poly-A, ARE and TATA motifs cause premature polyadenylation, mRNA instability and cryptic transcription |
| 5 | No repeats | `Repeated 10-mers = 0` | Repeats cause synthesis dropout and recombination. The grid's reference note says vendors tolerate repeated 9-mers but not 10-mers |
| 6 | No long homopolymers | `Longest homopolymer ≤ 5` | Runs of 6 or more cause polymerase slippage and synthesis errors |
| 7 | Local GC in range | `GC% min (50nt win) ≥ 30` and `GC% max (50nt win) ≤ 70` | Local GC above 70% fails synthesis, below 30% is unstable |
| 8 | Overall GC in a strict window | `GC% global` 45–55 | Project lead's option 1. CAI and GC rise together, so GC is fixed first and CAI second (Section 5) |
| 9 | Rank by CAI | Highest first, keep the top 10 | The objective of the optimization. Equal CAI keeps the grid's variant order, so the result is reproducible |

### Why step 2 comes before step 3

Step 2 filters on **how a variant was made** (a setting), whereas steps 4 to 8 test **the DNA itself**. Identical DNA can come from variants with different settings. If duplicates were removed first, the kept copy could be one made without the constraint, and a valid sequence would be dropped wrongly. Filtering first and removing duplicates second avoids this. On this data both orders give the same count (3,302 distinct sequences), so nothing was lost, but this is the safe order for other data.

---

## 4. The rare-codon rule, and how it is applied

The project lead's rule is "rare codons ≤ 5%". The grid output has no column that measures rare codons, but it does record the **setting** that controls them: `rare_codon_min_freq`.

- When `rare_codon_min_freq = 0.11`, the grid applied the `AvoidRareCodons` hard constraint. No codon with a frequency below 0.11 within its amino-acid family can appear. For human these are TCG, CTA, CGT, GCG, TTA, CCG and ACG.
- A sequence from such a variant has 0 of those codons, which satisfies "≤ 5%". The rule is applied as a constraint guarantee, not as a measured percentage.
- It is set in 1,296 of the 2,592 variants, so requiring it halves the pool: 7,330 distinct sequences become 3,302 across the 8 genes.

**Checked separately (not used to select):**
- All 80 selected sequences contain 0 codons below the 0.11 cut-off.
- The project scorer's own rare-codon measure (relative adaptiveness w < 0.3) is at most 0.98% on the selected sequences, with none above 5%.

---

## 5. Why CAI is ranked last, after the GC window

CAI and GC rise together because the best human codons are GC-rich. Within a GC 40–60% window, the correlation between CAI and GC was 0.86–0.92 across the genes, and 0.60–0.79 inside the 45–55% window (measured on the pool before the rare-codon step). Ranking by CAI alone would push every pick to the top GC edge.

Setting the GC window first and ranking by CAI second keeps GC in a controlled range, which is the project lead's option 1. Even so, the highest-CAI sequences in a window sit near its upper edge. All selected sequences have GC between 54.0 and 55.0% except two genes, CTNNB1 (low of 53.5%) and ALDH1A3 (low of 51.7%).

---

## 6. Funnel: human run

Each row shows the number of distinct sequences left after that step and all steps above it.

| Stage | CTNNB1 | IRF4 | FOXA1 | ERG | AURKB | AURKA | PREP | ALDH1A3 |
|---|---|---|---|---|---|---|---|---|
| Variants run | 2,592 | 2,592 | 2,592 | 2,592 | 2,592 | 2,592 | 2,592 | 2,592 |
| 1. Status OK | 2,520 | 2,589 | 2,412 | 2,580 | 2,584 | 2,538 | 2,352 | 2,189 |
| 2. Rare-codon constraint on | 1,260 | 1,295 | 1,116 | 1,284 | 1,296 | 1,276 | 1,164 | 1,067 |
| 3. Distinct DNA | 397 | 573 | 428 | 430 | 515 | 316 | 363 | 280 |
| 4. No forbidden motifs | 397 | 573 | 428 | 430 | 515 | 316 | 363 | 280 |
| 5. No repeated 10-mers | 237 | 401 | 186 | 298 | 364 | 242 | 226 | 160 |
| 6. Homopolymer ≤ 5 | 199 | 352 | 178 | 229 | 336 | 214 | 185 | 121 |
| 7. Window GC 30–70% | 199 | 352 | 178 | 229 | 336 | 214 | 185 | 121 |
| 8. GC global 45–55% | 62 | 190 | 93 | 113 | 209 | 126 | 121 | 69 |
| **9. Selected (top 10 by CAI)** | **10** | **10** | **10** | **10** | **10** | **10** | **10** | **10** |

The motif gate (step 4) and window-GC gate (step 7) remove nothing because the grid already enforces them in every variant. They remain as safety checks.

---

## 7. Which grid columns are used, and which are not

**Settings (how a sequence was made), as seen in the 2,592-variant grid:**

| Setting | Values in the grid | How it is handled |
|---|---|---|
| `rare_codon_min_freq` | 0.11 in 1,296 variants, none in 1,296 | **Used** (step 2) |
| `avoid_motifs` | True in all | Always on, so no filter is needed. The measured `Forbidden motifs left` is checked as well |
| `avoid_restriction_sites` | True in all (13 enzymes) | Always on, so no filter is needed |
| `gc_min`, `gc_max`, `gc_window` | 0.30/0.40/0.46 and 0.60/0.66/0.70, over 50 nt, in all | Always on. The measured window GC is checked (step 7) |
| `homopolymer_max` | 4, 5 or 6 | Judged on the measured `Longest homopolymer` instead |
| `kmer` | none, 9, 10 or 15 | Judged on the measured `Repeated 10-mers` instead |
| `codon_method` | none, `use_best_codon`, `match_codon_usage` | No restriction |
| `avoid_hairpins` | **False in all** | Never applied (see Section 8) |
| `tm_min`/`tm_max`, `terminal_gc_min`/`max` | on in 1,296 variants | Not required by the project lead, so not filtered |

**Measured columns used:** `Status`, `Forbidden motifs left`, `Repeated 10-mers`, `Longest homopolymer`, `GC% min/max (50nt win)`, `GC% global`, `CAI`, and `Optimized DNA Sequence` (duplicates).

**Not used:** `Objective score` (not comparable across codon methods), `Codon-match score`, `% identity` and `Nucleotides changed` (they measure distance from a reverse-translation seed, so they mean nothing here), `5' Tm` and `3' Tm`, `Runtime`, `Seed used`, `Host`, `Domain` and `Input kind`.

---

## 8. Result for the human run

| Gene | CAI | GC % | Codon methods in the 10 |
|---|---|---|---|
| CTNNB1 | 0.824–0.902 | 53.52–54.97 | `use_best_codon` 8, `match_codon_usage` 2 |
| IRF4 | 0.941–0.944 | 54.55–54.99 | `use_best_codon` 10 |
| FOXA1 | 0.898–0.917 | 54.59–54.94 | `use_best_codon` 10 |
| ERG | 0.942–0.958 | 54.36–54.95 | `use_best_codon` 10 |
| AURKB | 0.934–0.941 | 54.38–54.91 | `use_best_codon` 10 |
| AURKA | 0.932–0.952 | 54.01–55.00 | `use_best_codon` 10 |
| PREP | 0.955–0.958 | 54.27–54.88 | `use_best_codon` 10 |
| ALDH1A3 | 0.845–0.942 | 51.69–54.95 | `use_best_codon` 7, `match_codon_usage` 3 |

---

## 9. Limitations and open points

1. **Hairpins were never controlled.** `avoid_hairpins` is False in all 2,592 variants, so no hairpin constraint was applied and the filter has nothing to select on.
2. **Splice sites were never controlled.** The grid has no splice-site constraint.
3. **The rare-codon rule is applied through the constraint, not a measured percentage.** The 0.11 frequency cut-off and the scorer's w < 0.3 measure are different definitions. The check in Section 4 shows no sequence is above 1%.
4. **The strict 45–55% window caps CAI.** The GeneArt sequences have GC of about 56.5–57.8%, outside this window. A wider window such as 40–60% allows higher CAI, but then the top-CAI picks sit at GC 58–60%.
5. **Top-ranked sequences sit near the upper edge of the GC window** (Section 5).
6. **For CTNNB1 and ALDH1A3, the 10th-best CAI is much lower than the best.** The 62 and 69 sequences left in the window run out of high-CAI options (CTNNB1 0.824–0.902, ALDH1A3 0.845–0.942).
7. **The thresholds are the project lead's choices**, not universal rules. Edit the constants at the top of the script to change them.
8. **Duplicates keep the first variant only.** The names of the other variants that produced the same DNA are listed in the output.
9. **The grid's CAI is used for ranking.** It differs from the scorer's CAI by up to about 0.0025, which is negligible here.

---

## 10. Checks performed

- A separate row-by-row recomputation chose the same 10 variants in the same order for all 8 genes.
- All 80 selected rows satisfy every gate, none duplicate another, and each gene has 10.
- 0 codons below the 0.11 cut-off in the 80 selected sequences. The scorer's rare-codon measure is at most 0.98%.
- A gene with no survivors gives a warning and the other genes still complete. A missing input column stops with a clear message.
- Two runs give identical output.
- `codon_grid_search.py` is byte-identical to the original.

---

## 11. Next step

`Sequence_scoring/compare_top_vs_geneart.py` scores each selected sequence against the gene's GeneArt sequence, using `codon_compare.py` and `codon_features.py` unchanged:

```bash
cd ../Sequence_scoring
pip install seqfold          # needed for the three RNA-structure metrics
python compare_top_vs_geneart.py
```
