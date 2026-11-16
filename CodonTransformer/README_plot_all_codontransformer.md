# `plot_all_codontransformer.py` — README

Every plot and summary table of **GeneArt vs all CodonTransformer rows** (human, 8 genes). It reads the score table
written by `score_codontransformer.py` and draws distributions, clusters, win/draw/loss summaries, positional profiles,
liability maps and host-reference plots. CodonTransformer is never mixed with DNA Chisel or Codon Harmony here.

**Locked version (v2).** SHA-256 of the delivered file (v1 hash `9cc390f9…c414b47` is superseded by the dot-colour change in §1, row 4):

```
e19c41d50b8975b032b94aa968cc2a5cc54457b0dd75fd754b86de5494b81607
```

Check on your machine with `sha256sum plot_all_codontransformer.py` (Windows: `certutil -hashfile plot_all_codontransformer.py SHA256`).
Any edit changes the hash; treat a different hash as a different version.

---

## 1. Relation to the other plot scripts

Derived from `plot_all_codonharmony.py`, which is itself derived from `plot_all.py`. Plotting, clustering, win/draw/loss and
profile code is identical. Changes (verified by diff):

| # | Change |
|---|---|
| 1 | Labels / names: `Codon Harmony` / `CodonHarmony` → `CodonTransformer` (titles, legends, file names, column names, sheet `CodonTransformer`). |
| 2 | Defaults: `--scores outputs/scores/scores_codontransformer.xlsx`; grid `../CodonTransformer/outputs/human_codontransformer_grid.xlsx`; `--out-dir outputs/plots_codontransformer`. |
| 3 | **Only functional change:** pair-plot file tag. The Codon Harmony version used `Variant.split("_")[0]` (cell id, e.g. `C017`, unique per row). For CodonTransformer that would be `T0.2`, shared by many draws, so files would overwrite each other. The full `Variant` (`T0.2_p0.8_d1`) is used. |
| 4 | **Visibility (v2):** the CodonTransformer dots of the distribution plots were small, faint grey-blue (`#7A8B99`, size 3, alpha 0.22). They are now blue `#2B6CB0` (the CodonTransformer colour of the win/loss and scatter plots), size 16, alpha 0.60; the legend and the caption say "blue dot". Colours/sizes only: no data or metric changes. |

It imports `codon_features`, `codon_compare`, `score_all_sequences` and `run_comparison` (host-reference loaders only); none is modified.

---

## 2. Inputs (same folder layout as the other tools)

| File | Used for |
|---|---|
| `outputs/scores/scores_codontransformer.xlsx` (sheets `CodonTransformer`, `GeneArt`) | every number in the distribution / cluster / heatmap / win-rate / scatter plots |
| `human_codontransformer_grid.xlsx` (sheet `Results`) + original trials workbook | only the plots that need the DNA (profiles, motif map, codon usage, agreement) |
| `host_files/Human/` (`human_codon.csv`, `human_tgcn.csv`, `Human.txt`, dinucleotide `.svg`) | scorer references; host lines / plots (skipped with a warning if `Human.txt` is missing) |
| `codon_features.py`, `codon_compare.py`, `score_all_sequences.py`, `run_comparison.py` | imported |

`pip install numpy pandas matplotlib openpyxl scikit-learn seqfold tqdm`

---

## 2a. Why `score_all_sequences.py` is imported

`score_all_sequences.py` is the DNA Chisel scoring script, but here it is used only as a **shared toolbox**. No DNA Chisel grid or
sequence is read, so nothing from DNA Chisel is mixed into the CodonTransformer plots. This script uses:

| From `score_all_sequences.py` | Used for |
|---|---|
| `load_geneart` | reads the 8 GeneArt sequences and proteins from the trials workbook |
| `translate` | translation (protein) check of the DNA |
| `ROOT` | project folder (finds `host_files/Human` and the trials workbook) |
| `DNA` | the column name `Optimized DNA Sequence` |
| `CODON_CSV`, `TGCN_CSV` | paths of `human_codon.csv` / `human_tgcn.csv`, so the DNA that is re-scored for the profile plots uses the same reference files as the score table |

Why import instead of copying the functions:
- **One source of truth.** The scoring script and the plot scripts of every tool use identical definitions and reference paths; a copy could drift if one of them is edited later.
- **Importing does not run it.** Its `main()` is behind `if __name__ == "__main__"`, so only the functions and constants are loaded.
- **It is not modified**, and it is not part of the locked plot script.

Requirement: `score_all_sequences.py` must be in the same folder as `plot_all_codontransformer.py`, together with `codon_features.py`,
`codon_compare.py` and `run_comparison.py`.

---

## 3. Run

```bash
python plot_all_codontransformer.py
python plot_all_codontransformer.py --workers 40 --clusters 3 --pair-plots 10
python plot_all_codontransformer.py --no-profiles        # only the plots that need the scores file
```

| Option | Default | Meaning |
|---|---|---|
| `--scores` | `outputs/scores/scores_codontransformer.xlsx` | score table |
| `--grid` | `../CodonTransformer/outputs/human_codontransformer_grid.xlsx` | DNA source (same file as scoring) |
| `--trials` | auto | original trials workbook |
| `--out-dir` | `outputs/plots_codontransformer` | output folder |
| `--clusters` | 3 | K-means Low/Mid/High per metric per gene (2–5) |
| `--workers` | min(40, cores−1) | re-scoring of distinct DNA for profile plots |
| `--pair-plots N` | 10 | per gene, the N best-CAI distinct CodonTransformer sequences get a codon_compare overlay + scorecard + 9-panel figure (0 = none, negative = all) |
| `--no-profiles` | off | skip everything that needs DNA |
| `--no-host` | off | no host-reference lines / plots |
| `--dpi` | 130 | figure resolution |
| `--stop-policy` | `strip` | **must equal the value used for `score_codontransformer.py`** (default `strip` for both). All 640 CodonTransformer sequences end in a stop codon; `asis` would not match the score table. |
| `--max-rows-per-gene` | 0 | TEST ONLY |

---

## 4. What is (and is not) scored again

The 22 metrics in the distribution, cluster, heatmap, win-rate and scatter plots are **read from the score table**, not recomputed.
The profile / motif / codon plots need the DNA, so `codon_features.compute_features` is run once per distinct sequence
(488 for the delivered grid; cached in `profile_cache.pkl`). Its scalar scores are compared with the score table in the sheet
`Check_vs_scores_all` as a safeguard.

---

## 5. Outputs (`outputs/plots_codontransformer/`)

| Folder / file | Content |
|---|---|
| `distributions/` | one image per metric (22): every gene, every CodonTransformer row, GeneArt marked |
| `clusters/` | one image per metric (22): K-means Low/Mid/High per gene, GeneArt marked |
| `extras/` | GC1, GC2, ENC deviation, distance to host |
| `summary/` | heatmaps (GeneArt percentile / cluster / beats), win/draw/loss heatmap, win rate per metric, tally per gene, CAI/tAI-vs-GC scatter, ENC vs GC3 |
| `profiles/` | per gene, 12 positional panels GeneArt vs CodonTransformer |
| `geneart_features/` | codon_features 9-panel figure of each GeneArt sequence |
| `liabilities/` | motif map (restriction / splice / homopolymer positions), restriction sites per enzyme |
| `host_reference/` | codon usage, GC by position, dinucleotides vs the human host (needs `host_files/Human`) |
| `agreement/` | per-codon agreement of the CodonTransformer rows with GeneArt |
| `more_metrics/` | distributions + clusters for Kozak, GC min/std 60 bp, CpG count, longest repeat, AT%, other scalars |
| `pairs/`, `pair_features/` | with `--pair-plots N`; files named `<gene>_<Variant>_…png` |
| `summary_tables.xlsx`, `row_clusters.csv` | every number behind the plots (sheets: Read_me, GeneArt_vs_CodonTransformer, Cluster_ranges, Win_counts, Pivot_*, Extras_*, Restriction_per_enzyme, More_*, Check_vs_scores_all) |

---

## 6. How to read the results

- A **pair** is one CodonTransformer row vs the GeneArt sequence of the same gene (80 pairs per gene, 640 in total, 488 distinct sequences).
  Win / draw / loss follows `codon_compare` (higher / lower / GC band / not scored). Percentages describe CodonTransformer's
  distribution against **one** GeneArt sequence per gene (n = 8 genes); they are descriptive, not a significance test.
- CodonTransformer runs **without constraints** (only `match_protein`). Restriction, splice, homopolymer, repeat and GC-peak metrics
  are measured, not optimised, whereas GeneArt is a constrained vendor product, so these compare an unconstrained model with a
  vendor sequence. The 22 metrics mix codon-usage metrics with manufacturability metrics; read them by family.
- The metrics are in-silico proxies, not expression or function.
- **Blank (NA) cells in `summary_tables.xlsx` are expected**: see section 6a.
- CAI / rare-codon / tAI use your `human_codon.csv` and `human_tgcn.csv`; CodonTransformer never used them for generation.
- Do not overlay these plots on DNA Chisel's (top-10 by CAI) without saying the selection differs.

---

## 6a. Why some cells in `summary_tables.xlsx` are blank (NA)

**The metrics are calculated; only one column is blank.** The blank is `GeneArt_beats_pct`, which answers "in what % of the
CodonTransformer rows is GeneArt better". That needs a better/worse direction. In `codon_compare.COMPARE_METRICS` five metrics are
deliberately marked `na` (no universal direction: a higher or lower value is not automatically better):

`ENC`, `tAI_ramp`, `minmax_mean`, `codon_autocorr_lag1`, `GC3_pct`

For these the sheet still reports where GeneArt sits among the CodonTransformer rows (GeneArt value, min / median / max, percentile); only
the "beats" percentage is left empty. They are also not counted in the win / draw / loss tallies (same as in the DNA Chisel and Codon Harmony tables).

Example from the delivered run (gene CTNNB1):

| metric | better | GeneArt | CT min | CT median | CT max | GeneArt percentile | `GeneArt_beats_pct` |
|---|---|---|---|---|---|---|---|
| ENC | na | 34.96 | 33.53 | 34.15 | 42.01 | 68.75 | blank |
| GC3_pct | na | 69.78 | 72.09 | 79.77 | 80.41 | 0.00 | blank |
| CAI | higher | 0.9296 | 0.8680 | 0.9301 | 0.9360 | 45.62 | 45.62 |
| restriction_sites_total | lower | 9 | 11 | 13 | 17 | 0.00 | 100.00 |

Where blanks occur (all in `GeneArt_beats_pct` / `Pivot_beats_pct`):

| Sheet | Blank for | Count in the delivered run |
|---|---|---|
| `GeneArt_vs_CodonTransformer` | the 5 `na` metrics above | 40 (5 metrics x 8 genes) |
| `Pivot_beats_pct` | the same 5 metrics (columns) | 5 columns x 8 genes |
| `Extras_GeneArt_vs_CodonTransformer` | GC1, GC2, ENC deviation and the host-distance columns (context columns, no rule) | 64 (all rows) |
| `More_GeneArt_vs_CodonTransformer` | AT_pct, MFE_mean_window, minmax_min, minmax_pct_rare_windows, tAI_body, tAI_head, tm_mean (no rule) | 56 |

What would indicate a problem: a blank in any **value** column (GeneArt, CodonTransformer min / median / max, percentile) or a blank in a metric
that has a rule (higher / lower / band). None occurs in the delivered run, and `row_clusters.csv` has no blanks (a cluster named `Constant`
only means all rows of that gene share one value). If `seqfold` were missing, `TIS_accessibility_5p`, `MFE_min_window` and
`stable_structure_count` would be empty in the score table itself and a warning would be printed (see section 7).

---

## 7. Checks on the first real run

1. `Check_vs_scores_all` / console line `Check against scores_codontransformer.xlsx: 0 metric(s) with differences`. Any difference means a different grid file,
   different reference files, `seqfold` present in only one run, or a different `--stop-policy`.
2. Pair-plot files carry the full variant name and are not overwritten (`pairs/` has 2 files per selected sequence).
3. Without `seqfold`, `TIS_accessibility_5p`, `MFE_min_window`, `stable_structure_count` and the MFE profile panel are empty (warning printed).

---

## 8. Verification status

Run on the delivered files (score table, grid workbook, trials workbook, real `codon_features.py`):
- `--no-profiles`: 51 PNG + tables, no errors.
- Full run with profiles and `--pair-plots 2`: exit 0, 144 PNG, unique pair-file names.
- In that run all 13 nucleotide-level metrics (ENC, GC, GC max 60 bp, GC3, CpG, both hairpin counts, splice, homopolymer, restriction,
  both repeat counts, tm_spread) matched the score table exactly on all 648 sequences, which confirms the stop-strip handling is consistent.

**Not verified:** the test used stand-in reference files and no `seqfold`, so the 9 metrics that depend on `human_codon.csv`, `human_tgcn.csv`
or `seqfold` (CAI, rare codons, tAI, tAI ramp, MinMax, autocorrelation, TIS, MFE, stable structures) differed from your real score table, as
expected, and the host-reference plots were skipped. On your machine with the real files and `seqfold` these must show 0 differences.
