# plot_all.py: every plot of GeneArt vs all DNA Chisel rows

**Script:** `Sequence_scoring/plot_all.py` (replaces `plot_scores_all.py` and `plot_profiles_host.py`; delete those two)
**Needs next to it (unchanged):** `codon_features.py`, `codon_compare.py`, `score_all_sequences.py`, `run_comparison.py`
**Output:** `outputs/plots/`

Comparison unit everywhere: one **DNA Chisel row** (every Status-OK grid row) against the **GeneArt sequence of the same gene**. Genes are never mixed.

## 1. Run

```bash
cd Sequence_scoring
pip install numpy pandas matplotlib openpyxl scikit-learn seqfold tqdm
python plot_all.py                                  # everything, 10 pair figures per gene
python plot_all.py --workers 40 --clusters 3 --pair-plots 10
python plot_all.py --no-profiles                    # only plots that need scores_all.xlsx (fast, no DNA)
python plot_all.py --grid ../DNAChisel/outputs/human_grid_1.xlsx --scores outputs/scores/scores_all.xlsx
```

Use the **same grid file** that `score_all_sequences.py` used. The script stops with a message if a row of `scores_all.xlsx` has no DNA in the grid file.

| Option | Meaning |
|---|---|
| `--clusters 2..5` | K-means clusters per metric per gene (default 3 = Low / Mid / High) |
| `--pair-plots N` | per gene, the N best-CAI distinct DNA Chisel sequences get a `codon_compare` overlay + scorecard and the `codon_features` 9-panel figure. Default 10, `0` = none, `-1` = every row (tens of thousands of files) |
| `--no-profiles` | skip everything that needs the DNA |
| `--no-host` | no human-host lines or host plots |
| `--workers` | processes for the one-time scoring of the DNA (default min(40, cores-1)) |
| `--max-rows-per-gene` | TEST ONLY: random subset of rows for every plot |

## 2. What is read and what is computed

- **From `scores_all.xlsx` only (no scoring):** distributions, clusters, both heatmaps, win/draw/loss heatmap, win rate, tally, CAI/tAI scatter, ENC-vs-GC3. Better/worse directions and the GC band come from `codon_compare` (imported, not modified).
- **Needs the DNA:** positional profiles, motif positions, restriction sites per enzyme, codon / dinucleotide usage, agreement with GeneArt, host distances. `scores_all.xlsx` does not contain these, so `codon_features.compute_features` runs once per distinct sequence (same scorer, same `human_codon.csv` / `human_tgcn.csv`). Results are cached in `outputs/plots/profile_cache.pkl`; a second run does not score again.
- **Safeguard:** the scalar scores returned in that step are compared with `scores_all.xlsx` (sheet `Check_vs_scores_all`, differences must be 0; a warning is printed otherwise).
- **seqfold missing:** `TIS_accessibility_5p`, `MFE_min_window`, `stable_structure_count` have empty columns. They are skipped (and listed in the console) instead of being plotted, and the MFE profile panel is absent. Install seqfold and re-run `score_all_sequences.py` and this script to get them.
- `kozak_score` is plotted in `more_metrics/`; without `--utr5` it only sees the +4 base, so it is nearly binary.

## 3. Output files

| Folder / file | What it holds |
|---|---|
| `distributions/NN_<metric>.png` | one image per metric (22 with seqfold) |
| `clusters/NN_<metric>.png` | one image per metric |
| `more_metrics/distributions/`, `more_metrics/clusters/` | one image each for: `kozak_score`, `GC_min_60bp`, `GC_std_60bp`, `CpG_count`, `longest_repeat_bp`, `AT_pct`, `tAI_head`, `tAI_body`, `minmax_min`, `minmax_pct_rare_windows`, `MFE_mean_window` (needs seqfold), `tm_mean`. Same layout as the 22 main metrics; the audit ones are computed from the DNA with the baseline-audit definitions (60-nt non-overlapping windows, repeats >= 8 bp) |
| `extras/XNN_<name>.png` | GC1, GC2, ENC deviation, distances to the host (one image each) |
| `summary/` | `heatmap_geneart_cluster`, `heatmap_geneart_beats_dnachisel`, `win_heatmap_per_gene`, `win_rate_per_metric`, `tally_per_gene`, `scatter_cai_gc_tai_gc3`, `enc_vs_gc3` |
| `profiles/<gene>_profiles_GeneArt_vs_DNAChisel.png` | 12 positional panels per gene |
| `geneart_features/<gene>_GeneArt_profiles.png` | `codon_features` 9-panel figure of each GeneArt sequence |
| `liabilities/` | `motif_map_restriction_splice_homopolymer.png`, `restriction_sites_per_enzyme.png` |
| `host_reference/` | `codon_usage_per_1000_vs_host`, `gc_by_codon_position_vs_host`, `dinucleotides_vs_host` |
| `agreement/codon_agreement_with_geneart.png` | where along the CDS DNA Chisel uses the GeneArt codon |
| `pairs/`, `pair_features/` | per-pair overlay, scorecard and 9-panel figure (`--pair-plots`) |
| `summary_tables.xlsx`, `row_clusters.csv` | the numbers behind all plots (Section 5) |

## 4. How to read each plot

### Per-metric images (one image per score, all genes, all rows)

**`distributions/`**: one horizontal strip per gene, same x-scale for all genes.
- Grey dots: every DNA Chisel row of that gene (vertical position is random jitter only). `n` under the gene name is the row count.
- Box: quartiles (25-75%). Thick bar: median.
- Red diamond: the GeneArt sequence, value printed above it.
- Right column: GeneArt **percentile** among that gene's rows (0 = below all rows, 50 = median, 100 = above all rows, ties count half). `*` = GeneArt lies outside the whole DNA Chisel range.
- Title states the rule: higher / lower is better, closer to the GC-band midpoint (52.5%), or context only.
- Reference marks (dashed lines / shaded band, named in the legend): human GC band 40-65%, GC band midpoint, audit thresholds (CAI > 0.85, rare codons > 5%, GC max 60 bp > 70%, CpG > 6 per 100 bp, 5' hairpins >= 3, GC3 60-80%), and the **human host** value for GC_pct, GC3_pct and ENC (and GC1, GC2 in `extras/`).

**`clusters/`**: K-means (scikit-learn, 1-D, run separately inside each gene on its DNA Chisel rows).
- Dots coloured Low / Mid / High (colour is the **numeric range**, not good or bad; use the rule in the title). Dashed vertical line: boundary between clusters. White star: GeneArt.
- Right table: rows per cluster (and %) with the value range of each cluster, and **the cluster GeneArt falls into** (nearest cluster centre; `*` = outside the DNA Chisel range).
- Fewer than k distinct values: fewer clusters are made; one value (e.g. counts that are always 0): the panel says "single value" and GeneArt is reported as same / above / below.
- Identical DNA appears in several grid rows and keeps its weight in the cluster sizes (all rows are kept).

### `summary/`

- **`heatmap_geneart_cluster.png`** (genes x metrics): colour = cluster GeneArt falls into (Low / Mid / High range of that gene's DNA Chisel rows, grey = DNA Chisel has a single value); number = GeneArt percentile; `*` = outside the DNA Chisel range. Arrows after the metric name give the rule.
- **`heatmap_geneart_beats_dnachisel.png`**: cell = % of the gene's DNA Chisel rows that GeneArt beats on that metric (0 = worse than all, 50 = tie / middle, 100 = better than all). Last row = mean of genes. Only metrics with a better/worse rule.
- **`win_heatmap_per_gene.png`**: pair-wise result. Cell text, top to bottom: % of pairs DNA Chisel better / % draw / % GeneArt better. Colour = DNA Chisel % minus GeneArt % (blue = DNA Chisel ahead, red = GeneArt ahead). Last row = all genes pooled.
- **`win_rate_per_metric.png`**: per scored metric, share of all pairs won by DNA Chisel / draw / won by GeneArt.
- **`tally_per_gene.png`**: left, mean number of scored metrics per pair won by each side, per gene; right, spread over the pairs (one dot per pair, bar = median).
- **`scatter_cai_gc_tai_gc3.png`**: CAI against overall GC and tAI against GC3; dots = DNA Chisel rows, stars = GeneArt, colour = gene.
- **`enc_vs_gc3.png`**: ENC against GC3 with Wright's expected curve (Nc expected from GC3 alone) and the human host point. Points below the curve are more codon-biased than GC3 explains.

Winner rule (from `codon_compare`): higher / lower / band (GC% closer to 52.5%); ENC, tAI_ramp, minmax_mean, codon_autocorr_lag1 and GC3_pct are shown but not scored. All scored metrics count equally in the tallies; use the distributions for the size of a difference.

### `profiles/` and `geneart_features/`

Per gene, 12 panels along the CDS: GC (45 nt), MFE (40 nt, needs seqfold), %MinMax, rolling CAI, rolling tAI, rolling rare-codon %, GC3, windowed ENC, CpG density, hairpin density, autocorrelation of per-codon tAI, fragment Tm. Red = GeneArt; blue line = median of all DNA Chisel rows at that position; shaded bands = 25-75% and 5-95% of the rows. ACF and Tm panels are bars: GeneArt (red) and DNA Chisel median with quartiles (blue). Marks: human GC band, CAI 0.85, rare-codon 5%, stable-structure cutoff on the MFE panel. A red line leaving the bands marks a region where GeneArt differs from nearly all DNA Chisel rows. `geneart_features/` is the unmodified 9-panel `codon_features` figure of the GeneArt sequence.

### `liabilities/`

- **Motif map**: per gene, three rows (restriction, splice, homopolymer >= 6) along the CDS. Blue shading = % of the gene's DNA Chisel rows with a hit in that 30-nt bin; red diamonds = GeneArt hits (restriction hits labelled with the enzyme). The title gives the counts and the % of DNA Chisel rows with at least one hit.
- **Restriction sites per enzyme**: enzymes x genes, three panels: GeneArt sites, DNA Chisel mean sites per row, % of DNA Chisel rows with a site. Only enzymes with a hit are listed; the rest are named in the footer.

### `host_reference/`

Host values are the CoCoPUTs human data in `host_files/Human` (codon counts and GC1/2/3 from `Human.txt`, ENc 55.0, dinucleotides read from the SVG, normalised to 1000, +/-0.5). They are context, never a scored winner.
- **Codon usage per 1000**: host vs GeneArt (8 genes pooled) vs DNA Chisel (mean of the rows of each gene, genes pooled), grouped by amino acid.
- **GC by codon position**: overall, GC1, GC2, GC3 for host, GeneArt (dots = genes) and DNA Chisel (mean, whiskers = quartiles).
- **Dinucleotides**: total and codon-junction dinucleotides per 1000 for host, GeneArt (dots = genes) and DNA Chisel (mean, whiskers = quartiles).

### `agreement/`

Per gene, % of DNA Chisel rows that use exactly the GeneArt codon at each codon position (light = per codon, line = 9-codon average, dashed = gene mean). Shows where the two designs coincide (e.g. start region) and where they differ.

### `extras/`

Same image type as `distributions/`: GC1, GC2, ENC deviation (expected minus observed Nc, relative), and distances of each sequence to the host (L1 distance of codon usage, GC by position, ENC difference, dinucleotides). For the distances, lower = closer to the human host; context only.

## 5. `summary_tables.xlsx`

`GeneArt_vs_DNAChisel` (one row per metric x gene: GeneArt value, DNA Chisel min / median / max, percentile, % of rows GeneArt beats, cluster, outside-range flag), `Cluster_ranges`, `Win_counts`, `Pivot_percentile`, `Pivot_beats_pct`, `Pivot_cluster`, and with the DNA: `Extras_GeneArt_vs_DNAChisel`, `Extras_DNAChisel_rows` (per row: GC1, GC2, ENC deviation, host distances, % nucleotide identity and number of differing codons vs GeneArt), `Extras_GeneArt`, `Restriction_per_enzyme`, `Check_vs_scores_all`. `row_clusters.csv` gives the Low / Mid / High cluster of every DNA Chisel row for every metric.

## 6. Notes for interpretation

- The grid avoids all 34 restriction enzymes and the 4 splice motifs of the scorer (`DEFAULT_SITES` and `SPLICE_SITES` in `codon_grid_search_new.py` match the scorer's lists), so those two metrics are a fair comparison.
- DNA Chisel rows are the output of 2,592 settings per gene, not a random sample, and the codon metrics (CAI, rare codons, tAI) use the same `human_codon.csv` DNA Chisel optimised against.
- Percentiles and cluster sizes are over rows, so DNA that occurs in several grid rows counts several times.
- `hairpin_count_*` is a stem heuristic; compare between sequences, not against zero.
- Not included: PCA (separate suggestion) and sequence length (identical for all sequences of a gene).

## 7. Testing

The script was run end to end on synthetic data built with your scripts (3 genes, ~76 rows each, with and without a stand-in for seqfold, with and without the dinucleotide SVG, `--no-profiles`, `--clusters 4`, `--max-rows-per-gene`). It has not been run on your real files, so check the first `Check_vs_scores_all` sheet and the console line "Check against scores_all.xlsx: 0 metric(s) with differences".
