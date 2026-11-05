# Score plots from `scores_all.xlsx`

**Script:** `Sequence_scoring/plot_scores.py`
**Input:** `outputs/scores/scores_all.xlsx` (sheet `DNAChisel` = 19,788 Status-OK rows, sheet `GeneArt` = 8 rows)
**Output:** `outputs/score_plots/` (22 + 22 PNGs, 1 heatmap, 1 Excel)

Nothing is re-scored: every number is read from `scores_all.xlsx`. `codon_compare` is imported only for the metric order and the "better" direction. Metrics whose column is empty (e.g. seqfold not installed) are skipped and listed in the console.

## Run

```bash
cd Sequence_scoring
pip install numpy pandas matplotlib openpyxl scikit-learn
python plot_scores.py                      # 3 clusters per gene and metric
python plot_scores.py --k 4                # more clusters
python plot_scores.py --scores <path> --out-dir <path>
```

## Output files

| File | What it shows |
|---|---|
| `distributions/NN_<metric>.png` (22) | One image per metric. One horizontal strip per gene. |
| `clusters/NN_<metric>.png` (22) | One image per metric. One panel per gene, rows coloured by cluster, GeneArt overlaid. |
| `heatmap_geneart_position.png` | All 22 metrics x 8 genes: where GeneArt sits inside that gene's DNA Chisel distribution. |
| `cluster_summary.xlsx` | The numbers behind the cluster plots and heatmap (sheets below). |

## 1. Distribution plots (`distributions/`)

For one metric, each gene gets a horizontal strip (x axis = metric value, same scale for all genes so genes can be compared).

- **Blue dots:** every DNA Chisel row of that gene (`n` in the axis label). Vertical position is random jitter only.
- **Thick dark bar:** interquartile range (25-75%) of the DNA Chisel rows; **white tick:** median; **thin line:** min to max.
- **Yellow diamond:** the GeneArt sequence of that gene, with its value printed above.
- The title states the direction (higher / lower is better, or context only for the metrics `codon_compare` does not score).

Read it as: is the diamond inside, at the edge of, or outside the cloud of DNA Chisel results?

## 2. Cluster plots (`clusters/`)

For each metric and each gene separately, the DNA Chisel values are split by KMeans (scikit-learn, 1-D, `random_state=0`) into **Low, Mid, High** clusters, sorted by cluster centre.

- **Dots:** one per DNA Chisel row, coloured by cluster (blue = Low, yellow = Mid, red = High). Vertical position is random jitter.
- **Black dashed line:** the GeneArt value. The panel title says which cluster it falls into (the cluster whose centre is nearest), and flags "outside DNA Chisel range" when GeneArt is below the minimum or above the maximum of the rows.
- **Text under each panel:** value range and number of rows of each cluster.

Notes:
- Low/Mid/High refer to the **value**, not to good/bad. Use the direction in the title: for `CAI` High is good, for `CpG_density_per100bp` Low is good.
- Clustering is done per gene because genes have different baselines.
- If a gene has only 2 distinct values, 2 clusters (Low/High) are made; if only 1, the panel shows "Single value" (typical for counts that are always 0).
- Several grid rows can have identical DNA and therefore identical scores; all rows are kept, so such sequences weigh more in the cluster sizes.

## 3. Heatmap (`heatmap_geneart_position.png`)

Rows = metrics, columns = genes. Each cell answers: **"Compared with this gene's DNA Chisel rows, is the GeneArt value low, middle or high?"**

- **Colour and top number = percentile of GeneArt** among that gene's DNA Chisel rows (mid-rank for ties). 0% (dark blue) = GeneArt is below all DNA Chisel rows, 50% (white) = at the median, 100% (dark red) = above all rows.
- **Bottom number = the GeneArt value** itself.
- The label of each row shows the direction: `[higher=better]`, `[lower=better]`, `[closest to 52.5]` for GC_pct, `[context]` for metrics that are not scored.
- Colour is comparable only within a metric (across genes), not between different metrics. A 50% cell for a metric where all values are identical (e.g. repeat_count_20bp = 0) means "tie with every row".

## 4. `cluster_summary.xlsx`

| Sheet | Content |
|---|---|
| `Cluster_detail` | One row per gene x metric: GeneArt value, GeneArt cluster, percentile, outside-range flag, DNA Chisel min / median / max, range and row count of each cluster |
| `GeneArt_cluster_matrix` | Metric x gene table of the cluster GeneArt falls into |
| `GeneArt_percentile` | Metric x gene table of the percentile used in the heatmap |

## Not included

PCA (suggested separately) is not part of this script.
