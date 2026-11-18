# Drug Combination PowerPoint Report

`test_report_3.py` turns the output of `drug_synergy_2.py` into a PowerPoint file:
one summary slide plus one slide per drug pair, each with a table, a network picture and the classification.

It is the adapted version of `test_report_2.py` and follows the corrected pipeline (targets only, no neighbour expansion, seeded p-values).

## Files

| File | Role |
|---|---|
| `drug_synergy_2.py` | Run first. Produces `Drug_combination_results.xlsx` |
| `test_report_3.py` | This report script |
| `Drug_target.xlsx` | Input of `drug_synergy_2.py`. The report reads its `Module_genes` sheet |
| `Drug_combination_results.xlsx` | Output of `drug_synergy_2.py`. The report reads its `Drug_combination_results` sheet |
| `STRING_PPI_filtered.txt` | Tab-separated network with columns `source`, `target`, `Weight` |

`shortest_paths_cache.pkl` and the precompute script are **not** used by the report.

## Requirements

```
pip install pandas networkx plotly kaleido python-pptx openpyxl
```

`kaleido` is needed for `fig.write_image`. Newer kaleido versions also need Chrome installed. If image export fails, check this first.

## How to run

1. Run `drug_synergy_2.py` so `Drug_combination_results.xlsx` is up to date.
2. Put all files from the table above in one folder.
3. Run:

```
python test_report_3.py
```

Output:
- `Drug_Combination_Report.pptx`
- `report_images/`, one PNG per pair (can be deleted afterwards)

## What each slide shows

**Summary slide** (switch off with `SUMMARY_SLIDE = False`): all pairs with disease, drug names, Z (A), Z (B), Sab and classification.

**Pair slide:**
- **Table:** drug names, targets, z-score with p-value for each drug, and Sab.
  - Green z cell: z < 0 and p < 0.05 (close to the disease).
  - Amber z cell: z < 0 but p >= 0.05 (borderline; the same case `drug_synergy_2.py` warns about in its log).
  - p is shown as `<0.001` for the smallest value, since p cannot go below 1/1001.
- **Network picture:**
  - Red: disease module genes found in STRING.
  - Blue: Drug A targets. Orange: Drug B targets. Large and labelled. A red outline means the target is also a disease gene.
  - Grey: genes on one shortest path (in hops) from each target to its nearest disease gene.
- **Classification box:** green for Complementary Exposure, yellow for Overlapping Exposure, grey for anything else (for example an error).

## Settings (top of the script)

| Name | Default | Meaning |
|---|---|---|
| `NETWORK_FILE` | `STRING_PPI_filtered.txt` | Network file |
| `RESULTS_FILE` | `Drug_combination_results.xlsx` | Results from `drug_synergy_2.py` |
| `INPUT_FILE` | `Drug_target.xlsx` | Source of module genes |
| `OUT_PPTX` | `Drug_Combination_Report.pptx` | Output file |
| `ALPHA` | `0.05` | Cut-off for the green/amber colouring. Keep it equal to `alpha` in `drug_synergy_2.py` |
| `SUMMARY_SLIDE` | `True` | Add the summary slide |

## Notes and limitations

- **Hops, not weights:** distances are hop counts, as in the cache used by `drug_synergy_2.py`. `Weight` or `combined_score` is not used.
- **Ties:** if several shortest paths have equal length, only one is drawn. The distance is the same, the intermediate genes may differ.
- **Targets only:** the picture shows the same sets the numbers were calculated on, with no neighbour expansion.
- **Single-target drugs** (BRAF, EGFR) have d_AA or d_BB = 0 by convention. This is a limitation of the analysis, not of the report.
- **Missing genes:** targets not in STRING are listed on the slide. Disease genes not in STRING are not drawn (they are listed in the `drug_synergy_2.py` log).
- **Failed pairs:** rows whose classification starts with `Error` get a slide with a message instead of a picture.
- **Disease names must match:** the `Disease` value in the results must exactly match a column header in `Module_genes`, otherwise that slide is skipped with a warning.
- **Testing status:** the script was run on your real `STRING_PPI_filtered.txt` (14,301 nodes, 125,344 edges) and your real result and input files, with plotly replaced by a stub. Checked: the node counts (469 of 485 Melanoma and 532 of 559 CRC genes found in STRING, as in your review), the Sab values (0.5 and 2.0, recalculated by hand from the network), the hop distances behind the connector paths, and the slide layout. Not checked: the plotly/kaleido image export itself, which needs a run on your machine (about 2 s for the network and layout steps).
- **Targets that are disease genes:** in your current data every target (BRAF, MAP2K1, MAP2K2, EGFR) is also in its disease module, so distance is 0 and no grey connector path is drawn. Grey nodes appear only when a target is outside the module.

## Changes compared with `test_report_2.py`

All changes are marked `# [CHANGED n]` in the code. In short:

1. Module genes are read from `Drug_target.xlsx`, not from the results file.
2. Columns `Drug_A` and `Drug_B` replace `Drugs`. `P_value_separation` no longer exists.
3. The graph uses targets only plus shortest paths to the disease module, instead of targets plus all neighbours.
4. The table shows z with p, and Sab. The classification box is coloured by class.
5. Layout fixes: no empty title placeholder, and the picture no longer overlaps the classification box.
6. New summary slide, and failed pairs no longer crash the script.
7. Connector path fix: a target that is itself a disease gene gets no path (distance 0), consistent with your proximity calculation.
