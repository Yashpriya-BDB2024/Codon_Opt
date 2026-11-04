# GeneArt (WT) vs top-10 DNA Chisel: scoring, comparison and plots

**Script:** `Sequence_scoring/run_comparison.py`
**Input:** `human_strictGC_top10.xlsx` (80 sequences) + the GeneArt sequences in the trials workbook (8 genes)
**Output:** `Sequence_scoring/outputs/comparison/` with `comparison_results.xlsx` and 265 PNG plots

`codon_compare.py` and `codon_features.py` are **not modified**. The script imports them and only adds what is missing: the batch loop, a QC gate, the host-reference plots, summary plots and a few extra columns.

---

## 1. How to run

```bash
cd Sequence_scoring
pip install numpy pandas matplotlib openpyxl seqfold     # seqfold is needed for 3 RNA-structure metrics
python run_comparison.py                                  # all 80 pairs at once
python run_comparison.py --workers 4 --no-pair-plots      # faster: skips the 248 per-sequence PNGs
```

Defaults: top file `../DNAChisel/outputs/human_strictGC_top10.xlsx`, trials file found by glob in the parent folder, host files from `../host_files/Human/`. Add `--utr5 GCCACC` only if you want a Kozak −6..−1 context (it is applied to both sides).

**Run it with seqfold installed.** The numbers in this document were produced in an environment where seqfold could not be installed, so `TIS_accessibility_5p`, `MFE_min_window` and `stable_structure_count` are blank and counted as "na" (20 of the 23 scored metrics are therefore tallied here, and the MFE panel is missing from the overlay plots). Everything else is unaffected. Re-run on your machine to fill these three in; the Read_me sheet records whether seqfold was present.

---

## 2. What is compared with what

| | Label | Source |
|---|---|---|
| A | `WT` | GeneArt sequence, column "DNA optimized sequence" of the trials workbook (this is the "WT" in the scoring guide's usage example) |
| B | `DNAChisel` | each of the 10 rows per gene in the `Top` sheet (rank 1 = highest CAI) |

For every gene, **each of the 10 sequences is compared with that gene's own GeneArt sequence** (8 genes × 10 = 80 pairs, genes are never mixed). Both sides use the same `human_codon.csv` and `human_tgcn.csv`, the same scorer and the same enzyme list (`codon_features.RESTRICTION_SITES`, 34 enzymes).

**QC gate (stops the run if it fails):** every DNA Chisel sequence translates to the input protein, has the same length as GeneArt and has no internal stop. Passed for all 80. Nucleotide identity to GeneArt is 79.9–89.0% (mean 86.3%), i.e. 93–422 codons differ per pair.

**What is reused from the given scripts**

| Given function | Used for |
|---|---|
| `codon_compare.compare()` | 23 scored metrics, deltas and winners for every pair |
| `codon_compare.overlay_profiles()`, `scorecard()` | per-pair figures, exactly as described in the scoring guide |
| `codon_features.plot_profiles()` | the 9-panel single-sequence figure (adds autocorrelation and fragment Tm) |
| `codon_features.scan_positions()`, `codons_of()`, … | extra columns below |

---

## 3. Metric coverage

All 23 metrics in `COMPARE_METRICS` are compared, and **every scalar that `compute_features` returns** (about 70, including the per-enzyme `RS_*` counts, Kozak flags and Tm statistics) is written to the `All_scores` sheet, so nothing in the panel is dropped.

Winner rules are those of `codon_compare`: higher, lower, band (GC% closer to the human band midpoint 52.5) or na (shown with a delta but not scored: ENC, tAI_ramp, minmax_mean, codon_autocorr_lag1, GC3_pct).

Extra columns added (context, never counted as wins):

| Column | Why |
|---|---|
| `RS_in_grid_enzymes(13)` / `RS_in_other_enzymes(21)` | splits the restriction count by whether the grid told DNA Chisel to avoid that enzyme (Section 6) |
| `GC1`, `GC2` | GC by codon position, to compare with the host chart in `host_files/Human` |
| `codon_usage_L1_to_host`, `GC_pos_L1_to_host`, `ENC_abs_diff_to_host`, `dinuc_total_L1_to_host`, `dinuc_junction_L1_to_host` | distance to the human host reference (Human.txt counts, CoCoPUTs ENc 55.012, dinucleotide chart) |
| `pct_identity_nt`, `codons_differing` | how far each sequence is from GeneArt |

---

## 4. Output files

| File / folder | Content |
|---|---|
| `comparison_results.xlsx` | `Read_me`, `Overall_win_rate` (per metric over 80 pairs), `Gene_summary` (WT vs DNA Chisel mean/min/max and win counts per gene × metric), `Pair_summary` (wins, QC, identity, variant, codon method, per pair), `Long_all_pairs` (80 × 23 rows), `All_scores` (88 rows, every scalar) |
| `plots/pairs/` (160) | `<gene>_rankNN_overlay_profiles.png` and `_scorecard.png` from `codon_compare` |
| `plots/features/` (88) | `codon_features` 9-panel figure for each GeneArt and each DNA Chisel sequence |
| `plots/genes/` (8) | per gene: GeneArt vs all 10 DNA Chisel sequences on the 7 positional profiles |
| `plots/summary/` (5) | `metrics_by_gene` (every metric, every gene), `win_heatmap`, `win_rate_per_metric`, `tally_per_gene`, `cai_vs_gc_tai_vs_gc3` |
| `plots/host_reference/` (4) | codon usage per 1000 vs host, GC overall/GC1/GC2/GC3 vs host, dinucleotides (total and junction) vs host, ENC vs host ENc |

The host-reference plots follow the charts supplied in `host_files/Human` (codon frequencies per 1000, GC% by position, dinucleotide frequencies, effective number of codons). Codon and GC values come from the exact counts in `Human.txt`. Dinucleotide host values are read from the SVG bars and normalised to 1000, so treat them as ±0.5 per 1000.

---

## 5. Results (human, 80 pairs, 20 scored metrics without seqfold)

Per-pair wins: DNA Chisel 3.8–5.7, GeneArt 3.4–6.0 per gene, the rest draws. Neither side wins overall; they win different metrics.

| Metric | DNA Chisel wins | GeneArt wins | Draw | Reading |
|---|---|---|---|---|
| CAI | 0 | 80 | 0 | GeneArt 0.93–0.97 vs DNA Chisel mean 0.88–0.96 (Section 7) |
| tAI | 36 | 44 | 0 | close; PREP 10/10 to DNA Chisel, CTNNB1 and FOXA1 10/10 to GeneArt |
| rare_codon_pct | 0 | 8 | 72 | 8 sequences carry codons below w = 0.3 (CTNNB1 ranks 9–10, ALDH1A3 ranks 6–10, ERG rank 9), at most 0.98%; the rest are 0 for both sides |
| GC_pct (band) | 80 | 0 | 0 | by construction: DNA Chisel is held to 45–55%, GeneArt is 56.5–57.8%; both are inside 40–65% |
| GC_max_60bp | 55 | 0 | 25 | |
| CpG density | 80 | 0 | 0 | DNA Chisel 2.5–4.8 vs GeneArt 3.6–6.2 per 100 bp |
| hairpin_count_all / 5′ | 48 / 52 | 27 / 14 | 5 / 14 | heuristic count; probably helped by the 10-mer uniqueness constraint |
| restriction_sites_total | 18 | 41 | 21 | **explained by the enzyme list, Section 6** |
| splice_sites_total | 0 | 52 | 28 | GeneArt has 0–1 motifs; the grid has no splice constraint |
| tm_spread | 9 | 71 | 0 | fragment Tm is more uniform in GeneArt |
| homopolymer ≥6, repeats ≥10/≥20, Kozak | 0 | 0 | 80 each | all identical (Kozak sees only the +4 base without `--utr5`) |

GeneArt CAI per gene against the best and mean DNA Chisel value: CTNNB1 0.930 vs 0.898 / 0.876, IRF4 0.960 vs 0.942 / 0.940, FOXA1 0.928 vs 0.912 / 0.905, ERG 0.965 vs 0.955 / 0.947, AURKB 0.952 vs 0.938 / 0.934, AURKA 0.952 vs 0.952 / 0.943, PREP 0.968 vs 0.956 / 0.955, ALDH1A3 0.960 vs 0.941 / 0.903.

---

## 6. Is the restriction-enzyme list the reason GeneArt wins? Partly, and the evidence is in the data

- The scorer counts **34 enzymes** (`RESTRICTION_SITES`) in **both** sequences. The grid told DNA Chisel to avoid only **13** (`DEFAULT_SITES` in `codon_grid_search.py`: BsaI, BsmBI, BbsI, SapI, EcoRI, XhoI, NotI, XbaI, SpeI, KpnI, SalI, AscI, PacI). All 13 are inside the 34. The other 21 (AarI, AflII, AgeI, BamHI, BglII, ClaI, DraI, EcoRV, FseI, HindIII, HpaI, MluI, NcoI, NdeI, NheI, PmeI, PstI, SacI, SbfI, SmaI, SwaI) were never constrained.
- Result: both sides have **0 sites in the 13 grid enzymes** (80/80 DNA Chisel sequences and 8/8 GeneArt sequences). **Every counted site is in the other 21.** DNA Chisel averages 6.9 of them per sequence (range 1–14), GeneArt 3–11 depending on the gene (mean 5.6). The most frequent offenders in DNA Chisel are PstI, BglII, NcoI, SacI, BamHI, NheI and AarI.
- So the metric is not unfair in the sense of a hard-coded advantage, since the same list is applied to both. It does measure something DNA Chisel was never asked to do. It is unlikely by chance that GeneArt has zero sites for all 13 grid enzymes in all 8 genes (several are expected per gene), which suggests GeneArt avoids a standard list; I could not verify which list it uses.
- **Sensitivity check:** among the distinct sequences that pass the strict-GC filter (62–209 per gene), none has zero sites across all 34 enzymes, and the minimum is 1–4 per gene. Even with a 40–60% window only 2 sequences in the whole pool are clean. The grid therefore never explored this region, and no filter on the existing output can fix it.
- **Fix:** put the 34 enzyme names into `DEFAULT_SITES` (the grid reads `config.get("restriction_sites", DEFAULT_SITES)`; DNA Chisel resolves `<enzyme>_site` names through Biopython, so confirm each name loads) and re-run. Then the restriction comparison becomes a fair test.

Literature check (what the published descriptions imply):

- DNA Chisel removes the patterns you list; its documentation example shows `AvoidPattern("BsaI_site")` removing BsaI sites on both strands, and nothing else (Zulkower & Rosser 2020, *Bioinformatics* 36:4508).
- GeneArt's GeneOptimizer is a multi-parameter optimizer. The published description lists codon usage, GC content, exclusion of DNA motifs and repetitions, and avoidance of stable structures among its criteria, and names restriction sites and splice-site patterns as local properties it handles (Raab et al. 2010, *Syst Synth Biol* 4:215).
- This matches what the data show: GeneArt wins splice sites and restriction totals, which DNA Chisel was not set up to control, and DNA Chisel wins GC window and CpG, which the grid did control.

---

## 7. Why GeneArt wins CAI, and why it is not a tool failure

The lead's rule fixes overall GC at 45–55%, and GeneArt sits at 56.5–57.8%. The best human codons are GC-rich, so the window caps CAI. Using the grid output directly:

| Window | Best pool CAI per gene | GC of that sequence | vs GeneArt CAI |
|---|---|---|---|
| strict 45–55% | 0.898–0.956 | 54.7–55.0% | below GeneArt in all 8 genes (AURKA is a tie: 0.9515 vs 0.9518) |
| wide 40–60% | 0.947–0.991 | 58.6–60.0% | **above GeneArt in all 8 genes** |

So DNA Chisel can beat GeneArt on CAI, but only by going to GC 59–60%, which is further from the host than GeneArt. The comparison therefore reflects the strict-GC choice (option 1), not the optimizer.

---

## 8. Caveats to state when presenting these results

1. **Selection bias.** DNA Chisel rows are the best of 2,592 settings per gene chosen on CAI; GeneArt is a single vendor sequence.
2. **Circularity.** CAI, rare-codon % and tAI use the same `human_codon.csv` DNA Chisel optimised against.
3. **Different rare-codon definitions.** The grid cut-off is a codon frequency of 0.11; the scorer uses w = f/f_max < 0.3. A sequence can pass the first and still count some rare codons in the second (ALDH1A3 up to 0.98%, always below 5%).
4. **Hairpins and splice sites were never constrained** in the grid, so wins or losses there are incidental.
5. **GC_pct "wins" are a band-midpoint rule.** Both sides are inside the 40–65% human band.
6. **Both sets are far more codon-biased than the host.** ENC is about 31 (DNA Chisel) and 31.5 (GeneArt) against the host's 55.0, with GC3 65% and 71% (see `host_reference` plots).
7. **Kozak does not discriminate** without `--utr5`: only the +4 base is scored.
8. **Weights are equal.** The tally counts each metric once; a CAI loss counts the same as a CpG win. Use `Gene_summary` for the magnitude of each difference.
9. The three seqfold metrics are not in the numbers above (Section 1).

---

## 9. Checks performed

- `filter_top_strict_gc.py` re-run on `human_grid.xlsx` reproduces `human_strictGC_top10.xlsx` exactly (same 80 variants, same order, same DNA).
- All 80 DNA Chisel sequences and all 8 GeneArt sequences translate to the input protein (the trials sequences have no stop codon, so none is stripped).
- 8 genes × 10 pairs = 80 comparisons; every pair uses that gene's own GeneArt sequence.
- Host GC1/GC2/GC3 computed from `Human.txt` counts (55.8 / 42.4 / 55.4%) agree with the values read from the host GC% chart (about 55.5 / 42.2 / 55.1%).
- `codon_compare.py` and `codon_features.py` are imported as supplied; no line was changed.
