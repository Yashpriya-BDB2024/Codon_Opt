# plot_ncbi.py: plots of the NCBI reference CDS (NCBI only)

**Script:** `Sequence_scoring/plot_ncbi.py` (put it next to `plot_all.py`; it only borrows its drawing helpers, unmodified)
**Output:** `outputs/plots_ncbi/`

Only the NCBI sequences are used. No GeneArt, no DNA Chisel: nothing of them is read, drawn or used for axes, percentiles or clusters.

## 1. Run

```bash
cd Sequence_scoring
python plot_ncbi.py
python plot_ncbi.py --ncbi-scores outputs/scores/scores_ncbi_1.xlsx --trials "../Codon_optimization_trials_27102026.xlsx"
```

| Input | Used for | Default location |
|---|---|---|
| `scores_ncbi*.xlsx` (sheet `NCBI`) | the metric images | `outputs/scores/scores_ncbi*.xlsx` (else the project folder) |
| trials workbook, column "NCBI sequence" | positional profiles, motif map | newest `Codon*ptimization*trials*.xlsx` (`--trials` to choose) |

Optional switches (not needed for the normal run): `--no-profiles` (metric images only, no DNA needed), `--no-host` (no human-host lines and no `host_reference/` plots), `--dpi` (image resolution).

## 2. Output files and how to read them

| File | What it shows |
|---|---|
| `metrics/NN_<metric>.png` | **one image per score in `scores_ncbi`** (22; metrics with no data, e.g. the 3 seqfold metrics without seqfold, are skipped and listed). One row per gene, green square = NCBI value (printed above it), same x-scale for all genes. The title gives the rule (higher / lower is better, GC band). Reference marks: human GC band 40-65% and midpoint, audit thresholds (CAI > 0.85, rare codons > 5%, GC max 60 bp > 70%, CpG > 6 per 100 bp, 5' hairpins >= 3, GC3 60-80%), and the human host value for GC, GC3 and ENC. |
| `profiles/<gene>_profiles_NCBI.png` | **positional profiles, one figure per gene**, 12 panels along the CDS: GC (45 nt), MFE (40 nt, needs seqfold), %MinMax, rolling CAI, rolling tAI, rolling rare-codon %, GC3, windowed ENC, CpG density, hairpin density, autocorrelation of per-codon tAI, fragment Tm. Marks: human GC band, CAI 0.85, rare codons 5%, stable-structure cutoff on the MFE panel. A dip in rolling CAI / tAI or a spike in GC marks a region to inspect. |
| `liabilities/motif_map_NCBI.png` | where the sequence liabilities sit: restriction sites (enzyme named), splice motifs and homopolymers >= 6 along each CDS. This is the positional view of three of the 22 scores. |
| `host_reference/codon_usage_per_1000_vs_host.png` | codon usage per 1000 codons grouped by amino acid: human host (CoCoPUTs) vs NCBI (8 genes pooled). |
| `host_reference/gc_by_codon_position_vs_host.png` | overall GC, GC1, GC2, GC3 for the host and NCBI (black dots = each gene). |
| `host_reference/dinucleotides_vs_host.png` | total and codon-junction dinucleotides per 1000, host vs NCBI (dots = each gene). Needs the dinucleotide SVG in `host_files/Human`. |
| `summary/enc_vs_gc3.png` | ENC against GC3 with Wright's expected curve; a star per gene and the human host point. Points below the curve are more codon-biased than GC3 explains. |
| `ncbi_tables.xlsx` | `NCBI_scores`, `Check_vs_scores_ncbi`, `Restriction_per_enzyme`, `Motif_hits` (position and name of every hit), `GC_by_position` |

That is 22 + 8 + 1 + 3 + 1 images for 8 genes. Nothing else is needed to compare the three sources: the same metric images and the same profile / motif figure types exist for GeneArt and DNA Chisel in `plot_all.py`.

## 3. Notes

- **Why compare human sequences with the human host?** The host reference is the genome-wide average codon usage (CoCoPUTs), not these genes. A natural gene can sit away from that average (e.g. in GC3), so the plots show how typical the NCBI genes are, a baseline for how far GeneArt and DNA Chisel moved from average human usage.
- **Axes are scaled to the NCBI values only**, so they are not identical to the GeneArt / DNA Chisel images; read the printed values and the reference lines (same for every image) when comparing.
- **Stop codon.** The NCBI sequences in the trials workbook end with a stop codon (GeneArt and DNA Chisel do not). The script tries the sequence as written and with the stop removed, and uses the form that reproduces `scores_ncbi` exactly (printed in the console, sheet `Check_vs_scores_ncbi`). If neither matches, check that `host_files/Human/human_codon.csv`, `human_tgcn.csv` and seqfold are the ones used for `scores_ncbi`.
- Only the 8 human genes are used; E. coli / CHO rows of the trials workbook are ignored. The script warns if an NCBI sequence does not translate to its input protein.

## 4. Testing

Run end to end on synthetic data (3 genes) built with your scripts. It has not been run on your real files; check the console line "NCBI sequence ...: 0 metric(s) differ from the score table".
