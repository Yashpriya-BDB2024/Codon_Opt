# Score table: DNA Chisel sequences and GeneArt

**Script:** `Sequence_scoring/score_all_sequences.py`
**Input:** the `Results` sheet of the grid output (`human_grid.xlsx`) and the original trials workbook (GeneArt sequences)
**Output:** `Sequence_scoring/outputs/scores/scores_all.csv` and `scores_all.xlsx`

`codon_features.py` and `codon_compare.py` are imported and not modified. Every value in the table is read from `codon_features.compute_features()['scores']`, so there is one definition of each metric, the same one used by `codon_compare`.

## 1. How to run

```bash
cd Sequence_scoring
pip install numpy pandas openpyxl seqfold tqdm
python score_all_sequences.py --grid ../DNAChisel/outputs/human_grid_1.xlsx
python score_all_sequences.py --workers 40        # default: min(40, cores-1), maximum 61
```

If more than one trials workbook sits in the project folder, pass the original one with `--trials`. Only the first four columns of the trials sheet are read (name, protein, host, GeneArt DNA).

## 2. What is in the table

| | |
|---|---|
| Rows | every grid row with `Status = OK` (19,788) plus the GeneArt sequence of each gene (8). Failed and `Constraints not met` rows are left out because they have no usable sequence or did not meet the constraints. |
| Columns | `Name`, `Source` (`DNAChisel` or `GeneArt`), `Variant` (the variant name from the grid, or `GeneArt`), then the 22 metrics below |
| Excel sheets | `DNAChisel` (19,788 rows) and `GeneArt` (8 rows); the CSV has both together |

The 22 metrics are exactly the comparison panel of `codon_compare.py` without Kozak. Each distinct DNA is scored once and the result copied to every row with that DNA: identical DNA always gives identical scores, so this only saves run time (19,788 rows have 7,160 distinct sequences). Before scoring, the script checks that every sequence translates to its input protein.

## 3. The 22 metrics

All use the human codon table `human_codon.csv` and the tRNA table `human_tgcn.csv`. "Better" is the rule `codon_compare` uses to pick a winner; **na** means it is shown but not scored as a win.

### Codon usage

| Column | What it measures and how | Better |
|---|---|---|
| `CAI` | Codon Adaptation Index. For each codon, w = its frequency ÷ the frequency of the most-used codon for the same amino acid (from `human_codon.csv`). CAI = geometric mean of w over all codons, leaving out Met and Trp (one codon only). Range 0–1; 1 = every codon is the host's favourite | higher |
| `rare_codon_pct` | % of all codons whose w is below 0.3 | lower |
| `ENC` | Effective number of codons (Wright 1990), from the sequence's own codon counts. Nc = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6, where F is the probability that two codons of one amino acid are identical, averaged within each degeneracy class. 20 = strongest bias, 61 = no bias | na |
| `tAI` | tRNA adaptation index (dos Reis 2004). Each codon gets a weight from the tRNA gene copy numbers in `human_tgcn.csv`, including wobble pairing; tAI = geometric mean of the weights (Met left out) | higher |
| `tAI_ramp` | Mean tAI (geometric) of the codons after the first 30 minus that of the first 30 codons. Positive = slow start (a "ramp") | na |
| `minmax_mean` | %MinMax (Clarke & Clark 2008): in each sliding 18-codon window the summed codon frequencies are compared with what the average and the rarest/commonest codons would give; +100 = all most-frequent codons, −100 = all rarest. The column is the mean over all windows | na |
| `codon_autocorr_lag1` | Autocorrelation at lag 1 of the per-codon tAI weights: do fast or slow codons sit next to similar ones. −1 to +1 | na |

### Composition

| Column | What it measures and how | Better |
|---|---|---|
| `GC_pct` | (G + C) ÷ length × 100 | band: closer to the middle of the human band 40–65% (52.5%) wins |
| `GC_max_60bp` | The sequence is cut into non-overlapping 60-nt blocks; the GC% of the highest block | lower |
| `GC3_pct` | GC% at the third base of each codon | na |
| `CpG_density_per100bp` | Number of "CG" dinucleotides ÷ length × 100 | lower |

### RNA structure (needs `seqfold`)

| Column | What it measures and how | Better |
|---|---|---|
| `TIS_accessibility_5p` | Folding free energy ΔG (kcal/mol, 37 °C, `seqfold`) of the first 48 nt. Less negative = start region more open | higher |
| `MFE_min_window` | ΔG of every 40-nt window, stepping 12 nt; the lowest (most negative, most stable) value | higher |
| `stable_structure_count` | How many of those windows have ΔG below −15 kcal/mol | lower |
| `hairpin_count_all` | Heuristic. For each 6-nt stem along the sequence, count 1 if its reverse complement appears within the next 200 nt, after at least a 3-nt loop. Counts stem positions, not distinct hairpins | lower |
| `hairpin_count_5prime` | The same count on the first 100 nt only | lower |

### Sequence liabilities

| Column | What it measures and how | Better |
|---|---|---|
| `splice_sites_total` | Occurrences (forward strand) of GTAAG, GTGAG, TTTTTAG and TTTCAG | lower |
| `homopolymer_runs_ge6` | Number of runs of 6 or more identical bases | lower |
| `restriction_sites_total` | Occurrences of the recognition sites of the 34 enzymes in the scorer's list, both strands (a palindromic site is counted once) | lower |
| `repeat_count_10bp` | Number of distinct 10-nt pieces that occur at least twice | lower |
| `repeat_count_20bp` | The same for 20-nt pieces | lower |

### Synthesis

| Column | What it measures and how | Better |
|---|---|---|
| `tm_spread` | The sequence is cut into 60-nt fragments (a short tail is merged into the last one); each fragment's melting temperature is calculated with nearest-neighbour parameters (SantaLucia 1998, 50 mM Na⁺, 250 nM strand). Spread = highest − lowest Tm (°C); a small spread means the fragments are balanced | lower |

## 4. Notes

- The three structure metrics (`TIS_accessibility_5p`, `MFE_min_window`, `stable_structure_count`) are filled when `seqfold` is installed. The script prints a warning if it is not found; nothing is switched off in the code.
- The grid's enzyme list (34) and splice motifs (4) match the scorer's, so `restriction_sites_total` and `splice_sites_total` are a fair comparison.
- The scorer's hairpin count is a simple stem heuristic and is not the same definition as the grid's `AvoidHairpins` constraint; compare it between sequences, not against zero.
- The table's `CAI` is the scorer's. It differs slightly from the grid's own `CAI` column (up to 0.016) because the scorer leaves out Met and Trp.
