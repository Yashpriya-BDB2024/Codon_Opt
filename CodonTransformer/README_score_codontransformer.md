# `score_codontransformer.py` — README

Score table for the **CodonTransformer** grid: every `Status == OK` sequence of the grid workbook plus the
**GeneArt** sequence of each gene, scored on the same 22 metrics as the DNA Chisel and Codon Harmony tables.
CodonTransformer rows are never mixed with DNA Chisel or Codon Harmony rows.

---

## 1. Relation to the other scoring scripts

`score_codontransformer.py` is `score_codonharmony.py` with **only** these edits (verified by diff):

| What | Codon Harmony script | This script |
|---|---|---|
| `SOURCE` label | `CodonHarmony` | `CodonTransformer` |
| Default `--grid` | `CodonHarmony/outputs/human_codonharmony_grid.xlsx` | `CodonTransformer/outputs/human_codontransformer_grid.xlsx` |
| Output files | `scores_codonharmony.csv/.xlsx` | `scores_codontransformer.csv/.xlsx` |
| Sheet name | `CodonHarmony` | `CodonTransformer` |
| Messages / docstring | Codon Harmony | CodonTransformer |

Scoring is **unchanged**. `codon_features.py`, `codon_compare.py` and `score_all_sequences.py` are imported, not
modified; every metric comes from `score_all_sequences.score_one` → `codon_features.compute_features()`, with the
same `human_codon.csv` and `human_tgcn.csv` as the other tools.

---

## 2. Inputs

| File | Needed |
|---|---|
| `human_codontransformer_grid.xlsx` (from `run_human_parallel_codontransformer.py`) | sheet `Results` with `Name`, `Variant`, `Status`, `Optimized DNA Sequence` |
| Original trials workbook (`Codon*ptimization*trials*.xlsx`) | GeneArt sequences + input proteins (human rows) |
| `host_files/Human/human_codon.csv`, `human_tgcn.csv` | reference data, same as the other tools |
| `score_all_sequences.py`, `codon_features.py`, `codon_compare.py` | same folder |

`pip install numpy pandas openpyxl seqfold tqdm` (without `seqfold`, `TIS_accessibility_5p`, `MFE_min_window` and
`stable_structure_count` are empty and a warning is printed).

---

## 3. Run

```bash
python score_codontransformer.py --grid path/to/human_codontransformer_grid.xlsx
python score_codontransformer.py --grid ... --workers 40      # default min(40, cores-1), max 61
python score_codontransformer.py --grid ... --trials path/to/Codon_optimization_trials_27102026.xlsx
```

Do **not** use `--stop-policy asis` (see §4). The default `strip` is correct.

---

## 4. Terminal stop codon (important for CodonTransformer)

CodonTransformer appends a stop codon to its output: in the delivered grid, **640 of 640** sequences end in a stop
(length = 3 × aa + 3). The GeneArt sequences have none (0 of 8).

`compute_features` strips a trailing stop only for codon-level metrics (CAI, tAI, ENC …); nucleotide-level metrics
(GC, CpG, structure, hairpins, repeats, restriction/splice counts, Tm fragments) use the DNA as given. With the
default `--stop-policy strip`, one terminal stop codon is removed from **every** sequence of **both** sources before
scoring, so all sequences cover the same nucleotides. The script then checks that each sequence is exactly
3 × protein length and stops the run if not. `asis` would score CodonTransformer over 3 extra nt and is not
comparable with GeneArt.

---

## 5. What the script does

1. Reads `Results`, keeps `Status == OK`, cleans the DNA (`codon_features.clean_seq`).
2. Loads GeneArt DNA and proteins for the human genes from the trials workbook.
3. Reports terminal stops per source, applies the stop policy.
4. **Safeguards (the run exits on failure):** every gene in the grid exists in the trials workbook; every sequence
   translates to its input protein; every sequence is exactly 3 × protein length (strip policy).
5. Scores each **distinct** DNA once (identical DNA gives identical scores) and copies the result to every row
   with that DNA.
6. Writes the tables.

---

## 6. Output

`outputs/scores/scores_codontransformer.csv` (all rows) and `scores_codontransformer.xlsx` (sheets
`CodonTransformer`, `GeneArt`). Columns: `Name`, `Source`, `Variant` + 22 metrics:

`CAI, rare_codon_pct, ENC, tAI, tAI_ramp, minmax_mean, codon_autocorr_lag1, GC_pct, GC_max_60bp, GC3_pct,
CpG_density_per100bp, TIS_accessibility_5p, MFE_min_window, stable_structure_count, hairpin_count_all,
hairpin_count_5prime, splice_sites_total, homopolymer_runs_ge6, restriction_sites_total, repeat_count_10bp,
repeat_count_20bp, tm_spread` (Kozak excluded, as for the other tools).

`Variant` is `T<temperature>_p<top_p>_d<draw>`, e.g. `T0.2_p0.8_d1`; `Name` + `Variant` is unique.

---

## 7. Things to know when comparing

- **Rows per gene:** 80 (16 temperature × top_p cells × 5 draws), not one per cell. Draws can repeat, especially
  at low temperature: in the delivered grid there are 54–63 distinct DNA per gene (488 distinct overall), each
  scored once. Identical sequences carry identical scores; decide in the comparison step whether duplicates are
  counted or collapsed.
- **Status OK** means "translates to the input protein" (a validity check), not a constraint pass, as for Codon
  Harmony. CodonTransformer has no constraints for restriction/splice sites, homopolymers, GC windows etc.; these
  are only measured, in `restriction_sites_total`, `splice_sites_total`, `homopolymer_runs_ge6`, `GC_*`.
- **Grid-sheet columns are not for cross-tool use:** `GC% global`, `Length (nt)`, `CAI`, `Tm` in the grid workbook
  were computed on the raw sequence including the stop codon and with the grid script's own formulas. Compare
  tools only on the scoring table.
- **Species:** generation used the pretrained `Homo sapiens` organism token; `human_codon.csv` only affects
  scoring (CAI etc.), never generation.

---

## 8. Checks on the first real run

Expected console output:

```
Terminal stop codon found: CodonTransformer 640/640 | GeneArt 0/8
stop-policy strip: removed 640 (CodonTransformer) + 0 (GeneArt) terminal stop codon(s) ...
Translation check: all 648 sequences encode their input protein.
Span check: every sequence is exactly 3 x protein length ...
Status OK rows: 640 | GeneArt: 8 | distinct DNA to score: 488
```

If `seqfold` or the reference CSVs are missing, the three structure metrics and/or tAI will be empty: fix before using
the table.

---

## 9. Verification status

Verified on the delivered grid workbook, with the real `codon_features.py` / `score_all_sequences.py`: sheet and
column layout, uniqueness of `Name`+`Variant`, gene match, translation, span, stop handling, and that
`compute_features` accepts the stripped sequences. **Not** verified by the author of this README: a full run with
`human_codon.csv`, `human_tgcn.csv` and `seqfold`, and the 40-worker pool on your machine.
