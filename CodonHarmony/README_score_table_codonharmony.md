# Score table: Codon Harmony and GeneArt

**Script:** `Sequence_scoring/score_codonharmony.py`
**Input:** the `Results` sheet of the Codon Harmony grid output (`human_codonharmony_grid.xlsx`) and the original trials workbook (GeneArt sequences and input proteins)
**Output:** `Sequence_scoring/outputs/scores/scores_codonharmony.csv` and `scores_codonharmony.xlsx`

Purpose: benchmark Codon Harmony against GeneArt. DNA Chisel is **not** part of this table; the three sources are not mixed.

## 1. The scoring is unchanged

`codon_features.py`, `codon_compare.py` and `score_all_sequences.py` are imported and **not modified**. Every value is read from `codon_features.compute_features()['scores']` through `score_all_sequences.score_one`, with the same `human_codon.csv` and `human_tgcn.csv`. The 22 metrics, their definitions and their column order are those of `README_score_table.md`; read that file for what each metric means. The new script only does the following:

- reads the Codon Harmony workbook
- keeps the `Status = OK` rows
- adds the GeneArt rows
- checks the inputs
- hands each distinct DNA to the existing scorer once
- writes the table

## 2. How to run

```bash
cd Sequence_scoring
pip install numpy pandas openpyxl seqfold tqdm
python score_codonharmony.py --grid ../CodonHarmony/outputs/human_codonharmony_grid.xlsx
python score_codonharmony.py --workers 40            # default: min(40, cores-1), maximum 61
python score_codonharmony.py --stop-policy asis      # see section 4
```

If more than one trials workbook sits in the project folder, pass the original one with `--trials`. The trials workbook is **required** (it holds the GeneArt DNA and the input proteins); the script stops with a message if it cannot find it. Install `seqfold`, or `TIS_accessibility_5p`, `MFE_min_window` and `stable_structure_count` will be empty (the script prints a warning).

## 3. What is in the table

| | |
|---|---|
| Rows | every Codon Harmony grid row with `Status = OK` (350 in the supplied workbook) plus the GeneArt sequence of each gene (8) |
| Columns | `Name`, `Source` (`CodonHarmony` or `GeneArt`), `Variant` (the grid cell, e.g. `C017_ic5_c10_rel0.2_hp5`, or `GeneArt`), then the 22 metrics |
| Excel sheets | `CodonHarmony` and `GeneArt`; the CSV has both together |

Rows left out: `no solution within max_relax` (55) and `TRANSLATION MISMATCH` (27), because they have no usable sequence or do not encode the input protein.

`Status = OK` means only that the sequence translates back to the input protein. It does **not** mean that any liability constraint was met. Judge liabilities from the scored columns (`restriction_sites_total`, `splice_sites_total`, `homopolymer_runs_ge6`, and so on).

OK rows per gene in the supplied workbook (out of 54 cells each): ALDH1A3 48, AURKA 32, AURKB 48, CTNNB1 46, ERG 48, FOXA1 48, IRF4 32, PREP 48. Codon Harmony gives one sequence per cell, GeneArt gives one sequence per gene, and the numbers of Codon Harmony rows differ by gene, so summarise per gene before comparing.

## 4. Terminal stop codon (input handling only, not scoring)

`compute_features` strips a trailing stop codon for the codon-level metrics (CAI, rare codons, ENC, tAI and so on). It computes the nucleotide-level metrics (GC, CpG, structure, hairpins, repeats, restriction and splice counts, Tm fragments) on the DNA as given. Codon Harmony writes no stop codon; a GeneArt sequence may end in one. So that every sequence is scored over the same nucleotides, `--stop-policy` decides what is handed to the scorer:

| Policy | Behaviour |
|---|---|
| `strip` (default) | removes one terminal stop codon, if present, from every sequence of both sources. It does nothing when no sequence has one. Afterwards every sequence must be exactly 3 × its protein length, which is checked; the run stops if not |
| `asis` | passes the DNA unchanged, as `score_all_sequences.py` does. A warning is printed if the two sources differ in stop-codon handling |

The number of terminal stop codons found and removed per source is always printed. The default differs from `score_all_sequences.py`, which scores GeneArt as given, so confirm the policy with the project lead.

## 5. Safeguards in the script

- every sequence must translate to its input protein (the protein comes from the trials workbook)
- under `strip`, every sequence must be exactly 3 × protein length
- every gene in the grid must be in the trials workbook
- each distinct DNA is scored once and the result copied to every row with that DNA; identical DNA always gives identical scores

## 6. Comparing the two sources

- Use only the columns of this table. They come from one scorer for both sources.
- Do **not** compare these columns with the extra audit columns in the Codon Harmony workbook (`Splice sites remaining`, `Restriction sites remaining`, `GC% min/max (50nt win)`, `CAI` and so on). They are computed by the Codon Harmony script with different counting rules. In the supplied workbook, both count columns were reproduced exactly from the DNA, and the differences are fully explained by two rules:
  - **Splice sites:** the scorer counts the forward strand only; Codon Harmony counts both strands.
  - **Restriction sites:** the scorer counts a palindromic site once; Codon Harmony counts it twice. 29 of the 34 enzymes are palindromic (the exceptions are BsaI, BsmBI, BbsI, SapI and AarI).
  - Both use the same 34 enzymes and the same 4 splice motifs. Neither definition is wrong.
- The scorer's `CAI` is computed by the scorer for both sources; it differs from the `CAI` column in the Codon Harmony workbook, which uses the Codon Harmony script's own formula.
- The scorer's hairpin count is a simple stem heuristic; compare it between sequences, not against zero.

## 7. Observation in the supplied Codon Harmony output (to check)

In `human_codonharmony_grid.xlsx`, the 350 OK sequences still contain restriction sites. Only 3 of 350 are free of all 34 enzymes, and only 13 of 350 are free of the original 13 (BsaI, BsmBI, BbsI, SapI, EcoRI, XhoI, NotI, XbaI, SpeI, KpnI, SalI, AscI, PacI). The grid script passes the 34 enzymes to Codon Harmony, and the cause has not been established. This is a property of the Codon Harmony output, not of the scoring, and `restriction_sites_total` reports it correctly. Check it before drawing conclusions about restriction-site avoidance.

## 8. Verification status

Checked on the supplied files: the workbook layout, the 350 OK rows (all distinct, all ACGT, no empty DNA, none ending in a stop codon), exact agreement of `GC_pct` with the workbook's `GC% global`, and every branch of the stop-codon handling (no stops, stop under `strip`, stop under `asis`, wrong extra codon rejected). The scoring files were compared byte for byte with the uploaded copies.

Not yet checked: a run with the real `human_codon.csv`, `human_tgcn.csv`, `seqfold` and the real trials workbook. Check that the first run reports 350 and 8 rows, a clean stop-codon line, and non-empty structure columns.
