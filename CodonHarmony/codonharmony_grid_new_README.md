# `codonharmony_grid_new.py` + `run_human_parallel_codonharmony.py` — README & change log

Derived from the project lead's `codonharmony_grid.py` and its documentation. The edits exist for one reason:
**so Codon Harmony's results can be scored side by side with the DNA Chisel results** (`codon_grid_search_new.py`,
`run_human_parallel_new.py`). Codon Harmony's optimisation logic, grid defaults, host threshold, seeding, de-duplication
and translation check are **unchanged**. Every edited line is tagged `# [NEW]` / `# [CHANGED]` in the code.

---

## 1. Changes made (and why each is justified)

| # | Change | Justification |
|---|---|---|
| 1 | **Restriction enzymes: 13 → 34** (the exact `DEFAULT_SITES` list of the DNA Chisel script). `load_backend` also checks every name exists in `Bio.Restriction`. | Site avoidance is a *requirement* on the final sequence, not a tuning knob. If the two tools avoid different enzyme sets, a sequence from one tool can contain a site the other was forbidden to produce, and "restriction-site-free" would not be a fair comparison. The lead's own documentation says the list is meant to be "the same as the DNA Chisel grid"; DNA Chisel's list is now 34, so this keeps that intent. |
| 2 | **New measured columns:** `Splice sites remaining` (DNA Chisel's 4 splice motifs), `Restriction sites remaining`, `Nucleotides changed`, `% identity`, `Codon-match score`. | DNA Chisel's workbook reports these, so comparison needs them. They are *read off the finished sequence*; they do not alter optimisation. Codon Harmony removes splice sites with its own hard-coded patterns, so auditing with DNA Chisel's list gives one common yardstick. |
| 3 | **CAI kept exactly as the lead wrote it** (no extra CAI column). | Each tool's CAI is reported the way its own script computes it; the two are compared as reported (see §3). |
| 4 | **Output layout aligned with DNA Chisel:** `Variant`, `Host`, `Domain`, `Runtime (s)` columns; sheets `Summary` (per variant), `Results`, `Grid`, plus `Distinct per protein` (the lead's original summary). | One scoring script can then read both workbooks. `Variant` looks like `C017_ic5_c10_rel0.2_hp5`. Failed cells get an empty DNA, as in DNA Chisel. `Distinct per protein` now lists proteins with 0 valid sequences too (previously omitted). |
| 5 | **Bug fix:** scratch folder `"/tmp"` → `tempfile.gettempdir()`. | `/tmp` does not exist on Windows, so the profile-JSON write would fail there. |
| 6 | **New optional flag `--host-threshold`** (default = lead's 0.10). | Lets the lead test 0.11 without editing code. No default behaviour change. |
| 7 | Row/summary building moved into functions (`run_cell_row`, `build_summary`, `write_workbook`, …). | So the parallel runner and the sequential script produce identical output. Verified: same sequences, seeds and metrics. |

---

## 2. Deliberately NOT changed — and why

### Homopolymer axis stays `3, 4, 5`
- Tools are compared on the **measured** `Longest homopolymer` of the finished sequences, not on the setting that was passed in. Matching the setting values is therefore not needed for a valid comparison.
- `3,4,5` is the lead's agreed design; changing it breaks reproducibility with anything the lead has already run.
- The two sweeps already overlap (DNA Chisel: 4, 5, 6; Codon Harmony: 3, 4, 5), and 4 and 5 are common to both.
- If a matched range is wanted, no code edit is needed: `--homopolymer 4,5,6`.

### Host threshold stays `0.10`
- It is Codon Harmony's own default and part of the lead's "agreed design".
- The DNA Chisel value (`rare_codon_min_freq` 0.11) differs by 0.01, which is not a meaningful difference, and the two parameters belong to different algorithms (a profile filter vs. a hard constraint), so equal numbers would not make them equivalent anyway.
- Available as `--host-threshold 0.11` if ever wanted.

### Other things that cannot or should not be matched
- Grid axes (54 cells vs 2,592 variants): the tools have different parameters.
- Codon Harmony's hard-coded GC windows (not 0.44–0.68 / 50 nt), 9-nt repeat window, hairpin stem, and splice patterns. Compare through the measured `GC% min/max`, `Repeated 10-mers`, `Splice sites remaining` columns.
- `AWTAAA / ATTTA / TATAAA` cannot be avoided by Codon Harmony; they are only measured (`Forbidden motifs left`), as in the lead's script.
- Terminal Tm / GC constraints have no Codon Harmony equivalent (Tm still measured).
- `Status = OK` means "translates back to the input protein" here, but "all constraints passed" in DNA Chisel. For Codon Harmony judge compliance with `Restriction sites remaining`, `Splice sites remaining` and `Longest homopolymer`.

---

## 3. CAI

`CAI` is the lead's calculation, unchanged: `exp(mean log w)` over codons of multi-codon families (stop and single-codon amino acids M, W are skipped). The DNA Chisel script's `cai_value()` divides by *all* codons instead, so its values run slightly higher for the same sequence. The comparison uses each tool's CAI as reported; the formulas differ, so treat small absolute differences between tools with that in mind. If an exact like-for-like CAI is ever needed, rescore all sequences from both tools with one function.

---

## 3b. What one run produces ("Variant" = grid cell)

One Codon Harmony run returns **one** optimised sequence. The grid then runs every cell (inner_cycles × cycles × max_relax × homopolymer = 54 by default) for each gene, so each gene gets **up to 54 candidate sequences**, one per cell. `Variant` (e.g. `C017_ic5_c10_rel0.2_hp5`) is just the label of that cell, the equivalent of DNA Chisel's variant name. Different cells can converge to the same sequence; `Distinct per protein` reports how many unique valid sequences each gene actually ended up with.

---

## 4. Running it

**Sequential** (lead's way, any input workbook with `Name` + protein column):
```bash
python codonharmony_grid_new.py --input human_proteins.xlsx --output outputs/human_ch_grid.xlsx \
    --species human --table-dir codon_tables
```

**Parallel, for the 8 human sequences** (mirrors `run_human_parallel_new.py`; reads the same trials workbook and the same `host_files/Human/human_codon.csv`):
```bash
python run_human_parallel_codonharmony.py --test        # 2 seqs x 2 cells - smoke test first
python run_human_parallel_codonharmony.py               # full run, 8 x 54 = 432 runs
python run_human_parallel_codonharmony.py --resume      # continue after an interruption
```
- Each (sequence, cell) run is an independent task; seeds are `42 + cell index`, so results are **identical to a sequential run** whatever the worker count.
- Default 40 workers, capped at 61 (Windows limit). `--workers N` to change.
- Progress bar via `tqdm`. Finished runs are checkpointed to `<output>.partial.jsonl`; `--resume` reuses them, and the file is deleted after a successful run.
- All grid overrides work: `--inner-cycles`, `--cycles`, `--max-relax`, `--homopolymer`, `--host-threshold`, `--seed`.
- Keep `biopython_compat.py` in the same folder. `pip install dnachisel` in the same environment to fill `Nucleotides changed`, `% identity`, `Codon-match score` (blank otherwise).
- **Runtime:** with ~30 min per sequence's grid (as your lead says), 8 sequences in parallel take about as long as the slowest sequence; splitting by cell lets all cores be used.

---

## 5. Verify on your machine first

The sandbox used for writing this had no `codon_harmony`, Biopython or DNA Chisel. The grid, metrics, workbook writing, parallel runner, and resume were tested against a stubbed Codon Harmony only. Please:
1. Run `--test` of the parallel runner (real Codon Harmony): confirms the 34-enzyme list is accepted and workers start.
2. Confirm `Restriction sites remaining` is 0 on `OK` rows.
