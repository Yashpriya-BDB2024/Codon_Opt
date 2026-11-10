# `codonharmony_grid.py` — how it works

*Companion documentation for the Codon Harmony grid driver, written to parallel the DNA Chisel grid doc. Explains how Codon Harmony optimizes, why the grid (not repeated draws) supplies the variants, how species tables are wired, the environment constraint, and the output. Verified against codon-harmony 1.0.0 source, run end-to-end on Biopython 1.77.*

---

## 1. What it does, and how it differs from the other two grids

Codon Harmony is a harmonization **optimizer**. Each run reverse-translates the protein, then iteratively mutates codons to (a) match the host codon-usage *profile* and (b) remove sequence liabilities, and returns the **single** highest-CAI sequence that stayed within an allowed deviation from the host profile. It is neither a constraint enumerator (DNA Chisel) nor a multi-draw sampler (CodonTransformer): **one run → one sequence.**

So diversity has to come from the **grid itself** — every `(inner_cycles, cycles, max_relax, local_homopolymer_threshold)` cell is one run and one candidate. No `R`/repeats: the grid is the variant set. Because Codon Harmony converges toward a best-CAI sequence, cells whose parameters don't materially change that objective can land on the **same** sequence, so the script de-duplicates across cells and reports the number of **distinct** sequences per protein.

---

## 2. How a single run works

Per run (confirmed from source):

1. **Reverse-translate** the protein to a starting DNA.
2. **Outer loop — `cycles` independent trials.** Each trial resamples codons and uses a progressively relaxed profile-match target: `relax = 1 + max_relax · (trial / cycles)`, so later trials tolerate more deviation.
3. **Inner loop — `inner_cycles` cleanup passes per trial**, in order: GC-content scan → remove alternate start sites (if enabled) → remove repeats (9-nt window) → remove local homopolymers (> `local_homopolymer_threshold`) → remove hairpins (stem 10) → remove restriction sites.
4. After cleanup, the trial's deviation from the host profile is measured; if it exceeds `max_relax` the trial is discarded.
5. Across all surviving trials, the sequence with the **best CAI** is kept. Splice sites are removed at the end if enabled.

If no trial stays within `max_relax`, the run returns no DNA — the script records `no solution within max_relax` for that cell (seen with tight relax on hard proteins; it's informative, not an error).

---

## 3. The grid

**Swept axes** (Cartesian product; defaults, override on the CLI):

| Axis | Default | Effect |
|---|---|---|
| `inner_cycles` | `5, 10` | Refinement passes per trial. More = more thorough liability cleanup. |
| `cycles` | `5, 10, 20` | Independent trials; keeps the best. More = better optimization and a wider relax schedule. |
| `max_relax` | `0.1, 0.2, 0.3` | Max fractional deviation from the host profile allowed. The main "how much codon diversity is tolerated" dial; too tight → `no solution`. |
| `local_homopolymer_threshold` | `3, 4, 5` | Longest consecutive-NT run permitted. |

Defaults = 2×3×3×3 = 54 cells per protein.

**Fixed top-level config** (edit at the top of the script):
- `DEFAULT_RESTRICTION_ENZYMES` — the same 13-enzyme list as the DNA Chisel grid (`BsaI, BsmBI, BbsI, SapI, EcoRI, XhoI, NotI, XbaI, SpeI, KpnI, SalI, AscI, PacI`); all names validated against Codon Harmony's Biopython `RestrictionBatch`. Extend the list to add more.
- `host_threshold = 0.10` (codon-harmony default — lowest codon fraction per amino acid kept).
- Splice-site removal **on**, alternate-start-site removal **off**, set automatically by domain (eukaryote vs prokaryote). Note: `--remove-splice-sites` in codon-harmony is a boolean toggle, not a list — the splice patterns are hardcoded; there is no custom splice list to pass.

**Not CLI-tunable (hardcoded inside codon-harmony's `data` module)**, so the grid cannot sweep them without a monkeypatch: GC-content windows (IDT 20 nt 0.15–0.90; twist 50 nt 0.15–0.80; IDT_long 100 nt 0.28–0.68; twist_long ×3 0.30–0.65), the repeat-detection window (9 nt), hairpin stem length (10), and the homopolymer codon span (2). Codon Harmony also has **no generic motif list** — it removes restriction sites, splice sites, start sites, repeats, homopolymers and hairpins, nothing like DNA Chisel's `AWTAAA/ATTTA/TATAAA` screen (that screen is still *measured* post-hoc for comparability).

---

## 4. Species & codon usage table

Pass `--species {human, cho, ecoli}`. A hardcoded registry (`SPECIES_TABLES` at the top of the script) maps each to a CSV in `--table-dir` (default `codon_tables/`):

```
human -> human_codon.csv
cho   -> cho_codon.csv
ecoli -> ecoli_codon.csv
```

These are **your CoCoPUTs CSVs** — the same `amino_acid,codon,fraction` tables built for the DNA Chisel / scoring work. There is **no `--host` and no download**; put the three files in `codon_tables/` and the selected one is used. Add a species by extending `SPECIES_TABLES` and `DOMAIN`.

Codon Harmony's own profile format differs from the DNA Chisel table — it wants nested `{aa: {RNA_codon: frequency}}`. The script **converts on the fly** (reads your CSV, normalises per amino-acid family, emits the nested RNA-keyed JSON to a temp file) and feeds that to `--local-host-profile`. Validated against Codon Harmony's own loader. `--domain` is set from the species (human/cho = eukaryote → splice-site removal on; ecoli = prokaryote → start-site removal on, Shine-Dalgarno motif audit).

---

## 5. Reproducibility — the random seed

Codon Harmony's entire stochastic core is Python's standard-library `random` (every codon resampling is a `random.choice/choices/randint`), and it never seeds itself. The script therefore calls `random.seed(--seed + cell_index)` before each run. This makes the **whole grid reproducible** — re-running the identical command reproduces byte-identical sequences (verified, and it holds across separate processes without needing `PYTHONHASHSEED`) — while giving each cell its own seed so cells remain distinct draws. `--seed` (default 42) sets the base.

---

## 6. A sequence-integrity caveat (why there's a `Translates to protein` check)

Codon Harmony's cleanup routines do **not** always preserve the exact protein — in testing, one cell emitted a sequence one codon too long whose C-terminus translated differently from the input. This is a real quirk of codon-harmony 1.0.0. The script guards against it: every output is translated and checked against the input protein, flagged `TRANSLATION MISMATCH` in `Status` if it fails, and **excluded from the distinct-sequence count**. Always filter on `Status == OK` before using sequences; the QC column exists precisely because this tool can occasionally alter the protein.

---

## 7. Output

Three sheets.

- **`Sequences`** — one row per cell: `Name`, `Species`, `Input kind`, `Cell`, the four axis values, `Seed`, `Status` (`OK` / `TRANSLATION MISMATCH` / `no solution within max_relax` / `FAILED: …`), then the shared metrics — translates-to-protein, length, global GC, worst-case 50-nt-window GC (min/max), longest homopolymer, repeated 10-mers, CAI (from the selected table), 5′/3′ terminal Tm, forbidden motifs (domain-specific) — and the DNA.
- **`Summary`** — one row per protein: `Grid cells` and **`Distinct sequences`** (valid, de-duplicated). This is the headline — how many usable variants the grid actually produced.
- **`Grid`** — the cells.

---

## 8. Environment (no downgrade needed) & licensing

Codon Harmony 1.0.0 imports several APIs removed in Biopython ≥ 1.78 — `Bio.Alphabet` (IUPAC alphabets), `Bio.SeqUtils.CodonUsage` (`SynonymousCodons`, `CodonsDict`, `CodonAdaptationIndex`), `Bio.SeqUtils.GC`, and `Seq.tomutable()` / `MutableSeq.toseq()`. Rather than pin an old Biopython or spin up a separate env, this driver ships with **`biopython_compat.py`**, a shim that reinstates exactly those pieces and is imported **first** (before codon-harmony). So Codon Harmony runs unchanged on the **same modern Biopython** as your other tools — no downgrade, no venv.

```bash
pip install codon-harmony pandas numpy openpyxl       # your normal environment
python codonharmony_grid.py --input proteins.xlsx --output cho_ch.xlsx \
    --species cho --table-dir codon_tables
```

Keep `biopython_compat.py` next to `codonharmony_grid.py`. The shim is **behaviour-equivalent**: with the same seed and parameters it was verified to produce **byte-for-byte identical** sequences to a real Biopython-1.77 run (the key detail is that the shim vendors Biopython's exact `SynonymousCodons` ordering, which codon-harmony's `random.choices` and dict iteration depend on — a different order would change the sampling stream). If you ever prefer the isolated route instead, a `biopython==1.77` venv also works and needs no shim.

Other flags: `--input-type`/`--seq-col`, `--limit`, `--seed`, and grid overrides `--inner-cycles` / `--cycles` / `--max-relax` / `--homopolymer` (comma-separated).

Codon Harmony is MIT-licensed (clear for commercial use); it depends on Biopython (BSD-style) and python_codon_tables (MIT). Since you supply the codon table file and never use `--host`, there is no runtime download and no third-party table-licensing question in the optimization path — the provenance of the CoCoPUTs CSVs you supply is whatever you've cleared.
