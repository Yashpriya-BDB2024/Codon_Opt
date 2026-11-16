# `codontransformer_grid.py` — how it works

*Companion documentation for the CodonTransformer grid driver, written to parallel the DNA Chisel grid doc. Explains how CodonTransformer generates sequences, what the grid sweeps and why it's small, how the model (not a codon table) carries species, and what the output reports. Verified against CodonTransformer 1.6.7 source.*

---

## 1. What it does, and how it differs from the DNA Chisel grid

CodonTransformer is a **generative sampler**, not a constraint solver. Where DNA Chisel searches a combinatorial space of hard constraints, CodonTransformer is a BigBird masked-language model that, given a protein and an organism, emits a probability distribution over codon tokens and *samples* DNA sequences from it. Diversity therefore comes from stochastic decoding — temperature and nucleus (top_p) sampling — not from a grid of constraints. That makes the tunable surface small: the whole grid is `temperature × top_p`, and for each combination the model draws `num_sequences` sequences per protein.

For each generated sequence the script records the same sequence metrics as the DNA Chisel grid (so the tools are comparable) plus **per-combination diversity** — how many of the draws are unique and their mean pairwise identity.

---

## 2. The decoding settings

All live on `predict_dna_sequence`. Three are **fixed** by design; two form the **grid**.

**Fixed (the agreed spec):**

| Setting | Value | Why |
|---|---|---|
| `deterministic` | `False` | `True` is greedy decoding — one most-likely sequence, zero diversity. Sampling must be on for the grid to mean anything. |
| `match_protein` | `True` | Restricts sampling to codons of the correct amino acid, so every output back-translates to the input protein even at high temperature (see §3). |
| `num_sequences` | `5` | Draws per protein per (temperature, top_p) cell. Configurable with `--num-sequences`. |

**Grid axes (swept, Cartesian product):**

| Axis | Default values | Effect |
|---|---|---|
| `temperature` | `0.2, 0.4, 0.6, 0.8` | The primary diversity dial. The model's own guidance: 0.2 conservative, 0.5 medium, 0.8 high. Higher = more diverse, and (without `match_protein`) rising risk of invalid translation. |
| `top_p` | `0.8, 0.9, 0.95, 1.0` | Nucleus-sampling cutoff. Lower restricts to the highest-probability codons (conservative); closer to 1.0 admits the tail (more diverse). Trades off against temperature. |

Defaults give 16 cells × 5 draws = 80 sequences per protein. Set the two lists with `--temperatures` and `--top-ps`.

`attention_type` (`original_full`, default, vs `block_sparse`) and `device` are not diversity knobs — `original_full` is correct for antibody-length CDS; device is auto (`cuda` if available, else CPU).

---

## 3. How generation works (and why there is no codon table in it)

The path, confirmed from source:

1. The protein and the `organism` (converted to an `organism_id`) are tokenised into one input.
2. One forward pass produces `logits` over codon tokens for each position. **All `num_sequences` draws are sampled from this single forward pass**, so more draws are cheap.
3. If `match_protein=True`, a mask built from the genetic code (`AMINO_ACID_TO_INDEX` — which codons encode each residue) sets the logits of every non-matching codon to `-inf`. This is a hard validity mask, not a weighting.
4. Sampling applies temperature scaling then nucleus (top_p) truncation to the (masked) logits and draws a codon per position.

The crucial consequence: **CodonTransformer never consumes a codon usage table at generation time.** Organism-specific codon preference was learned during training and lives in the model weights, selected by the `organism` token (or by a fine-tuned checkpoint). The only lookup table in the generation path is the plain genetic code, used for the `match_protein` mask. So the `--codon-table` flag in this script is **post-hoc scoring only** — it computes CAI on the already-generated sequences and has zero influence on what the model produces. It exists so CodonTransformer output can be measured on the same CAI yardstick as the DNA Chisel runs.

This is the opposite of DNA Chisel, where the codon table *is* the optimization target.

---

## 4. Model, organism, and CHO

`--model` / `--tokenizer` default to the public `adibvafa/CodonTransformer` weights, loaded once (the exact `AutoTokenizer` + `BigBirdForMaskedLM` pattern) and reused across the whole grid. `--organism` is passed through unchanged (name or ID).

**CHO is not in the 164 pretrained organisms** — verified against the full `ORGANISM2ID`; the only mammalian proxies are `Homo sapiens` and `Mus musculus` (no *Cricetulus griseus*, no hamster). Since CodonTransformer carries species in its weights, the *only* way to get CHO codon usage into generation is to point `--model`/`--tokenizer` at a **fine-tuned CHO checkpoint** (your Stage-1 model). A CoCoPUTs CHO table cannot steer generation — it can only score the output afterward via `--codon-table`. Without a checkpoint, your options are a proxy organism (human/mouse) for generation.

---

## 5. Input

Excel with a `Name` column plus a protein (or DNA) column. `--input-type auto` (default) detects per sequence from the alphabet; DNA is translated to protein first, since CodonTransformer takes a protein. Column auto-detected from common names or set with `--seq-col`.

---

## 6. Output

Three sheets.

- **`Sequences`** — one row per draw: `Name`, `Organism`, `Input kind`, `Temperature`, `top_p`, `Draw #`, `Status` (`OK`, or `TRANSLATION MISMATCH` if a draw doesn't back-translate — a QC guard, rare under `match_protein`), then the shared metrics — translates-to-protein, length, global GC, worst-case GC over 50-nt windows (min/max), longest homopolymer, repeated 10-mers, CAI (if `--codon-table` given), 5′/3′ terminal Tm, forbidden motifs (`--domain` eukaryote/prokaryote) — and the DNA.
- **`Summary`** — one row per protein × (temperature, top_p): `n_draws`, `n_unique`, `n_valid`, `Unique fraction`, **`Mean pairwise identity %`** (lower = more diverse), and the mean of each metric. This is the diversity readout — at low temperature you'll typically see `n_unique` drop below 5 (the model repeating itself), which tells you where the useful part of the range begins.
- **`Grid`** — the (temperature, top_p) cells.

A generation failure for a cell (e.g. CUDA OOM) is caught and recorded in `Status`, not allowed to crash the run. `--seed` seeds Torch before each call for reproducible draws.

---

## 7. Command reference

```bash
# default grid, public model, human
python codontransformer_grid.py --input proteins.xlsx --output ct_grid.xlsx \
    --organism "Homo sapiens"

# fine-tuned CHO checkpoint, CAI scored against your CHO table
python codontransformer_grid.py --input proteins.xlsx --output cho_ct.xlsx \
    --model /path/to/cho_finetuned --tokenizer /path/to/cho_finetuned \
    --organism "Cricetulus griseus" --codon-table cho_codon.csv

# custom decoding grid
python codontransformer_grid.py --input proteins.xlsx --output ct.xlsx \
    --organism "Homo sapiens" --temperatures 0.4,0.6,0.8 --top-ps 0.9,0.95,1.0
```

`--num-sequences` (draws per cell, default 5), `--input-type`/`--seq-col`, `--domain` (motif screen), `--limit`, `--seed`, `--device` round out the flags.

---

## 8. Environment & licensing

Needs `torch`, `transformers`, `CodonTransformer`, plus `pandas`/`numpy`/`openpyxl`/`biopython` (Tm). Runs on GPU if present, else CPU (slow). The base weights auto-download from HuggingFace on first use; a local `--model` path skips that. CodonTransformer is MIT-licensed; the public weights carry their own HuggingFace model license — clear the checkpoint you actually ship. This grid shares the Biopython/pandas stack with the DNA Chisel grid, so the two can live in one environment; only Codon Harmony needs isolation.
