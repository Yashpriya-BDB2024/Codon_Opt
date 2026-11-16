"""
codontransformer_grid.py

Grid / sweep driver for CodonTransformer diverse-sequence generation - the
generative counterpart to codon_grid_search.py (which drives DNA Chisel).

CodonTransformer is a sampler, not a constraint solver, so the "grid" here is a
small product of the two decoding knobs that control diversity - temperature and
top_p - and for each (temperature, top_p) combination the model draws
`num_sequences` sequences per protein. The other decoding settings are fixed, as
specified:
      deterministic = False     (sampling on; required for any diversity)
      match_protein = True      (sampling restricted to codons of the correct
                                 amino acid -> every output back-translates to
                                 the input protein, even at high temperature)
      num_sequences = 5         (draws per protein per combination; configurable)

For each generated sequence it records the same sequence metrics as the DNA
Chisel grid (so the two tools are directly comparable), plus per-combination
diversity metrics (how many of the N draws are unique, and their mean pairwise
identity).

Input
-----
Excel with a "Name" column and a sequence column, protein or DNA:
  --input-type protein | dna | auto (default: detect). DNA is translated to the
  protein first (CodonTransformer takes a protein and an organism).

Model / organism
----------------
--model / --tokenizer default to the public "adibvafa/CodonTransformer" weights;
point --model at a fine-tuned checkpoint (e.g. your CHO-fine-tuned model) to
generate for a host the base model doesn't cover. CHO (C. griseus) is NOT in the
164 pretrained organisms - use a fine-tuned checkpoint or a proxy organism
(Homo sapiens / Mus musculus). --organism is passed through unchanged.

Usage
-----
    python codontransformer_grid.py --input proteins.xlsx --output ct_grid.xlsx \
        --organism "Homo sapiens" \
        --temperatures 0.2,0.4,0.6,0.8 --top-ps 0.8,0.9,0.95,1.0

    # fine-tuned CHO checkpoint, with CAI scored against your CHO codon table
    python codontransformer_grid.py --input proteins.xlsx --output cho_ct.xlsx \
        --model /path/to/cho_finetuned --tokenizer /path/to/cho_finetuned \
        --organism "Cricetulus griseus" --codon-table cho_cocoputs.csv

Requires: torch, transformers, CodonTransformer, pandas, numpy, openpyxl,
biopython (Tm). Runs on GPU if available, else CPU (slow). The model is loaded
once and reused across the whole grid.
"""

import argparse
import itertools
import json
import math
import os
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Standard genetic code (translation + codon-table loading)
# ---------------------------------------------------------------------------
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R',
    'CGA': 'R', 'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}
STANDARD_AAS = set("ACDEFGHIKLMNPQRSTVWY")

NAME_COL = "Name"
DNA_COL_CANDIDATES = ["Generated DNA Sequence", "DNA Sequence", "DNA", "CDS", "CDS Sequence"]
PROT_COL_CANDIDATES = ["Protein Sequence", "Protein", "Amino Acid Sequence", "AA Sequence", "Peptide"]
SHARED_COL_CANDIDATES = ["Sequence"]

AUDIT_KMER = 10
GC_WINDOW_AUDIT = 50
TM_AUDIT_WINDOW = 30

EUKARYOTE_MOTIFS = ["AWTAAA", "ATTTA", "TATAAA"]
PROKARYOTE_MOTIFS = ["AGGAGG"]

# Fixed generation settings (per the agreed spec; exposed as flags but default here)
FIXED_DETERMINISTIC = False
FIXED_MATCH_PROTEIN = True
DEFAULT_NUM_SEQUENCES = 5


# ===========================================================================
# Sequence helpers
# ===========================================================================
def translate(dna):
    dna = dna.upper().replace("U", "T")
    return "".join(GENETIC_CODE.get(dna[i:i + 3], "X") for i in range(0, len(dna) - 2, 3))


def protein_from_cell(raw, input_type):
    """Normalise one input cell to a clean protein sequence (+ detected kind)."""
    seq = "".join(str(raw).split()).upper()
    kind = input_type
    if input_type == "auto":
        kind = "dna" if (set(seq) <= set("ACGTUN") and len(seq) % 3 == 0) else "protein"
    if kind == "dna":
        prot = translate(seq).rstrip("*")
        if "X" in prot or "*" in prot:
            raise ValueError("DNA does not translate cleanly to a protein (internal stop / bad codon)")
        return prot, "dna->protein"
    prot = seq.rstrip("*")
    bad = set(prot) - STANDARD_AAS
    if bad:
        raise ValueError(f"invalid amino acid(s): {sorted(bad)}")
    return prot, "protein"


def pick_seq_column(df, seq_col, input_type):
    if seq_col:
        if seq_col not in df.columns:
            raise ValueError(f"--seq-col '{seq_col}' not in columns: {list(df.columns)}")
        return seq_col
    candidates = []
    if input_type in ("auto", "protein"):
        candidates += PROT_COL_CANDIDATES
    if input_type in ("auto", "dna"):
        candidates += DNA_COL_CANDIDATES
    candidates += SHARED_COL_CANDIDATES
    for col in candidates:
        if col in df.columns:
            return col
    others = [c for c in df.columns if c != NAME_COL]
    if len(others) == 1:
        return others[0]
    raise ValueError(f"Could not auto-detect the sequence column; pass --seq-col. Columns: {list(df.columns)}")


# ===========================================================================
# Codon usage table loader (file only; same format as codon_grid_search.py) -
# used purely to score CAI on the generated sequences. Optional.
# ===========================================================================
def _normalise_table(raw):
    tbl = {}
    for aa, codons in raw.items():
        codons = {c.upper(): float(v) for c, v in codons.items() if v is not None}
        total = sum(codons.values())
        if total > 0:
            tbl[aa] = {c: v / total for c, v in codons.items()}
    return tbl


def _first_present(cols, candidates):
    for c in candidates:
        if c in cols:
            return c
    return None


def load_codon_table(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Codon table file not found: {path}")
    low = path.lower()
    if low.endswith((".csv", ".tsv")):
        df = pd.read_csv(path, sep="\t" if low.endswith(".tsv") else ",")
        df.columns = [c.strip().lower() for c in df.columns]
        aa_col = _first_present(df.columns, ["amino_acid", "aa", "aminoacid", "residue"])
        codon_col = _first_present(df.columns, ["codon", "triplet"])
        val_col = _first_present(df.columns, ["fraction", "frequency", "freq", "count", "number", "proportion"])
        if codon_col is None or val_col is None:
            raise ValueError(f"Codon-table CSV needs a 'codon' and a value column. Got: {list(df.columns)}")
        raw = defaultdict(dict)
        for _, r in df.iterrows():
            codon = str(r[codon_col]).strip().upper().replace("U", "T")
            if codon not in GENETIC_CODE:
                continue
            aa = str(r[aa_col]).strip().upper() if aa_col else GENETIC_CODE[codon]
            if len(aa) != 1:
                aa = GENETIC_CODE[codon]
            try:
                raw[aa][codon] = float(r[val_col])
            except (ValueError, TypeError):
                continue
        return _normalise_table(raw)
    if low.endswith(".json"):
        data = json.loads(Path(path).read_text())
        if all(isinstance(v, dict) for v in data.values()):
            return _normalise_table(data)
        raw = defaultdict(dict)
        for codon, val in data.items():
            c = str(codon).strip().upper().replace("U", "T")
            if c in GENETIC_CODE:
                raw[GENETIC_CODE[c]][c] = float(val)
        return _normalise_table(raw)
    raise ValueError(f"Unrecognised codon-table file type (use .csv/.tsv/.json): {path}")


def cai(dna, table):
    """Codon Adaptation Index (Sharp & Li 1987): geometric mean of w=f/f_max
    within each amino-acid family, over coding codons (M, W and stops excluded)."""
    if table is None:
        return None
    wmax = {aa: max(codons.values()) for aa, codons in table.items()}
    logs = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3]
        aa = GENETIC_CODE.get(codon)
        if aa in (None, "*"):
            continue
        fam = table.get(aa)
        if not fam or len(fam) < 2:    # skip single-codon families (M, W) and missing
            continue
        f = fam.get(codon)
        if f and wmax[aa] > 0:
            logs.append(math.log(f / wmax[aa]))
    return round(math.exp(sum(logs) / len(logs)), 4) if logs else None


# ===========================================================================
# Sequence metrics (self-contained; mirror codon_grid_search.py)
# ===========================================================================
def gc_fraction(seq):
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


def window_gc_extremes(seq, window=GC_WINDOW_AUDIT):
    if len(seq) <= window:
        g = gc_fraction(seq) * 100
        return g, g
    gcs = [gc_fraction(seq[i:i + window]) * 100 for i in range(len(seq) - window + 1)]
    return min(gcs), max(gcs)


def longest_homopolymer(seq):
    return max((len(m.group(0)) for m in re.finditer(r"(.)\1*", seq)), default=0)


def repeated_kmers(seq, k=AUDIT_KMER):
    if len(seq) < k:
        return 0
    counts = {}
    for i in range(len(seq) - k + 1):
        counts[seq[i:i + k]] = counts.get(seq[i:i + k], 0) + 1
    return sum(1 for c in counts.values() if c > 1)


def terminal_tm(seq, window=TM_AUDIT_WINDOW):
    try:
        from Bio.SeqUtils import MeltingTemp as mt
        w = min(window, len(seq) // 2) or len(seq)
        return round(mt.Tm_NN(seq[:w]), 1), round(mt.Tm_NN(seq[-w:]), 1)
    except Exception:
        return None, None


def motifs_remaining(seq, motifs):
    iupac = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]",
             "S": "[GC]", "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]",
             "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}
    rc = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    total = 0
    for motif in motifs:
        pat = "".join(iupac.get(b, b) for b in motif)
        total += len(re.findall(f"(?={pat})", seq)) + len(re.findall(f"(?={pat})", rc))
    return total


def seq_metrics(dna, protein, table, motifs):
    dna = dna.upper().replace("U", "T")
    back = translate(dna).rstrip("*")
    valid = (back == protein)
    wmin, wmax = window_gc_extremes(dna)
    tm5, tm3 = terminal_tm(dna)
    return {
        "Translates to protein": valid,
        "Length (nt)": len(dna),
        "GC% global": round(gc_fraction(dna) * 100, 2),
        f"GC% min ({GC_WINDOW_AUDIT}nt win)": round(wmin, 2),
        f"GC% max ({GC_WINDOW_AUDIT}nt win)": round(wmax, 2),
        "Longest homopolymer": longest_homopolymer(dna),
        f"Repeated {AUDIT_KMER}-mers": repeated_kmers(dna),
        "CAI": cai(dna, table),
        f"5' Tm ({TM_AUDIT_WINDOW}nt)": tm5,
        f"3' Tm ({TM_AUDIT_WINDOW}nt)": tm3,
        "Forbidden motifs left": motifs_remaining(dna, motifs),
    }


def mean_pairwise_identity(seqs):
    """Mean % identity over all pairs of equal-length draws (lower = more diverse)."""
    seqs = [s for s in seqs if s]
    if len(seqs) < 2:
        return None
    idents = []
    for a, b in itertools.combinations(seqs, 2):
        n = min(len(a), len(b))
        if n == 0:
            continue
        same = sum(x == y for x, y in zip(a[:n], b[:n]))
        idents.append(100 * same / n)
    return round(sum(idents) / len(idents), 2) if idents else None


# ===========================================================================
# Backend (isolated so the rest of the script is testable without torch)
# ===========================================================================
def build_predict_fn(model_path, tokenizer_path, organism, device_str,
                     attention_type, num_sequences, seed):
    """Load the model/tokenizer once and return predict_fn(protein, temperature,
    top_p) -> list[DNASequencePrediction]. All fixed settings are bound here."""
    import torch
    from transformers import AutoTokenizer, BigBirdForMaskedLM
    from CodonTransformer.CodonPrediction import predict_dna_sequence

    device = torch.device(device_str if device_str
                          else ("cuda" if torch.cuda.is_available() else "cpu"))
    print(f"Loading tokenizer: {tokenizer_path}")
    tokenizer = AutoTokenizer.from_pretrained(tokenizer_path)
    print(f"Loading model: {model_path}  -> device {device}")
    model = BigBirdForMaskedLM.from_pretrained(model_path).to(device)
    model.eval()

    def predict_fn(protein, temperature, top_p):
        if seed is not None:
            torch.manual_seed(seed)          # reproducible draws per call
        out = predict_dna_sequence(
            protein=protein,
            organism=organism,
            device=device,
            tokenizer=tokenizer,
            model=model,
            attention_type=attention_type,
            deterministic=FIXED_DETERMINISTIC,
            match_protein=FIXED_MATCH_PROTEIN,
            num_sequences=num_sequences,
            temperature=temperature,
            top_p=top_p,
        )
        return out if isinstance(out, list) else [out]

    return predict_fn, str(device)


# ===========================================================================
# Grid run (pure; takes predict_fn so it can be tested with a stub)
# ===========================================================================
def run_grid(records, predict_fn, combos, table, motifs, organism):
    """records: list of (name, protein, kind). Returns (per-draw rows, per-combo rows)."""
    seq_rows, combo_rows = [], []
    for name, protein, kind in records:
        for temperature, top_p in combos:
            try:
                preds = predict_fn(protein, temperature, top_p)
                dnas = [p.predicted_dna.upper().replace("U", "T") for p in preds]
                err = None
            except Exception as e:
                dnas, err = [], f"GENERATION FAILED: {str(e)[:120]}"

            draw_metrics = []
            for i, dna in enumerate(dnas, 1):
                m = seq_metrics(dna, protein, table, motifs)
                draw_metrics.append(m)
                seq_rows.append({
                    NAME_COL: name, "Organism": organism, "Input kind": kind,
                    "Temperature": temperature, "top_p": top_p, "Draw": i,
                    "Status": "OK" if m["Translates to protein"] else "TRANSLATION MISMATCH",
                    **m, "Generated DNA Sequence": dna,
                })
            if err:
                seq_rows.append({NAME_COL: name, "Organism": organism, "Input kind": kind,
                                 "Temperature": temperature, "top_p": top_p, "Draw": 0,
                                 "Status": err})

            # per-combination diversity + mean metrics
            n = len(dnas)
            n_unique = len(set(dnas))
            n_valid = sum(1 for m in draw_metrics if m["Translates to protein"])
            metric_keys = [k for k in (draw_metrics[0] if draw_metrics else {})
                           if isinstance(draw_metrics[0][k], (int, float))]
            means = {f"mean {k}": round(float(np.mean([m[k] for m in draw_metrics])), 2)
                     for k in metric_keys} if draw_metrics else {}
            combo_rows.append({
                NAME_COL: name, "Organism": organism,
                "Temperature": temperature, "top_p": top_p,
                "n_draws": n, "n_unique": n_unique, "n_valid": n_valid,
                "Unique fraction": round(n_unique / n, 2) if n else None,
                "Mean pairwise identity %": mean_pairwise_identity(dnas),
                **means,
                "Status": err or "OK",
            })
            print(f"{str(name):<12} T={temperature:<4} top_p={top_p:<5} "
                  f"-> {n_unique}/{n} unique" + (f"  [{err}]" if err else ""))
    return seq_rows, combo_rows


def parse_floats(text):
    return [float(x) for x in str(text).replace(" ", "").split(",") if x != ""]


def main():
    ap = argparse.ArgumentParser(description="CodonTransformer temperature x top_p grid generator.")
    ap.add_argument("--input", required=True, help="Excel with 'Name' + a protein/DNA column")
    ap.add_argument("--output", required=True, help="Output .xlsx")
    ap.add_argument("--model", default="adibvafa/CodonTransformer",
                    help="HF name or path to a (fine-tuned) checkpoint")
    ap.add_argument("--tokenizer", default="adibvafa/CodonTransformer")
    ap.add_argument("--organism", default="Homo sapiens",
                    help="Organism name/ID (CHO not in pretrained vocab; use a fine-tuned model or a proxy)")
    ap.add_argument("--temperatures", default="0.2,0.4,0.6,0.8", help="Comma-separated")
    ap.add_argument("--top-ps", default="0.8,0.9,0.95,1.0", help="Comma-separated")
    ap.add_argument("--num-sequences", type=int, default=DEFAULT_NUM_SEQUENCES,
                    help="Draws per protein per (temperature, top_p) combination")
    ap.add_argument("--attention-type", default="original_full",
                    choices=["original_full", "block_sparse"])
    ap.add_argument("--domain", choices=["eukaryote", "prokaryote"], default="eukaryote")
    ap.add_argument("--codon-table", default=None, help="Optional codon table file for CAI scoring")
    ap.add_argument("--input-type", choices=["auto", "protein", "dna"], default="auto")
    ap.add_argument("--seq-col", default=None)
    ap.add_argument("--device", default=None, help="cuda / cpu (default: auto)")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    combos = [(t, p) for t in parse_floats(args.temperatures) for p in parse_floats(args.top_ps)]
    motifs = PROKARYOTE_MOTIFS if args.domain == "prokaryote" else EUKARYOTE_MOTIFS
    table = load_codon_table(args.codon_table) if args.codon_table else None

    data = pd.read_excel(args.input)
    if NAME_COL not in data.columns:
        raise ValueError(f"Missing column: {NAME_COL}")
    seq_col = pick_seq_column(data, args.seq_col, args.input_type)
    if args.limit:
        data = data.head(args.limit)

    records = []
    for _, r in data.iterrows():
        try:
            prot, kind = protein_from_cell(r[seq_col], args.input_type)
            records.append((r[NAME_COL], prot, kind))
        except Exception as e:
            print(f"{str(r[NAME_COL]):<12} INPUT ERROR: {e}")

    print(f"Proteins={len(records)}  combos={len(combos)} (temp x top_p)  "
          f"num_sequences={args.num_sequences}  organism={args.organism}  domain={args.domain}")

    predict_fn, device_used = build_predict_fn(
        args.model, args.tokenizer, args.organism, args.device,
        args.attention_type, args.num_sequences, args.seed)
    print(f"Device: {device_used}")

    seq_rows, combo_rows = run_grid(records, predict_fn, combos, table, motifs, args.organism)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(out) as writer:
        pd.DataFrame(combo_rows).to_excel(writer, sheet_name="Summary", index=False)
        pd.DataFrame(seq_rows).to_excel(writer, sheet_name="Sequences", index=False)
        pd.DataFrame([{"Temperature": t, "top_p": p} for t, p in combos]).to_excel(
            writer, sheet_name="Grid", index=False)
    print(f"\nSaved {len(seq_rows)} sequences across {len(combos)} combos to {out}")


if __name__ == "__main__":
    main()
