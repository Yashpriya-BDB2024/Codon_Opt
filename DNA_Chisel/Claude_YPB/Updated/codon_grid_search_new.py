"""
codon_grid_search_new.py

Grid / sweep driver for DNA Chisel codon-optimization variants.

Where the original dnachisel_variants.py ran 10 hand-written variants, this
script generates variants automatically by sweeping *ranges* of the relevant
DNA Chisel constraints/objectives (codon strategy, GC bounds, k-mer uniqueness,
rare-codon cutoff, homopolymer length, motif/site avoidance, hairpins, terminal
GC, and terminal melting temperature), runs each variant on every input
sequence, and records a full metrics table so you can see which cutoffs help.

Codon usage table (REQUIRED, --codon-table)
--------------------------------------------
The host codon usage table is supplied as a FILE. There is no species name,
no TaxID, no download, and no cached/bundled table - you pass the exact table
you have validated. The format matches the loader used by the codon-scoring
scripts (codon_score.py / codon_features.py):
  * CSV/TSV with columns  amino_acid, codon, fraction
    (the value column may instead be named frequency / freq / count / number /
    proportion; amino_acid may be omitted and is then derived from the codon).
  * JSON, either nested  {aa: {codon: value}}  or flat  {codon: value}.
Values are normalised per amino-acid family to sum to 1.0, which is exactly the
dict DNA Chisel expects, and the table is injected into every codon-aware spec
(CodonOptimize, AvoidRareCodons, MaximizeCAI, MatchTargetCodonUsage) via their
``codon_usage_table=`` argument. Build a CHO / C. griseus table from CoCoPUTs
(or a high-expression CHO transcriptome subset) once and pass it here.

Species handling
----------------
The grid and all numeric parameters are identical for every species - species
differences enter only through (a) the codon table you supply and (b) the
--domain flag, which selects the species-/domain-specific motif screen:
  * eukaryote (default): cryptic poly-A signal (AWTAAA), AU-rich element
    (ATTTA), TATA box (TATAAA).
  * prokaryote: internal Shine-Dalgarno / ribosome-binding core (AGGAGG).
Restriction sites, GC, k-mer, homopolymer, Tm and codon axes are shared.

Input sequences
---------------
Excel with a "Name" column and a sequence column, DNA or protein:
  --input-type dna | protein | auto (default: detect per sequence).
Protein is reverse-translated to a DNA seed (deterministic) then optimized.
Column auto-detected from common names, or set with --seq-col.

Two sweep modes
---------------
  * cartesian - full product of the grid, pruned + de-duplicated. Runs the
                whole grid by default; pass --max-variants N for a seeded
                random sample of N.
  * ofat      - one-factor-at-a-time: baseline, then vary one axis at a time.

Usage
-----
    python codon_grid_search_new.py --input genes.xlsx --output grid.xlsx \
        --codon-table cho_cocoputs.csv --domain eukaryote

    python codon_grid_search_new.py --input proteins.xlsx --output ofat.xlsx \
        --codon-table ecoli.csv --domain prokaryote \
        --mode ofat --limit 1 --input-type protein

    python codon_grid_search_new.py --dump-config default_grid.json   # edit the grid

Requires: dnachisel, pandas, numpy, openpyxl  (tested on dnachisel 3.2.16).
"""

import argparse
import itertools
import json
import math
import os
import re
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from dnachisel import (
    AvoidHairpins,
    AvoidPattern,
    AvoidRareCodons,
    CodonOptimize,
    DnaOptimizationProblem,
    EnforceGCContent,
    EnforceMeltingTemperature,
    EnforceTerminalGCContent,
    EnforceTranslation,
    MatchTargetCodonUsage,
    MaximizeCAI,
    UniquifyAllKmers,
    reverse_translate,
)

# ---------------------------------------------------------------------------
# Standard genetic code (used to map codons -> amino acids when loading tables)
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

# ---------------------------------------------------------------------------
# Input columns 
# ---------------------------------------------------------------------------
NAME_COL = "Name"
DNA_COL_CANDIDATES = ["Generated DNA Sequence", "DNA Sequence", "DNA",
                      "CDS Sequence", "CDS", "Nucleotide Sequence"]
PROT_COL_CANDIDATES = ["Protein Sequence", "Protein", "Amino Acid Sequence",
                       "AA Sequence", "Peptide"]
SHARED_COL_CANDIDATES = ["Sequence"]

AA_ALPHABET = STANDARD_AAS
DNA_ALPHABET = set("ACGTN")

# Audit windows for metrics
AUDIT_KMER = 10
GC_WINDOW_AUDIT = 50
TM_AUDIT_WINDOW = 30

# ---------------------------------------------------------------------------
# Species-/domain-specific motif screens. The grid is identical for all
# species; only this motif set switches with --domain.
#   eukaryote:  AWTAAA = poly-A signal (W=A/T), ATTTA = AU-rich element,
#               TATAAA = TATA box
#   prokaryote: AGGAGG = internal Shine-Dalgarno / ribosome-binding core
#               (removing it avoids spurious internal translation initiation /
#               ribosome pausing in E. coli-type hosts)
# ---------------------------------------------------------------------------
EUKARYOTE_MOTIFS = ["AWTAAA", "ATTTA", "TATAAA"]
PROKARYOTE_MOTIFS = ["AGGAGG"]

# Type-IIS + common cloning sites (shared across species)
DEFAULT_SITES = [
    "BsaI", "BsmBI", "BbsI", "SapI", "AarI",
    "EcoRI", "BamHI", "HindIII", "XhoI", "XbaI",
    "SalI", "NheI", "NcoI", "NdeI", "KpnI",
    "SacI", "PstI", "EcoRV", "HpaI", "SmaI",
    "AgeI", "MluI", "SpeI", "ClaI", "AflII",
    "BglII", "DraI", "PacI", "NotI", "AscI",
    "FseI", "PmeI", "SbfI", "SwaI",
]

# Cryptic splice-site consensus motifs (humans)
# Note:
# GTAAG and GTGAG are possible short subsequences of the human 5' donor consensus.
# TTTTTTAG and TTTCAG represent particular pyrimidine-rich/acceptor-like sequence patterns, although they do not cover the full diversity of human 3′ acceptors.
# "branch_point": "TACTAAC" removed as it is characteristic of budding yeast, not the human consensus.
SPLICE_SITES = {
    "5prime_donor_GTAAG": "GTAAG",
    "5prime_donor_GTGAG": "GTGAG",
    "3prime_acceptor_poly": "TTTTTAG",
    "3prime_acceptor_CAG": "TTTCAG",         
}

# ---------------------------------------------------------------------------
# Default search space (fixed for all species). "tm" and "terminal_gc" are
# grouped axes: None (off) or [window, lo, hi].
# ---------------------------------------------------------------------------
DEFAULT_CONFIG = {
    "baseline": {
        "codon_method": "match_codon_usage",
        "gc_min": 0.44,
        "gc_max": 0.68,
        "gc_window": 50,
        "gc_target": None,
        "kmer": 10,
        "rare_codon_min_freq": 0.11,
        "homopolymer_max": 5,
        "avoid_motifs": True,
        "avoid_restriction_sites": True,
	"avoid_splice_sites": True,
        "avoid_hairpins": True,
        "tm": None,
        "terminal_gc": None,
    },
    "grid": {
        "codon_method": [None, "use_best_codon", "match_codon_usage"],
        "gc_min": [0.30, 0.40, 0.46],
        "gc_max": [0.60, 0.66, 0.70],
        "gc_window": [50],
        "kmer": [None, 9, 10, 15],
        "rare_codon_min_freq": [None, 0.11],
        "homopolymer_max": [4, 5, 6],
        "avoid_motifs": [True],
        "avoid_restriction_sites": [True],
	"avoid_splice_sites": [True],
        "avoid_hairpins": [True],
        "tm": [None, [30, 55, 65]],
        "terminal_gc": [None, [30, 0.40, 0.60]],
    },
}

SCALAR_AXES = [
    "codon_method", "gc_min", "gc_max", "gc_window", "gc_target",
    "kmer", "rare_codon_min_freq", "homopolymer_max",
    "avoid_motifs", "avoid_restriction_sites",
    "avoid_splice_sites", "avoid_hairpins",
]
 
GROUPED_AXES = ["tm", "terminal_gc"]
AXIS_KEYS = SCALAR_AXES + GROUPED_AXES


# ===========================================================================
# Codon usage table loading (FILE ONLY - no species names, no downloads)
# Ported from the codon-scoring scripts' loader; the resolve-by-name /
# python_codon_tables branches are deliberately removed.
# ===========================================================================
def _normalise_table(raw):
    """Make every amino-acid family's fractions sum to 1.0; drop empty families.
    This is exactly the dict shape DNA Chisel's ``codon_usage_table`` expects."""
    tbl = {}
    for aa, codons in raw.items():
        codons = {c.upper(): float(v) for c, v in codons.items() if v is not None}
        total = sum(codons.values())
        if total <= 0:
            continue
        tbl[aa] = {c: v / total for c, v in codons.items()}
    return tbl


def _first_present(cols, candidates):
    for c in candidates:
        if c in cols:
            return c
    return None


def _load_table_csv(path):
    sep = "\t" if path.lower().endswith(".tsv") else ","
    df = pd.read_csv(path, sep=sep)
    df.columns = [c.strip().lower() for c in df.columns]
    aa_col = _first_present(df.columns, ["amino_acid", "aa", "aminoacid", "residue"])
    codon_col = _first_present(df.columns, ["codon", "triplet"])
    val_col = _first_present(df.columns,
                             ["fraction", "frequency", "freq", "count", "number", "proportion"])
    if codon_col is None or val_col is None:
        raise ValueError(
            "Codon-table CSV needs a 'codon' column and a value column "
            "(fraction/frequency/count). Got: {}".format(list(df.columns)))
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


def _load_table_json(path):
    with open(path) as fh:
        data = json.load(fh)
    if all(isinstance(v, dict) for v in data.values()):   # nested {aa:{codon:val}}
        return _normalise_table(data)
    raw = defaultdict(dict)                                # flat {codon:val}
    for codon, val in data.items():
        c = str(codon).strip().upper().replace("U", "T")
        if c in GENETIC_CODE:
            raw[GENETIC_CODE[c]][c] = float(val)
    return _normalise_table(raw)


def load_codon_table(path):
    """Load a host codon usage table from a FILE only. Accepts CSV/TSV
    (amino_acid,codon,fraction|frequency|count) or JSON (nested or flat).
    No species names, no TaxIDs, no downloads, no bundled/cached tables.
    Returns a normalised {aa: {codon: fraction}} dict for DNA Chisel."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Codon table file not found: {path}")
    low = path.lower()
    if low.endswith((".csv", ".tsv")):
        return _load_table_csv(path)
    if low.endswith(".json"):
        return _load_table_json(path)
    raise ValueError(f"Unrecognised codon-table file type (use .csv/.tsv/.json): {path}")


def validate_table(table):
    """Return the standard amino acids missing from the table (should be empty)."""
    return sorted(STANDARD_AAS - (set(table) - {"*"}))


# ===========================================================================
# Variant generation
# ===========================================================================
def _full_params(partial, baseline):
    params = {k: baseline.get(k) for k in AXIS_KEYS}
    params.update({k: v for k, v in partial.items() if k in AXIS_KEYS})
    return params


def _dedup_key(params):
    return tuple(tuple(params[k]) if isinstance(params[k], list) else params[k]
                 for k in AXIS_KEYS)


def _is_valid(params):
    gmin, gmax = params.get("gc_min"), params.get("gc_max")
    if gmin is not None and gmax is not None and gmin >= gmax:
        return False
    tm = params.get("tm")
    if tm is not None and tm[1] >= tm[2]:
        return False
    tg = params.get("terminal_gc")
    if tg is not None and tg[1] >= tg[2]:
        return False
    return True


def generate_cartesian(config, max_variants, seed):
    grid, baseline = config["grid"], config["baseline"]
    axes = list(grid.keys())
    combos, seen = [], set()
    for values in itertools.product(*(grid[a] for a in axes)):
        params = _full_params(dict(zip(axes, values)), baseline)
        if not _is_valid(params):
            continue
        key = _dedup_key(params)
        if key in seen:
            continue
        seen.add(key)
        combos.append(params)
    if max_variants and len(combos) > max_variants:
        rng = np.random.RandomState(seed)
        idx = sorted(rng.choice(len(combos), size=max_variants, replace=False))
        combos = [combos[i] for i in idx]
    return combos


def generate_ofat(config):
    grid, baseline = config["grid"], config["baseline"]
    combos = [_full_params({}, baseline)]
    seen = {_dedup_key(combos[0])}
    for axis, values in grid.items():
        for v in values:
            params = _full_params({axis: v}, baseline)
            if not _is_valid(params):
                continue
            key = _dedup_key(params)
            if key in seen:
                continue
            seen.add(key)
            combos.append(params)
    return combos


def name_variant(idx, params, baseline):
    short = {
        "codon_method": "codon", "gc_min": "gcmin", "gc_max": "gcmax",
        "gc_window": "gcwin", "gc_target": "gctgt", "kmer": "k",
        "rare_codon_min_freq": "rare", "homopolymer_max": "hp",
        "avoid_motifs": "motif", "avoid_restriction_sites": "site",
	"avoid_splice_sites": "splice",
        "avoid_hairpins": "hair", "tm": "tm", "terminal_gc": "tgc",
    }
    diffs = []
    for k in AXIS_KEYS:
        v, b = params.get(k), baseline.get(k)
        if v == b:
            continue
        if v is None:
            tag = f"{short[k]}off"
        elif k == "tm":
            tag = f"tm{v[1]:g}-{v[2]:g}"
        elif k == "terminal_gc":
            tag = f"tgc{v[1] * 100:g}-{v[2] * 100:g}"
        elif isinstance(v, bool):
            tag = f"{short[k]}{'Y' if v else 'N'}"
        elif isinstance(v, float):
            tag = f"{short[k]}{v:g}"
        else:
            tag = f"{short[k]}{v}"
        diffs.append(tag)
    label = "baseline" if not diffs else "_".join(diffs)
    return f"V{idx:03d}_{label}"


def flatten_params(params):
    out = {k: params.get(k) for k in SCALAR_AXES}
    tm = params.get("tm")
    out["tm_window"], out["tm_min"], out["tm_max"] = tm if tm else (None, None, None)
    tg = params.get("terminal_gc")
    (out["terminal_gc_window"], out["terminal_gc_min"],
     out["terminal_gc_max"]) = tg if tg else (None, None, None)
    return out


# ===========================================================================
# Problem construction
# ===========================================================================
def build_problem(sequence, params, table, motifs, sites):
    """Translate a variant's parameter dict into a DnaOptimizationProblem.

    `table` is a normalised {aa:{codon:fraction}} dict, passed to every
    codon-aware spec via codon_usage_table=.
    """
    L = len(sequence)
    constraints = [EnforceTranslation()]
    objectives = []

    if params.get("codon_method"):
        objectives.append(CodonOptimize(codon_usage_table=table,
                                        method=params["codon_method"]))

    gmin, gmax = params.get("gc_min"), params.get("gc_max")
    if gmin is not None or gmax is not None:
        constraints.append(EnforceGCContent(
            mini=gmin if gmin is not None else 0.0,
            maxi=gmax if gmax is not None else 1.0,
            window=params.get("gc_window")))

    if params.get("gc_target") is not None:
        objectives.append(EnforceGCContent(
            target=params["gc_target"], window=params.get("gc_window")))

    hp = params.get("homopolymer_max")
    if hp is not None:
        constraints += [AvoidPattern(f"{hp + 1}x{base}") for base in "ACGT"]

    if params.get("kmer"):
        constraints.append(UniquifyAllKmers(int(params["kmer"])))

    if params.get("rare_codon_min_freq") is not None:
        constraints.append(AvoidRareCodons(
            min_frequency=params["rare_codon_min_freq"], codon_usage_table=table))

    if params.get("avoid_motifs"):
    	constraints += [AvoidPattern(m) for m in motifs]
    
    if params.get("avoid_restriction_sites"):
    	constraints += [AvoidPattern(f"{enzyme}_site") for enzyme in sites]

    if params.get("avoid_splice_sites"):
    	constraints += [AvoidPattern(motif) for motif in SPLICE_SITES.values()]

    if params.get("avoid_hairpins"):
    	constraints.append(AvoidHairpins())

    tg = params.get("terminal_gc")
    if tg:
        w, lo, hi = tg
        w = min(int(w), L // 2) or 1
        constraints.append(EnforceTerminalGCContent(window_size=w, mini=lo, maxi=hi))

    tm = params.get("tm")
    if tm:
        w, lo, hi = tm
        w = min(int(w), L // 2) or 1
        constraints.append(EnforceMeltingTemperature(mini=lo, maxi=hi, location=(0, w)))
        constraints.append(EnforceMeltingTemperature(mini=lo, maxi=hi, location=(L - w, L)))

    return DnaOptimizationProblem(
        sequence=sequence, constraints=constraints, objectives=objectives, logger=None)


def run_variant(sequence, params, table, motifs, sites, seeds=(42,)):
    """Resolve + optimize one variant, retrying across seeds on solver failure.
    Returns (sequence, status, obj_score, seed_used)."""
    last_error = None
    for seed in seeds:
        np.random.seed(seed)
        try:
            problem = build_problem(sequence, params, table, motifs, sites)
            problem.resolve_constraints()
            problem.optimize()
            ok = problem.all_constraints_pass()
            status = "OK" if ok else "Constraints not met"
            try:
                score = round(problem.objectives_evaluations().scores_sum(), 3)
            except Exception:
                score = None
            return problem.sequence, status, score, seed
        except Exception as error:
            last_error = f"FAILED: {str(error)[:120]}"
    return sequence, last_error, None, None


# ===========================================================================
# Input loading (DNA or protein)
# ===========================================================================
def pick_seq_column(df, seq_col, input_type):
    if seq_col:
        if seq_col not in df.columns:
            raise ValueError(f"--seq-col '{seq_col}' not in columns: {list(df.columns)}")
        return seq_col
    candidates = []
    if input_type in ("auto", "dna"):
        candidates += DNA_COL_CANDIDATES
    if input_type in ("auto", "protein"):
        candidates += PROT_COL_CANDIDATES
    candidates += SHARED_COL_CANDIDATES
    for col in candidates:
        if col in df.columns:
            return col
    others = [c for c in df.columns if c != NAME_COL]
    if len(others) == 1:
        return others[0]
    raise ValueError("Could not auto-detect the sequence column; pass --seq-col. "
                     f"Columns present: {list(df.columns)}")


def classify_sequence(seq):
    chars = set(seq)
    if chars <= DNA_ALPHABET and len(seq) % 3 == 0:
        return "dna"
    if chars <= (AA_ALPHABET | {"*"}):
        return "protein"
    return "protein" if (chars - DNA_ALPHABET) else "dna"


def to_start_dna(raw, input_type):
    seq = "".join(str(raw).split()).upper()
    kind = classify_sequence(seq) if input_type == "auto" else input_type
    if kind == "protein":
        seq = seq.rstrip("*")
        bad = set(seq) - AA_ALPHABET
        if bad:
            raise ValueError(f"invalid amino acid(s): {sorted(bad)}")
        return reverse_translate(seq), "protein"
    seq = seq.replace("U", "T")
    bad = set(seq) - DNA_ALPHABET
    if bad:
        raise ValueError(f"invalid DNA base(s): {sorted(bad)}")
    if len(seq) % 3 != 0:
        raise ValueError(f"DNA length {len(seq)} is not a multiple of 3")
    return seq, "dna"


# ===========================================================================
# Metrics
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


def count_changes(before, after):
    return sum(a != b for a, b in zip(before, after))


def terminal_tm(seq, window=TM_AUDIT_WINDOW):
    try:
        w = min(window, len(seq) // 2) or len(seq)
        s5 = EnforceMeltingTemperature(mini=0, maxi=100, location=(0, w))
        s3 = EnforceMeltingTemperature(mini=0, maxi=100, location=(len(seq) - w, len(seq)))
        p = DnaOptimizationProblem(sequence=seq, constraints=[], objectives=[s5], logger=None)
        tm5 = float(re.search(r"Tm = ([\d.]+)", s5.evaluate(p).message).group(1))
        tm3 = float(re.search(r"Tm = ([\d.]+)", s3.evaluate(p).message).group(1))
        return round(tm5, 1), round(tm3, 1)
    except Exception:
        return None, None


def motifs_remaining(seq, motifs):
    iupac = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]",
             "S": "[GC]", "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]",
             "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}
    comp = str.maketrans("ACGT", "TGCA")
    rc = seq.translate(comp)[::-1]
    total = 0
    for motif in motifs:
        pat = "".join(iupac.get(b, b) for b in motif)
        total += len(re.findall(f"(?={pat})", seq)) + len(re.findall(f"(?={pat})", rc))
    return total

def splice_sites_remaining(seq):
    return motifs_remaining(seq, list(SPLICE_SITES.values()))

def cai_value(seq, table):
    try:
        spec = MaximizeCAI(codon_usage_table=table)
        p = DnaOptimizationProblem(sequence=seq, constraints=[EnforceTranslation()],
                                   objectives=[spec], logger=None)
        ev = p.objectives_evaluations().evaluations[0]
        n = len(seq) // 3
        return round(math.exp(ev.score / n), 4) if n else None
    except Exception:
        return None


def codon_match_score(seq, table):
    try:
        spec = MatchTargetCodonUsage(codon_usage_table=table)
        p = DnaOptimizationProblem(sequence=seq, constraints=[EnforceTranslation()],
                                   objectives=[spec], logger=None)
        return round(p.objectives_evaluations().evaluations[0].score, 3)
    except Exception:
        return None


def metrics_row(before, after, table, motifs):
    gc = gc_fraction(after) * 100
    wmin, wmax = window_gc_extremes(after)
    tm5, tm3 = terminal_tm(after)
    n_change = count_changes(before, after)
    return {
        "Nucleotides changed": n_change,
        "% identity": round(100 * (1 - n_change / len(before)), 2) if before else None,
        "GC% global": round(gc, 2),
        f"GC% min ({GC_WINDOW_AUDIT}nt win)": round(wmin, 2),
        f"GC% max ({GC_WINDOW_AUDIT}nt win)": round(wmax, 2),
        "Longest homopolymer": longest_homopolymer(after),
        f"Repeated {AUDIT_KMER}-mers": repeated_kmers(after),
        "CAI": cai_value(after, table),
        "Codon-match score": codon_match_score(after, table),
        f"5' Tm ({TM_AUDIT_WINDOW}nt)": tm5,
        f"3' Tm ({TM_AUDIT_WINDOW}nt)": tm3,
        "Forbidden motifs left": motifs_remaining(after, motifs),
	"Splice sites remaining": splice_sites_remaining(after),
    }


# ===========================================================================
# Main
# ===========================================================================
def load_config(path):
    cfg = json.loads(Path(path).read_text())
    merged = json.loads(json.dumps(DEFAULT_CONFIG))
    merged["baseline"].update(cfg.get("baseline", {}))
    if "grid" in cfg:
        merged["grid"] = cfg["grid"]
    for extra in ("motifs", "restriction_sites"):
        if extra in cfg:
            merged[extra] = cfg[extra]
    return merged


def main():
    parser = argparse.ArgumentParser(
        description="Grid/sweep driver for DNA Chisel codon-optimization variants.")
    parser.add_argument("--input", help="Excel with a 'Name' column + a sequence column")
    parser.add_argument("--output", help="Output .xlsx")
    parser.add_argument("--codon-table",
                        help="REQUIRED. Path to the host codon usage table FILE "
                             "(CSV/TSV: amino_acid,codon,fraction | JSON: nested or "
                             "flat). No species names or downloads are accepted.")
    parser.add_argument("--host-label", default=None,
                        help="Label for the host in the output (default: table file name).")
    parser.add_argument("--domain", choices=["eukaryote", "prokaryote"], default="eukaryote",
                        help="Selects the domain-specific motif screen: eukaryote = "
                             "poly-A/ARE/TATA (default); prokaryote = internal "
                             "Shine-Dalgarno. The grid is identical either way.")
    parser.add_argument("--config", help="JSON search-space config (see --dump-config)")
    parser.add_argument("--input-type", choices=["auto", "dna", "protein"], default="auto",
                        help="How to read the sequence column (default: auto-detect).")
    parser.add_argument("--seq-col", default=None,
                        help="Name of the sequence column (default: auto-detect).")
    parser.add_argument("--dump-config", metavar="PATH",
                        help="Write the default config to PATH and exit")
    parser.add_argument("--mode", choices=["cartesian", "ofat"], default="cartesian")
    parser.add_argument("--max-variants", type=int, default=0,
                        help="Cap for cartesian mode; 0 (default) = no cap, run the "
                             "full grid. A positive N takes a seeded random sample.")
    parser.add_argument("--limit", type=int, default=None,
                        help="Only first N input sequences (quick test)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--retries", type=int, default=3,
                        help="Extra random seeds to try if the solver dead-ends "
                             "on a variant (default 3; 0 disables retrying).")
    args = parser.parse_args()

    if args.dump_config:
        Path(args.dump_config).write_text(json.dumps(DEFAULT_CONFIG, indent=2))
        print(f"Wrote default config to {args.dump_config}")
        return
    if not args.input or not args.output:
        parser.error("--input and --output are required (unless --dump-config)")
    if not args.codon_table:
        parser.error("--codon-table is required: pass a CSV/TSV/JSON codon usage table file")

    # Load the host codon usage table (file only; injected into every codon spec)
    table = load_codon_table(args.codon_table)
    missing = validate_table(table)
    if missing:
        print(f"WARNING: codon table is missing amino-acid families {missing}; "
              "sequences using them will fail. Supply a complete table.")
    host_label = args.host_label or os.path.splitext(os.path.basename(args.codon_table))[0]

    config = load_config(args.config) if args.config else json.loads(json.dumps(DEFAULT_CONFIG))
    baseline = config["baseline"]
    # Motif screen: explicit config override, else domain default
    if "motifs" in config:
        motifs = config["motifs"]
    else:
        motifs = PROKARYOTE_MOTIFS if args.domain == "prokaryote" else EUKARYOTE_MOTIFS
    sites = config.get("restriction_sites", DEFAULT_SITES)

    if args.mode == "cartesian":
        combos = generate_cartesian(config, args.max_variants, args.seed)
    else:
        combos = generate_ofat(config)
    variants = [(name_variant(i + 1, p, baseline), p) for i, p in enumerate(combos)]
    print(f"Mode={args.mode}  host={host_label}  domain={args.domain}  "
          f"table_families={len(set(table) - {'*'})}  input-type={args.input_type}  "
          f"variants={len(variants)}")

    rng = np.random.RandomState(args.seed)
    fallback_seeds = [int(x) for x in rng.randint(1, 10 ** 6, size=max(args.retries, 0))]
    seeds = [args.seed] + fallback_seeds

    data = pd.read_excel(args.input)
    if NAME_COL not in data.columns:
        raise ValueError(f"Missing column: {NAME_COL}")
    seq_col = pick_seq_column(data, args.seq_col, args.input_type)
    print(f"Using sequence column: '{seq_col}'")
    if args.limit:
        data = data.head(args.limit)

    results = []
    for _, record in data.iterrows():
        name = record[NAME_COL]
        try:
            start, kind = to_start_dna(record[seq_col], args.input_type)
        except Exception as err:
            results.append({NAME_COL: name, "Variant": "-", "Host": host_label,
                            "Domain": args.domain, "Input kind": "?",
                            "Status": f"INPUT ERROR: {err}"})
            print(f"{str(name):<12} INPUT ERROR: {err}")
            continue

        for vname, params in variants:
            t0 = time.time()
            new_seq, status, obj_score, seed_used = run_variant(
                start, params, table, motifs, sites, seeds=seeds)
            row = {
                NAME_COL: name,
                "Variant": vname,
                "Host": host_label,
                "Domain": args.domain,
                "Input kind": kind,
                "Status": status,
                "Seed used": seed_used,
                "Objective score": obj_score,
                "Runtime (s)": round(time.time() - t0, 2),
            }
            row.update(flatten_params(params))
            if status.startswith("FAILED"):
                row["Optimized DNA Sequence"] = ""
            else:
                row.update(metrics_row(start, new_seq, table, motifs))
                row["Optimized DNA Sequence"] = new_seq
            results.append(row)
            print(f"{str(name):<12} {vname:<42} {status}")

    res = pd.DataFrame(results)

    graded = res[res["Variant"] != "-"].copy()
    metric_cols = [c for c in (
        "Nucleotides changed", "% identity", "GC% global",
        f"GC% min ({GC_WINDOW_AUDIT}nt win)", f"GC% max ({GC_WINDOW_AUDIT}nt win)",
        "Longest homopolymer", f"Repeated {AUDIT_KMER}-mers", "CAI",
        "Codon-match score", f"5' Tm ({TM_AUDIT_WINDOW}nt)", f"3' Tm ({TM_AUDIT_WINDOW}nt)",
        "Forbidden motifs left", "Runtime (s)") if c in graded.columns]
    if not graded.empty:
        summary = (graded.assign(_ok=(graded["Status"] == "OK").astype(int))
                   .groupby("Variant")
                   .agg(n_sequences=("Status", "size"),
                        n_passed=("_ok", "sum"),
                        **{c: (c, "mean") for c in metric_cols})
                   .reset_index())
        summary[metric_cols] = summary[metric_cols].round(2)
    else:
        summary = pd.DataFrame()

    grid_df = pd.DataFrame([{"Variant": v, **flatten_params(p)} for v, p in variants])

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(out) as writer:
        if not summary.empty:
            summary.to_excel(writer, sheet_name="Summary", index=False)
        res.to_excel(writer, sheet_name="Results", index=False)
        grid_df.to_excel(writer, sheet_name="Grid", index=False)
    print(f"\nSaved {len(res)} rows ({len(variants)} variants) to {out}")


if __name__ == "__main__":
    main()
