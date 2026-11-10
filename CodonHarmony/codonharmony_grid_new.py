"""
codonharmony_grid_new.py   (derived from the project lead's codonharmony_grid.py)

>>> CHANGES vs codonharmony_grid.py - made ONLY for comparability with the DNA Chisel
>>> grid (codon_grid_search_new.py + run_human_parallel_new.py). Codon Harmony's
>>> optimisation logic, grid defaults, seeding, de-duplication and translation QC
>>> are UNTOUCHED. Every edited spot is tagged  # [NEW] / # [CHANGED]  in the code.
>>>   1. Restriction enzymes: 13 -> the same 34 as DNA Chisel DEFAULT_SITES.
>>>   2. Post-hoc audit columns added: Splice sites remaining (DNA Chisel's SPLICE_SITES),
>>>      Restriction sites remaining, Nucleotides changed, % identity, Codon-match score.
>>>      CAI is the lead's original calculation, unchanged.
>>>   3. Output schema aligned with DNA Chisel: Variant / Host / Domain / Runtime (s)
>>>      columns; sheets Summary (per variant) + Results + Grid (+ "Distinct per protein").
>>>   4. New optional flag --host-threshold (default = lead's 0.10). Homopolymer values
>>>      stay 3,4,5 by default; use --homopolymer 4,5,6 to match DNA Chisel's range.
>>>   5. Bug fix: scratch folder now tempfile.gettempdir() (the old "/tmp" fallback
>>>      fails on Windows).
>>>   6. Row/summary building moved into reusable functions (run_cell_row, build_summary,
>>>      distinct_per_protein) so run_human_parallel_codonharmony.py produces
>>>      byte-identical results to a sequential run.

Grid / sweep driver for Codon Harmony - the third codon-optimization backend,
alongside codon_grid_search.py (DNA Chisel) and codontransformer_grid.py.

Codon Harmony is a harmonization OPTIMIZER: each run returns a single best-CAI
sequence that matches the host codon-usage profile while removing sequence
liabilities. It is not a multi-output generator, so diversity comes from the
GRID itself - every (inner_cycles, cycles, max_relax, homopolymer) cell is one
run and one candidate. Its sampling uses Python's stdlib `random`, so we seed it
per cell (seed = --seed + cell_index) to make the whole grid reproducible while
keeping each cell a distinct draw. Candidates are de-duplicated across cells and
the number of DISTINCT sequences per protein is reported.

Species / codon usage table
---------------------------
Pass --species {human, cho, ecoli}. Each maps (hardcoded registry below) to a
CSV in --table-dir, which must hold the three CoCoPUTs tables you built for the
DNA Chisel / scoring work (amino_acid,codon,fraction - the same files). No
--host, no download. The selected CSV is converted on the fly to Codon Harmony's
own --local-host-profile JSON format (nested {aa:{RNA_codon:freq}}) and fed in.
Add more species by extending SPECIES_TABLES and DOMAIN.

Fixed settings (per the agreed design)
--------------------------------------
    host_threshold   = 0.10           (codon-harmony default)
    remove splice    = True           (mammalian; human/cho)
    remove start     = False          (bacterial-only; off here)
    restriction enz  = DEFAULT_RESTRICTION_ENZYMES  (your DNA Chisel list; editable)
    no R / repeats   - the grid provides the variants

Environment
-----------
Codon Harmony 1.0.0 imports APIs removed in Biopython >= 1.78 (Bio.Alphabet,
Bio.SeqUtils.CodonUsage, GC, Seq.tomutable). Rather than pin an old Biopython or
use a separate env, this script imports `biopython_compat` (shipped alongside it)
FIRST, which reinstates exactly those pieces so codon-harmony runs unchanged on a
MODERN Biopython - verified byte-for-byte identical to a real Biopython-1.77 run.
Keep biopython_compat.py next to this file. Just:
    pip install codon-harmony pandas numpy openpyxl       # on your normal env
    python codonharmony_grid.py --input proteins.xlsx \
        --output cho_ch.xlsx --species cho --table-dir codon_tables

Metrics per sequence mirror the DNA Chisel / CodonTransformer grids (CAI, GC,
GC-window, homopolymer, repeated 10-mers, 5'/3' Tm, forbidden motifs) so all
three tools land on one comparable yardstick.
"""

import argparse
import itertools
import json
import math
import os
import random
import re
import tempfile
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

# Compatibility shim FIRST: lets Codon Harmony run on a modern Biopython (it uses
# APIs removed in >=1.78). Must be imported before codon_harmony, which this file
# imports lazily inside the backend functions. Keep biopython_compat.py alongside.
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import biopython_compat  # noqa: F401,E402

# ===========================================================================
# TOP-LEVEL CONFIG (edit here)
# ===========================================================================
# species -> CSV filename (in --table-dir). The CSVs are your CoCoPUTs tables
# (amino_acid,codon,fraction), the same ones used for DNA Chisel.
SPECIES_TABLES = {
    "human": "human_codon.csv",
    "cho":   "cho_codon.csv",
    "ecoli": "ecoli_codon.csv",
}
# eukaryote vs prokaryote decides the motif screen (and splice/start defaults)
DOMAIN = {"human": "eukaryote", "cho": "eukaryote", "ecoli": "prokaryote"}

# Restriction sites to remove.
# [CHANGED] was 13 enzymes; now the SAME 34-enzyme list as DEFAULT_SITES in
# codon_grid_search_new.py (DNA Chisel) so both tools avoid identical sites.
# Names must exist in Biopython's Bio.Restriction (load_backend validates them).
DEFAULT_RESTRICTION_ENZYMES = [
    "BsaI", "BsmBI", "BbsI", "SapI", "AarI",
    "EcoRI", "BamHI", "HindIII", "XhoI", "XbaI",
    "SalI", "NheI", "NcoI", "NdeI", "KpnI",
    "SacI", "PstI", "EcoRV", "HpaI", "SmaI",
    "AgeI", "MluI", "SpeI", "ClaI", "AflII",
    "BglII", "DraI", "PacI", "NotI", "AscI",
    "FseI", "PmeI", "SbfI", "SwaI",
]

# [NEW] DNA Chisel's cryptic splice-site motifs, copied verbatim from
# codon_grid_search_new.py. Codon Harmony removes splice sites with its OWN
# hard-coded patterns (not configurable), so these are used only to AUDIT the
# output ("Splice sites remaining") on the same yardstick as DNA Chisel.
SPLICE_SITES = {
    "5prime_donor_GTAAG": "GTAAG",
    "5prime_donor_GTGAG": "GTGAG",
    "3prime_acceptor_poly": "TTTTTAG",
    "3prime_acceptor_CAG": "TTTCAG",
}

# Grid axes (swept as a Cartesian product). Override on the CLI if needed.
GRID = {
    "inner_cycles": [5, 10],
    "cycles": [5, 10, 20],
    "max_relax": [0.1, 0.2, 0.3],
    "local_homopolymer_threshold": [3, 4, 5],
}

# Fixed settings
HOST_THRESHOLD = 0.10                   # codon-harmony default (override with --host-threshold)
REMOVE_SPLICE_SITES_DEFAULT = True      # eukaryote/mammalian
REMOVE_START_SITES_DEFAULT = False      # bacterial only

# Motif audit (post-hoc, for cross-tool comparability; codon-harmony itself has
# no generic motif list - it removes restriction/splice/start/repeat/homopolymer/hairpin)
EUKARYOTE_MOTIFS = ["AWTAAA", "ATTTA", "TATAAA"]
PROKARYOTE_MOTIFS = ["AGGAGG"]

AUDIT_KMER = 10
GC_WINDOW_AUDIT = 50
TM_AUDIT_WINDOW = 30

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
PROT_COL_CANDIDATES = ["Protein Sequence", "Protein", "Amino Acid Sequence", "AA Sequence", "Peptide"]
DNA_COL_CANDIDATES = ["Generated DNA Sequence", "DNA Sequence", "DNA", "CDS", "Sequence"]


# ===========================================================================
# Table loading + conversion to Codon Harmony's profile format
# ===========================================================================
def load_dnachisel_table(path):
    """Read a DNA Chisel codon table -> nested {aa: {DNA_codon: fraction}}."""
    low = path.lower()
    if low.endswith(".json"):
        data = json.loads(Path(path).read_text())
        flat = ({c.upper().replace("U", "T"): f for cod in data.values() for c, f in cod.items()}
                if all(isinstance(v, dict) for v in data.values())
                else {c.upper().replace("U", "T"): f for c, f in data.items()})
    else:
        df = pd.read_csv(path, sep="\t" if low.endswith(".tsv") else ",")
        df.columns = [c.strip().lower() for c in df.columns]
        cc = next((c for c in ("codon", "triplet") if c in df.columns), None)
        vc = next((c for c in ("fraction", "frequency", "freq", "count", "proportion") if c in df.columns), None)
        if not cc or not vc:
            raise ValueError(f"Codon-table CSV needs 'codon' + a value column. Got {list(df.columns)}")
        flat = {}
        for _, r in df.iterrows():
            c = str(r[cc]).strip().upper().replace("U", "T")
            if c in GENETIC_CODE:
                try:
                    flat[c] = float(r[vc])
                except (ValueError, TypeError):
                    pass
    nested = defaultdict(dict)
    for codon, freq in flat.items():
        nested[GENETIC_CODE[codon]][codon] = freq
    # per-family normalise (CAI + profile both want within-family fractions)
    out = {}
    for aa, codons in nested.items():
        tot = sum(codons.values())
        if tot > 0:
            out[aa] = {c: v / tot for c, v in codons.items()}
    return out


def to_codon_harmony_profile(nested_dna_table, out_path):
    """Write Codon Harmony's --local-host-profile: nested {aa:{RNA_codon:freq}}."""
    prof = {aa: {c.replace("T", "U"): f for c, f in codons.items()}
            for aa, codons in nested_dna_table.items()}
    Path(out_path).write_text(json.dumps(prof, indent=1))
    return out_path


# ===========================================================================
# Sequence helpers / metrics (self-contained; mirror the other two grids)
# ===========================================================================
def translate(dna):
    dna = dna.upper().replace("U", "T")
    return "".join(GENETIC_CODE.get(dna[i:i + 3], "X") for i in range(0, len(dna) - 2, 3))


def protein_from_cell(raw, input_type):
    seq = "".join(str(raw).split()).upper()
    kind = input_type
    if input_type == "auto":
        kind = "dna" if (set(seq) <= set("ACGTUN") and len(seq) % 3 == 0) else "protein"
    if kind == "dna":
        prot = translate(seq).rstrip("*")
        if "X" in prot or "*" in prot:
            raise ValueError("DNA does not translate cleanly")
        return prot, "dna->protein"
    prot = seq.rstrip("*")
    bad = set(prot) - STANDARD_AAS
    if bad:
        raise ValueError(f"invalid amino acid(s): {sorted(bad)}")
    return prot, "protein"


def pick_seq_column(df, seq_col, input_type):
    if seq_col:
        if seq_col not in df.columns:
            raise ValueError(f"--seq-col '{seq_col}' not in {list(df.columns)}")
        return seq_col
    for col in (PROT_COL_CANDIDATES + DNA_COL_CANDIDATES):
        if col in df.columns:
            return col
    others = [c for c in df.columns if c != NAME_COL]
    if len(others) == 1:
        return others[0]
    raise ValueError(f"Could not auto-detect sequence column; pass --seq-col. Columns: {list(df.columns)}")


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


def cai(dna, nested_table):
    wmax = {aa: max(c.values()) for aa, c in nested_table.items()}
    logs = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3]
        aa = GENETIC_CODE.get(codon)
        if aa in (None, "*"):
            continue
        fam = nested_table.get(aa)
        if not fam or len(fam) < 2:
            continue
        f = fam.get(codon)
        if f and wmax[aa] > 0:
            logs.append(math.log(f / wmax[aa]))
    return round(math.exp(sum(logs) / len(logs)), 4) if logs else None


# [NEW] ----------------------------------------------------------------------
def splice_sites_remaining(seq):
    return motifs_remaining(seq, list(SPLICE_SITES.values()))


def restriction_sites_remaining(seq, enzymes=None):
    """Recognition sites (both strands) of the avoided enzymes still present."""
    try:
        import Bio.Restriction as BR
        sites = [getattr(BR, e).site for e in (enzymes or DEFAULT_RESTRICTION_ENZYMES)]
        return motifs_remaining(seq, sites)
    except Exception:
        return None


def start_dna_from_protein(protein):
    """Same deterministic starting DNA DNA Chisel uses, so 'Nucleotides changed' /
    '% identity' mean the same thing in both tools. Needs dnachisel (optional)."""
    try:
        from dnachisel import reverse_translate
        return reverse_translate(protein)
    except Exception:
        return None


def codon_match_score(dna, nested_table):
    """DNA Chisel MatchTargetCodonUsage score (optional - needs dnachisel)."""
    try:
        from dnachisel import DnaOptimizationProblem, EnforceTranslation, MatchTargetCodonUsage
        spec = MatchTargetCodonUsage(codon_usage_table=nested_table)
        p = DnaOptimizationProblem(sequence=dna, constraints=[EnforceTranslation()],
                                   objectives=[spec], logger=None)
        return round(p.objectives_evaluations().evaluations[0].score, 3)
    except Exception:
        return None


def count_changes(before, after):
    return sum(a != b for a, b in zip(before, after))
# ---------------------------------------------------------------------------


def seq_metrics(dna, protein, nested_table, motifs):
    dna = dna.upper().replace("U", "T")
    valid = translate(dna).rstrip("*") == protein
    wmin, wmax = window_gc_extremes(dna)
    tm5, tm3 = terminal_tm(dna)
    before = start_dna_from_protein(protein)          # [NEW]
    n_change = count_changes(before, dna) if before else None
    return {
        "Translates to protein": valid,
        "Length (nt)": len(dna),
        "Nucleotides changed": n_change,                                         # [NEW]
        "% identity": round(100 * (1 - n_change / len(before)), 2) if before else None,  # [NEW]
        "GC% global": round(gc_fraction(dna) * 100, 2),
        f"GC% min ({GC_WINDOW_AUDIT}nt win)": round(wmin, 2),
        f"GC% max ({GC_WINDOW_AUDIT}nt win)": round(wmax, 2),
        "Longest homopolymer": longest_homopolymer(dna),
        f"Repeated {AUDIT_KMER}-mers": repeated_kmers(dna),
        "CAI": cai(dna, nested_table),                                           # lead's original, unchanged
        "Codon-match score": codon_match_score(dna, nested_table),               # [NEW]
        f"5' Tm ({TM_AUDIT_WINDOW}nt)": tm5,
        f"3' Tm ({TM_AUDIT_WINDOW}nt)": tm3,
        "Forbidden motifs left": motifs_remaining(dna, motifs),
        "Splice sites remaining": splice_sites_remaining(dna),                   # [NEW]
        "Restriction sites remaining": restriction_sites_remaining(dna),         # [NEW]
    }


# ===========================================================================
# Codon Harmony backend (isolated; imports happen here so the rest is testable)
# ===========================================================================
def load_backend(profile_json, restriction_enzymes, host_threshold=HOST_THRESHOLD):
    """Return (cut, host_profile, cra, rest_enz) from Codon Harmony."""
    from codon_harmony.util import codon_use
    from codon_harmony.data import RestrictionEnzymes
    cut, host_profile, cra = codon_use.host_codon_usage(
        host="0", threshold=host_threshold, table_path=profile_json)
    import Bio.Restriction as BR                                                  # [NEW] fail early on a bad name
    bad = [e for e in restriction_enzymes if not hasattr(BR, e)]
    if bad:
        raise ValueError(f"Enzyme name(s) not found in Bio.Restriction: {bad}")
    rest_enz = RestrictionEnzymes(restriction_enzymes)
    return cut, host_profile, cra, rest_enz


def harmonize_one(protein, cell, cut, host_profile, cra, rest_enz,
                  splice_sites, start_sites, seed):
    """One Codon Harmony run (one cell, one seed) -> (dna, best_cai) or (None, err)."""
    import contextlib, io
    import codon_harmony  # applies Seq.back_translate monkeypatch
    from codon_harmony.codon_harmony import _harmonize_sequence
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    from types import SimpleNamespace

    random.seed(seed)                       # reproducible, per-cell
    args = SimpleNamespace(
        cycles=cell["cycles"], inner_cycles=cell["inner_cycles"],
        max_relax=cell["max_relax"],
        local_homopolymer_threshold=cell["local_homopolymer_threshold"],
        splice_sites=splice_sites, start_sites=start_sites)
    rec = SeqRecord(Seq(protein, IUPAC.protein), id="seq", name="seq", description="")
    try:
        with contextlib.redirect_stderr(io.StringIO()):
            out = _harmonize_sequence(rec, args, cut, host_profile, cra, rest_enz)
        if "dna" not in out:
            return None, "no solution within max_relax"
        fasta = out["dna"]
        dna = "".join(l for l in fasta.splitlines() if not l.startswith(">")).upper()
        return dna, None
    except Exception as e:
        return None, f"FAILED: {str(e)[:120]}"


# ===========================================================================
# Grid run
# ===========================================================================
def build_cells(grid):
    axes = list(grid.keys())
    cells = []
    for values in itertools.product(*(grid[a] for a in axes)):
        cells.append(dict(zip(axes, values)))
    return cells


def name_cell(ci, cell):
    """[NEW] DNA Chisel-style 'Variant' label so both tools share a join key."""
    return (f"C{ci + 1:03d}_ic{cell['inner_cycles']}_c{cell['cycles']}"
            f"_rel{cell['max_relax']:g}_hp{cell['local_homopolymer_threshold']}")


def run_cell_row(name, protein, kind, ci, cell, backend, nested_table, motifs,
                 splice_sites, start_sites, species, base_seed, domain=None):
    """[NEW] One (protein, grid cell) -> one result row + the DNA (or None).
    Seed = base_seed + cell index, so the result does not depend on run order or on
    which process runs it (this is what makes the parallel runner reproducible)."""
    cut, host_profile, cra, rest_enz = backend
    t0 = time.time()
    dna, err = harmonize_one(protein, cell, cut, host_profile, cra, rest_enz,
                             splice_sites, start_sites, seed=base_seed + ci)
    row = {NAME_COL: name, "Variant": name_cell(ci, cell),
           "Host": f"{species}_codon", "Domain": domain,
           "Species": species, "Input kind": kind,
           "Cell": ci + 1, **cell, "Seed": base_seed + ci,
           "Runtime (s)": round(time.time() - t0, 2)}
    if err:
        row["Status"] = err
        row["Optimized DNA Sequence"] = ""                  # same as DNA Chisel FAILED rows
        return row, None
    m = seq_metrics(dna, protein, nested_table, motifs)
    row["Status"] = "OK" if m["Translates to protein"] else "TRANSLATION MISMATCH"
    row.update(m)
    row["Optimized DNA Sequence"] = dna
    return row, (dna if m["Translates to protein"] else None)


def distinct_per_protein(seq_rows, n_cells, species):
    """Lead's per-protein summary: distinct VALID sequences across the grid."""
    distinct = defaultdict(set)
    names = []
    for r in seq_rows:
        if r[NAME_COL] not in names:
            names.append(r[NAME_COL])
        if r["Status"] == "OK":
            distinct[r[NAME_COL]].add(r["Optimized DNA Sequence"])
    return [{NAME_COL: n, "Species": species, "Grid cells": n_cells,
             "Distinct sequences": len(distinct[n])} for n in names]


def build_summary(res):
    """[NEW] Per-Variant summary, same metrics / order / rounding as the DNA Chisel workbook."""
    ok = res["Status"] == "OK"
    metric_cols = [c for c in (
        "Nucleotides changed", "% identity", "GC% global",
        f"GC% min ({GC_WINDOW_AUDIT}nt win)", f"GC% max ({GC_WINDOW_AUDIT}nt win)",
        "Longest homopolymer", f"Repeated {AUDIT_KMER}-mers", "CAI",
        "Codon-match score", f"5' Tm ({TM_AUDIT_WINDOW}nt)", f"3' Tm ({TM_AUDIT_WINDOW}nt)",
        "Forbidden motifs left", "Splice sites remaining", "Restriction sites remaining",
        "Runtime (s)") if c in res.columns]
    summary = (res.assign(_ok=ok.astype(int)).groupby("Variant", sort=False)
               .agg(n_sequences=("Status", "size"), n_passed=("_ok", "sum"),
                    **{c: (c, "mean") for c in metric_cols})
               .reset_index())
    summary[metric_cols] = summary[metric_cols].round(2)
    return summary


def write_workbook(out, seq_rows, cells, species):
    """[NEW] Summary | Results | Grid | Distinct per protein (shared by sequential + parallel runs)."""
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    res = pd.DataFrame(seq_rows)
    grid_df = pd.DataFrame([{"Variant": name_cell(i, c), **c} for i, c in enumerate(cells)])
    per_protein = distinct_per_protein(seq_rows, len(cells), species)
    with pd.ExcelWriter(out) as writer:
        build_summary(res).to_excel(writer, sheet_name="Summary", index=False)
        res.to_excel(writer, sheet_name="Results", index=False)
        grid_df.to_excel(writer, sheet_name="Grid", index=False)
        pd.DataFrame(per_protein).to_excel(writer, sheet_name="Distinct per protein", index=False)
    return per_protein


def run_grid(records, cells, backend, nested_table, motifs,
             splice_sites, start_sites, species, base_seed, domain=None):
    seq_rows = []
    for name, protein, kind in records:
        for ci, cell in enumerate(cells):
            row, _ = run_cell_row(name, protein, kind, ci, cell, backend, nested_table, motifs,
                                  splice_sites, start_sites, species, base_seed, domain)
            seq_rows.append(row)
            print(f"{str(name):<12} cell {ci+1:>2}/{len(cells)} "
                  f"ic={cell['inner_cycles']} c={cell['cycles']} relax={cell['max_relax']} "
                  f"hp={cell['local_homopolymer_threshold']} -> " +
                  ("OK" if row["Status"] == "OK" else row["Status"]))
    return seq_rows, distinct_per_protein(seq_rows, len(cells), species)


def main():
    ap = argparse.ArgumentParser(description="Codon Harmony inner_cycles x cycles x max_relax x homopolymer grid.")
    ap.add_argument("--input", required=True, help="Excel with 'Name' + a protein/DNA column")
    ap.add_argument("--output", required=True)
    ap.add_argument("--species", required=True, choices=list(SPECIES_TABLES),
                    help="Selects the hardcoded codon table (human/cho/ecoli)")
    ap.add_argument("--table-dir", default="codon_tables",
                    help="Folder holding the per-species CoCoPUTs CSVs (default: ./codon_tables)")
    ap.add_argument("--input-type", choices=["auto", "protein", "dna"], default="auto")
    ap.add_argument("--seq-col", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--seed", type=int, default=42, help="Base seed; cell i uses seed+i")
    # optional grid overrides
    ap.add_argument("--inner-cycles", default=None, help="Comma-separated override")
    ap.add_argument("--cycles", default=None)
    ap.add_argument("--max-relax", default=None)
    ap.add_argument("--homopolymer", default=None)
    ap.add_argument("--host-threshold", type=float, default=HOST_THRESHOLD,   # [NEW]
                    help="Lowest codon fraction per amino acid kept (default 0.10 = lead's value)")
    args = ap.parse_args()

    # resolve codon table
    table_path = os.path.join(args.table_dir, SPECIES_TABLES[args.species])
    if not os.path.exists(table_path):
        raise FileNotFoundError(
            f"Codon table for '{args.species}' not found: {table_path}\n"
            f"Place your CoCoPUTs CSVs in '{args.table_dir}/' "
            f"({', '.join(SPECIES_TABLES.values())}).")
    nested_table = load_dnachisel_table(table_path)
    missing = sorted(STANDARD_AAS - set(nested_table))
    if missing:
        print(f"WARNING: table missing amino-acid families {missing}")
    scratch = os.environ.get("TMPDIR") or tempfile.gettempdir()              # [CHANGED] "/tmp" fails on Windows
    profile_json = to_codon_harmony_profile(nested_table, os.path.join(scratch, f"_ch_profile_{args.species}.json"))

    domain = DOMAIN[args.species]
    motifs = PROKARYOTE_MOTIFS if domain == "prokaryote" else EUKARYOTE_MOTIFS
    splice_sites = domain == "eukaryote"
    start_sites = domain == "prokaryote"

    grid = {k: list(v) for k, v in GRID.items()}
    if args.inner_cycles:
        grid["inner_cycles"] = [int(x) for x in args.inner_cycles.split(",")]
    if args.cycles:
        grid["cycles"] = [int(x) for x in args.cycles.split(",")]
    if args.max_relax:
        grid["max_relax"] = [float(x) for x in args.max_relax.split(",")]
    if args.homopolymer:
        grid["local_homopolymer_threshold"] = [int(x) for x in args.homopolymer.split(",")]
    cells = build_cells(grid)

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

    print(f"species={args.species} ({domain})  proteins={len(records)}  grid cells={len(cells)}  "
          f"splice={splice_sites} start={start_sites}  enzymes={len(DEFAULT_RESTRICTION_ENZYMES)}")

    backend = load_backend(profile_json, DEFAULT_RESTRICTION_ENZYMES, args.host_threshold)
    seq_rows, per_protein = run_grid(records, cells, backend, nested_table, motifs,
                                     splice_sites, start_sites, args.species, args.seed,
                                     domain=domain)

    out = Path(args.output)
    per_protein = write_workbook(out, seq_rows, cells, args.species)     # [CHANGED] see write_workbook
    print(f"\nSaved {len(seq_rows)} rows ({len(cells)} cells) to {out}")
    for p in per_protein:
        print(f"  {p[NAME_COL]}: {p['Distinct sequences']} distinct of {p['Grid cells']} cells")


if __name__ == "__main__":
    main()
