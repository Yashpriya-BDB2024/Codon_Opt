"""
biopython_compat.py — run Codon Harmony 1.0.0 on modern Biopython (>=1.78)
without downgrading or a separate environment.

Codon Harmony imports three things that newer Biopython removed/changed:
  * Bio.Alphabet (IUPAC, ProteinAlphabet, _get_base_alphabet)   — removed in 1.78
  * Bio.SeqUtils.CodonUsage (SynonymousCodons, CodonsDict,
    CodonAdaptationIndex)                                        — removed
  * Bio.SeqUtils.GC                                              — removed (-> gc_fraction)
  * Seq(data, alphabet) with a .alphabet attribute               — alphabet dropped

Importing THIS module *before* `codon_harmony` reinstates exactly those pieces
(as thin, behaviour-equivalent shims) so codon-harmony imports and runs
unchanged on the same Biopython your other tools use. Import order matters:

    import biopython_compat   # noqa: F401  (must come first)
    import codon_harmony
    from codon_harmony.codon_harmony import _harmonize_sequence

The CAI math and the SynonymousCodons / CodonsDict tables reproduce Biopython's
own definitions (BSD-licensed), rebuilt from the standard genetic code.
"""

import math
import sys
import types
from collections import defaultdict

import Bio
import Bio.SeqUtils as _sequtils
from Bio.Seq import Seq, MutableSeq
from Bio.SeqUtils import seq3

try:
    from Bio.SeqUtils import gc_fraction as _gc_fraction
except Exception:  # very old Biopython still has GC
    _gc_fraction = None

_GENETIC_CODE = {
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


# ---------------------------------------------------------------------------
# 1. Fake Bio.Alphabet
# ---------------------------------------------------------------------------
class _BaseAlphabet:
    pass


class ProteinAlphabet(_BaseAlphabet):
    pass


class DNAAlphabet(_BaseAlphabet):
    pass


class RNAAlphabet(_BaseAlphabet):
    pass


_PROT, _DNA, _RNA = ProteinAlphabet(), DNAAlphabet(), RNAAlphabet()

if "Bio.Alphabet" not in sys.modules:
    _alpha_mod = types.ModuleType("Bio.Alphabet")
    _alpha_mod.IUPAC = types.SimpleNamespace(
        protein=_PROT, unambiguous_dna=_DNA, unambiguous_rna=_RNA,
        ambiguous_dna=_DNA, ambiguous_rna=_RNA)
    _alpha_mod.ProteinAlphabet = ProteinAlphabet
    _alpha_mod.DNAAlphabet = DNAAlphabet
    _alpha_mod.RNAAlphabet = RNAAlphabet
    _alpha_mod.Alphabet = _BaseAlphabet
    _alpha_mod._get_base_alphabet = lambda a: a
    sys.modules["Bio.Alphabet"] = _alpha_mod
    Bio.Alphabet = _alpha_mod


# ---------------------------------------------------------------------------
# 2. Make Seq tolerate a second (alphabet) argument and carry .alphabet
# ---------------------------------------------------------------------------
if not getattr(Seq, "_compat_alphabet_patched", False):
    _orig_seq_init = Seq.__init__

    def _seq_init(self, data, alphabet=None, length=None, **kw):
        if length is None:
            _orig_seq_init(self, data)
        else:
            _orig_seq_init(self, data, length)
        self.__dict__["alphabet"] = alphabet if alphabet is not None else _DNA

    Seq.__init__ = _seq_init
    # default for any Seq created elsewhere without going through __init__
    if not hasattr(Seq, "alphabet"):
        Seq.alphabet = _DNA
    Seq._compat_alphabet_patched = True

# Seq.tomutable() (removed 1.79) and MutableSeq.toseq() (removed) — Codon Harmony
# round-trips through these constantly in its cleanup routines.
if not hasattr(Seq, "tomutable"):
    Seq.tomutable = lambda self: MutableSeq(str(self))
if not hasattr(MutableSeq, "toseq"):
    MutableSeq.toseq = lambda self: Seq(str(self))


# ---------------------------------------------------------------------------
# 3. Reconstruct Bio.SeqUtils.CodonUsage
# ---------------------------------------------------------------------------
# Biopython's exact SynonymousCodons (Bio.SeqUtils.CodonUsageIndices), reproduced
# verbatim — family order AND per-family codon order must match, because Codon
# Harmony's random.choices() draws and its dict-iteration order depend on them.
# Getting this exactly right makes the shim byte-identical to real Biopython.
_SYNONYMOUS_CODONS = {
    "CYS": ["TGT", "TGC"],
    "ASP": ["GAT", "GAC"],
    "SER": ["TCT", "TCG", "TCA", "TCC", "AGC", "AGT"],
    "GLN": ["CAA", "CAG"],
    "MET": ["ATG"],
    "ASN": ["AAC", "AAT"],
    "PRO": ["CCT", "CCG", "CCA", "CCC"],
    "LYS": ["AAG", "AAA"],
    "STOP": ["TAG", "TGA", "TAA"],
    "THR": ["ACC", "ACA", "ACG", "ACT"],
    "PHE": ["TTT", "TTC"],
    "ALA": ["GCA", "GCC", "GCG", "GCT"],
    "GLY": ["GGT", "GGG", "GGA", "GGC"],
    "ILE": ["ATC", "ATA", "ATT"],
    "LEU": ["TTA", "TTG", "CTC", "CTT", "CTG", "CTA"],
    "HIS": ["CAT", "CAC"],
    "ARG": ["CGA", "CGC", "CGG", "CGT", "AGG", "AGA"],
    "TRP": ["TGG"],
    "VAL": ["GTA", "GTC", "GTG", "GTT"],
    "GLU": ["GAG", "GAA"],
    "TYR": ["TAT", "TAC"],
}


def _build_synonymous_codons():
    return {k: list(v) for k, v in _SYNONYMOUS_CODONS.items()}


class CodonAdaptationIndex:
    """Minimal drop-in for Bio.SeqUtils.CodonUsage.CodonAdaptationIndex.
    Codon Harmony only calls set_cai_index(index) and cai_for_gene(dna)."""

    def __init__(self):
        self.index = {}

    def set_cai_index(self, index):
        self.index = dict(index)

    def cai_for_gene(self, dna_sequence):
        dna_sequence = str(dna_sequence).upper().replace("U", "T")
        total, n = 0.0, 0
        for i in range(0, len(dna_sequence) - 2, 3):
            codon = dna_sequence[i:i + 3]
            # skip non-degenerate residues (Met/Trp) and anything without a
            # positive relative-adaptiveness weight, matching Biopython
            if codon in ("ATG", "TGG"):
                continue
            w = self.index.get(codon, 0.0)
            if w > 0:
                total += math.log(w)
                n += 1
        return math.exp(total / n) if n else 0.0


if not hasattr(_sequtils, "CodonUsage"):
    _cu = types.ModuleType("Bio.SeqUtils.CodonUsage")
    _cu.SynonymousCodons = _build_synonymous_codons()
    _cu.CodonsDict = {codon: 0 for codon in _GENETIC_CODE}
    _cu.CodonAdaptationIndex = CodonAdaptationIndex
    sys.modules["Bio.SeqUtils.CodonUsage"] = _cu
    _sequtils.CodonUsage = _cu


# ---------------------------------------------------------------------------
# 4. Bio.SeqUtils.GC  (removed; reimplement via gc_fraction)
# ---------------------------------------------------------------------------
if not hasattr(_sequtils, "GC"):
    def GC(seq):
        s = str(seq)
        if _gc_fraction is not None:
            return _gc_fraction(s) * 100
        s = s.upper()
        return 100.0 * (s.count("G") + s.count("C")) / len(s) if s else 0.0

    _sequtils.GC = GC
