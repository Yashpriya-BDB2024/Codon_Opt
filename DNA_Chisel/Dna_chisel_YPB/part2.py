"""
analyze_sequences_part2.py

Extra sequence checks (local GC, codon usage, rare codons, homopolymers,
tandem repeats, unwanted motifs, restriction sites, cryptic splice sites,
hairpins) for the Homo sapiens DNA sequences.

Everything is pure Python (pandas, python_codon_tables), so no web
tools or non-commercial software are needed.

Usage:
    python analyze_sequences_part2.py --input human_codon_sequences.xlsx --output part2_report.xlsx
"""

import argparse
import math
import re
from collections import Counter
from pathlib import Path

import pandas as pd
from python_codon_tables import get_codons_table

# ----------------------------------------------------------------------------
# SETTINGS - edit these to change cutoffs / motif lists
# ----------------------------------------------------------------------------
GC_WINDOW = 50           # sliding window size (bp) for local GC
GC_MIN, GC_MAX = 40, 60  # acceptable local GC range (%)
KMER_SIZE = 12           # k-mer size used to find repeated segments
HAIRPIN_STEM = 8         # minimum stem length (bp) for a hairpin
HAIRPIN_LOOP = (3, 10)   # allowed loop size (min, max)
FIVE_PRIME_LEN = 60      # 5' region checked separately (affects translation start)
FIVE_PRIME_STEM = 6      # stem length used in the 5' region

STOP_CODONS = {"TAA", "TAG", "TGA"}

# Unwanted sequence motifs (sense strand only)
UNWANTED_MOTIFS = {
    "PolyA signal (AATAAA/ATTAAA)": "A[AT]TAAA",
    "AU-rich element (ATTTA)": "ATTTA",
    "TATA box (TATAAA)": "TATAAA",
    "Pol III terminator (TTTTT)": "TTTTT",
    "G-quadruplex (GGGxGGGxGGGxGGG)": "GGG.{1,7}GGG.{1,7}GGG.{1,7}GGG",
}

# Restriction sites commonly used in cloning (edit to match the company's list)
RESTRICTION_SITES = {
    "BsaI": "GGTCTC", "BsmBI": "CGTCTC", "BbsI": "GAAGAC", "SapI": "GCTCTTC",
    "EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT",
    "NdeI": "CATATG", "XhoI": "CTCGAG", "NotI": "GCGGCCGC",
    "XbaI": "TCTAGA", "SpeI": "ACTAGT", "PstI": "CTGCAG", "KpnI": "GGTACC",
}

# Cryptic splice site patterns (rough consensus screen, NOT a predictor)
SPLICE_DONOR = "[AC]AGGT[AG]AGT"             # MAG|GTRAGT
SPLICE_ACCEPTOR = "[CT]{8,}[ACGT]{0,2}[CT]AGG"  # polypyrimidine tract + YAG|G

# A sequence is flagged for a metric when its value is ABOVE this limit.
# (None = informational only, never flagged.)
FLAG_IF_ABOVE = {
    "Windows outside GC range (%)": 0,
    "Rare codons <10% (count)": 0,
    "Homopolymer runs >=6 (count)": 0,
    "Tandem repeats (count)": 0,
    "Repeated k-mers (count)": 0,
    "Unwanted motifs (total)": 0,
    "Restriction sites (total)": 0,
    "Cryptic donor sites": 0,
    "Cryptic acceptor sites": 0,
    "Hairpins (count)": 0,
    "5' hairpins (count)": 0,
}


# ----------------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------------
def clean_dna(seq):
    """Return an upper-case DNA string without spaces ('' if missing)."""
    return "" if pd.isna(seq) else "".join(str(seq).split()).upper()


def count_matches(pattern, dna):
    """Count (possibly overlapping) matches of a regex pattern in dna."""
    return len(re.findall(f"(?={pattern})", dna))


def reverse_complement(dna):
    return dna.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def get_codons(dna):
    """Split DNA into codons, dropping a final stop codon if present."""
    codons = [dna[i:i + 3] for i in range(0, len(dna) - len(dna) % 3, 3)]
    if codons and codons[-1] in STOP_CODONS:
        codons = codons[:-1]
    return codons


def build_codon_tables(species="h_sapiens"):
    """Return (codon -> amino acid, codon -> frequency within its amino acid)."""
    table = get_codons_table(species)
    codon_to_aa, codon_freq = {}, {}
    for aa, codons in table.items():
        if aa == "*":
            continue
        for codon, freq in codons.items():
            codon_to_aa[codon] = aa
            codon_freq[codon] = freq
    return codon_to_aa, codon_freq


# ----------------------------------------------------------------------------
# Analysis functions (each returns a small dict of results)
# ----------------------------------------------------------------------------
def analyze_local_gc(dna):
    """GC content in sliding windows, plus GC at the 3rd codon position."""
    gc = [1 if b in "GC" else 0 for b in dna]
    n_windows = len(dna) - GC_WINDOW + 1

    if n_windows < 1:
        return {"GC3 (%)": None, "Min window GC (%)": None,
                "Max window GC (%)": None, "Windows outside GC range (%)": None}

    # Running sum so each window costs O(1)
    window_gc = []
    current = sum(gc[:GC_WINDOW])
    window_gc.append(current / GC_WINDOW * 100)
    for i in range(GC_WINDOW, len(dna)):
        current += gc[i] - gc[i - GC_WINDOW]
        window_gc.append(current / GC_WINDOW * 100)

    outside = sum(1 for w in window_gc if w < GC_MIN or w > GC_MAX)
    third = dna[2::3]
    gc3 = sum(1 for b in third if b in "GC") / len(third) * 100

    return {
        "GC3 (%)": round(gc3, 2),
        "Min window GC (%)": round(min(window_gc), 2),
        "Max window GC (%)": round(max(window_gc), 2),
        "Windows outside GC range (%)": round(outside / len(window_gc) * 100, 2),
    }


def analyze_codon_usage(dna, codon_to_aa, codon_freq):
    """CAI, rare codon counts and deviation from the human codon usage table."""
    codons = [c for c in get_codons(dna) if c in codon_to_aa]
    if not codons:
        return {"CAI": None, "Rare codons <10% (count)": None,
                "Rare codons <20% (count)": None, "Codon usage deviation": None}

    # Relative adaptiveness w = freq / freq of the best codon for that amino acid
    best = {}
    for codon, freq in codon_freq.items():
        aa = codon_to_aa[codon]
        best[aa] = max(best.get(aa, 0), freq)

    # CAI = geometric mean of w (Met and Trp have one codon, so skip them)
    logs = [math.log(max(codon_freq[c] / best[codon_to_aa[c]], 0.01))
            for c in codons if codon_to_aa[c] not in "MW"]
    cai = math.exp(sum(logs) / len(logs)) if logs else None

    rare10 = sum(1 for c in codons if codon_freq[c] < 0.10)
    rare20 = sum(1 for c in codons if codon_freq[c] < 0.20)

    # Deviation: how far the codon choice per amino acid is from the human table
    counts = Counter(codons)
    aa_totals = Counter(codon_to_aa[c] for c in codons)
    deviation = 0.0
    for codon, freq in codon_freq.items():
        aa = codon_to_aa[codon]
        if aa_totals[aa] == 0:
            continue
        observed = counts[codon] / aa_totals[aa]
        deviation += abs(observed - freq) * aa_totals[aa]
    deviation /= len(codons)

    return {
        "CAI": round(cai, 3) if cai else None,
        "Rare codons <10% (count)": rare10,
        "Rare codons <20% (count)": rare20,
        "Codon usage deviation": round(deviation, 3),
    }


def analyze_homopolymers(dna):
    """Count runs of 6+ identical bases (A/T runs and G/C runs separately)."""
    runs = [m.group() for m in re.finditer(r"(.)\1{5,}", dna)]
    return {
        "Homopolymer runs >=6 (count)": len(runs),
        "Longest homopolymer (bp)": max((len(r) for r in runs), default=0),
    }


def analyze_repeats(dna):
    """Tandem repeats (unit 2-6 bp, 4+ copies) and repeated k-mers."""
    tandem = [m.group() for m in re.finditer(r"(.{2,6}?)\1{3,}", dna)]

    kmers = Counter(dna[i:i + KMER_SIZE] for i in range(len(dna) - KMER_SIZE + 1))
    repeated = sum(1 for count in kmers.values() if count > 1)

    return {
        "Tandem repeats (count)": len(tandem),
        "Longest tandem repeat (bp)": max((len(t) for t in tandem), default=0),
        "Repeated k-mers (count)": repeated,
    }


def analyze_motifs(dna):
    """Unwanted motifs (sense strand) and CpG dinucleotides."""
    result = {}
    total = 0
    for name, pattern in UNWANTED_MOTIFS.items():
        n = count_matches(pattern, dna)
        result[name] = n
        total += n
    result["Unwanted motifs (total)"] = total
    result["CpG count"] = count_matches("CG", dna)  # informational
    return result


def analyze_restriction_sites(dna):
    """Count restriction sites on both strands."""
    result = {}
    total = 0
    for enzyme, site in RESTRICTION_SITES.items():
        sites = {site, reverse_complement(site)}  # palindromes collapse to 1
        n = sum(count_matches(s, dna) for s in sites)
        result[f"Site {enzyme}"] = n
        total += n
    result["Restriction sites (total)"] = total
    return result


def analyze_splice_sites(dna):
    """Rough count of cryptic donor / acceptor splice-like motifs."""
    return {
        "Cryptic donor sites": count_matches(SPLICE_DONOR, dna),
        "Cryptic acceptor sites": count_matches(SPLICE_ACCEPTOR, dna),
    }


def count_hairpins(dna, stem, loop_min, loop_max):
    """Count simple hairpins: a stem-bp segment followed (after a short loop)
    by its reverse complement. Only a proxy for real RNA folding."""
    count = 0
    for i in range(len(dna) - 2 * stem - loop_min + 1):
        arm = dna[i:i + stem]
        if set(arm) <= {"A", "T"}:      # AT-only stems are weak, skip them
            continue
        target = reverse_complement(arm)
        for loop in range(loop_min, loop_max + 1):
            start = i + stem + loop
            if dna[start:start + stem] == target:
                count += 1
                break
    return count


def analyze_hairpins(dna):
    """Hairpins in the whole sequence and in the 5' region."""
    return {
        "Hairpins (count)": count_hairpins(dna, HAIRPIN_STEM, *HAIRPIN_LOOP),
        "5' hairpins (count)": count_hairpins(
            dna[:FIVE_PRIME_LEN], FIVE_PRIME_STEM, *HAIRPIN_LOOP),
        "5' GC (first 60 bp, %)": round(
            sum(1 for b in dna[:FIVE_PRIME_LEN] if b in "GC")
            / max(len(dna[:FIVE_PRIME_LEN]), 1) * 100, 2),
    }


def analyze_one_sequence(dna, codon_to_aa, codon_freq):
    """Run every check on one DNA sequence and merge the results."""
    results = {}
    results.update(analyze_local_gc(dna))
    results.update(analyze_codon_usage(dna, codon_to_aa, codon_freq))
    results.update(analyze_homopolymers(dna))
    results.update(analyze_repeats(dna))
    results.update(analyze_motifs(dna))
    results.update(analyze_restriction_sites(dna))
    results.update(analyze_splice_sites(dna))
    results.update(analyze_hairpins(dna))
    return results


# ----------------------------------------------------------------------------
# Flags + summary
# ----------------------------------------------------------------------------
def make_flags(row):
    """List the metrics for which this sequence is above its limit."""
    flagged = [m for m, limit in FLAG_IF_ABOVE.items()
               if limit is not None and pd.notna(row.get(m)) and row[m] > limit]
    return "; ".join(flagged) if flagged else "None"


def make_summary(report):
    """One row per metric: min / mean / max and how many sequences are flagged."""
    rows = []
    for metric, limit in FLAG_IF_ABOVE.items():
        values = report[metric].dropna()
        rows.append({
            "Metric": metric,
            "Flag if above": limit,
            "Min": values.min(), "Mean": round(values.mean(), 2),
            "Max": values.max(),
            "Sequences flagged": int((values > limit).sum()),
            "Total sequences": len(values),
        })
    # Informational metrics (no flagging)
    for metric in ["CAI", "GC3 (%)", "Min window GC (%)", "Max window GC (%)",
                   "Codon usage deviation", "CpG count"]:
        values = report[metric].dropna()
        rows.append({"Metric": metric, "Flag if above": "info only",
                     "Min": values.min(), "Mean": round(values.mean(), 3),
                     "Max": values.max(), "Sequences flagged": "",
                     "Total sequences": len(values)})
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Part 2 sequence analysis.")
    parser.add_argument("--input", required=True, help="Excel with generated DNA")
    parser.add_argument("--output", required=True, help="Output Excel report")
    args = parser.parse_args()

    input_path = Path(args.input).resolve()
    output_path = Path(args.output).resolve()
    if input_path == output_path:
        raise ValueError("Input and output files must be different.")

    dna_column = "Generated DNA Sequence"
    host_column = "Codon optimization- Host"
    name_column = "Name"

    data = pd.read_excel(input_path)
    if dna_column not in data.columns:
        raise ValueError(f"Missing column: {dna_column}")

    # Keep only Homo sapiens (same rule as part 1)
    if host_column in data.columns:
        hosts = data[host_column].astype("string").str.replace(r"\s+", "", regex=True)
        data = data[hosts.str.casefold() == "homosapiens"].copy()

    codon_to_aa, codon_freq = build_codon_tables("h_sapiens")

    rows = [analyze_one_sequence(clean_dna(s), codon_to_aa, codon_freq)
            for s in data[dna_column]]
    results = pd.DataFrame(rows, index=data.index)

    # Keep the report light: name (if present) + results only
    id_cols = [c for c in (name_column,) if c in data.columns]
    report = pd.concat([data[id_cols], results], axis=1)
    report["Flags"] = report.apply(make_flags, axis=1)

    summary = make_summary(report)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output_path) as writer:
        report.to_excel(writer, sheet_name="Per sequence", index=False)
        summary.to_excel(writer, sheet_name="Summary", index=False)

    print(f"Sequences analysed: {len(report)}")
    print(summary[["Metric", "Sequences flagged", "Total sequences"]].to_string(index=False))
    print(f"Report saved to: {output_path}")


if __name__ == "__main__":
    main()