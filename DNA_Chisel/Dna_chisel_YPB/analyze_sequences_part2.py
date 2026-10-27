"""
Reads human_codon_sequences.xlsx (the file made by main.py) and measures properties NOT already covered by
analyze_sequences.py (that script already reports length, base counts, total GC, longest homopolymer, ATG/stop
and translation check).

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


NAME_COL = "Name"
GENEART_COL = "DNA optimized sequence"      # reference
GENERATED_COL = "Generated DNA Sequence"    # our main.py output

GC_WINDOWS = (50, 100)        # window sizes (bp) for local GC
RARE_GRID = (0.05, 0.10, 0.15, 0.20)   # codon-frequency thresholds to inspect
KMER_SIZES = range(6, 31)     # k values tried when looking for repeats
TANDEM_UNIT = (2, 6)          # repeat unit length range (bp)
STEM_LOOP = (3, 10)           # hairpin loop size range (nt)
STEM_MIN_PAIRS = 4            # shortest stem considered a hairpin
FIVE_PRIME_LEN = 50           # 5' region checked separately for hairpins
MIN_EXPECTED_USES = 5         # only list "unused" codons if >=5 uses were expected

# Sense-strand motifs. CpG is informational (is GeneArt avoiding it or not?)
MOTIFS = {
    "polyA signal": "A[AT]TAAA",
    "AU-rich ATTTA": "ATTTA",
    "TATA box": "TATAAA",
    "G-quadruplex": "GGG.{1,7}GGG.{1,7}GGG.{1,7}GGG",
    "CpG": "CG",
}

# Edit to the enzymes that matter for your cloning vector / assembly method.
RESTRICTION_SITES = {
    "BsaI": "GGTCTC", "BsmBI": "CGTCTC", "BbsI": "GAAGAC", "SapI": "GCTCTTC",
    "EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT", "NdeI": "CATATG",
    "XhoI": "CTCGAG", "NotI": "GCGGCCGC", "XbaI": "TCTAGA", "SpeI": "ACTAGT",
    "PstI": "CTGCAG", "KpnI": "GGTACC", "NcoI": "CCATGG", "SalI": "GTCGAC",
    "NheI": "GCTAGC", "EcoRV": "GATATC", "AscI": "GGCGCGCC", "PacI": "TTAATTAA",
}

# Rough consensus screens (NOT a splice predictor, see notes in the explanation)
SPLICE_DONOR = "AGGT[AG]AG"            # exon end ...AG|GTRAG...
SPLICE_ACCEPTOR = "[CT]{8,}[ACGT]CAG"  # polypyrimidine tract + NCAG

PAIR = {"A": "T", "T": "A", "G": "C", "C": "G"}


# ----------------------------------------------------------------------------
# Small helpers
# ----------------------------------------------------------------------------
def clean_dna(seq):
    """Upper-case DNA without spaces/line breaks."""
    return "".join(str(seq).split()).upper()


def reverse_complement(dna):
    return dna.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def count_regex(pattern, dna):
    """Count (overlapping) matches of a regex."""
    return len(re.findall(f"(?={pattern})", dna))


def count_site(site, dna):
    """Count a restriction site on both strands (palindromes counted once)."""
    return sum(count_regex(s, dna) for s in {site, reverse_complement(site)})


def split_codons(dna):
    """Split into codons; a final stop codon (if any) is dropped."""
    codons = [dna[i:i + 3] for i in range(0, len(dna) - len(dna) % 3, 3)]
    if codons and codons[-1] in {"TAA", "TAG", "TGA"}:
        codons = codons[:-1]
    return codons


def load_codon_tables(species="h_sapiens"):
    """Human codon table in the same format DNA Chisel uses.
    'freq' = frequency of a codon among the codons of its own amino acid."""
    table = get_codons_table(species)
    aa_of, freq, best = {}, {}, {}
    for aa, codons in table.items():
        if aa == "*":
            continue
        best[aa] = max(codons, key=codons.get)
        for codon, f in codons.items():
            aa_of[codon] = aa
            freq[codon] = f
    return {"aa": aa_of, "freq": freq, "best": best,
            "n_syn": Counter(aa_of.values())}


# ----------------------------------------------------------------------------
# Property 1: local GC and GC3
# ----------------------------------------------------------------------------
def window_gc(dna, window):
    """GC% of every sliding window (step 1), using a running sum."""
    is_gc = [1 if b in "GC" else 0 for b in dna]
    current = sum(is_gc[:window])
    values = [current]
    for i in range(window, len(dna)):
        current += is_gc[i] - is_gc[i - window]
        values.append(current)
    return [100 * v / window for v in values]


def gc_stats(dna):
    """Lowest/highest local GC for each window size, and GC at codon position 3."""
    result = {}
    for w in GC_WINDOWS:
        values = window_gc(dna, w)
        result[f"GC min {w}bp (%)"] = round(min(values), 1)
        result[f"GC max {w}bp (%)"] = round(max(values), 1)
    third = dna[2::3]
    result["GC3 (%)"] = round(100 * sum(b in "GC" for b in third) / len(third), 1)
    return result


# ----------------------------------------------------------------------------
# Property 2: codon usage and rare codons
# ----------------------------------------------------------------------------
def codon_usage_stats(codons, tables):
    """CAI, share of best codons, and distance from the human codon table."""
    aa_of, freq, best = tables["aa"], tables["freq"], tables["best"]
    multi = [c for c in codons if tables["n_syn"][aa_of[c]] > 1]  # skip Met/Trp

    # CAI: geometric mean of (codon freq / freq of best codon of that aa)
    logs = [math.log(max(freq[c] / freq[best[aa_of[c]]], 0.001)) for c in multi]
    cai = math.exp(sum(logs) / len(logs))

    best_share = 100 * sum(c == best[aa_of[c]] for c in multi) / len(multi)

    # Distance: % of codons that would have to change to match the table exactly
    counts = Counter(codons)
    aa_counts = Counter(aa_of[c] for c in codons)
    diff = sum(abs(counts[c] / aa_counts[aa] - freq[c]) * aa_counts[aa]
               for c, aa in aa_of.items() if aa_counts[aa])
    distance = 100 * diff / (2 * len(codons))

    return {"CAI": round(cai, 3),
            "Best-codon share (%)": round(best_share, 1),
            "Distance from human table (%)": round(distance, 1)}


def rare_codon_stats(codons, tables):
    """Lowest codon frequency used, and % of codons under each threshold."""
    freq = tables["freq"]
    result = {"Lowest codon freq used": round(min(freq[c] for c in codons), 3)}
    for t in RARE_GRID:
        pct = 100 * sum(freq[c] < t for c in codons) / len(codons)
        result[f"Codons <{t:.2f} (%)"] = round(pct, 2)
    return result


def find_unused_codons(geneart_seqs, tables):
    """Codons GeneArt never used although the table says it should have (>=5 expected uses).
    Helps to see whether GeneArt avoids more than just 'rare' codons."""
    codons = [c for s in geneart_seqs for c in split_codons(s)]
    used = Counter(codons)
    aa_totals = Counter(tables["aa"][c] for c in codons)
    rows = []
    for codon, aa in tables["aa"].items():
        expected = aa_totals[aa] * tables["freq"][codon]
        if used[codon] == 0 and expected >= MIN_EXPECTED_USES:
            rows.append({"Codon": codon, "Amino acid": aa,
                         "Human table freq": round(tables["freq"][codon], 3),
                         "Expected uses": round(expected, 1)})
    rows.sort(key=lambda r: -r["Human table freq"])
    return pd.DataFrame(rows, columns=["Codon", "Amino acid",
                                       "Human table freq", "Expected uses"])


# ----------------------------------------------------------------------------
# Property 3: homopolymers (limits per base are learned from GeneArt)
# ----------------------------------------------------------------------------
def longest_run_per_base(dna):
    longest = {b: 0 for b in "ACGT"}
    for m in re.finditer(r"(.)\1*", dna):
        longest[m.group(1)] = max(longest[m.group(1)], len(m.group()))
    return longest


def count_long_runs(dna, run_limits):
    """Number of runs longer than the GeneArt maximum for that base."""
    return sum(1 for m in re.finditer(r"(.)\1*", dna)
               if len(m.group()) > run_limits[m.group(1)])


# ----------------------------------------------------------------------------
# Property 4: tandem repeats and repeated k-mers
# ----------------------------------------------------------------------------
def tandem_tracts(dna):
    """Lengths of tandem repeats: a 2-6 bp unit repeated at least 3 times.
    Single-base units are skipped (those are homopolymers)."""
    lo, hi = TANDEM_UNIT
    pattern = rf"(.{{{lo},{hi}}}?)\1{{2,}}"
    return [len(m.group()) for m in re.finditer(pattern, dna)
            if len(set(m.group(1))) > 1]


def count_repeated_kmers(dna, k):
    """Number of k-mers that occur more than once (either strand)."""
    seen = Counter()
    for i in range(len(dna) - k + 1):
        kmer = dna[i:i + k]
        seen[min(kmer, reverse_complement(kmer))] += 1
    return sum(1 for n in seen.values() if n > 1)


def smallest_unique_k(seqs):
    """Smallest k for which none of the sequences has a repeated k-mer."""
    for k in KMER_SIZES:
        if all(count_repeated_kmers(s, k) == 0 for s in seqs):
            return k
    return None


# ----------------------------------------------------------------------------
# Property 5: hairpins (simple stem-loop finder, hydrogen-bond score)
# ----------------------------------------------------------------------------
def find_stem_scores(dna):
    """Score every perfect stem (>= STEM_MIN_PAIRS) around a 3-10 nt loop.
    Score = 3 per G-C pair + 2 per A-T pair (number of hydrogen bonds)."""
    n = len(dna)
    scores = []
    for loop_start in range(1, n):
        for loop in range(STEM_LOOP[0], STEM_LOOP[1] + 1):
            left, right = loop_start - 1, loop_start + loop
            pairs = score = 0
            while left >= 0 and right < n and PAIR[dna[left]] == dna[right]:
                score += 3 if dna[left] in "GC" else 2
                pairs += 1
                left -= 1
                right += 1
            if pairs >= STEM_MIN_PAIRS:
                scores.append(score)
    return scores


# ----------------------------------------------------------------------------
# Limits learned from GeneArt
# ----------------------------------------------------------------------------
def derive_limits(geneart_seqs):
    """Read cutoffs from the GeneArt sequences."""
    runs = [longest_run_per_base(s) for s in geneart_seqs]
    return {
        "runs": {b: max(r[b] for r in runs) for b in "ACGT"},
        "tandem": max(max(tandem_tracts(s), default=0) for s in geneart_seqs),
        "k": smallest_unique_k(geneart_seqs),
        "stem": max(max(find_stem_scores(s), default=0) for s in geneart_seqs),
        "stem_5p": max(max(find_stem_scores(s[:FIVE_PRIME_LEN]), default=0)
                       for s in geneart_seqs),
    }


# ----------------------------------------------------------------------------
# Analyse one sequence
# ----------------------------------------------------------------------------
def analyze_sequence(dna, tables, limits):
    """All part-2 properties for one DNA sequence."""
    codons = split_codons(dna)
    result = {}
    result.update(gc_stats(dna))
    result.update(codon_usage_stats(codons, tables))
    result.update(rare_codon_stats(codons, tables))

    result["Runs above GeneArt limit (count)"] = count_long_runs(dna, limits["runs"])

    tracts = tandem_tracts(dna)
    result["Longest tandem repeat (bp)"] = max(tracts, default=0)
    result["Tandem repeats above GeneArt max (count)"] = sum(
        t > limits["tandem"] for t in tracts)
    k = limits["k"]
    result["Repeated k-mers (count)"] = count_repeated_kmers(dna, k) if k else None

    for name, pattern in MOTIFS.items():
        result[f"Motif {name} (count)"] = count_regex(pattern, dna)

    sites = {name: count_site(site, dna) for name, site in RESTRICTION_SITES.items()}
    for name, n in sites.items():
        result[f"Site {name} (count)"] = n
    result["Restriction sites total (count)"] = sum(sites.values())

    result["Cryptic donor-like motifs (count)"] = count_regex(SPLICE_DONOR, dna)
    result["Cryptic acceptor-like motifs (count)"] = count_regex(SPLICE_ACCEPTOR, dna)

    stems = find_stem_scores(dna)
    stems_5p = find_stem_scores(dna[:FIVE_PRIME_LEN])
    result["Strongest stem score"] = max(stems, default=0)
    result["Stems above GeneArt max (count)"] = sum(s > limits["stem"] for s in stems)
    result["5' strongest stem score"] = max(stems_5p, default=0)
    result["5' stems above GeneArt max (count)"] = sum(
        s > limits["stem_5p"] for s in stems_5p)
    return result


# ----------------------------------------------------------------------------
# Compare Generated vs GeneArt, gene by gene
# ----------------------------------------------------------------------------
def worse_direction(metric):
    """Which direction means 'worse than GeneArt' for a metric (None = compare only)."""
    if metric.startswith("GC min") or metric.startswith("Lowest codon freq"):
        return "lower"
    if (metric.startswith("GC max") or metric.endswith("(count)")
            or metric.startswith("Codons <") or "stem score" in metric
            or metric.startswith("Longest tandem")):
        return "higher"
    return None


def compare_sources(report):
    """One row per metric: GeneArt vs Generated, and in how many genes
    Generated is worse than GeneArt for the same protein."""
    geneart = report[report["Source"] == "GeneArt"].set_index(NAME_COL)
    generated = report[report["Source"] == "Generated"].set_index(NAME_COL)
    generated = generated.loc[geneart.index]

    rows = []
    for metric in geneart.columns.drop("Source"):
        direction = worse_direction(metric)
        if direction == "lower":
            worse = int((generated[metric] < geneart[metric]).sum())
        elif direction == "higher":
            worse = int((generated[metric] > geneart[metric]).sum())
        else:
            worse = ""
        rows.append({
            "Metric": metric,
            "Worse if Generated is": direction or "compare only",
            "GeneArt min": geneart[metric].min(),
            "GeneArt mean": round(geneart[metric].mean(), 3),
            "GeneArt max": geneart[metric].max(),
            "Generated min": generated[metric].min(),
            "Generated mean": round(generated[metric].mean(), 3),
            "Generated max": generated[metric].max(),
            f"Genes where Generated is worse (of {len(geneart)})": worse,
        })
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# DNA Chisel hints (plain-text summary of what the data suggests)
# ----------------------------------------------------------------------------
def names_with_zero_in_geneart(report, prefix):
    """Metric names (without prefix/suffix) that GeneArt never contains."""
    ga = report[report["Source"] == "GeneArt"]
    cols = [c for c in ga.columns if c.startswith(prefix) and c.endswith("(count)")]
    zero = [c[len(prefix):-len(" (count)")] for c in cols if ga[c].sum() == 0]
    return [z for z in zero if z != "total"]


def make_hints(report, limits, unused_codons):
    ga = report[report["Source"] == "GeneArt"]
    gen = report[report["Source"] == "Generated"]
    hints = []

    for w in GC_WINDOWS:
        lo = math.floor(ga[f"GC min {w}bp (%)"].min())
        hi = math.ceil(ga[f"GC max {w}bp (%)"].max())
        gen_lo = gen[f"GC min {w}bp (%)"].min()
        gen_hi = gen[f"GC max {w}bp (%)"].max()
        hints.append(f"GC window {w} bp: GeneArt stays within {lo}-{hi}% "
                     f"(Generated {gen_lo}-{gen_hi}%)  ->  "
                     f"EnforceGCContent(mini={lo / 100:.2f}, maxi={hi / 100:.2f}, window={w})")

    runs = limits["runs"]
    patterns = ", ".join(f'"{runs[b] + 1}x{b}"' for b in "ACGT")
    hints.append(f"Homopolymers: GeneArt max run per base = {runs}  ->  AvoidPattern({patterns})")

    if limits["k"]:
        hints.append(f"Repeats: GeneArt has no repeated {limits['k']}-mer (either strand)  ->  "
                     f"UniquifyAllKmers({limits['k']})")
    hints.append(f"Tandem repeats: longest tract in GeneArt = {limits['tandem']} bp "
                 f"(Generated max = {gen['Longest tandem repeat (bp)'].max()} bp)")

    absent_motifs = names_with_zero_in_geneart(report, "Motif ")
    hints.append(f"Motifs never present in GeneArt: {absent_motifs or 'none'}  ->  AvoidPattern")
    absent_sites = names_with_zero_in_geneart(report, "Site ")
    hints.append(f"Restriction sites never present in GeneArt: {absent_sites or 'none'}  ->  "
                 f"AvoidPattern('<enzyme>_site')  (confirm with your cloning plan)")

    hints.append(f"Rare codons: lowest codon frequency GeneArt ever used = "
                 f"{ga['Lowest codon freq used'].min()}  ->  AvoidRareCodons(min_frequency ~ that value)")
    if not unused_codons.empty:
        top = ", ".join(f"{row['Codon']}({row['Human table freq']})"
                        for _, row in unused_codons.head(8).iterrows())
        hints.append(f"Codons GeneArt never used (table freq): {top}")

    hints.append(f"Codon strategy: best-codon share GeneArt {ga['Best-codon share (%)'].mean():.0f}% "
                 f"vs Generated {gen['Best-codon share (%)'].mean():.0f}%; "
                 f"distance from human table GeneArt {ga['Distance from human table (%)'].mean():.0f}% "
                 f"vs Generated {gen['Distance from human table (%)'].mean():.0f}%  "
                 f"->  decides use_best_codon vs match_codon_usage")

    gen_stem = gen["Strongest stem score"].max()
    gen_stem_5p = gen["5' strongest stem score"].max()
    hints.append(f"Hairpins: strongest stem score GeneArt max = {limits['stem']} "
                 f"(5' region: {limits['stem_5p']}); Generated max = {gen_stem} "
                 f"(5' region: {gen_stem_5p})")
    return hints


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Part 2 sequence analysis.")
    parser.add_argument("--input", required=True, help="human_codon_sequences.xlsx")
    parser.add_argument("--output", required=True, help="Output Excel report")
    args = parser.parse_args()

    input_path = Path(args.input).resolve()
    output_path = Path(args.output).resolve()
    if input_path == output_path:
        raise ValueError("Input and output files must be different.")

    data = pd.read_excel(input_path)
    for col in (NAME_COL, GENEART_COL, GENERATED_COL):
        if col not in data.columns:
            raise ValueError(f"Missing column: {col}. Found: {data.columns.tolist()}")

    # Clean and validate (part 1 already checks translation; here we only need ACGT)
    for col in (GENEART_COL, GENERATED_COL):
        data[col] = data[col].map(clean_dna)
        bad = data[~data[col].str.fullmatch(r"[ACGT]+")]
        if not bad.empty:
            raise ValueError(f"Invalid/empty DNA in '{col}' for: {bad[NAME_COL].tolist()}")

    tables = load_codon_tables("h_sapiens")
    limits = derive_limits(data[GENEART_COL].tolist())

    rows = []
    for _, r in data.iterrows():
        for source, col in (("GeneArt", GENEART_COL), ("Generated", GENERATED_COL)):
            row = {NAME_COL: r[NAME_COL], "Source": source}
            row.update(analyze_sequence(r[col], tables, limits))
            rows.append(row)
    report = pd.DataFrame(rows)

    comparison = compare_sources(report)
    unused_codons = find_unused_codons(data[GENEART_COL].tolist(), tables)
    hints = make_hints(report, limits, unused_codons)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output_path) as writer:
        report.to_excel(writer, sheet_name="Per sequence", index=False)
        comparison.to_excel(writer, sheet_name="GeneArt vs Generated", index=False)
        pd.DataFrame({"What the data suggests": hints}).to_excel(
            writer, sheet_name="DNA Chisel hints", index=False)
        unused_codons.to_excel(writer, sheet_name="Unused codons", index=False)

    print("\n".join(hints))
    print(f"\nReport saved to: {output_path}")


if __name__ == "__main__":
    main()