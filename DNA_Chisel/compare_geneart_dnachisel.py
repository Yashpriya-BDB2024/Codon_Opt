# -*- coding: utf-8 -*-
"""
DNA Chisel vs GeneArt - Full Comparison Script with Publication-Quality Plots

Run with:
    source ~/codon_optimization/dnachisel/venv/bin/activate
    python3 compare_geneart_dnachisel.py
"""

import pandas as pd
import sys
import os
import math
import numpy as np

from dnachisel import (
    DnaOptimizationProblem,
    CodonOptimize,
    EnforceGCContent,
    AvoidPattern,
    AvoidHairpins,
    EnforceTranslation,
    UniquifyAllKmers,
)

INPUT_FILE   = "codon_optimization_trials.xlsx"
OUTPUT_FILE  = "comparison_results.xlsx"
PLOTS_FOLDER = "comparison_plots"

HOST_MAP = {
    "homo sapiens":     "h_sapiens_9606",
    "homosapiens":      "h_sapiens_9606",
    "human":            "h_sapiens_9606",
    "e.coli":           "e_coli_316407",
    "e. coli":          "e_coli_316407",
    "ecoli":            "e_coli_316407",
    "escherichia coli": "e_coli_316407",
}

RESTRICTION_SITES = {
    "BsaI": "GGTCTC", "BsmBI": "CGTCTC", "EcoRI": "GAATTC",
    "NotI": "GCGGCCGC", "XhoI": "CTCGAG", "BamHI": "GGATCC",
    "NheI": "GCTAGC", "HindIII": "AAGCTT", "XbaI": "TCTAGA", "SalI": "GTCGAC",
}

SPLICE_SITES = {
    "5prime_donor": "GTAAG", "3prime_acceptor": "TTTTTAG", "branch_point": "TACTAAC",
}

HUMAN_CODON_TABLE = {
    'A': {'GCT': 0.27, 'GCC': 0.40, 'GCA': 0.23, 'GCG': 0.11},
    'R': {'CGT': 0.08, 'CGC': 0.19, 'CGA': 0.11, 'CGG': 0.21, 'AGA': 0.20, 'AGG': 0.20},
    'N': {'AAT': 0.46, 'AAC': 0.54},
    'D': {'GAT': 0.46, 'GAC': 0.54},
    'C': {'TGT': 0.45, 'TGC': 0.55},
    'Q': {'CAA': 0.27, 'CAG': 0.73},
    'E': {'GAA': 0.42, 'GAG': 0.58},
    'G': {'GGT': 0.16, 'GGC': 0.34, 'GGA': 0.25, 'GGG': 0.25},
    'H': {'CAT': 0.41, 'CAC': 0.59},
    'I': {'ATT': 0.36, 'ATC': 0.48, 'ATA': 0.16},
    'L': {'TTA': 0.07, 'TTG': 0.13, 'CTT': 0.13, 'CTC': 0.20, 'CTA': 0.07, 'CTG': 0.41},
    'K': {'AAA': 0.42, 'AAG': 0.58},
    'M': {'ATG': 1.00},
    'F': {'TTT': 0.45, 'TTC': 0.55},
    'P': {'CCT': 0.28, 'CCC': 0.33, 'CCA': 0.27, 'CCG': 0.11},
    'S': {'TCT': 0.15, 'TCC': 0.22, 'TCA': 0.15, 'TCG': 0.06, 'AGT': 0.15, 'AGC': 0.24},
    'T': {'ACT': 0.25, 'ACC': 0.36, 'ACA': 0.28, 'ACG': 0.11},
    'W': {'TGG': 1.00},
    'Y': {'TAT': 0.43, 'TAC': 0.57},
    'V': {'GTT': 0.18, 'GTC': 0.24, 'GTA': 0.11, 'GTG': 0.47},
    '*': {'TAA': 0.30, 'TAG': 0.24, 'TGA': 0.47},
}

ECOLI_CODON_TABLE = {
    'A': {'GCT': 0.17, 'GCC': 0.26, 'GCA': 0.22, 'GCG': 0.35},
    'R': {'CGT': 0.36, 'CGC': 0.36, 'CGA': 0.07, 'CGG': 0.11, 'AGA': 0.07, 'AGG': 0.04},
    'N': {'AAT': 0.49, 'AAC': 0.51},
    'D': {'GAT': 0.63, 'GAC': 0.37},
    'C': {'TGT': 0.46, 'TGC': 0.54},
    'Q': {'CAA': 0.34, 'CAG': 0.66},
    'E': {'GAA': 0.68, 'GAG': 0.32},
    'G': {'GGT': 0.35, 'GGC': 0.37, 'GGA': 0.13, 'GGG': 0.15},
    'H': {'CAT': 0.57, 'CAC': 0.43},
    'I': {'ATT': 0.51, 'ATC': 0.42, 'ATA': 0.07},
    'L': {'TTA': 0.14, 'TTG': 0.13, 'CTT': 0.12, 'CTC': 0.10, 'CTA': 0.04, 'CTG': 0.47},
    'K': {'AAA': 0.74, 'AAG': 0.26},
    'M': {'ATG': 1.00},
    'F': {'TTT': 0.58, 'TTC': 0.42},
    'P': {'CCT': 0.18, 'CCC': 0.13, 'CCA': 0.20, 'CCG': 0.49},
    'S': {'TCT': 0.17, 'TCC': 0.15, 'TCA': 0.14, 'TCG': 0.14, 'AGT': 0.16, 'AGC': 0.25},
    'T': {'ACT': 0.19, 'ACC': 0.40, 'ACA': 0.17, 'ACG': 0.25},
    'W': {'TGG': 1.00},
    'Y': {'TAT': 0.59, 'TAC': 0.41},
    'V': {'GTT': 0.28, 'GTC': 0.20, 'GTA': 0.17, 'GTG': 0.35},
    '*': {'TAA': 0.61, 'TAG': 0.09, 'TGA': 0.30},
}

CODON_TABLES   = {'h_sapiens_9606': HUMAN_CODON_TABLE, 'e_coli_316407': ECOLI_CODON_TABLE}
GC_CONSTRAINTS = {'h_sapiens_9606': (0.40, 0.65), 'e_coli_316407': (0.40, 0.55)}

def get_codon_table(species_id):
    return CODON_TABLES.get(species_id, HUMAN_CODON_TABLE)

def get_gc_constraints(species_id):
    return GC_CONSTRAINTS.get(species_id, (0.40, 0.65))

def calculate_cai(seq, codon_table):
    if not seq or not codon_table:
        return None
    seq = seq.upper().replace(" ", "")
    if seq[-3:] in ["TAA", "TAG", "TGA"]:
        seq = seq[:-3]
    codon_to_aa = {c: aa for aa, codons in codon_table.items() for c in codons}
    max_freq    = {aa: max(codons.values()) for aa, codons in codon_table.items() if codons}
    log_sum, count = 0, 0
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if len(codon) < 3: continue
        aa = codon_to_aa.get(codon)
        if not aa or aa not in codon_table: continue
        freq  = codon_table[aa].get(codon, 0)
        max_f = max_freq.get(aa, 1)
        if max_f > 0 and freq > 0:
            log_sum += math.log(freq / max_f)
            count   += 1
    return round(math.exp(log_sum / count), 4) if count > 0 else None

def rare_codon_percent(seq, codon_table, threshold=0.3):
    if not seq or not codon_table: return None
    seq = seq.upper().replace(" ", "")
    if seq[-3:] in ["TAA", "TAG", "TGA"]: seq = seq[:-3]
    codon_to_aa = {c: aa for aa, codons in codon_table.items() for c in codons}
    max_freq    = {aa: max(codons.values()) for aa, codons in codon_table.items() if codons}
    rare, total = 0, 0
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if len(codon) < 3: continue
        aa = codon_to_aa.get(codon)
        if not aa or aa not in codon_table: continue
        freq  = codon_table[aa].get(codon, 0)
        max_f = max_freq.get(aa, 1)
        if max_f > 0 and (freq / max_f) < threshold: rare += 1
        total += 1
    return round(rare / total * 100, 2) if total > 0 else None

def gc_content(seq):
    if not seq or not isinstance(seq, str): return None
    seq = seq.upper().replace(" ", "")
    return round(sum(1 for b in seq if b in "GC") / len(seq) * 100, 2) if seq else None

def gc_windows(seq, window=60):
    if not seq or len(seq) < window: return None, None
    seq  = seq.upper().replace(" ", "")
    vals = [gc_content(seq[i:i+window]) for i in range(0, len(seq) - window, window)]
    return round(min(vals), 2), round(max(vals), 2)

def gc_std(seq, window=60):
    if not seq or len(seq) < window * 2: return None
    seq  = seq.upper().replace(" ", "")
    vals = [gc_content(seq[i:i+window]) for i in range(0, len(seq) - window, window)]
    if len(vals) < 2: return None
    mean = sum(vals) / len(vals)
    return round(math.sqrt(sum((v - mean)**2 for v in vals) / len(vals)), 2)

def gc3_content(seq):
    if not seq or len(seq) < 3: return None
    seq    = seq.upper().replace(" ", "")
    thirds = [seq[i+2] for i in range(0, len(seq)-2, 3)]
    return round(sum(1 for b in thirds if b in "GC") / len(thirds) * 100, 2) if thirds else None

def scan_restriction_sites(seq):
    if not seq or not isinstance(seq, str): return {k: 0 for k in RESTRICTION_SITES}, 0
    seq = seq.upper().replace(" ", "")
    complement = str.maketrans("ATGC", "TACG")
    counts, total = {}, 0
    for name, pattern in RESTRICTION_SITES.items():
        rc = pattern[::-1].translate(complement)
        fwd = seq.count(pattern)
        rev = seq.count(rc) if rc != pattern else 0
        counts[name] = fwd + rev
        total += fwd + rev
    return counts, total

def scan_splice_sites(seq):
    if not seq: return 0
    seq = seq.upper().replace(" ", "")
    return sum(seq.count(p) for p in SPLICE_SITES.values())

def cpg_count(seq):
    if not seq: return None
    return seq.upper().replace(" ", "").count("CG")

def cpg_density(seq):
    if not seq: return None
    s = seq.replace(" ", "")
    return round(cpg_count(seq) / len(s) * 100, 2) if s else None

def find_repetitions(seq, min_repeat_len=10, min_copies=2):
    if not seq or len(seq) < min_repeat_len * 2: return 0
    seq = seq.upper().replace(" ", "")
    count, seen = 0, set()
    for i in range(len(seq) - min_repeat_len):
        kmer = seq[i:i+min_repeat_len]
        if kmer in seen: continue
        seen.add(kmer)
        if seq.count(kmer) >= min_copies: count += 1
    return count

def longest_repeat(seq, min_len=8):
    if not seq or len(seq) < min_len * 2: return 0
    seq = seq.upper().replace(" ", "")
    for length in range(len(seq)//2, min_len-1, -1):
        seen = set()
        for i in range(len(seq) - length + 1):
            kmer = seq[i:i+length]
            if kmer in seen: return length
            seen.add(kmer)
    return 0

def homopolymer_count(seq, min_run=6):
    if not seq: return 0
    seq = seq.upper().replace(" ", "")
    return sum(seq.count(b * min_run) for b in "ATGC")

def count_hairpins(seq, stem_size=6, loop_min=3, window=200):
    if not seq or len(seq) < stem_size * 2 + loop_min: return 0
    seq = seq.upper().replace(" ", "")
    complement = str.maketrans("ATGC", "TACG")
    count = 0
    for i in range(len(seq) - stem_size):
        stem         = seq[i:i+stem_size]
        rc_stem      = stem[::-1].translate(complement)
        search_start = i + stem_size + loop_min
        search_end   = min(i + window, len(seq) - stem_size)
        if search_start >= search_end: continue
        if rc_stem in seq[search_start:search_end]: count += 1
    return count

def hairpin_count_5prime(seq, stem_size=6, region_size=100):
    if not seq or len(seq) < stem_size * 2: return 0
    return count_hairpins(seq[:region_size], stem_size=stem_size)

def percent_identity(seq1, seq2):
    if not seq1 or not seq2: return None
    s1 = seq1.upper().replace(" ", "")
    s2 = seq2.upper().replace(" ", "")
    if len(s1) != len(s2): return None
    return round(sum(a == b for a, b in zip(s1, s2)) / len(s1) * 100, 2)

def count_differences(seq1, seq2):
    if not seq1 or not seq2: return None
    s1 = seq1.upper().replace(" ", "")
    s2 = seq2.upper().replace(" ", "")
    if len(s1) != len(s2): return "Different lengths ({} vs {})".format(len(s1), len(s2))
    diff = sum(a != b for a, b in zip(s1, s2))
    return "{} bp ({}%)".format(diff, round(diff / len(s1) * 100, 1))

def backtranslate(protein, codon_table):
    dna = ""
    for aa in protein:
        if aa not in codon_table: continue
        codons = codon_table[aa]
        dna   += max(codons, key=lambda c: codons[c])
    return dna

def normalize_host(host_str):
    if not host_str or not isinstance(host_str, str): return None, "Unknown host"
    key     = host_str.strip().lower()
    species = HOST_MAP.get(key)
    if not species:
        for k, v in HOST_MAP.items():
            if k in key or key in k:
                species = v
                break
    return species, host_str.strip()

def get_all_metrics(seq, codon_table, label):
    if not seq: return {}
    rs_counts, rs_total = scan_restriction_sites(seq)
    gc60_min, gc60_max  = gc_windows(seq, 60)
    gc30_min, gc30_max  = gc_windows(seq, 30)
    gc_ov               = gc_content(seq)
    m = {}
    m["{} CAI".format(label)]                       = calculate_cai(seq, codon_table)
    m["{} rare codon % (w<0.3)".format(label)]      = rare_codon_percent(seq, codon_table, 0.3)
    m["{} GC% overall".format(label)]               = gc_ov
    m["{} GC% min (60bp)".format(label)]            = gc60_min
    m["{} GC% max (60bp)".format(label)]            = gc60_max
    m["{} GC% min (30bp)".format(label)]            = gc30_min
    m["{} GC% max (30bp)".format(label)]            = gc30_max
    m["{} GC std dev".format(label)]                = gc_std(seq, 60)
    m["{} GC3 content %".format(label)]             = gc3_content(seq)
    m["{} AT% overall".format(label)]               = round(100 - gc_ov, 2) if gc_ov else None
    m["{} total restriction sites".format(label)]   = rs_total
    for site, cnt in rs_counts.items():
        m["{} {}".format(label, site)]              = cnt
    m["{} cryptic splice sites".format(label)]      = scan_splice_sites(seq)
    m["{} CpG count".format(label)]                 = cpg_count(seq)
    m["{} CpG density (per 100bp)".format(label)]   = cpg_density(seq)
    m["{} repeat count (>=10bp)".format(label)]     = find_repetitions(seq, 10, 2)
    m["{} repeat count (>=20bp)".format(label)]     = find_repetitions(seq, 20, 2)
    m["{} longest repeat (bp)".format(label)]       = longest_repeat(seq)
    m["{} homopolymer runs (>=6)".format(label)]    = homopolymer_count(seq, 6)
    m["{} hairpin count (all)".format(label)]       = count_hairpins(seq, stem_size=6)
    m["{} hairpin count (5' region)".format(label)] = hairpin_count_5prime(seq)
    m["{} length (bp)".format(label)]               = len(seq.replace(" ", ""))
    return m

def run_dnachisel(protein, species_id, codon_table):
    gc_min, gc_max = get_gc_constraints(species_id)
    start_seq      = backtranslate(protein, codon_table)
    for attempt, constraints in enumerate([
        [EnforceTranslation(), EnforceGCContent(mini=gc_min, maxi=gc_max, window=60),
         AvoidPattern("BsaI_site"), AvoidPattern("BsmBI_site"),
         AvoidHairpins(stem_size=8, hairpin_window=100), UniquifyAllKmers(k=9)],
        [EnforceTranslation(), EnforceGCContent(mini=gc_min, maxi=gc_max, window=60),
         AvoidPattern("BsaI_site"), AvoidPattern("BsmBI_site"), UniquifyAllKmers(k=9)],
        [EnforceTranslation(), EnforceGCContent(mini=gc_min, maxi=gc_max, window=60),
         AvoidPattern("BsaI_site"), AvoidPattern("BsmBI_site")],
    ]):
        try:
            labels  = ["Success", "Success (no hairpin constraint)", "Success (minimal constraints)"]
            errors  = [None, "Hairpin constraint removed", "Hairpin and repeat constraints removed"]
            problem = DnaOptimizationProblem(
                sequence=start_seq, constraints=constraints,
                objectives=[CodonOptimize(species=species_id, method="use_best_codon")]
            )
            problem.max_random_iters = 10000
            problem.resolve_constraints()
            problem.optimize()
            return {"status": labels[attempt], "optimized_seq": problem.sequence,
                    "constraints_pass": problem.all_constraints_pass(),
                    "cai_score": problem.objectives_text_summary(), "error": errors[attempt]}
        except Exception:
            continue
    return {"status": "Failed", "optimized_seq": None, "constraints_pass": None,
            "cai_score": None, "error": "All constraint levels failed"}

# ── Main ──────────────────────────────────────────────────────
print("=" * 60)
print("  DNA Chisel vs GeneArt -- Full Parameter Comparison")
print("=" * 60)

if not os.path.exists(INPUT_FILE):
    print("ERROR: {} not found.".format(INPUT_FILE))
    sys.exit(1)

df = pd.read_excel(INPUT_FILE)
df.columns = df.columns.str.strip()
print("  Found {} sequences\n".format(len(df)))

NAME_COL = "Name"; PROTEIN_COL = "Amino acid sequence"
HOST_COL = "Codon optimization- Host"; GENEART_COL = "DNA optimized sequence"

missing = [c for c in [NAME_COL, PROTEIN_COL, HOST_COL, GENEART_COL] if c not in df.columns]
if missing:
    print("ERROR: Missing columns: {}".format(missing)); sys.exit(1)

results = []
for idx, row in df.iterrows():
    name        = row[NAME_COL]
    protein     = str(row[PROTEIN_COL]).strip() if pd.notna(row[PROTEIN_COL]) else ""
    host_raw    = row[HOST_COL]
    geneart_seq = str(row[GENEART_COL]).strip() if pd.notna(row[GENEART_COL]) else ""
    print("  [{}/{}] {} -- Host: {}".format(idx+1, len(df), name, host_raw))
    species_id, host_clean = normalize_host(host_raw)
    if not protein or not species_id:
        reason = "no protein" if not protein else "unknown host"
        print("    SKIP: {}".format(reason))
        results.append({"Name": name, "Host": host_clean, "Status": "Skipped -- " + reason})
        continue
    codon_table    = get_codon_table(species_id)
    gc_min, gc_max = get_gc_constraints(species_id)
    dc             = run_dnachisel(protein, species_id, codon_table)
    dc_seq         = dc["optimized_seq"] or ""
    print("    Status: {} | GC: {}-{}%".format(dc["status"], int(gc_min*100), int(gc_max*100)))
    result = {
        "Name": name, "Host": host_clean, "Species ID": species_id,
        "Protein length (aa)": len(protein), "DNA Chisel Status": dc["status"],
        "Constraints passed": dc["constraints_pass"], "CAI score (text)": dc["cai_score"],
        "Error": dc["error"],
        "% identity (GeneArt vs DNAChisel)": percent_identity(geneart_seq, dc_seq),
        "Positions differing": count_differences(geneart_seq, dc_seq),
    }
    result.update(get_all_metrics(geneart_seq, codon_table, "GeneArt"))
    result.update(get_all_metrics(dc_seq,      codon_table, "DNAChisel"))
    delta_pairs = [
        ("CAI", "GeneArt CAI", "DNAChisel CAI"),
        ("rare codon % (w<0.3)", "GeneArt rare codon % (w<0.3)", "DNAChisel rare codon % (w<0.3)"),
        ("GC% overall", "GeneArt GC% overall", "DNAChisel GC% overall"),
        ("GC% max (60bp)", "GeneArt GC% max (60bp)", "DNAChisel GC% max (60bp)"),
        ("GC std dev", "GeneArt GC std dev", "DNAChisel GC std dev"),
        ("total restriction sites", "GeneArt total restriction sites", "DNAChisel total restriction sites"),
        ("cryptic splice sites", "GeneArt cryptic splice sites", "DNAChisel cryptic splice sites"),
        ("CpG density", "GeneArt CpG density (per 100bp)", "DNAChisel CpG density (per 100bp)"),
        ("repeat count (>=10bp)", "GeneArt repeat count (>=10bp)", "DNAChisel repeat count (>=10bp)"),
        ("longest repeat", "GeneArt longest repeat (bp)", "DNAChisel longest repeat (bp)"),
        ("hairpin count (all)", "GeneArt hairpin count (all)", "DNAChisel hairpin count (all)"),
        ("hairpin 5' region", "GeneArt hairpin count (5' region)", "DNAChisel hairpin count (5' region)"),
    ]
    for label, ga_col, dc_col in delta_pairs:
        ga_val = result.get(ga_col); dc_val = result.get(dc_col)
        if ga_val is not None and dc_val is not None:
            try: result["Delta {}".format(label)] = round(float(dc_val) - float(ga_val), 4)
            except: result["Delta {}".format(label)] = None
        else: result["Delta {}".format(label)] = None
    results.append(result)
    ga_cai = result.get("GeneArt CAI"); dc_cai = result.get("DNAChisel CAI")
    if ga_cai is not None and dc_cai is not None:
        print("    CAI: GeneArt={} | DNAChisel={}".format(ga_cai, dc_cai))

# ── Write Excel ───────────────────────────────────────────────
print("\nWriting {}...".format(OUTPUT_FILE))
results_df  = pd.DataFrame(results)
summary_cols = [
    "Name", "Host", "Protein length (aa)", "DNA Chisel Status", "Constraints passed",
    "% identity (GeneArt vs DNAChisel)", "Positions differing",
    "GeneArt CAI", "DNAChisel CAI", "Delta CAI",
    "GeneArt rare codon % (w<0.3)", "DNAChisel rare codon % (w<0.3)", "Delta rare codon % (w<0.3)",
    "GeneArt GC% overall", "DNAChisel GC% overall", "Delta GC% overall",
    "GeneArt GC% max (60bp)", "DNAChisel GC% max (60bp)", "Delta GC% max (60bp)",
    "GeneArt GC% min (60bp)", "DNAChisel GC% min (60bp)",
    "GeneArt GC std dev", "DNAChisel GC std dev", "Delta GC std dev",
    "GeneArt GC3 content %", "DNAChisel GC3 content %",
    "GeneArt total restriction sites", "DNAChisel total restriction sites", "Delta total restriction sites",
    "GeneArt cryptic splice sites", "DNAChisel cryptic splice sites", "Delta cryptic splice sites",
    "GeneArt CpG density (per 100bp)", "DNAChisel CpG density (per 100bp)", "Delta CpG density",
    "GeneArt repeat count (>=10bp)", "DNAChisel repeat count (>=10bp)", "Delta repeat count (>=10bp)",
    "GeneArt longest repeat (bp)", "DNAChisel longest repeat (bp)", "Delta longest repeat",
    "GeneArt homopolymer runs (>=6)", "DNAChisel homopolymer runs (>=6)",
    "GeneArt hairpin count (all)", "DNAChisel hairpin count (all)", "Delta hairpin count (all)",
    "GeneArt hairpin count (5' region)", "DNAChisel hairpin count (5' region)", "Delta hairpin 5' region",
]
with pd.ExcelWriter(OUTPUT_FILE, engine="openpyxl") as writer:
    results_df[[c for c in summary_cols if c in results_df.columns]].to_excel(writer, sheet_name="Summary", index=False)
    results_df.to_excel(writer, sheet_name="Full Results", index=False)
    results_df[["Name","Host","GeneArt sequence","DNAChisel sequence"]].to_excel(writer, sheet_name="Sequences", index=False) if all(c in results_df.columns for c in ["GeneArt sequence","DNAChisel sequence"]) else None
    from openpyxl.styles import Font, PatternFill, Alignment
    for sn in writer.sheets:
        ws = writer.sheets[sn]
        for col in ws.columns:
            ws.column_dimensions[col[0].column_letter].width = min(max((len(str(c.value)) for c in col if c.value), default=10) + 4, 50)
        for cell in ws[1]:
            cell.font = Font(bold=True, name="Arial", size=9)
            cell.fill = PatternFill("solid", start_color="D9E1F2")
            cell.alignment = Alignment(horizontal="center", wrap_text=True)
        ws.row_dimensions[1].height = 45
print("Done! Saved to {}".format(OUTPUT_FILE))

# ══════════════════════════════════════════════════════════════
# PUBLICATION-QUALITY SCATTER PLOTS
# ══════════════════════════════════════════════════════════════
print("\nGenerating plots...")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from matplotlib import rcParams

os.makedirs(PLOTS_FOLDER, exist_ok=True)

# ── Design system ─────────────────────────────────────────────
BG_COLOR    = "#FAFAFA"
PANEL_COLOR = "#FFFFFF"
GRID_COLOR  = "#E0E4EF"
TEXT_COLOR  = "#1A1D2E"
MUTED_COLOR = "#6B7280"

COL_HUMAN   = "#1565C0"   # dark blue -- human sequences
COL_ECOLI   = "#E65100"   # dark amber -- E. coli sequences
COL_BETTER  = "#2E7D32"   # dark green -- DNA Chisel wins
COL_WORSE   = "#C62828"   # dark red -- GeneArt wins
COL_DRAW    = "#78909C"   # grey -- draw

REF_GREEN   = "#2E7D32"
REF_RED     = "#C62828"
REF_AMBER   = "#E65100"

rcParams.update({
    "font.family":       "DejaVu Sans",
    "text.color":        TEXT_COLOR,
    "axes.facecolor":    PANEL_COLOR,
    "axes.edgecolor":    GRID_COLOR,
    "axes.labelcolor":   TEXT_COLOR,
    "axes.grid":         True,
    "grid.color":        GRID_COLOR,
    "grid.linewidth":    0.6,
    "grid.alpha":        0.6,
    "xtick.color":       MUTED_COLOR,
    "ytick.color":       MUTED_COLOR,
    "figure.facecolor":  BG_COLOR,
    "savefig.facecolor": BG_COLOR,
    "savefig.edgecolor": "none",
    "legend.facecolor":  PANEL_COLOR,
    "legend.edgecolor":  GRID_COLOR,
})

def is_ecoli(row_name, df_res):
    row = df_res[df_res["Name"] == row_name]
    if row.empty: return False
    return "coli" in str(row.iloc[0].get("Species ID", "")).lower()

def make_slopegraph(df_res, ga_col, dc_col, title,
                    better="lower", ref_bands=None, ref_lines=None,
                    ref_annotation="", filename="plot.png"):
    """
    Slopegraph: each sequence is a labelled line from GeneArt to DNA Chisel.
    - Rising line  = DNA Chisel higher
    - Falling line = GeneArt higher
    - Line color   = good (green) or bad (red) depending on 'better' direction
    - Dot shape    = circle (human) or triangle (E. coli)
    - Reference band shown on Y axis
    """
    valid = df_res.dropna(subset=[ga_col, dc_col]).copy()
    valid = valid[valid[ga_col].apply(lambda x: isinstance(x, (int, float)))]
    valid = valid[valid[dc_col].apply(lambda x: isinstance(x, (int, float)))]
    if valid.empty:
        print("  Skipping {} -- no valid data".format(title))
        return

    ga_vals  = valid[ga_col].astype(float).values
    dc_vals  = valid[dc_col].astype(float).values
    names    = valid["Name"].values
    species  = valid["Species ID"].values if "Species ID" in valid.columns else [""] * len(names)

    fig, ax = plt.subplots(figsize=(7, 8))
    fig.patch.set_facecolor(BG_COLOR)
    ax.set_facecolor(PANEL_COLOR)

    all_vals = list(ga_vals) + list(dc_vals)
    val_min  = min(all_vals)
    val_max  = max(all_vals)
    pad      = (val_max - val_min) * 0.18 if val_max != val_min else 1.0
    y_lo     = val_min - pad
    y_hi     = val_max + pad

    # ── Reference bands on Y axis ─────────────────────────────
    if ref_bands:
        for (blo, bhi, bcolor, balpha, blabel) in ref_bands:
            ax.axhspan(max(blo, y_lo), min(bhi, y_hi),
                       color=bcolor, alpha=balpha, zorder=1, label=blabel)

    # ── Reference lines ───────────────────────────────────────
    if ref_lines:
        for (rval, rcolor, rstyle, rlabel) in ref_lines:
            ax.axhline(rval, color=rcolor, linewidth=1.5,
                       linestyle=rstyle, alpha=0.9, zorder=2, label=rlabel)

    # ── Draw each sequence as a slope line ────────────────────
    x_ga = 0.0
    x_dc = 1.0

    for i, (ga_v, dc_v, name, sp) in enumerate(zip(ga_vals, dc_vals, names, species)):
        ecoli  = "coli" in str(sp).lower()
        marker = "^" if ecoli else "o"
        dot_color = COL_ECOLI if ecoli else COL_HUMAN

        # Line color: green if improvement, red if worse, grey if equal
        delta = dc_v - ga_v
        if better == "lower":
            if delta < -0.001:   lcolor = COL_BETTER
            elif delta > 0.001:  lcolor = COL_WORSE
            else:                lcolor = COL_DRAW
        elif better == "higher":
            if delta > 0.001:    lcolor = COL_BETTER
            elif delta < -0.001: lcolor = COL_WORSE
            else:                lcolor = COL_DRAW
        else:
            lcolor = COL_DRAW

        # Draw connecting line
        ax.plot([x_ga, x_dc], [ga_v, dc_v],
                color=lcolor, linewidth=1.8, alpha=0.75,
                zorder=3, solid_capstyle="round")

        # Draw dots
        ax.scatter(x_ga, ga_v, color=dot_color, marker=marker,
                   s=90, edgecolors=lcolor, linewidths=1.8,
                   zorder=5, alpha=0.95)
        ax.scatter(x_dc, dc_v, color=dot_color, marker=marker,
                   s=90, edgecolors=lcolor, linewidths=1.8,
                   zorder=5, alpha=0.95)

        # Labels -- GeneArt side (left)
        offset = (y_hi - y_lo) * 0.012
        ax.text(x_ga - 0.04, ga_v, str(name)[:9],
                ha="right", va="center", fontsize=7.5,
                color=TEXT_COLOR, alpha=0.85)

        # Labels -- DNA Chisel side (right) -- show delta
        delta_str = "{:+.2f}".format(delta)
        ax.text(x_dc + 0.04, dc_v, delta_str,
                ha="left", va="center", fontsize=7.5,
                color=lcolor, fontweight="bold", alpha=0.9)

    # ── X axis styling ────────────────────────────────────────
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(y_lo, y_hi)
    ax.set_xticks([x_ga, x_dc])
    ax.set_xticklabels(["GeneArt", "DNA Chisel"],
                       fontsize=12, fontweight="bold", color=TEXT_COLOR)
    ax.tick_params(axis="x", length=0, pad=10)
    ax.tick_params(axis="y", labelsize=9, labelcolor=MUTED_COLOR)

    # Remove top/right/bottom spines, keep left
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_edgecolor(GRID_COLOR)

    ax.set_title(title, fontsize=13, fontweight="bold",
                 color=TEXT_COLOR, pad=16)
    ax.set_ylabel("Value", fontsize=10, color=MUTED_COLOR, labelpad=8)
    ax.grid(axis="y", color=GRID_COLOR, linewidth=0.6, alpha=0.6, zorder=0)
    ax.set_axisbelow(True)

    # ── Better direction label ────────────────────────────────
    direction_map = {"lower": "Lower is better", "higher": "Higher is better",
                     "middle": "Middle range is best"}
    direction_color = {"lower": COL_BETTER, "higher": COL_BETTER, "middle": REF_AMBER}
    ax.text(0.99, 0.98, direction_map.get(better, ""),
            transform=ax.transAxes, fontsize=8.5,
            color=direction_color.get(better, MUTED_COLOR),
            va="top", ha="right", fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor=PANEL_COLOR,
                      edgecolor=GRID_COLOR, alpha=0.85))

    # ── Reference annotation ──────────────────────────────────
    if ref_annotation:
        ax.text(0.99, 0.02, ref_annotation, transform=ax.transAxes,
                fontsize=7.5, color=MUTED_COLOR, va="bottom", ha="right",
                style="italic",
                bbox=dict(boxstyle="round,pad=0.3", facecolor=PANEL_COLOR,
                          edgecolor=GRID_COLOR, alpha=0.75))

    # ── Delta column header ───────────────────────────────────
    ax.text(x_dc + 0.04, y_hi, "delta",
            ha="left", va="top", fontsize=7.5,
            color=MUTED_COLOR, style="italic")

    # ── Legend ────────────────────────────────────────────────
    legend_elements = [
        plt.Line2D([0],[0], color=COL_BETTER, linewidth=2,
                   label="DNA Chisel better" if better != "higher" else "DNA Chisel higher"),
        plt.Line2D([0],[0], color=COL_WORSE,  linewidth=2,
                   label="GeneArt better" if better != "higher" else "GeneArt higher"),
        plt.Line2D([0],[0], color=COL_DRAW,   linewidth=2, label="Equal"),
        plt.Line2D([0],[0], marker="o", color="none",
                   markerfacecolor=COL_HUMAN, markersize=8, label="Human"),
        plt.Line2D([0],[0], marker="^", color="none",
                   markerfacecolor=COL_ECOLI, markersize=8, label="E. coli"),
    ]
    ax.legend(handles=legend_elements, loc="lower left",
              fontsize=7.5, framealpha=0.9, edgecolor=GRID_COLOR,
              ncol=2)

    plt.tight_layout(pad=1.8)
    out_path = os.path.join(PLOTS_FOLDER, filename)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved: {}".format(out_path))


# ── Plot definitions ──────────────────────────────────────────
rdf = results_df.copy()
for col in rdf.columns:
    try: rdf[col] = pd.to_numeric(rdf[col], errors="ignore")
    except: pass

plots = [
    dict(
        ga_col="GeneArt CAI", dc_col="DNAChisel CAI",
        title="1. Codon Adaptation Index (CAI)",
        better="higher",
        ref_bands=[(0.85, 1.01, REF_GREEN, 0.08, "Excellent (>0.85)")],
        ref_lines=[(0.85, REF_GREEN, "--", "CAI = 0.85 (excellent threshold)")],
        ref_annotation="CAI > 0.85 = well optimized",
        filename="01_CAI.png",
    ),
    dict(
        ga_col="GeneArt rare codon % (w<0.3)", dc_col="DNAChisel rare codon % (w<0.3)",
        title="2. Rare Codon % (w < 0.3)",
        better="lower",
        ref_bands=[(0, 5, REF_GREEN, 0.08, "Safe (<5%)")],
        ref_lines=[(5, REF_RED, "--", ">5% = problematic")],
        ref_annotation="< 5% rare codons recommended",
        filename="02_rare_codons.png",
    ),
    dict(
        ga_col="GeneArt GC% overall", dc_col="DNAChisel GC% overall",
        title="3. GC Content Overall (%)",
        better="middle",
        ref_bands=[(40, 65, REF_GREEN, 0.08, "Ideal: 40-65% (human)")],
        ref_lines=[(40, REF_AMBER, ":", "min 40%"), (65, REF_AMBER, ":", "max 65%")],
        ref_annotation="Ideal: 40-65% (human), 40-55% (E. coli)",
        filename="03_GC_overall.png",
    ),
    dict(
        ga_col="GeneArt GC% max (60bp)", dc_col="DNAChisel GC% max (60bp)",
        title="4. GC% Max per 60bp Window",
        better="lower",
        ref_bands=[(0, 70, REF_GREEN, 0.06, "Safe (<70%)")],
        ref_lines=[(70, REF_RED, "--", "70% = synthesis failure risk")],
        ref_annotation="Spikes > 70% cause local synthesis failure",
        filename="04_GC_max.png",
    ),
    dict(
        ga_col="GeneArt GC% min (60bp)", dc_col="DNAChisel GC% min (60bp)",
        title="5. GC% Min per 60bp Window",
        better="higher",
        ref_bands=[(30, 100, REF_GREEN, 0.06, "Safe (>30%)")],
        ref_lines=[(30, REF_RED, "--", "30% = instability threshold")],
        ref_annotation="Drops < 30% GC cause DNA instability",
        filename="05_GC_min.png",
    ),
    dict(
        ga_col="GeneArt GC std dev", dc_col="DNAChisel GC std dev",
        title="6. GC Standard Deviation (Smoothness)",
        better="lower",
        ref_bands=[(0, 3, REF_GREEN, 0.08, "Smooth (<3%)")],
        ref_lines=[(3, REF_AMBER, "--", "std dev = 3% (smooth threshold)")],
        ref_annotation="Lower = smoother GC distribution",
        filename="06_GC_std.png",
    ),
    dict(
        ga_col="GeneArt GC3 content %", dc_col="DNAChisel GC3 content %",
        title="7. GC3 Content % (Wobble Position)",
        better="middle",
        ref_bands=[(60, 80, REF_GREEN, 0.07, "Balanced: 60-80%")],
        ref_lines=[(80, REF_AMBER, "--", "80% = aggressive CAI bias")],
        ref_annotation="High GC3 = aggressive CAI maximization",
        filename="07_GC3.png",
    ),
    dict(
        ga_col="GeneArt total restriction sites", dc_col="DNAChisel total restriction sites",
        title="8. Total Restriction Sites",
        better="lower",
        ref_bands=[(0, 0.5, REF_GREEN, 0.15, "Ideal: 0 sites")],
        ref_lines=None,
        ref_annotation="0 restriction sites = safe for cloning",
        filename="08_restriction_sites.png",
    ),
    dict(
        ga_col="GeneArt cryptic splice sites", dc_col="DNAChisel cryptic splice sites",
        title="9. Cryptic Splice Sites",
        better="lower",
        ref_bands=[(0, 0.5, REF_GREEN, 0.15, "Ideal: 0 sites")],
        ref_lines=None,
        ref_annotation="0 = no accidental mRNA splicing",
        filename="09_splice_sites.png",
    ),
    dict(
        ga_col="GeneArt CpG density (per 100bp)", dc_col="DNAChisel CpG density (per 100bp)",
        title="10. CpG Density (per 100bp)",
        better="lower",
        ref_bands=[(0, 6, REF_GREEN, 0.07, "Acceptable: <6")],
        ref_lines=[(6, REF_RED, "--", ">6 = methylation silencing risk")],
        ref_annotation="High CpG triggers epigenetic silencing",
        filename="10_CpG_density.png",
    ),
    dict(
        ga_col="GeneArt repeat count (>=10bp)", dc_col="DNAChisel repeat count (>=10bp)",
        title="11. Repeat Count (>=10bp)",
        better="lower",
        ref_bands=[(0, 0.5, REF_GREEN, 0.15, "Ideal: 0 repeats")],
        ref_lines=None,
        ref_annotation="Repeats cause genetic instability",
        filename="11_repeat_count.png",
    ),
    dict(
        ga_col="GeneArt longest repeat (bp)", dc_col="DNAChisel longest repeat (bp)",
        title="12. Longest Repeat (bp)",
        better="lower",
        ref_bands=[(0, 15, REF_GREEN, 0.07, "Safe: <15bp")],
        ref_lines=[(15, REF_RED, "--", "15bp = recombination risk")],
        ref_annotation="Repeats > 15bp trigger recombination",
        filename="12_longest_repeat.png",
    ),
    dict(
        ga_col="GeneArt homopolymer runs (>=6)", dc_col="DNAChisel homopolymer runs (>=6)",
        title="13. Homopolymer Runs (>=6 bases)",
        better="lower",
        ref_bands=[(0, 0.5, REF_GREEN, 0.15, "Ideal: 0 runs")],
        ref_lines=None,
        ref_annotation="Homopolymers cause insertions/deletions",
        filename="13_homopolymer.png",
    ),
    dict(
        ga_col="GeneArt hairpin count (all)", dc_col="DNAChisel hairpin count (all)",
        title="14. Total Hairpin Count",
        better="lower",
        ref_bands=None,
        ref_lines=None,
        ref_annotation="Fewer hairpins = smoother translation",
        filename="14_hairpins_all.png",
    ),
    dict(
        ga_col="GeneArt hairpin count (5' region)", dc_col="DNAChisel hairpin count (5' region)",
        title="15. Hairpins in 5' Region (first 100bp)",
        better="lower",
        ref_bands=[(0, 3, REF_GREEN, 0.1, "Safe: <3")],
        ref_lines=[(3, REF_RED, "--", ">=3 blocks ribosome binding")],
        ref_annotation="5' hairpins block ribosome access",
        filename="15_hairpins_5prime.png",
    ),
]

for p in plots:
    make_slopegraph(rdf, **p)

# ── Summary scorecard ─────────────────────────────────────────
print("  Generating summary scorecard...")
metric_labels = [
    "CAI", "Rare Codon %", "GC Overall", "GC Max (60bp)",
    "GC Min (60bp)", "GC Std Dev", "GC3 %", "Restriction Sites",
    "Splice Sites", "CpG Density", "Repeats (>=10bp)", "Longest Repeat",
    "Homopolymers", "Hairpins (all)", "Hairpins (5')"
]
ga_cols_sum = [
    "GeneArt CAI", "GeneArt rare codon % (w<0.3)", "GeneArt GC% overall",
    "GeneArt GC% max (60bp)", "GeneArt GC% min (60bp)", "GeneArt GC std dev",
    "GeneArt GC3 content %", "GeneArt total restriction sites",
    "GeneArt cryptic splice sites", "GeneArt CpG density (per 100bp)",
    "GeneArt repeat count (>=10bp)", "GeneArt longest repeat (bp)",
    "GeneArt homopolymer runs (>=6)", "GeneArt hairpin count (all)",
    "GeneArt hairpin count (5' region)"
]
dc_cols_sum  = [c.replace("GeneArt", "DNAChisel") for c in ga_cols_sum]
betters_sum  = ["higher","lower","middle","lower","higher","lower",
                "middle","lower","lower","lower","lower","lower",
                "lower","lower","lower"]

dc_wins = ga_wins = draws = 0
bar_colors = []
win_labels = []

for ga_c, dc_c, bet in zip(ga_cols_sum, dc_cols_sum, betters_sum):
    valid = rdf.dropna(subset=[ga_c, dc_c])
    if valid.empty:
        bar_colors.append(COL_DRAW); win_labels.append("N/A"); draws += 1; continue
    ga_v = pd.to_numeric(valid[ga_c], errors="coerce").mean()
    dc_v = pd.to_numeric(valid[dc_c], errors="coerce").mean()
    if bet == "higher":
        winner = "dc" if dc_v > ga_v + 0.001 else ("ga" if ga_v > dc_v + 0.001 else "draw")
    elif bet == "lower":
        winner = "dc" if dc_v < ga_v - 0.001 else ("ga" if ga_v < dc_v - 0.001 else "draw")
    else:
        winner = "draw"
    if winner == "dc":
        dc_wins += 1; bar_colors.append(COL_BETTER); win_labels.append("DNA Chisel")
    elif winner == "ga":
        ga_wins += 1; bar_colors.append(COL_WORSE);  win_labels.append("GeneArt")
    else:
        draws   += 1; bar_colors.append(COL_DRAW);   win_labels.append("Draw")

fig, (ax_bar, ax_pie) = plt.subplots(1, 2, figsize=(14, 7),
                                       gridspec_kw={"width_ratios": [2.5, 1]})
fig.patch.set_facecolor(BG_COLOR)
for i, (label, wl, bc) in enumerate(zip(metric_labels, win_labels, bar_colors)):
    ax_bar.barh(i, 1, color=bc, height=0.65, edgecolor=BG_COLOR, linewidth=0.5)
    ax_bar.text(-0.02, i, label, va="center", ha="right",
                fontsize=10, color=TEXT_COLOR, fontweight="500")
    ax_bar.text(0.5, i, wl, va="center", ha="center",
                fontsize=9, color="white" if bc != COL_DRAW else TEXT_COLOR,
                fontweight="bold")
ax_bar.set_xlim(-1.8, 1.2)
ax_bar.set_ylim(-0.8, len(metric_labels) - 0.2)
ax_bar.set_facecolor(PANEL_COLOR)
ax_bar.axis("off")
ax_bar.set_title("Parameter-by-Parameter Winner", fontsize=13,
                  fontweight="bold", color=TEXT_COLOR, pad=16)

pie_vals   = [v for v in [dc_wins, ga_wins, draws] if v > 0]
pie_labels_full = []
pie_colors_full = []
if dc_wins > 0: pie_labels_full.append("DNA Chisel ({})".format(dc_wins)); pie_colors_full.append(COL_BETTER)
if ga_wins > 0: pie_labels_full.append("GeneArt ({})".format(ga_wins)); pie_colors_full.append(COL_WORSE)
if draws   > 0: pie_labels_full.append("Draw ({})".format(draws)); pie_colors_full.append(COL_DRAW)

wedges, texts, autotexts = ax_pie.pie(
    pie_vals, labels=pie_labels_full, colors=pie_colors_full,
    autopct="%1.0f%%", startangle=90,
    textprops={"color": TEXT_COLOR, "fontsize": 10},
    wedgeprops={"edgecolor": BG_COLOR, "linewidth": 2}
)
for at in autotexts:
    at.set_color("white"); at.set_fontweight("bold"); at.set_fontsize(10)
ax_pie.set_facecolor(PANEL_COLOR)
ax_pie.set_title("Overall Summary -- {} parameters compared".format(len(metric_labels)),
                  fontsize=12, fontweight="bold", color=TEXT_COLOR, pad=12)
fig.suptitle("GeneArt vs DNA Chisel -- All Parameters Scorecard",
             fontsize=15, fontweight="bold", color=TEXT_COLOR, y=1.01)
plt.tight_layout(pad=2.0)
scorecard_path = os.path.join(PLOTS_FOLDER, "00_scorecard_summary.png")
plt.savefig(scorecard_path, dpi=150, bbox_inches="tight")
plt.close()
print("  Saved: {}".format(scorecard_path))

print("\n" + "=" * 60)
print("  All done!")
print("  Excel:  {}".format(OUTPUT_FILE))
print("  Plots:  {}/  ({} files)".format(PLOTS_FOLDER, len(plots) + 1))
print("=" * 60)