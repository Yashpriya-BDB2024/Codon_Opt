# -*- coding: utf-8 -*-
"""
run_comparison.py  -  GeneArt (WT) vs the top-10 DNA Chisel sequences, all genes at once
========================================================================================

Drop this file next to codon_compare.py / codon_features.py (Sequence_scoring/).
It IMPORTS both and does not modify them:

    codon_compare.compare()          -> every WT-vs-Chisel pair, the 23 scored metrics + winners
    codon_compare.overlay_profiles() -> the per-pair overlay figure (as in the scoring guide)
    codon_compare.scorecard()        -> the per-pair scorecard figure (as in the scoring guide)
    codon_features.plot_profiles()   -> the 9-panel single-sequence figure, for WT and every Chisel sequence
    codon_features.*                 -> scan_positions / codon helpers for the few extras below

Label A = "WT"  = the GeneArt sequence in the trials workbook ("DNA optimized sequence").
Label B = "DNAChisel" = a row of the 'Top' sheet written by filter_top_strict_gc.py.
Both go through the SAME scorer, the SAME human_codon.csv and the SAME human_tgcn.csv,
exactly like the usage example in the scoring guide (--label-a WT --label-b Optimized).

Run (from Sequence_scoring/, after pip install seqfold):
    python run_comparison.py
    python run_comparison.py --workers 4 --no-pair-plots        # faster, skips the 160 per-pair PNGs

What it adds on top of the given scripts (nothing existing is re-implemented):
  * batch loop over 8 genes x 10 sequences (80 pairs) + one results workbook
  * QC gate: translation == protein, same length, no internal stop
  * restriction sites split into the 13 enzymes the grid forced DNA Chisel to avoid vs the
    21 other enzymes in the scorer's list (answers "was the enzyme list fair?")
  * host-reference plots from host_files/Human (codon usage /1000, GC1/GC2/GC3, dinucleotides, ENc)
  * summary plots: metric-by-gene, win heat-map, win rate per metric, tally, CAI-vs-GC, gene profiles
"""

from __future__ import annotations

import argparse
import ast
import re
import sys
import time
import warnings
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

import codon_compare as cc            # given, unchanged
import codon_features as cf           # given, unchanged

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
LABEL_A, LABEL_B = "WT", "DNAChisel"             # WT = GeneArt sequence from the trials file
HOST_LABEL = "human"
HOST_DIR = ROOT / "host_files" / "Human"
CODON_CSV = HOST_DIR / "human_codon.csv"
TGCN_CSV = HOST_DIR / "human_tgcn.csv"
HOST_COUNTS_TXT = HOST_DIR / "Human.txt"          # CoCoPUTs codon counts (exact per-1000, GC1/2/3)
HOST_DINUC_SVG = HOST_DIR / "Homo sapiens (9606) Dinucleotide Frequencies (per 1000 dinucleotides).svg"
HOST_ENC = 55.012                                  # CoCoPUTs ENc, screenshot in host_files/Human
HOST_ENC_GCCORR = 55.407                           # CoCoPUTs ENc (GC corrected), same screenshot
GRID_SCRIPT = ROOT / "DNAChisel" / "codon_grid_search.py"
# Enzymes the grid told DNA Chisel to avoid (codon_grid_search.DEFAULT_SITES). Re-read from the
# grid script below if it is present; this list is the fallback.
GRID_ENZYMES_FALLBACK = ["BsaI", "BsmBI", "BbsI", "SapI", "EcoRI", "XhoI", "NotI",
                         "XbaI", "SpeI", "KpnI", "SalI", "AscI", "PacI"]
SENSE = [c for c, a in cf.GENETIC_CODE.items() if a != "*"]
DINUC = [a + b for a in "TCAG" for b in "TCAG"]
C_A, C_B = cc.C_A, cc.C_B


# ===========================================================================
# loading
# ===========================================================================
def grid_enzymes():
    try:
        tree = ast.parse(GRID_SCRIPT.read_text(encoding="utf-8"))
        for node in tree.body:
            if isinstance(node, ast.Assign) and any(getattr(t, "id", "") == "DEFAULT_SITES" for t in node.targets):
                return list(ast.literal_eval(node.value))
    except Exception:
        pass
    return GRID_ENZYMES_FALLBACK


def load_wt():
    """Trials workbook: headers on row 3, data from column B (same as run_human_parallel.py)."""
    trials = next(ROOT.glob("Codon*ptimization*trials*.xlsx"))
    d = pd.read_excel(trials, header=2).iloc[:, 1:5]
    d.columns = ["Name", "Protein", "Host", "WT_DNA"]
    d = d.dropna(subset=["Name", "Protein"])
    d = d[d["Host"].astype(str).str.lower().str.contains("homo")].reset_index(drop=True)
    d["Protein"] = d["Protein"].astype(str).str.strip().str.upper()
    d["WT_DNA"] = d["WT_DNA"].map(cf.clean_seq)
    return d


def translate(seq):
    return "".join(cf.GENETIC_CODE.get(seq[i:i + 3], "X") for i in range(0, len(seq) - len(seq) % 3, 3))


# ===========================================================================
# one WT-vs-Chisel pair  (runs in a worker process)
# ===========================================================================
def _pair_job(a):
    gene, rank, wt, opt, utr5 = a
    ra, rb, rows = cc.compare(wt, opt, LABEL_A, LABEL_B, host=HOST_LABEL, host_label=HOST_LABEL,
                              codon_table=str(CODON_CSV), tgcn=str(TGCN_CSV), utr5=utr5)
    return gene, rank, ra, rb, rows


# ===========================================================================
# extras computed with codon_features helpers (no new scoring definitions of existing metrics)
# ===========================================================================
def restriction_split(seq, grid_set):
    grid = {k: v for k, v in cf.RESTRICTION_SITES.items() if k in grid_set}
    other = {k: v for k, v in cf.RESTRICTION_SITES.items() if k not in grid_set}
    n_grid = sum(cf.scan_positions(seq, grid, True)[0].values())
    n_other = sum(cf.scan_positions(seq, other, True)[0].values())
    return n_grid, n_other


def pos_gc(seq):
    cod = cf.codons_of(seq)
    n = len(cod) or 1
    return [100.0 * sum(c[i] in "GC" for c in cod) / n for i in range(3)]


def codon_per1000(seq):
    cnt = Counter(cf.coding_codons(seq))
    n = sum(cnt.values()) or 1
    return {c: 1000.0 * cnt.get(c, 0) / n for c in SENSE}


def dinuc_per1000(seq):
    seq = cf.strip_trailing_stop(cf.clean_seq(seq))
    tot = Counter(seq[i:i + 2] for i in range(len(seq) - 1))
    jun = Counter(seq[i:i + 2] for i in range(2, len(seq) - 1, 3))
    nt, nj = sum(tot.values()) or 1, sum(jun.values()) or 1
    return ({d: 1000.0 * tot.get(d, 0) / nt for d in DINUC},
            {d: 1000.0 * jun.get(d, 0) / nj for d in DINUC})


def load_host_counts():
    """Human.txt (CoCoPUTs): codon, per-1000, (count). Returns per-1000 over SENSE codons + GC1/2/3."""
    if not HOST_COUNTS_TXT.exists():
        return None
    txt = HOST_COUNTS_TXT.read_text(encoding="utf-8", errors="ignore")
    cnt = {m[0]: int(m[2]) for m in re.findall(r"([ACGT]{3})\s+([\d.]+)\s+\(\s*(\d+)\)", txt)}
    if len(cnt) != 64:
        return None
    sense_tot = sum(cnt[c] for c in SENSE)
    per1000 = {c: 1000.0 * cnt[c] / sense_tot for c in SENSE}
    gc = [100.0 * sum(cnt[c] for c in SENSE if c[i] in "GC") / sense_tot for i in range(3)]
    return {"per1000": per1000, "gc1": gc[0], "gc2": gc[1], "gc3": gc[2], "gc": float(np.mean(gc))}


def load_host_dinuc():
    """Read the two bar series (total / junction) out of the CoCoPUTs SVG and normalise each to 1000."""
    if not HOST_DINUC_SVG.exists():
        return None
    t = HOST_DINUC_SVG.read_text(encoding="utf-8", errors="ignore")
    bars = re.findall(r'<rect x="(\d+)" y="\d+" width="\d+" height="(\d+)" stroke="none" '
                      r'stroke-width="0" fill="(#7a29cc|#7acc29)"', t)
    # bars live inside the plot area (x 90..1190); the legend swatches sit at x ~ 1197 and are skipped
    tot = [int(h) for x, h, f in bars if f == "#7a29cc" and 90 <= int(x) < 1190]
    jun = [int(h) for x, h, f in bars if f == "#7acc29" and 90 <= int(x) < 1190]
    if len(tot) != 16 or len(jun) != 16:
        return None
    order = re.findall(r"<text[^>]*>([ACGT]{2})</text>", t)[:16]       # label order in the chart
    if len(order) != 16:
        return None
    nt, nj = float(sum(tot)), float(sum(jun))
    return ({d: 1000 * h / nt for d, h in zip(order, tot)}, {d: 1000 * h / nj for d, h in zip(order, jun)})


# ===========================================================================
# plots
# ===========================================================================
def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.spines.top": False, "axes.spines.right": False})
    return plt


GENE_COL = ["#1565C0", "#E65100", "#2E7D32", "#6A1B9A", "#C62828", "#00838F", "#AD1457", "#5D4037"]


def plot_metric_by_gene(long, genes, path, scored_only=False):
    plt = _plt()
    metrics = [m for m, d in cc.COMPARE_METRICS]
    metrics = [m for m in metrics if long[long.metric == m][[LABEL_A, LABEL_B]].notna().any().any()]
    ncol = 4
    nrow = int(np.ceil(len(metrics) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.0 * nrow))
    axs = np.atleast_1d(axs).ravel()
    direction = dict(cc.COMPARE_METRICS)
    rng = np.random.RandomState(1)
    for ax, m in zip(axs, metrics):
        sub = long[long.metric == m]
        for gi, g in enumerate(genes):
            s = sub[sub.Gene == g]
            vb = pd.to_numeric(s[LABEL_B], errors="coerce").dropna().values
            va = pd.to_numeric(s[LABEL_A], errors="coerce").dropna().values
            if len(vb):
                ax.scatter(gi + rng.uniform(-0.18, 0.18, len(vb)), vb, s=14, color=C_B, alpha=0.55, zorder=2)
                ax.hlines(vb.mean(), gi - 0.3, gi + 0.3, color=C_B, lw=2, zorder=3)
            if len(va):
                ax.scatter([gi], [va[0]], marker="D", s=46, color=C_A, edgecolor="white", zorder=4)
        ax.set_xticks(range(len(genes)))
        ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=7.5)
        ax.set_title("{}  [{}]".format(m, {"higher": "higher better", "lower": "lower better",
                                           "band": "closer to GC band mid", "na": "not scored"}[direction[m]]),
                     fontsize=9, fontweight="bold", loc="left")
        ax.grid(axis="y", alpha=0.3)
    for ax in axs[len(metrics):]:
        ax.axis("off")
    h = [plt.Line2D([], [], marker="D", ls="", color=C_A, label="WT (GeneArt)"),
         plt.Line2D([], [], marker="o", ls="", color=C_B, alpha=.6, label="DNAChisel top-10 (dots) / mean (bar)")]
    fig.legend(handles=h, loc="upper right", fontsize=10, ncol=2)
    fig.suptitle("Every compared metric, per gene: GeneArt WT vs DNA Chisel top-10", fontsize=15, fontweight="bold", y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_win_heatmap(long, genes, path):
    plt = _plt()
    sc = long[long.w.isin(["A", "B", "draw"])]
    metrics = [m for m, d in cc.COMPARE_METRICS if d != "na" and m in set(sc.metric)]
    M = np.full((len(genes) + 1, len(metrics)), np.nan)
    T = [[""] * len(metrics) for _ in range(len(genes) + 1)]
    for j, m in enumerate(metrics):
        for i, g in enumerate(genes + ["ALL"]):
            s = sc[(sc.metric == m) & ((sc.Gene == g) if g != "ALL" else True)]
            if len(s):
                b, a, d = (s.w == "B").sum(), (s.w == "A").sum(), (s.w == "draw").sum()
                M[i, j] = (b - a) / len(s) * 100          # net result: draws sit at 0 (neutral)
                T[i][j] = "{}-{}-{}".format(b, a, d)
    fig, ax = plt.subplots(figsize=(0.82 * len(metrics) + 2.5, 0.62 * (len(genes) + 1) + 2.2))
    im = ax.imshow(np.ma.masked_invalid(M), cmap="RdYlGn", vmin=-100, vmax=100, aspect="auto")
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels(metrics, rotation=60, ha="right", fontsize=8.5)
    ax.set_yticks(range(len(genes) + 1))
    ax.set_yticklabels(genes + ["ALL"], fontsize=9)
    for i in range(len(genes) + 1):
        for j in range(len(metrics)):
            ax.text(j, i, T[i][j], ha="center", va="center", fontsize=7)
    ax.axhline(len(genes) - 0.5, color="k", lw=1.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.025)
    cb.set_label("net result: green = DNAChisel ahead, red = GeneArt WT ahead, yellow = draws / even")
    ax.set_title("Who wins each metric?  cell text = DNAChisel wins - WT wins - draws (out of 10 pairs per gene, 80 in ALL)",
                 fontsize=11, fontweight="bold", loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_win_rate(long, path):
    plt = _plt()
    sc = long[long.w.isin(["A", "B", "draw"])]
    rows = []
    for m, d in cc.COMPARE_METRICS:
        s = sc[sc.metric == m]
        if len(s):
            rows.append((m, (s.w == "B").mean() * 100, (s.w == "draw").mean() * 100, (s.w == "A").mean() * 100))
    rows.sort(key=lambda r: r[1])
    fig, ax = plt.subplots(figsize=(10, 0.36 * len(rows) + 1.6))
    y = np.arange(len(rows))
    b = np.array([r[1] for r in rows]); dr = np.array([r[2] for r in rows]); a = np.array([r[3] for r in rows])
    ax.barh(y, b, color=cc.GREEN, label="DNAChisel wins")
    ax.barh(y, dr, left=b, color=cc.DRAW, label="draw")
    ax.barh(y, a, left=b + dr, color=cc.RED, label="GeneArt WT wins")
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=9)
    ax.set_xlabel("% of the 80 pairs"); ax.legend(loc="lower right", fontsize=9)
    ax.set_title("Win rate per scored metric over all 80 pairs", fontsize=12, fontweight="bold", loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_tally(pairs, genes, path):
    plt = _plt()
    g = pairs.groupby("Gene")[["wins_" + LABEL_B, "draws", "wins_" + LABEL_A]].mean().reindex(genes)
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(genes))
    ax.bar(x, g["wins_" + LABEL_B], color=cc.GREEN, label="DNAChisel wins")
    ax.bar(x, g["draws"], bottom=g["wins_" + LABEL_B], color=cc.DRAW, label="draws")
    ax.bar(x, g["wins_" + LABEL_A], bottom=g["wins_" + LABEL_B] + g["draws"], color=cc.RED, label="GeneArt WT wins")
    ax.set_xticks(x); ax.set_xticklabels(genes)
    ax.set_ylabel("mean number of scored metrics (per pair)")
    ax.legend(fontsize=9)
    ax.set_title("Scorecard tally per gene (mean over the 10 pairs)", fontsize=12, fontweight="bold", loc="left")
    fig.tight_layout(); fig.savefig(path, dpi=140, bbox_inches="tight"); plt.close(fig)


def plot_cai_gc(allsc, genes, path):
    plt = _plt()
    fig, axs = plt.subplots(1, 2, figsize=(13, 5))
    for ax, (xm, ym, ttl) in zip(axs, [("GC_pct", "CAI", "CAI vs overall GC%"), ("GC3_pct", "tAI", "tAI vs GC3%")]):
        for gi, g in enumerate(genes):
            s = allsc[(allsc.Gene == g) & (allsc.Role == LABEL_B)]
            w = allsc[(allsc.Gene == g) & (allsc.Role == LABEL_A)]
            ax.scatter(s[xm], s[ym], s=22, color=GENE_COL[gi % 8], alpha=.7, label=g)
            ax.scatter(w[xm], w[ym], s=130, marker="*", color=GENE_COL[gi % 8], edgecolor="k", zorder=5)
        ax.set_xlabel(xm); ax.set_ylabel(ym); ax.grid(alpha=.3)
        ax.set_title(ttl + "  (dots = DNAChisel, stars = GeneArt WT)", fontsize=10.5, fontweight="bold", loc="left")
    axs[0].axvspan(45, 55, color="#FFF3E0", zorder=0); axs[0].legend(fontsize=7.5, ncol=2)
    fig.tight_layout(); fig.savefig(path, dpi=140, bbox_inches="tight"); plt.close(fig)


def plot_gene_profiles(gene, ra_wt, rbs, path):
    plt = _plt()
    A = ra_wt["profiles"]
    keys = [("GC % (nt window)", lambda P: P["nt_window"]["gc"], "nt"),
            ("MFE dG (nt window)", lambda P: P["mfe_window"], "nt"),
            ("%MinMax (codon)", lambda P: P["codon_window"]["minmax"], "codon"),
            ("Rolling CAI", lambda P: P["codon_window"]["rolling_cai"], "codon"),
            ("Rolling tAI", lambda P: P["codon_window"]["rolling_tai"], "codon"),
            ("Windowed ENC", lambda P: P["codon_window"]["enc_window"], "codon"),
            ("CpG density /100bp", lambda P: P["nt_window"]["cpg_density"], "nt")]
    keys = [k for k in keys if k[1](A).get("value")]
    fig, axs = plt.subplots(len(keys), 1, figsize=(11, 2.1 * len(keys)))
    axs = np.atleast_1d(axs)
    for ax, (title, get, xl) in zip(axs, keys):
        for k, rb in enumerate(rbs):
            p = get(rb["profiles"])
            if p.get("value"):
                ax.plot(p["index"], p["value"], color=C_B, lw=0.8, alpha=0.45, label=LABEL_B + " (10)" if k == 0 else None)
        p = get(A)
        ax.plot(p["index"], p["value"], color=C_A, lw=1.8, label="WT (GeneArt)")
        ax.set_title(title, fontsize=10, loc="left", fontweight="bold"); ax.set_xlabel(xl)
        ax.grid(alpha=.3); ax.legend(fontsize=8, loc="best")
    fig.suptitle("{}: WT vs the 10 DNA Chisel sequences".format(gene), fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.99]); fig.savefig(path, dpi=130, bbox_inches="tight"); plt.close(fig)


def plot_host_codon(host, wt_seqs, ch_seqs, path):
    plt = _plt()
    wt = Counter(); ch = Counter()
    for s in wt_seqs: wt.update(cf.coding_codons(s))
    for s in ch_seqs: ch.update(cf.coding_codons(s))
    order = [c for c in [b1 + b2 + b3 for b1 in "TCAG" for b2 in "TCAG" for b3 in "TCAG"] if c in SENSE]
    order.sort(key=lambda c: (cf.GENETIC_CODE[c], c))
    nw, nc = sum(wt.values()), sum(ch.values())
    x = np.arange(len(order))
    fig, ax = plt.subplots(figsize=(19, 5.2))
    ax.bar(x - 0.28, [host["per1000"][c] for c in order], 0.27, color="#9E9E9E", label="Human host (CoCoPUTs)")
    ax.bar(x, [1000 * wt[c] / nw for c in order], 0.27, color=C_A, label="GeneArt WT (8 genes pooled)")
    ax.bar(x + 0.28, [1000 * ch[c] / nc for c in order], 0.27, color=C_B, label="DNAChisel top-10 (80 pooled)")
    ax.set_xticks(x); ax.set_xticklabels(["{}\n{}".format(cf.GENETIC_CODE[c], c) for c in order], fontsize=6.5)
    ax.set_ylabel("per 1000 sense codons"); ax.legend(); ax.grid(axis="y", alpha=.3)
    ax.set_title("Codon usage per 1000 codons, grouped by amino acid", fontsize=12, fontweight="bold", loc="left")
    fig.tight_layout(); fig.savefig(path, dpi=140, bbox_inches="tight"); plt.close(fig)


def plot_host_gc(host, allsc, genes, path):
    plt = _plt()
    cols = ["GC_pct", "GC1", "GC2", "GC3_pct"]
    names = ["GC overall", "GC1", "GC2", "GC3"]
    hv = [host["gc"], host["gc1"], host["gc2"], host["gc3"]]
    fig, ax = plt.subplots(figsize=(10, 5.2))
    x = np.arange(4)
    ax.bar(x - 0.28, hv, 0.27, color="#9E9E9E", label="Human host")
    w = allsc[allsc.Role == LABEL_A]; c = allsc[allsc.Role == LABEL_B]
    ax.bar(x, [w[k].mean() for k in cols], 0.27, color=C_A, label="GeneArt WT (mean of 8)")
    ax.bar(x + 0.28, [c[k].mean() for k in cols], 0.27, color=C_B, label="DNAChisel (mean of 80)")
    for i, k in enumerate(cols):
        ax.scatter(np.full(len(w), i), w[k], s=14, color="k", zorder=4)
        ax.scatter(np.full(len(c), i + 0.28) + np.random.RandomState(i).uniform(-.06, .06, len(c)), c[k], s=8, color="k", alpha=.5, zorder=4)
    ax.axhspan(45, 55, xmin=0, xmax=0.25, color="#FFF3E0", zorder=0)
    ax.set_xticks(x); ax.set_xticklabels(names); ax.set_ylabel("GC %"); ax.legend(); ax.grid(axis="y", alpha=.3)
    ax.set_title("GC by codon position vs host (black dots = individual sequences)", fontsize=12, fontweight="bold", loc="left")
    fig.tight_layout(); fig.savefig(path, dpi=140, bbox_inches="tight"); plt.close(fig)


def plot_host_dinuc(hd, wt_d, ch_d, path):
    plt = _plt()
    fig, axs = plt.subplots(2, 1, figsize=(13, 8))
    for ax, (k, ttl) in zip(axs, [(0, "Total dinucleotides"), (1, "Junction dinucleotides (codon 3rd base + next codon 1st base)")]):
        x = np.arange(16)
        ax.bar(x - 0.28, [hd[k][d] for d in DINUC], 0.27, color="#9E9E9E", label="Human host")
        ax.bar(x, [np.mean([s[k][d] for s in wt_d]) for d in DINUC], 0.27, color=C_A, label="GeneArt WT")
        ax.bar(x + 0.28, [np.mean([s[k][d] for s in ch_d]) for d in DINUC], 0.27, color=C_B, label="DNAChisel")
        ax.set_xticks(x); ax.set_xticklabels(DINUC); ax.set_ylabel("per 1000"); ax.grid(axis="y", alpha=.3)
        ax.set_title(ttl, fontsize=11, fontweight="bold", loc="left"); ax.legend(fontsize=8)
    fig.suptitle("Dinucleotide frequencies per 1000 (host values read from the CoCoPUTs chart, normalised to 1000)", fontsize=11)
    fig.tight_layout(); fig.savefig(path, dpi=140, bbox_inches="tight"); plt.close(fig)


def plot_enc(allsc, genes, path):
    plt = _plt()
    fig, ax = plt.subplots(figsize=(10, 5))
    rng = np.random.RandomState(2)
    for gi, g in enumerate(genes):
        s = allsc[(allsc.Gene == g) & (allsc.Role == LABEL_B)]["ENC"]
        w = allsc[(allsc.Gene == g) & (allsc.Role == LABEL_A)]["ENC"]
        ax.scatter(gi + rng.uniform(-.15, .15, len(s)), s, s=16, color=C_B, alpha=.6)
        ax.scatter([gi], w, marker="D", s=60, color=C_A, edgecolor="white", zorder=4)
    ax.axhline(HOST_ENC, color="#616161", ls="--"); ax.text(len(genes) - .5, HOST_ENC + .3, "host ENc 55.0", ha="right", fontsize=9)
    ax.axhline(HOST_ENC_GCCORR, color="#9E9E9E", ls=":")
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes); ax.set_ylabel("ENC (Nc)")
    ax.grid(axis="y", alpha=.3)
    ax.set_title("Effective number of codons vs the human host value (20 = max bias, 61 = none)", fontsize=11, fontweight="bold", loc="left")
    fig.tight_layout(); fig.savefig(path, dpi=140, bbox_inches="tight"); plt.close(fig)


# ===========================================================================
# main
# ===========================================================================
def main():
    ap = argparse.ArgumentParser(description="GeneArt (WT) vs top-10 DNA Chisel, all genes at once.")
    ap.add_argument("--top", default=str(ROOT / "DNAChisel" / "outputs" / "human_strictGC_top10.xlsx"))
    ap.add_argument("--out-dir", default=str(HERE / "outputs" / "comparison"))
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--utr5", default="", help="optional upstream bases for the Kozak -6..-1 context (same for both)")
    ap.add_argument("--no-pair-plots", action="store_true", help="skip the 2 PNGs per pair (160 files)")
    args = ap.parse_args()
    t0 = time.time()
    out = Path(args.out_dir)
    (out / "plots" / "pairs").mkdir(parents=True, exist_ok=True)
    (out / "plots" / "genes").mkdir(parents=True, exist_ok=True)
    (out / "plots" / "summary").mkdir(parents=True, exist_ok=True)
    (out / "plots" / "host_reference").mkdir(parents=True, exist_ok=True)

    try:
        import seqfold  # noqa: F401
        have_seqfold = True
    except Exception:
        have_seqfold = False
        warnings.warn("seqfold is NOT installed: TIS_accessibility_5p, MFE_min_window and "
                      "stable_structure_count will be blank (pip install seqfold, then re-run).")

    wt = load_wt()
    top = pd.read_excel(args.top, sheet_name="Top")
    genes = [g for g in wt["Name"] if g in set(top["Name"])]
    print("genes:", genes, "| pairs:", int(top["Name"].isin(genes).sum()), "| seqfold:", have_seqfold)
    grid_set = set(grid_enzymes())

    # ---- QC + job list -------------------------------------------------------------------------
    jobs, qc = [], {}
    for _, w in wt.iterrows():
        g = w["Name"]
        t = top[top["Name"] == g].sort_values("Rank")
        for _, r in t.iterrows():
            opt = cf.clean_seq(r["Optimized DNA Sequence"])
            wd = w["WT_DNA"]
            p_wt, p_op = translate(cf.strip_trailing_stop(wd)), translate(cf.strip_trailing_stop(opt))
            same_len = len(wd) == len(opt)
            ident = (100.0 * sum(a == b for a, b in zip(wd, opt)) / len(wd)) if same_len else None
            ncod = sum(wd[i:i + 3] != opt[i:i + 3] for i in range(0, min(len(wd), len(opt)) - 2, 3))
            qc[(g, int(r["Rank"]))] = {
                "same_length": same_len, "protein_matches_input": p_op == w["Protein"],
                "WT_protein_matches_input": p_wt == w["Protein"], "internal_stop": "*" in p_op,
                "pct_identity_nt": None if ident is None else round(ident, 2),
                "codons_differing": ncod, "codons_total": len(wd) // 3,
                "Variant": r["Variant"], "codon_method": r.get("codon_method"), "grid_CAI": r["CAI"]}
            jobs.append((g, int(r["Rank"]), wd, opt, args.utr5))
    bad = [k for k, v in qc.items() if not (v["protein_matches_input"] and v["same_length"] and not v["internal_stop"])]
    if bad:
        sys.exit("QC FAILED (protein differs / length differs / internal stop): {}".format(bad[:10]))
    print("QC ok: all {} DNAChisel sequences translate to the input protein, same length as GeneArt, no internal stop."
          .format(len(jobs)))

    # ---- score every pair with the given compare() ---------------------------------------------
    if args.workers > 1:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            res = list(ex.map(_pair_job, jobs, chunksize=1))
    else:
        res = [_pair_job(j) for j in jobs]
    print("scored {} pairs in {:.0f}s".format(len(res), time.time() - t0))

    long_rows, pair_rows, score_rows = [], [], []
    wt_scores_done = set()
    by_gene = {g: [] for g in genes}
    ra_gene = {}
    seq_of = {(j[0], j[1]): j[3] for j in jobs}
    wt_of = {j[0]: j[2] for j in jobs}
    for gene, rank, ra, rb, rows in res:
        q = qc[(gene, rank)]
        for r in rows:
            long_rows.append({"Gene": gene, "Rank": rank, "metric": r["metric"], LABEL_A: r[LABEL_A], LABEL_B: r[LABEL_B],
                              "delta": r["delta"], "better": r["better"], "winner": r["winner"], "w": r["_w"]})
        sc = [r for r in rows if r["_w"] in ("A", "B", "draw")]
        pair_rows.append({"Gene": gene, "Rank": rank, "Variant": q["Variant"], "codon_method": q["codon_method"],
                          "scored_metrics": len(sc), "wins_" + LABEL_B: sum(r["_w"] == "B" for r in sc),
                          "wins_" + LABEL_A: sum(r["_w"] == "A" for r in sc), "draws": sum(r["_w"] == "draw" for r in sc),
                          **{k: v for k, v in q.items() if k not in ("Variant", "codon_method")}})
        by_gene[gene].append(rb)
        ra_gene[gene] = ra
        # per-pair figures straight from codon_compare
        if not args.no_pair_plots:
            cc._init_mpl()
            tag = "{}_rank{:02d}".format(gene, rank)
            cc.overlay_profiles(ra, rb, LABEL_A, LABEL_B, str(out / "plots" / "pairs" / (tag + "_overlay_profiles.png")))
            cc.scorecard(rows, LABEL_A, LABEL_B, str(out / "plots" / "pairs" / (tag + "_scorecard.png")))
        # all scalar scores of the Chisel side
        score_rows.append({"Gene": gene, "Role": LABEL_B, "Rank": rank, **rb["scores"]})
        if gene not in wt_scores_done:
            score_rows.append({"Gene": gene, "Role": LABEL_A, "Rank": 0, **ra["scores"]})
            wt_scores_done.add(gene)

    long = pd.DataFrame(long_rows)
    pairs = pd.DataFrame(pair_rows)
    allsc = pd.DataFrame(score_rows)

    # ---- extras (restriction split, per-position GC, host distance) ------------------------------
    host = load_host_counts()
    hdin = load_host_dinuc()
    ex_rows = []
    for _, r in allsc.iterrows():
        seq = wt_of[r.Gene] if r.Role == LABEL_A else seq_of[(r.Gene, int(r.Rank))]
        ng, no = restriction_split(seq, grid_set)
        g1, g2, g3 = pos_gc(seq)
        row = {"RS_in_grid_enzymes(13)": ng, "RS_in_other_enzymes(21)": no, "GC1": g1, "GC2": g2}
        if host:
            p = codon_per1000(seq)
            row["codon_usage_L1_to_host"] = sum(abs(p[c] - host["per1000"][c]) for c in SENSE) / 10.0   # % of 1000
            row["GC_pos_L1_to_host"] = abs(g1 - host["gc1"]) + abs(g2 - host["gc2"]) + abs(g3 - host["gc3"])
        row["ENC_abs_diff_to_host"] = abs(r["ENC"] - HOST_ENC) if pd.notna(r["ENC"]) else None
        if hdin:
            td, jd = dinuc_per1000(seq)
            row["dinuc_total_L1_to_host"] = sum(abs(td[d] - hdin[0][d]) for d in DINUC)
            row["dinuc_junction_L1_to_host"] = sum(abs(jd[d] - hdin[1][d]) for d in DINUC)
        ex_rows.append(row)
    allsc = pd.concat([allsc.reset_index(drop=True), pd.DataFrame(ex_rows)], axis=1)
    pairs = pairs.merge(allsc[allsc.Role == LABEL_B][["Gene", "Rank", "RS_in_grid_enzymes(13)", "RS_in_other_enzymes(21)"]],
                        on=["Gene", "Rank"])

    # ---- gene + overall summaries ----------------------------------------------------------------
    gs = []
    for (g, m), s in long.groupby(["Gene", "metric"], sort=False):
        a = pd.to_numeric(s[LABEL_A], errors="coerce"); b = pd.to_numeric(s[LABEL_B], errors="coerce")
        gs.append({"Gene": g, "metric": m, "better": s["better"].iloc[0], "WT (GeneArt)": a.iloc[0] if a.notna().any() else None,
                   "DNAChisel mean": b.mean(), "DNAChisel min": b.min(), "DNAChisel max": b.max(),
                   "n DNAChisel wins": int((s.w == "B").sum()), "n WT wins": int((s.w == "A").sum()),
                   "n draws": int((s.w == "draw").sum())})
    gene_summary = pd.DataFrame(gs)
    ov = []
    for m, s in long.groupby("metric", sort=False):
        sc = s[s.w.isin(["A", "B", "draw"])]
        ov.append({"metric": m, "better": s["better"].iloc[0], "pairs scored": len(sc),
                   "DNAChisel wins": int((sc.w == "B").sum()), "WT wins": int((sc.w == "A").sum()),
                   "draws": int((sc.w == "draw").sum()),
                   "mean delta (Chisel - WT)": pd.to_numeric(s["delta"], errors="coerce").mean()})
    overall = pd.DataFrame(ov)

    # ---- workbook ---------------------------------------------------------------------------------
    notes = pd.DataFrame({"Notes": [
        "WT = GeneArt sequence from the trials workbook. DNAChisel = a row of the Top sheet (strict-GC top 10 by CAI).",
        "Winners come from codon_compare.compare() unchanged: higher / lower / band (GC_pct vs human band 40-65, mid 52.5) / na (shown, not scored).",
        "seqfold installed when scored: {}. If False, TIS_accessibility_5p, MFE_min_window, stable_structure_count are blank (counted as 'na').".format(have_seqfold),
        "Restriction sites: the scorer counts 34 enzymes. The grid made DNA Chisel avoid only {} of them; GeneArt was never constrained by the grid. "
        "See columns RS_in_grid_enzymes(13) / RS_in_other_enzymes(21).".format(len(grid_set)),
        "CAI, rare_codon_pct and tAI use the same human_codon.csv that DNA Chisel optimised against (circular advantage for DNA Chisel).",
        "DNAChisel rows are the best of 2,592 settings per gene, selected on CAI (selection bias); WT is one vendor sequence.",
        "Hairpins and splice sites were never constrained in the grid, so any win there is incidental.",
        "Host reference (Human.txt exact counts; dinucleotides read from the CoCoPUTs SVG, +/-0.5 per 1000; ENc from the screenshot) is context, not a scored winner."]})
    xl = out / "comparison_results.xlsx"
    with pd.ExcelWriter(xl, engine="openpyxl") as w:
        notes.to_excel(w, sheet_name="Read_me", index=False)
        overall.to_excel(w, sheet_name="Overall_win_rate", index=False)
        gene_summary.to_excel(w, sheet_name="Gene_summary", index=False)
        pairs.to_excel(w, sheet_name="Pair_summary", index=False)
        long.drop(columns="w").to_excel(w, sheet_name="Long_all_pairs", index=False)
        allsc.to_excel(w, sheet_name="All_scores", index=False)
        for ws in w.book.worksheets:
            for col in ws.columns:
                ws.column_dimensions[col[0].column_letter].width = min(34, max(10, max(len(str(c.value or "")) for c in col[:60]) + 2))
            ws.freeze_panes = "A2"
    print("wrote", xl)

    # ---- codon_features' own 9-panel single-sequence figure, for WT and every DNA Chisel sequence ----
    if not args.no_pair_plots:
        F = out / "plots" / "features"
        F.mkdir(parents=True, exist_ok=True)
        for g in genes:
            cf.plot_profiles(dict(ra_gene[g], name="{} WT (GeneArt)".format(g)), str(F / "{}_WT_profiles.png".format(g)))
            for k, rb in enumerate(by_gene[g], start=1):
                cf.plot_profiles(dict(rb, name="{} DNAChisel rank {}".format(g, k)),
                                 str(F / "{}_DNAChisel_rank{:02d}_profiles.png".format(g, k)))

    # ---- plots ------------------------------------------------------------------------------------
    S = out / "plots" / "summary"
    plot_metric_by_gene(long, genes, S / "metrics_by_gene.png")
    plot_win_heatmap(long, genes, S / "win_heatmap.png")
    plot_win_rate(long, S / "win_rate_per_metric.png")
    plot_tally(pairs, genes, S / "tally_per_gene.png")
    plot_cai_gc(allsc, genes, S / "cai_vs_gc_tai_vs_gc3.png")
    for g in genes:
        plot_gene_profiles(g, ra_gene[g], by_gene[g], out / "plots" / "genes" / "{}_profiles_WT_vs_top10.png".format(g))
    H = out / "plots" / "host_reference"
    plot_enc(allsc, genes, H / "enc_vs_host.png")
    if host:
        plot_host_codon(host, [wt_of[g] for g in genes], [seq_of[k] for k in seq_of], H / "codon_usage_per_1000_vs_host.png")
        plot_host_gc(host, allsc, genes, H / "gc_by_position_vs_host.png")
    if hdin:
        wd = [dinuc_per1000(wt_of[g]) for g in genes]
        cd = [dinuc_per1000(seq_of[k]) for k in seq_of]
        plot_host_dinuc(hdin, wd, cd, H / "dinucleotides_vs_host.png")

    # ---- console summary ----------------------------------------------------------------------------
    pd.set_option("display.width", 200)
    print("\nOverall (80 pairs):")
    print(overall.round(3).to_string(index=False))
    print("\nTally per gene (mean wins per pair):")
    print(pairs.groupby("Gene")[["wins_" + LABEL_B, "wins_" + LABEL_A, "draws"]].mean().reindex(genes).round(2))
    print("\nDone in {:.0f}s. Outputs in {}".format(time.time() - t0, out))


if __name__ == "__main__":
    main()
