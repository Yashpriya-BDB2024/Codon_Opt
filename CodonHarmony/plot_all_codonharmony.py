# -*- coding: utf-8 -*-
"""
plot_all_codonharmony.py  -  every plot of GeneArt vs ALL Codon Harmony rows
============================================================================

Derived from plot_all.py (the DNA Chisel version, left untouched). The plotting, clustering, win/draw/loss and
profile code is IDENTICAL; only what is required to point it at Codon Harmony was changed:

  1. Labels / sheet names: "DNA Chisel" -> "Codon Harmony", "DNAChisel" -> "CodonHarmony" everywhere (titles,
     legends, file names, table column names). The score table's sheet is therefore "CodonHarmony" (as written by
     score_codonharmony.py), not "DNAChisel".
  2. Defaults: --scores outputs/scores/scores_codonharmony.xlsx, grid = ../CodonHarmony/outputs/
     human_codonharmony_grid.xlsx, --out-dir outputs/plots_codonharmony (so the DNA Chisel plots are never overwritten).
  3. --stop-policy (default strip), the same option as score_codonharmony.py: the DNA that is re-scored for the
     profile / motif / codon plots must be the DNA that was scored for the table, otherwise Check_vs_scores_all
     reports differences. Use the same value you used for score_codonharmony.py.

Inputs   scores_codonharmony.xlsx (score_codonharmony.py; sheets CodonHarmony and GeneArt, 22 metrics)
         the Codon Harmony grid workbook (sheet Results) and the original trials workbook (GeneArt DNA):
         only for the plots that need the DNA itself
Imported, NOT modified: codon_features, codon_compare, run_comparison, score_all_sequences

Run from Sequence_scoring/   (pip install numpy pandas matplotlib openpyxl scikit-learn seqfold tqdm)
    python plot_all_codonharmony.py
    python plot_all_codonharmony.py --workers 40 --clusters 3 --pair-plots 10
    python plot_all_codonharmony.py --no-profiles          # only the plots that need the scores file
Outputs: outputs/plots_codonharmony/   (same folders and files as plot_all.py)
"""


from __future__ import annotations

import argparse
import os
import pickle
import sys
import textwrap
import warnings
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap, Normalize, to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from sklearn.cluster import KMeans

import codon_compare as cc            # given, unchanged: metric list, better-direction, GC band
import codon_features as cf           # given, unchanged
import run_comparison as rc           # given, unchanged: host-reference loaders, codon / dinucleotide helpers
import score_all_sequences as sas     # given, unchanged: GeneArt loader, reference file paths, translate

HERE = Path(__file__).resolve().parent
ROOT = sas.ROOT
DNA = sas.DNA
LA, LB = "GeneArt", "CodonHarmony"
C_HOST = "#9E9E9E"
# run_comparison looks for the host files next to the parent folder; point it at the folder sas found
rc.HOST_COUNTS_TXT = ROOT / "host_files" / "Human" / "Human.txt"
rc.HOST_DINUC_SVG = ROOT / "host_files" / "Human" / "Homo sapiens (9606) Dinucleotide Frequencies (per 1000 dinucleotides).svg"


DIRECTION = dict(cc.COMPARE_METRICS)
BAND = cc.gc_band_for("human")                      # (40, 65)
BAND_MID = sum(BAND) / 2.0                          # 52.5

FAMILIES = [
    ("Codon usage", ["CAI", "rare_codon_pct", "ENC", "tAI", "tAI_ramp", "minmax_mean", "codon_autocorr_lag1"]),
    ("Composition", ["GC_pct", "GC_max_60bp", "GC3_pct", "CpG_density_per100bp"]),
    ("RNA structure", ["TIS_accessibility_5p", "MFE_min_window", "stable_structure_count",
                       "hairpin_count_all", "hairpin_count_5prime"]),
    ("Sequence liabilities", ["splice_sites_total", "homopolymer_runs_ge6", "restriction_sites_total",
                              "repeat_count_10bp", "repeat_count_20bp"]),
    ("Synthesis", ["tm_spread"]),
]
ARROW = {"higher": "\u2191", "lower": "\u2193", "band": "\u25CE", "na": "\u2013", "dist": "\u2193"}
RULE_TEXT = {
    "higher": "higher is better",
    "lower": "lower is better",
    "band": "closer to the middle of the human GC band ({:g}\u2013{:g}%, midpoint {:g}%) is better".format(*BAND, BAND_MID),
    "na": "no better/worse rule (shown for context, not scored)",
    "dist": "distance to the human host reference: lower = closer to host (context, not scored)",
}

# Reference marks drawn on the distribution / cluster plots. Thresholds are the ones listed in the baseline
# scoring audit ("Reference thresholds the plots apply"); the human GC band is the one codon_compare uses.
# plot_profiles_host.py adds the human host values (GC1/GC2/GC3, ENc) at run time.
REF_SPANS = {"GC_pct": [(BAND[0], BAND[1], "human GC band {:g}\u2013{:g}%".format(*BAND))],
             "GC3_pct": [(60.0, 80.0, "audit: 60\u201380% balanced")]}
REF_LINES = {"GC_pct": [(BAND_MID, "GC band midpoint {:g}%".format(BAND_MID), "#2E7D32")],
             "CAI": [(0.85, "audit: CAI > 0.85 excellent", "#6A1B9A")],
             "rare_codon_pct": [(5.0, "audit: > 5% problematic", "#6A1B9A")],
             "GC_max_60bp": [(70.0, "audit: > 70% local synthesis failure", "#6A1B9A")],
             "CpG_density_per100bp": [(6.0, "audit: > 6 methylation-silencing risk", "#6A1B9A")],
             "hairpin_count_5prime": [(3.0, "audit: \u2265 3 blocks ribosome binding", "#6A1B9A")]}


def draw_refs(ax, metric):
    for x0, x1, _ in REF_SPANS.get(metric, []):
        ax.axvspan(x0, x1, color="#2E7D32", alpha=0.07, zorder=0)
    for x, _, col in REF_LINES.get(metric, []):
        ax.axvline(x, color=col, ls="--", lw=1.1, zorder=0)


def ref_handles(metric):
    h = [Patch(color="#2E7D32", alpha=0.2, label=l) for _, _, l in REF_SPANS.get(metric, [])]
    return h + [Line2D([], [], color=c, ls="--", lw=1.1, label=l) for _, l, c in REF_LINES.get(metric, [])]
NAMES = {1: ["Constant"], 2: ["Low", "High"], 3: ["Low", "Mid", "High"],
         4: ["Low", "Mid-low", "Mid-high", "High"], 5: ["Very low", "Low", "Mid", "High", "Very high"]}
C_DNA, C_GA, C_NOSPREAD = "#7A8B99", "#C62828", "#BDBDBD"
EPS = 1e-9

C_DC = "#2B6CB0"          # Codon Harmony in the win/loss and scatter plots (GeneArt = C_GA)
GENE_COL = ["#1565C0", "#E65100", "#2E7D32", "#6A1B9A", "#C62828", "#00838F", "#AD1457", "#5D4037"]


def mlabel(m):
    return "{} {}".format(m, ARROW[DIRECTION[m]])


def fmt(v):
    return "{:.4g}".format(v)


def palette(k):
    cm = LinearSegmentedColormap.from_list("lmh", ["#0072B2", "#E69F00", "#D55E00"])
    return [cm(x) for x in np.linspace(0, 1, k)] if k > 1 else [cm(0.0)]


def pct_below(x, v):
    """Percentile of v among x: % of x below v, ties count half (0 = below all, 100 = above all)."""
    return 100.0 * ((x < v - EPS).sum() + 0.5 * (np.abs(x - v) <= EPS).sum()) / len(x)


def beat_share(x, v, direction):
    """% of Codon Harmony rows that GeneArt is better than (ties count half). NaN for 'na' metrics."""
    if direction == "higher":
        return pct_below(x, v)
    if direction == "lower":
        return 100.0 - pct_below(x, v)
    if direction == "band":
        d, dv = np.abs(x - BAND_MID), abs(v - BAND_MID)
        return 100.0 * ((d > dv + EPS).sum() + 0.5 * (np.abs(d - dv) <= EPS).sum()) / len(x)
    return float("nan")


def cluster_1d(x, gv, k):
    """K-means on one metric of one gene. Returns (rank label per row, sorted centres, GeneArt rank)."""
    kk = min(k, len(np.unique(x)))
    if kk == 1:
        return np.zeros(len(x), int), np.array([float(x[0])]), 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        km = KMeans(n_clusters=kk, n_init=10, random_state=42).fit(x.reshape(-1, 1))
    order = np.argsort(km.cluster_centers_.ravel())
    rank = np.empty(kk, int)
    rank[order] = np.arange(kk)
    return rank[km.labels_], km.cluster_centers_.ravel()[order], int(rank[km.predict([[gv]])[0]])


def heat_code(rank, kk, k):
    """Colour code for the cluster heatmap: position on the k-step Low..High scale; k = 'no spread'."""
    if kk == 1:
        return k
    return int(round(rank * (k - 1) / (kk - 1)))


def _frame(metric, genes, subtitle, ratio, headers):
    """Figure with a plot axis (left) and a text-table axis (right) that share the gene rows."""
    n = len(genes)
    W, H = 15.0, 0.95 * n + 3.2
    fig = plt.figure(figsize=(W, H))
    gs = fig.add_gridspec(1, 2, width_ratios=ratio, wspace=0.02, left=1.1 / W, right=0.99,
                          top=1 - 1.55 / H, bottom=1.35 / H)
    ax, axt = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1], sharey=None)
    d = DIRECTION[metric]
    fig.text(0.01, 1 - 0.25 / H, "{}   ({})".format(metric, RULE_TEXT[d]), fontsize=13.5, fontweight="bold", va="center")
    fig.text(0.01, 1 - 0.62 / H, textwrap.fill(subtitle, 170), fontsize=9, color="#444", va="top")
    for a_ in (ax, axt):
        a_.set_ylim(n - 0.5, -0.5)
    ax.set_yticks(range(n))
    ax.grid(axis="x", alpha=0.3)
    ax.set_xlabel(metric)
    axt.axis("off")
    axt.set_xlim(0, 1)
    for x, txt in headers:
        axt.text(x, -0.62, txt, fontsize=8.5, fontweight="bold", va="bottom")
    draw_refs(ax, metric)
    return fig, ax, axt


def plot_distribution(metric, dna, ga, genes, st, path, dpi):
    sub = st[st["metric"] == metric].set_index("gene")
    fig, ax, axt = _frame(metric, genes,
                          "Each grey dot = one Codon Harmony row of that gene (all Status-OK grid rows). Box = quartiles, thick bar = "
                          "median. Red diamond = GeneArt sequence of that gene. Same x-scale for every gene. Right column = GeneArt "
                          "percentile among that gene's Codon Harmony rows (0 = below all, 50 = median, 100 = above all).",
                          [4.2, 0.8], [(0.05, "GeneArt\npercentile")])
    rng = np.random.RandomState(0)
    for i, g in enumerate(genes):
        x = dna.loc[dna["Name"] == g, metric].to_numpy(float)
        ax.scatter(x, i + rng.uniform(-0.3, 0.3, len(x)), s=3, color=C_DNA, alpha=0.22, linewidths=0,
                   rasterized=True, zorder=1)
        q1, med, q3 = np.percentile(x, [25, 50, 75])
        ax.add_patch(plt.Rectangle((q1, i - 0.25), max(q3 - q1, 1e-9), 0.5, fill=False, ec="#1f2d3a", lw=1.3, zorder=3))
        ax.vlines(med, i - 0.25, i + 0.25, colors="#1f2d3a", linewidth=2.6, zorder=4)
        gv = float(sub.loc[g, "GeneArt"])
        ax.scatter([gv], [i], marker="D", s=75, color=C_GA, edgecolor="white", linewidths=1, zorder=6)
        ax.annotate(fmt(gv), (gv, i), xytext=(0, 12), textcoords="offset points", ha="center", fontsize=8,
                    color=C_GA, fontweight="bold", zorder=7)
        star = "*" if sub.loc[g, "GeneArt_outside_CodonHarmony_range"] else ""
        axt.text(0.05, i, "{:.0f}{}".format(sub.loc[g, "GeneArt_percentile"], star), va="center", fontsize=11,
                 fontweight="bold")
    ax.set_yticklabels(["{}\nn={:,}".format(g, int((dna["Name"] == g).sum())) for g in genes], fontsize=9)
    handles = [Line2D([], [], marker="o", ls="", color=C_DNA, alpha=0.6, markersize=5, label="Codon Harmony row"),
               Line2D([], [], color="#1f2d3a", lw=2.6, label="median (box = quartiles)"),
               Line2D([], [], marker="D", ls="", color=C_GA, markeredgecolor="white", markersize=9, label="GeneArt")]
    handles += ref_handles(metric)
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.45, 0.035), ncol=4, fontsize=8.5, frameon=False)
    fig.text(0.01, 0.008, "* GeneArt value lies outside the whole Codon Harmony range of that gene", fontsize=8, color="#444")
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def plot_clusters(metric, dna, ga, genes, st, cl, k, path, dpi):
    sub = st[st["metric"] == metric].set_index("gene")
    csub = cl[cl["metric"] == metric]
    cols = palette(k)
    ext = cols + [C_NOSPREAD]
    step = 0.78 / k
    heads = [(j * step, NAMES[k][j] + " cluster\nrows (%), [range]") for j in range(k)] + [(0.84, "GeneArt\nfalls in")]
    fig, ax, axt = _frame(metric, genes,
                          "K-means (k={}) on this metric, run separately inside each gene on its Codon Harmony rows. Dots are coloured by "
                          "cluster, dashed line = boundary between clusters, white star = GeneArt value of that gene. Right: rows per "
                          "cluster and value range of each cluster, and the cluster GeneArt falls in.".format(k),
                          [3.6, 3.4], heads)
    rng = np.random.RandomState(0)
    for i, g in enumerate(genes):
        x = dna.loc[dna["Name"] == g, metric].to_numpy(float)
        r = sub.loc[g]
        kk = int(r["n_clusters_used"])
        names = NAMES[kk]
        lab = dna.loc[dna["Name"] == g, "_rank_" + metric].to_numpy(int)
        cc_ = [ext[heat_code(j, kk, k)] for j in range(kk)]
        ax.scatter(x, i + rng.uniform(-0.3, 0.3, len(x)), s=3, c=[cc_[j] for j in lab], alpha=0.4,
                   linewidths=0, rasterized=True, zorder=1)
        for b in r["_bounds"]:
            ax.vlines(b, i - 0.42, i + 0.42, colors="#222", linestyles="--", linewidth=1, zorder=4)
        gv = float(r["GeneArt"])
        gcol = cc_[0] if kk == 1 else cc_[names.index(r["GeneArt_cluster"])]
        ax.scatter([gv], [i], marker="*", s=330, color="white", edgecolor="black", linewidths=1.4, zorder=6)
        if kk == k:
            for j in range(k):
                rr = csub[(csub["gene"] == g) & (csub["cluster"] == NAMES[k][j])].iloc[0]
                axt.text(j * step, i, "{:,} ({:.0f}%)\n[{} \u2013 {}]".format(rr["n_rows"], rr["pct_of_gene_rows"],
                         fmt(rr["min"]), fmt(rr["max"])), va="center", fontsize=8, color=cols[j])
        else:
            axt.text(0.0, i, ("single value {} in all {:,} rows (no clusters possible)".format(fmt(x[0]), len(x)) if kk == 1
                              else "only {} distinct values:\n".format(len(np.unique(x))) + "\n".join(
                "{} {:,} rows [{} \u2013 {}]".format(names[j], int((lab == j).sum()), fmt(x[lab == j].min()),
                                                    fmt(x[lab == j].max())) for j in range(kk))),
                     va="center", fontsize=8, color="#333")
        star = " *" if r["GeneArt_outside_CodonHarmony_range"] else ""
        axt.text(0.84, i, r["GeneArt_cluster"] + star, va="center", fontsize=11, fontweight="bold",
                 color="#555" if kk == 1 else gcol)
    ax.set_yticklabels(genes, fontsize=9.5)
    handles = [Line2D([], [], marker="o", ls="", color=cols[j], markersize=6, label=NAMES[k][j] + " cluster")
               for j in range(k)]
    handles += [Line2D([], [], color="#222", ls="--", label="cluster boundary"),
                Line2D([], [], marker="*", ls="", color="white", markeredgecolor="black", markersize=14, label="GeneArt")]
    handles += ref_handles(metric)
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.035), ncol=k + 2, fontsize=8.5, frameon=False)
    fig.text(0.01, 0.008, "* GeneArt value lies outside the whole Codon Harmony range of that gene.  Low/Mid/High describe the numeric "
             "range of the metric only, not good or bad: see the better/worse rule in the title.", fontsize=8, color="#444")
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def _groups(ax, metrics, n_rows):
    pos = 0
    for fam, ms in FAMILIES:
        cnt = sum(m in metrics for m in ms)
        if not cnt:
            continue
        x0, x1 = pos - 0.45, pos + cnt - 0.55
        ax.plot([x0, x1], [-0.62, -0.62], color="#333", lw=1.5, clip_on=False)
        ax.text((x0 + x1) / 2, -0.75, fam, ha="center", va="bottom", fontsize=9.5, fontweight="bold")
        if pos:
            ax.axvline(pos - 0.5, color="white", lw=3)
        pos += cnt


def heatmap_cluster(st, genes, metrics, k, path, dpi):
    code = np.array([[int(st[(st.metric == m) & (st.gene == g)]["_code"].iloc[0]) for m in metrics] for g in genes])
    pct = np.array([[st[(st.metric == m) & (st.gene == g)]["GeneArt_percentile"].iloc[0] for m in metrics] for g in genes])
    out = np.array([[bool(st[(st.metric == m) & (st.gene == g)]["GeneArt_outside_CodonHarmony_range"].iloc[0])
                     for m in metrics] for g in genes])
    colors = palette(k) + [C_NOSPREAD]
    fig, ax = plt.subplots(figsize=(0.78 * len(metrics) + 3, 0.62 * len(genes) + 3.6))
    ax.imshow(code, cmap=ListedColormap(colors), norm=BoundaryNorm(np.arange(-0.5, k + 1.5), k + 1), aspect="auto")
    for i in range(len(genes)):
        for j in range(len(metrics)):
            rgb = to_rgb(colors[code[i, j]])
            dark = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2] < 0.55
            ax.text(j, i, "{:.0f}{}".format(pct[i, j], "*" if out[i, j] else ""), ha="center", va="center", fontsize=8.5,
                    color="white" if dark else "black", fontweight="bold" if out[i, j] else "normal")
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels([mlabel(m) for m in metrics], rotation=50, ha="right", fontsize=8.5)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=10)
    ax.set_xticks(np.arange(-0.5, len(metrics), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(genes), 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.2)
    ax.tick_params(which="minor", length=0)
    _groups(ax, metrics, len(genes))
    ax.set_title("Which cluster does the GeneArt sequence fall into?\n"
                 "Colour = cluster (K-means inside each gene, on the Codon Harmony rows of that gene) of the metric's value range; "
                 "number = GeneArt percentile among those rows",
                 fontsize=11, fontweight="bold", loc="left", pad=48)
    handles = [Patch(color=colors[j], label="{} range".format(NAMES[k][j])) for j in range(k)]
    handles += [Patch(color=C_NOSPREAD, label="Codon Harmony has a single value (no spread)")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=k + 1, frameon=False, fontsize=9,
              title="Cell colour", title_fontsize=9)
    fig.text(0.01, 0.01, "Number in cell: percentile of GeneArt among the Codon Harmony rows of the same gene "
             "(0 = below all, 50 = median, 100 = above all; ties count half).  * = GeneArt is outside the Codon Harmony range of that gene.\n"
             "Arrows after the metric name: \u2191 higher is better, \u2193 lower is better, \u25CE closer to GC-band midpoint is better, "
             "\u2013 no better/worse rule.  Low/Mid/High describe the numeric range only, not good or bad.",
             fontsize=8, color="#333", va="bottom")
    fig.tight_layout(rect=[0, 0.07, 1, 1])
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def heatmap_beats(st, genes, metrics, path, dpi):
    ms = [m for m in metrics if DIRECTION[m] != "na"]
    mat = np.array([[st[(st.metric == m) & (st.gene == g)]["GeneArt_beats_pct"].iloc[0] for m in ms] for g in genes])
    mat = np.vstack([mat, mat.mean(axis=0)])
    rows = list(genes) + ["Mean of genes"]
    fig, ax = plt.subplots(figsize=(0.78 * len(ms) + 3.6, 0.62 * len(rows) + 3.8))
    im = ax.imshow(mat, cmap="RdYlGn", vmin=0, vmax=100, aspect="auto")
    for i in range(len(rows)):
        for j in range(len(ms)):
            ax.text(j, i, "{:.0f}".format(mat[i, j]), ha="center", va="center", fontsize=8.5,
                    fontweight="bold" if i == len(rows) - 1 else "normal")
    ax.set_xticks(range(len(ms)))
    ax.set_xticklabels([mlabel(m) for m in ms], rotation=50, ha="right", fontsize=8.5)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontsize=10)
    ax.set_xticks(np.arange(-0.5, len(ms), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.2)
    ax.tick_params(which="minor", length=0)
    ax.axhline(len(genes) - 0.5, color="black", lw=2)
    _groups(ax, ms, len(rows))
    ax.set_title("How often is GeneArt better than Codon Harmony?\n"
                 "Cell = % of the Codon Harmony rows of the same gene that the GeneArt sequence beats on that metric",
                 fontsize=11, fontweight="bold", loc="left", pad=48)
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.05, pad=0.30, aspect=40)
    cb.set_ticks([0, 25, 50, 75, 100])
    cb.set_ticklabels(["0\nGeneArt worse than\nall Codon Harmony rows", "25", "50\ntie / middle", "75",
                       "100\nGeneArt better than\nall Codon Harmony rows"], fontsize=8.5)
    note = ("Better/worse follows codon_compare: \u2191 higher wins, \u2193 lower wins, \u25CE GC% closer to {:g}% wins. Ties count "
            "half, so a metric where both sides have the same value in every row shows 50. Metrics with no better/worse rule "
            "(ENC, tAI_ramp, minmax_mean, codon_autocorr_lag1, GC3_pct) are not shown here (see the cluster heatmap)."
            ).format(BAND_MID)
    fig.text(0.01, 0.01, textwrap.fill(note, 200), fontsize=8, color="#333", va="bottom")
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def row_winners(x, v, direction):
    """+1 = Codon Harmony row wins, -1 = GeneArt wins, 0 = draw; same rule as codon_compare.winner_of."""
    if direction == "higher":
        d = x - v
    elif direction == "lower":
        d = v - x
    else:                                                    # band
        d = abs(v - BAND_MID) - np.abs(x - BAND_MID)
    return np.where(d > 1e-9, 1, np.where(d < -1e-9, -1, 0))


def plot_win_rate(wc, scored, path, dpi):
    tot = wc.groupby("metric")[["CodonHarmony_wins", "draws", "GeneArt_wins"]].sum().loc[scored]
    pct = 100 * tot.div(tot.sum(axis=1), axis=0).sort_values("CodonHarmony_wins")
    fig, ax = plt.subplots(figsize=(11, 0.42 * len(pct) + 2.4))
    y = np.arange(len(pct))
    ax.barh(y, pct["CodonHarmony_wins"], color=C_DC, label="Codon Harmony row better")
    ax.barh(y, pct["draws"], left=pct["CodonHarmony_wins"], color="#BDBDBD", label="draw")
    ax.barh(y, pct["GeneArt_wins"], left=pct["CodonHarmony_wins"] + pct["draws"], color=C_GA, label="GeneArt better")
    ax.set_yticks(y)
    ax.set_yticklabels([mlabel(m) for m in pct.index], fontsize=9)
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of all pairs  (pair = one Codon Harmony row vs the GeneArt sequence of the same gene; n = {:,})"
                  .format(int(tot.iloc[0].sum())))
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False)
    ax.set_title("Win rate per scored metric over all pairs\nWinner rule of codon_compare: \u2191 higher, \u2193 lower, "
                 "\u25CE GC% closer to {:g}%; ties are draws".format(BAND_MID), fontsize=11.5, fontweight="bold", loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_tally(wins, genes, path, dpi):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.2), gridspec_kw={"width_ratios": [1, 1.3]})
    x = np.arange(len(genes))
    mean = {k: [float((w == v).sum(axis=1).mean()) for w in (wins[g] for g in genes)] for k, v in
            (("B", 1), ("D", 0), ("A", -1))}
    ax1.bar(x, mean["B"], color=C_DC, label="Codon Harmony better")
    ax1.bar(x, mean["D"], bottom=mean["B"], color="#BDBDBD", label="draw")
    ax1.bar(x, mean["A"], bottom=np.array(mean["B"]) + np.array(mean["D"]), color=C_GA, label="GeneArt better")
    ax1.set_xticks(x)
    ax1.set_xticklabels(genes, rotation=40, ha="right")
    ax1.set_ylabel("mean number of scored metrics per pair")
    ax1.legend(fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3, frameon=False)
    ax1.set_title("Mean scorecard per pair", fontsize=11, fontweight="bold", loc="left")
    rng = np.random.RandomState(3)
    for i, g in enumerate(genes):
        v = (wins[g] == 1).sum(axis=1).to_numpy()
        ax2.scatter(v + rng.uniform(-0.3, 0.3, len(v)), i + rng.uniform(-0.3, 0.3, len(v)), s=3, color=C_DC, alpha=0.2,
                    linewidths=0, rasterized=True)
        ax2.vlines(np.median(v), i - 0.35, i + 0.35, colors="black", linewidth=2.5)
    ax2.set_yticks(range(len(genes)))
    ax2.set_yticklabels(genes)
    ax2.set_ylim(len(genes) - 0.5, -0.5)
    ax2.set_xlabel("number of scored metrics (out of {}) on which the Codon Harmony row beats GeneArt".format(wins[genes[0]].shape[1]))
    ax2.grid(axis="x", alpha=0.3)
    ax2.set_title("Spread over Codon Harmony rows (one dot = one pair, bar = median)", fontsize=11, fontweight="bold", loc="left")
    fig.suptitle("Scorecard tally per gene (all pairs)", fontsize=13, fontweight="bold", x=0.01, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_scatter(dna, ga, genes, path, dpi):
    fig, axs = plt.subplots(1, 2, figsize=(14, 5.6))
    for ax, (xm, ym) in zip(axs, [("GC_pct", "CAI"), ("GC3_pct", "tAI")]):
        draw_refs_x = REF_SPANS.get(xm, [])
        for x0, x1, _ in draw_refs_x:
            ax.axvspan(x0, x1, color="#2E7D32", alpha=0.07, zorder=0)
        for gi, g in enumerate(genes):
            d = dna[dna["Name"] == g]
            ax.scatter(d[xm], d[ym], s=5, color=GENE_COL[gi % 8], alpha=0.25, linewidths=0, rasterized=True)
            r = ga[ga["Name"] == g].iloc[0]
            ax.scatter([r[xm]], [r[ym]], marker="*", s=230, color=GENE_COL[gi % 8], edgecolor="black", linewidths=1.2,
                       zorder=5)
        for y, lab, col in [(l[0], l[1], l[2]) for l in REF_LINES.get(ym, [])]:
            ax.axhline(y, color=col, ls="--", lw=1)
        ax.set_xlabel(xm)
        ax.set_ylabel(ym)
        ax.grid(alpha=0.3)
        ax.set_title("{} vs {}".format(ym, xm), fontsize=11, fontweight="bold", loc="left")
    h = [Line2D([], [], marker="o", ls="", color=GENE_COL[i % 8], label=g) for i, g in enumerate(genes)]
    h += [Line2D([], [], marker="o", ls="", color="#888", alpha=0.5, markersize=4, label="dot = Codon Harmony row"),
          Line2D([], [], marker="*", ls="", color="white", markeredgecolor="black", markersize=13,
                 label="star = GeneArt (gene colour)")]
    fig.legend(handles=h, loc="lower center", ncol=5, frameon=False, fontsize=9)
    fig.suptitle("Codon optimality against GC: every Codon Harmony row and the GeneArt sequence of each gene",
                 fontsize=12.5, fontweight="bold", x=0.01, ha="left")
    fig.tight_layout(rect=[0, 0.1, 1, 0.95])
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def analyse(dna, ga, metrics, genes, k):
    stats, clus = [], []
    rowcl = dna[["Name", "Variant"]].copy()
    for m in metrics:
        rowcl[m + "_cluster"] = ""
        for g in genes:
            idx = dna.index[dna["Name"] == g]
            x = dna.loc[idx, m].to_numpy(float)
            gv = float(ga.loc[ga["Name"] == g, m].iloc[0])
            lab, cen, gr = cluster_1d(x, gv, k)
            kk = len(cen)
            names = NAMES[kk]
            rowcl.loc[idx, m + "_cluster"] = [names[i] for i in lab]
            bounds = [(cen[i] + cen[i + 1]) / 2 for i in range(kk - 1)]
            gname = names[gr]
            if kk == 1:      # no spread in Codon Harmony: say where GeneArt is relative to that single value
                gname = "Same value" if abs(gv - x[0]) <= EPS else ("Above range" if gv > x[0] else "Below range")
            for i in range(kk):
                xs = x[lab == i]
                clus.append({"metric": m, "gene": g, "cluster": names[i], "rank": i, "n_clusters_used": kk,
                             "n_rows": len(xs), "pct_of_gene_rows": round(100 * len(xs) / len(x), 2),
                             "min": xs.min(), "max": xs.max(), "centre": cen[i],
                             "GeneArt_here": i == gr})
            stats.append({"metric": m, "gene": g, "better": DIRECTION[m], "GeneArt": gv,
                          "CodonHarmony_n": len(x), "CodonHarmony_min": x.min(), "CodonHarmony_median": float(np.median(x)),
                          "CodonHarmony_max": x.max(), "GeneArt_percentile": round(pct_below(x, gv), 2),
                          "GeneArt_beats_pct": round(beat_share(x, gv, DIRECTION[m]), 2),
                          "GeneArt_cluster": gname, "n_clusters_used": kk,
                          "GeneArt_outside_CodonHarmony_range": bool(gv < x.min() - EPS or gv > x.max() + EPS),
                          "_code": heat_code(gr, kk, k), "_bounds": bounds})
    return pd.DataFrame(stats), pd.DataFrame(clus), rowcl


PANELS = [
    ("gc", "GC % (45-nt window)", "nt", ("nt_window", "gc")),
    ("mfe", "MFE dG, 40-nt window (kcal/mol)", "nt", ("mfe_window",)),
    ("minmax", "%MinMax (18-codon window)", "codon", ("codon_window", "minmax")),
    ("rolling_cai", "Rolling CAI (18 codons)", "codon", ("codon_window", "rolling_cai")),
    ("rolling_tai", "Rolling tAI (18 codons)", "codon", ("codon_window", "rolling_tai")),
    ("rolling_rare", "Rolling rare-codon % (18 codons)", "codon", ("codon_window", "rolling_rare_pct")),
    ("gc3", "GC3 % (20-codon window)", "codon", ("codon_window", "gc3_window")),
    ("enc", "Windowed ENC (60 codons)", "codon", ("codon_window", "enc_window")),
    ("cpg", "CpG density /100 bp (50-nt window)", "nt", ("nt_window", "cpg_density")),
    ("hairpin", "Hairpin density (60-nt window)", "nt", ("nt_window", "hairpin_density")),
]
EXTRA_DIR = {"GC1": "na", "GC2": "na", "codon_usage_L1_to_host": "dist", "GC_pos_L1_to_host": "dist",
             "ENC_abs_diff_to_host": "dist", "dinuc_total_L1_to_host": "dist", "dinuc_junction_L1_to_host": "dist"}


EXTRA_DIR = {"GC1": "na", "GC2": "na", "ENC_deviation": "na",
             "codon_usage_L1_to_host": "dist", "GC_pos_L1_to_host": "dist", "ENC_abs_diff_to_host": "dist",
             "dinuc_total_L1_to_host": "dist", "dinuc_junction_L1_to_host": "dist"}


# remaining metrics: Kozak + the baseline-audit metrics that scores_codonharmony.xlsx does not hold + other scalars of the scorer
MORE_DIR = {"kozak_score": "higher", "GC_min_60bp": "higher", "GC_std_60bp": "lower", "CpG_count": "lower",
            "longest_repeat_bp": "lower", "AT_pct": "na", "tAI_head": "na", "tAI_body": "na", "minmax_min": "na",
            "minmax_pct_rare_windows": "na", "MFE_mean_window": "na", "tm_mean": "na"}
REF_LINES.update({"GC_min_60bp": [(30.0, "audit: < 30% DNA instability", "#6A1B9A")],
                  "GC_std_60bp": [(3.0, "audit: < 3% smooth", "#6A1B9A")],
                  "longest_repeat_bp": [(15.0, "audit: > 15 bp recombination risk", "#6A1B9A")]})


def longest_repeat(seq, kmin=8):
    """Length of the longest exact repeat (a substring occurring at least twice, overlaps allowed); 0 if none >= kmin."""
    def has(k):
        seen = set()
        for i in range(len(seq) - k + 1):
            w = seq[i:i + k]
            if w in seen:
                return True
            seen.add(w)
        return False
    if len(seq) < 2 * kmin or not has(kmin):
        return 0
    lo, hi = kmin, len(seq) - 1            # has(lo) is true; find the largest k with has(k)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if has(mid):
            lo = mid
        else:
            hi = mid - 1
    return lo


def more_metrics(seq, scores):
    s = cf.clean_seq(seq)
    w = np.array([cf.gc_content(s[i:i + 60]) for i in range(0, len(s) - 60 + 1, 60)], float)
    gc = cf.gc_content(s)
    out = {"GC_min_60bp": float(w.min()) if len(w) else np.nan, "GC_std_60bp": float(w.std()) if len(w) > 1 else np.nan,
           "CpG_count": s.count("CG"), "longest_repeat_bp": longest_repeat(s), "AT_pct": 100.0 - gc}
    for k in ("kozak_score", "tAI_head", "tAI_body", "minmax_min", "minmax_pct_rare_windows", "MFE_mean_window", "tm_mean"):
        v = scores.get(k)
        out[k] = np.nan if v is None else float(v)
    return out


# --------------------------------------------------------------------------------------------------
# scoring of the DNA (once per distinct sequence, in worker processes) + cache
# --------------------------------------------------------------------------------------------------
def _job(seq):
    r = cf.compute_features(seq, host="human", host_label="human", codon_table=str(sas.CODON_CSV),
                            tgcn=str(sas.TGCN_CSV))
    P = r["profiles"]
    prof = {}
    for key, _t, _x, path in PANELS:
        d = P
        for k in path:
            d = d[k]
        prof[key] = (np.asarray(d["index"], np.float32), np.asarray(d["value"], np.float32)) if d.get("value") else None
    a = P["autocorrelation"]
    prof["acf"] = (np.asarray(a["lag"], np.float32), np.asarray(a["acf"], np.float32)) if a.get("acf") else None
    t = P["tm_fragments"]
    prof["tm"] = (np.asarray(t["index"], np.float32), np.asarray(t["tm"], np.float32)) if t and t.get("tm") else None
    cod = Counter(cf.coding_codons(seq))
    td, jd = rc.dinuc_per1000(seq)
    extra = {"cod": np.array([cod.get(c, 0) for c in rc.SENSE], np.float32),
             "gc": rc.pos_gc(seq),
             "dt": np.array([td[d] for d in rc.DINUC], np.float32),
             "dj": np.array([jd[d] for d in rc.DINUC], np.float32)}
    return seq, {"scores": r["scores"], "prof": prof, "pos": r["positions"], "extra": extra, "length": len(seq)}


def fingerprint():
    return [(f.name, f.stat().st_size, int(f.stat().st_mtime)) for f in (sas.CODON_CSV, sas.TGCN_CSV, Path(cf.__file__))
            if Path(f).exists()]


def score_all(seqs, workers, cache_path):
    data, fp = {}, fingerprint()
    if cache_path.exists():
        try:
            with open(cache_path, "rb") as fh:
                c = pickle.load(fh)
            if c.get("fp") == fp:
                data = c["data"]
            else:
                print("cache is from different reference/scorer files: scoring again")
        except Exception:
            pass
    todo = [s for s in seqs if s not in data]
    print("distinct sequences: {} | cached: {} | to score: {}".format(len(seqs), len(seqs) - len(todo), len(todo)))
    if todo:
        try:
            from tqdm import tqdm
        except Exception:
            tqdm = lambda x, **k: x                                    # noqa: E731
        if workers <= 1:
            for s in tqdm(todo, desc="Scoring", unit="seq"):
                data[s] = _job(s)[1]
        else:
            with ProcessPoolExecutor(max_workers=workers) as ex:
                for s, r in tqdm(ex.map(_job, todo, chunksize=4), total=len(todo), desc="Scoring", unit="seq"):
                    data[s] = r
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, "wb") as fh:
            pickle.dump({"fp": fp, "data": data}, fh, protocol=4)
    return data


def check_against_table(dna, ga, data, metrics):
    """Scalar scores returned by the scorer here must equal the scores_codonharmony.xlsx values."""
    rows = []
    for m in metrics:
        diffs = []
        for seq, v in list(zip(dna["_seq"], dna[m])) + list(zip(ga["_seq"], ga[m])):
            x = data[seq]["scores"].get(m)
            if x is None or pd.isna(v):
                diffs.append(0.0 if (x is None and pd.isna(v)) else np.inf)
            else:
                diffs.append(abs(float(x) - float(v)))
        d = np.array(diffs)
        rows.append({"metric": m, "sequences_compared": len(d), "differences_gt_1e-6": int((d > 1e-6).sum()),
                     "max_abs_difference": float(np.max(d))})
    return pd.DataFrame(rows)


def expand(res, name):
    """compact cache entry -> the structure codon_compare / codon_features plot functions expect."""
    P = {"codon_window": {}, "nt_window": {}}
    for key, _t, _x, path in PANELS:
        v = res["prof"][key]
        d = {"index": [], "value": []} if v is None else {"index": v[0].tolist(), "value": v[1].tolist()}
        if len(path) == 1:
            P[path[0]] = d
        else:
            P[path[0]][path[1]] = d
    a, t = res["prof"]["acf"], res["prof"]["tm"]
    P["autocorrelation"] = {"lag": [] if a is None else a[0].tolist(), "acf": [] if a is None else a[1].tolist()}
    P["tm_fragments"] = {"index": [] if t is None else t[0].tolist(), "tm": [] if t is None else t[1].tolist()}
    return {"name": name, "profiles": P, "scores": res["scores"], "positions": res["pos"]}


def pair_rows(sa, sb):
    """Same rows as codon_compare.compare() builds, from scores that are already computed."""
    band = cc.gc_band_for("human")
    rows = []
    for metric, direction in cc.COMPARE_METRICS:
        av, bv = sa.get(metric), sb.get(metric)
        try:
            delta = round(float(bv) - float(av), 4)
        except (TypeError, ValueError):
            delta = None
        w = cc.winner_of(metric, av, bv, direction, band)
        rows.append({"metric": metric, LA: av, LB: bv, "delta": delta, "better": direction,
                     "winner": {"A": LA, "B": LB, "draw": "draw", "na": "-"}[w], "_w": w})
    return rows


def _aligned(arrs, ref_idx):
    out = np.full((len(arrs), len(ref_idx)), np.nan, np.float32)
    for i, (idx, val) in enumerate(arrs):
        if len(idx) == len(ref_idx) and np.array_equal(idx, ref_idx):
            out[i] = val
        else:
            _, a, b = np.intersect1d(ref_idx, idx, return_indices=True)
            out[i, a] = val[b]
    return out


def plot_motif_map(genes, ga_res, rows, data, path, dpi, bin_nt=30):
    cats = [("restriction", "Restriction site"), ("splice", "Splice motif"), ("homopolymer", "Homopolymer \u22656")]
    fig, axs = plt.subplots(len(genes), 1, figsize=(14, 2.0 * len(genes) + 1.6))
    axs = np.atleast_1d(axs)
    for ax, g in zip(axs, genes):
        seqs = rows[g]
        L = ga_res[g]["length"]
        nb = int(np.ceil(L / bin_nt))
        M = np.zeros((3, nb))
        anyhit = [0, 0, 0]
        for s in seqs:
            for ci, (k, _) in enumerate(cats):
                hits = data[s]["pos"][k]
                if hits:
                    anyhit[ci] += 1
                    for b in {int(h["pos"]) // bin_nt for h in hits}:
                        M[ci, min(b, nb - 1)] += 1
        M = 100.0 * M / len(seqs)
        ax.imshow(M, cmap="Blues", vmin=0, vmax=100, aspect="auto", extent=[0, nb * bin_nt, 2.5, -0.5])
        gp = ga_res[g]["pos"]
        counts = []
        for ci, (k, _) in enumerate(cats):
            hits = gp[k]
            counts.append(len(hits))
            if hits:
                ax.scatter([h["pos"] for h in hits], [ci] * len(hits), marker="D", s=46, color=C_GA, edgecolor="white",
                           linewidths=0.8, zorder=5)
            if k == "restriction":
                for h in hits:
                    ax.annotate(h["name"], (h["pos"], ci), xytext=(0, 8), textcoords="offset points", rotation=60,
                                fontsize=6.5, color=C_GA, ha="left", va="bottom")
        ax.set_yticks(range(3))
        ax.set_yticklabels([c[1] for c in cats], fontsize=8.5)
        ax.set_xlim(0, nb * bin_nt)
        ax.set_title("{}  ({:,} nt):  GeneArt has {} restriction / {} splice / {} homopolymer hits;  Codon Harmony rows with at "
                     "least one: {:.0f}% / {:.0f}% / {:.0f}%  (n = {:,})".format(
                         g, L, counts[0], counts[1], counts[2], *[100.0 * a / len(seqs) for a in anyhit], len(seqs)),
                     fontsize=9.5, loc="left", fontweight="bold", pad=14)
        ax.set_xlabel("nt position" if g == genes[-1] else "")
    fig.tight_layout(rect=[0, 0.06, 0.93, 0.96])
    cax = fig.add_axes([0.945, 0.2, 0.012, 0.6])
    cb = fig.colorbar(ScalarMappable(Normalize(0, 100), "Blues"), cax=cax)
    cb.set_label("% of the gene's Codon Harmony rows with a hit in that {}-nt bin".format(bin_nt), fontsize=8.5)
    fig.legend(handles=[Line2D([], [], marker="D", ls="", color=C_GA, markeredgecolor="white", markersize=8,
                               label="GeneArt hit (restriction: enzyme name; either strand, start position)")],
               loc="lower center", frameon=False, fontsize=9)
    fig.suptitle("Where the sequence liabilities sit along each CDS", fontsize=13, fontweight="bold", x=0.01, ha="left")
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_host_codon(host, ga_vec, dc_vec, path, dpi):
    order = [c for c in [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"] if c in rc.SENSE]
    order.sort(key=lambda c: (cf.GENETIC_CODE[c], c))
    pos = {c: i for i, c in enumerate(rc.SENSE)}
    per = lambda v: 1000.0 * v / v.sum()                                # noqa: E731
    g, d = per(ga_vec), per(dc_vec)
    x = np.arange(len(order))
    fig, ax = plt.subplots(figsize=(19, 5.4))
    ax.bar(x - 0.28, [host["per1000"][c] for c in order], 0.27, color=C_HOST, label="Human host (CoCoPUTs, Human.txt)")
    ax.bar(x, [g[pos[c]] for c in order], 0.27, color=C_GA, label="GeneArt (8 genes pooled)")
    ax.bar(x + 0.28, [d[pos[c]] for c in order], 0.27, color=C_DC,
           label="Codon Harmony (mean over the rows of each gene, genes pooled)")
    ax.set_xticks(x)
    ax.set_xticklabels(["{}\n{}".format(cf.GENETIC_CODE[c], c) for c in order], fontsize=6.5)
    ax.set_ylabel("per 1000 sense codons")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.3)
    ax.set_title("Codon usage per 1000 codons, grouped by amino acid", fontsize=12, fontweight="bold", loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_host_dinuc(hd, ga_arr, dc_arr, path, dpi):
    fig, axs = plt.subplots(2, 1, figsize=(13, 8.5))
    for ax, k, ttl in zip(axs, (0, 1), ("Total dinucleotides", "Junction dinucleotides (3rd base of a codon + 1st base of the next)")):
        x = np.arange(16)
        ga, dc = ga_arr[k], dc_arr[k]
        ax.bar(x - 0.28, [hd[k][d] for d in rc.DINUC], 0.27, color=C_HOST, label="Human host")
        ax.bar(x, ga.mean(axis=0), 0.27, color=C_GA, label="GeneArt (mean of the genes; dots = each gene)")
        ax.scatter(np.repeat(x, len(ga)), ga.T.ravel(), s=9, color="black", zorder=4)
        med = np.median(dc, axis=0)
        q1, q3 = np.percentile(dc, [25, 75], axis=0)
        ax.bar(x + 0.28, dc.mean(axis=0), 0.27, color=C_DC, label="Codon Harmony (mean over all rows; whiskers = quartiles)")
        ax.errorbar(x + 0.28, med, yerr=[np.maximum(med - q1, 0), np.maximum(q3 - med, 0)], fmt="none", ecolor="#222",
                    capsize=2, lw=0.9, zorder=4)
        ax.set_xticks(x)
        ax.set_xticklabels(rc.DINUC)
        ax.set_ylabel("per 1000")
        ax.grid(axis="y", alpha=0.3)
        ax.set_title(ttl, fontsize=11, fontweight="bold", loc="left")
        ax.legend(fontsize=8.5, frameon=False)
    fig.suptitle("Dinucleotide frequencies per 1000 (host values read from the CoCoPUTs chart, normalised to 1000, +/-0.5)",
                 fontsize=11, fontweight="bold", x=0.01, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------------------------------
# positional profiles (one figure per gene)
# --------------------------------------------------------------------------------------------------
def plot_gene_profiles(gene, ga_res, seqs, data, path, dpi):
    panels = [(k, t, x) for k, t, x, _ in PANELS] + [("acf", "Autocorrelation of per-codon tAI (ACF)", "lag"),
                                                    ("tm", "Fragment Tm (60-nt fragments, \u00b0C)", "fragment")]
    keep = [p for p in panels if ga_res["prof"][p[0]] is not None]
    ncol = 2
    nrow = int(np.ceil(len(keep) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(16, 2.7 * nrow + 1.8))
    axs = np.atleast_1d(axs).ravel()
    refs = set()
    for ax, (key, title, xl) in zip(axs, keep):
        gi, gv = ga_res["prof"][key]
        arrs = [data[s]["prof"][key] for s in seqs if data[s]["prof"][key] is not None]
        if not arrs:
            ax.axis("off")
            continue
        st = _aligned(arrs, gi)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p5, p25, p50, p75, p95 = np.nanpercentile(st, [5, 25, 50, 75, 95], axis=0)
        if key in ("acf", "tm"):
            ax.bar(gi - 0.2, gv, 0.4, color=C_GA)
            ax.bar(gi + 0.2, p50, 0.4, color=C_DC, yerr=[p50 - p25, p75 - p50], ecolor="#222", capsize=2)
        else:
            ax.fill_between(gi, p5, p95, color=C_DC, alpha=0.15, linewidth=0)
            ax.fill_between(gi, p25, p75, color=C_DC, alpha=0.30, linewidth=0)
            ax.plot(gi, p50, color=C_DC, lw=1.2)
            ax.plot(gi, gv, color=C_GA, lw=1.7)
        if key == "gc":                      # human GC band used by codon_compare
            ax.axhspan(BAND[0], BAND[1], color="#2E7D32", alpha=0.07, zorder=0)
            refs.add("band")
        if key == "rolling_cai":
            ax.axhline(0.85, color="#6A1B9A", ls="--", lw=1)
            refs.add("cai")
        if key == "rolling_rare":
            ax.axhline(5.0, color="#6A1B9A", ls="--", lw=1)
            refs.add("rare")
        if key == "mfe":
            ax.axhline(cf.DEFAULTS["stable_dg_cutoff"], color="#555", ls="--", lw=1)
            ax.text(0.995, cf.DEFAULTS["stable_dg_cutoff"], "stable-structure cutoff", transform=ax.get_yaxis_transform(),
                    ha="right", va="bottom", fontsize=7, color="#555")
        ax.set_title(title, fontsize=10, loc="left", fontweight="bold")
        ax.set_xlabel(xl)
        ax.grid(alpha=0.3)
    for ax in axs[len(keep):]:
        ax.axis("off")
    h = [Line2D([], [], color=C_GA, lw=2, label="GeneArt"),
         Line2D([], [], color=C_DC, lw=1.5, label="Codon Harmony median"),
         Patch(color=C_DC, alpha=0.30, label="Codon Harmony 25\u201375%"),
         Patch(color=C_DC, alpha=0.15, label="Codon Harmony 5\u201395%")]
    if "band" in refs:
        h.append(Patch(color="#2E7D32", alpha=0.2, label="human GC band {:g}\u2013{:g}%".format(*BAND)))
    if "cai" in refs:
        h.append(Line2D([], [], color="#6A1B9A", ls="--", label="audit: CAI > 0.85"))
    if "rare" in refs:
        h.append(Line2D([], [], color="#6A1B9A", ls="--", label="audit: rare codons > 5% problematic"))
    top = 0.935
    fig.suptitle("{}: positional profiles, GeneArt vs all {:,} Codon Harmony rows".format(gene, len(seqs)),
                 fontsize=13, fontweight="bold", x=0.01, ha="left", y=0.995)
    fig.text(0.01, 0.962, "Red line = GeneArt. Blue line = median of the Codon Harmony rows, bands = 25-75% and 5-95% of the rows "
             "at each position. Bars in the ACF and Tm panels: GeneArt (red) and Codon Harmony median with quartiles (blue).",
             fontsize=8.5, color="#444", va="center")
    fig.legend(handles=h, loc="upper left", ncol=7, frameon=False, fontsize=8.5, bbox_to_anchor=(0.005, top + 0.003))
    fig.tight_layout(rect=[0, 0, 1, top - 0.012])
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------------------------------
# win / draw / loss heatmap (pair = one Codon Harmony row vs the GeneArt sequence of the same gene)
# --------------------------------------------------------------------------------------------------
def plot_win_heatmap(wc, genes, scored, path, dpi):
    rows = list(genes) + ["All genes"]
    M = np.zeros((len(rows), len(scored)))
    T = [[""] * len(scored) for _ in rows]
    for j, m in enumerate(scored):
        sub = wc[wc["metric"] == m].set_index("gene")
        for i, g in enumerate(rows):
            s = sub.loc[g] if g != "All genes" else sub[["pairs", "CodonHarmony_wins", "GeneArt_wins", "draws"]].sum()
            n = float(s["pairs"])
            dc, ga_, dr = (100.0 * float(s[k]) / n for k in ("CodonHarmony_wins", "GeneArt_wins", "draws"))
            M[i, j] = dc - ga_
            T[i][j] = "{:.0f}\n{:.0f}\n{:.0f}".format(dc, dr, ga_)
    cmap = plt.get_cmap("RdBu")                  # red = GeneArt ahead, blue = Codon Harmony ahead
    norm = Normalize(-100, 100)
    fig, ax = plt.subplots(figsize=(0.82 * len(scored) + 3.4, 0.95 * len(rows) + 4.2))
    im = ax.imshow(M, cmap=cmap, norm=norm, aspect="auto")
    for i in range(len(rows)):
        for j in range(len(scored)):
            rgb = cmap(norm(M[i, j]))[:3]
            dark = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2] < 0.5
            ax.text(j, i, T[i][j], ha="center", va="center", fontsize=7.5, color="white" if dark else "black")
    ax.set_xticks(range(len(scored)))
    ax.set_xticklabels([mlabel(m) for m in scored], rotation=50, ha="right", fontsize=8.5)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontsize=10)
    ax.set_xticks(np.arange(-0.5, len(scored), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.2)
    ax.tick_params(which="minor", length=0)
    ax.axhline(len(genes) - 0.5, color="black", lw=2)
    _groups(ax, scored, len(rows))
    ax.set_title("Who wins each metric?  Pair = one Codon Harmony row against the GeneArt sequence of the same gene\n"
                 "Cell text, top to bottom: % of pairs Codon Harmony better / % draw / % GeneArt better.  "
                 "Colour = (Codon Harmony % - GeneArt %)", fontsize=11, fontweight="bold", loc="left", pad=48)
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.05, pad=0.32, aspect=40)
    cb.set_ticks([-100, -50, 0, 50, 100])
    cb.set_ticklabels(["-100\nGeneArt better\nin every pair", "-50", "0\neven or draws", "50",
                       "100\nCodon Harmony better\nin every pair"], fontsize=8.5)
    note = ("Winner rule of codon_compare: \u2191 higher wins, \u2193 lower wins, \u25CE GC% closer to {:g}% wins; ties are draws. "
            "Metrics without a better/worse rule (ENC, tAI_ramp, minmax_mean, codon_autocorr_lag1, GC3_pct) are not scored.").format(BAND_MID)
    fig.text(0.01, 0.01, textwrap.fill(note, 190), fontsize=8, color="#333", va="bottom")
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------------------------------
# ENC against GC3 (Wright's expected curve) and GC by codon position
# --------------------------------------------------------------------------------------------------
def plot_enc_gc3(dna, ga, genes, host, path, dpi):
    fig, ax = plt.subplots(figsize=(10, 6.4))
    s = np.linspace(0.02, 0.98, 200)
    ax.plot(100 * s, [cf.enc_expected_from_gc3(v) for v in s], color="#333", lw=1.6, zorder=1)
    for gi, g in enumerate(genes):
        d = dna[dna["Name"] == g]
        ax.scatter(d["GC3_pct"], d["ENC"], s=5, color=GENE_COL[gi % 8], alpha=0.25, linewidths=0, rasterized=True)
        r = ga[ga["Name"] == g].iloc[0]
        ax.scatter([r["GC3_pct"]], [r["ENC"]], marker="*", s=240, color=GENE_COL[gi % 8], edgecolor="black", linewidths=1.2, zorder=5)
    if host:
        ax.scatter([host["gc3"]], [rc.HOST_ENC], marker="X", s=150, color="black", zorder=6)
        ax.annotate("human host (GC3 {:.0f}%, ENc {:.1f})".format(host["gc3"], rc.HOST_ENC), (host["gc3"], rc.HOST_ENC),
                    xytext=(8, 8), textcoords="offset points", fontsize=9)
    ax.set_xlabel("GC3 (%)")
    ax.set_ylabel("ENC (Nc): 20 = strongest codon bias, 61 = none")
    ax.set_ylim(15, 63)
    ax.grid(alpha=0.3)
    h = [Line2D([], [], color="#333", lw=1.6, label="Nc expected from GC3 alone (Wright 1990): points below = more bias than GC3 explains")]
    h += [Line2D([], [], marker="o", ls="", color=GENE_COL[i % 8], label=g) for i, g in enumerate(genes)]
    h += [Line2D([], [], marker="o", ls="", color="#888", alpha=0.5, markersize=4, label="dot = Codon Harmony row"),
          Line2D([], [], marker="*", ls="", color="white", markeredgecolor="black", markersize=13, label="star = GeneArt (gene colour)")]
    if host:
        h.append(Line2D([], [], marker="X", ls="", color="black", markersize=9, label="human host"))
    fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, fontsize=8.5)
    fig.suptitle("Codon bias (ENC) against GC3: every Codon Harmony row and the GeneArt sequence of each gene",
                 fontsize=12, fontweight="bold", x=0.01, ha="left")
    fig.tight_layout(rect=[0, 0.13, 1, 0.96])
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_host_gc(host, dna, ga, ex_dna, ex_ga, path, dpi):
    names = ["GC overall", "GC1", "GC2", "GC3"]
    dcv = [dna["GC_pct"].to_numpy(float), ex_dna["GC1"].to_numpy(float), ex_dna["GC2"].to_numpy(float), dna["GC3_pct"].to_numpy(float)]
    gav = [ga["GC_pct"].to_numpy(float), ex_ga["GC1"].to_numpy(float), ex_ga["GC2"].to_numpy(float), ga["GC3_pct"].to_numpy(float)]
    hv = [host["gc"], host["gc1"], host["gc2"], host["gc3"]]
    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    x = np.arange(4)
    ax.bar(x - 0.28, hv, 0.27, color=C_HOST, label="Human host (Human.txt)")
    ax.bar(x, [v.mean() for v in gav], 0.27, color=C_GA, label="GeneArt (mean of the genes; black dots = each gene)")
    for i, v in enumerate(gav):
        ax.scatter(np.full(len(v), i), v, s=14, color="black", zorder=4)
    med = [np.median(v) for v in dcv]
    q = [np.percentile(v, [25, 75]) for v in dcv]
    ax.bar(x + 0.28, [v.mean() for v in dcv], 0.27, color=C_DC, label="Codon Harmony (mean of all rows; whiskers = quartiles)")
    ax.errorbar(x + 0.28, med, yerr=[[m - a for m, (a, b) in zip(med, q)], [b - m for m, (a, b) in zip(med, q)]],
                fmt="none", ecolor="#222", capsize=3, lw=1, zorder=4)
    ax.axhspan(BAND[0], BAND[1], xmin=0, xmax=0.25, color="#2E7D32", alpha=0.08, zorder=0)
    ax.text(-0.45, BAND[1] + 0.6, "human GC band {:g}\u2013{:g}% (overall GC only)".format(*BAND), fontsize=8, color="#2E7D32")
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("GC %")
    ax.set_ylim(0, max(85, max(float(np.max(v)) for v in gav + dcv) + 6))
    ax.legend(frameon=False, fontsize=9, loc="upper left", bbox_to_anchor=(0, -0.08), ncol=1)
    ax.grid(axis="y", alpha=0.3)
    ax.set_title("GC by codon position against the human host", fontsize=12, fontweight="bold", loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------------------------------
# liabilities: restriction sites per enzyme
# --------------------------------------------------------------------------------------------------
def restriction_table(genes, ga_res, rows, data):
    out = []
    for e in cf.RESTRICTION_SITES:
        for g in genes:
            cnt = np.array([data[s]["scores"].get("RS_" + e, 0) or 0 for s in rows[g]], float)
            out.append({"enzyme": e, "site": cf.RESTRICTION_SITES[e], "gene": g,
                        "GeneArt_sites": int(ga_res[g]["scores"].get("RS_" + e, 0) or 0),
                        "CodonHarmony_mean_sites_per_row": float(cnt.mean()),
                        "CodonHarmony_pct_rows_with_site": float(100.0 * (cnt > 0).mean()), "CodonHarmony_rows": len(cnt)})
    return pd.DataFrame(out)


def plot_restriction_enzymes(rt, genes, path, dpi):
    enz = list(cf.RESTRICTION_SITES)
    pv = lambda col: rt.pivot(index="enzyme", columns="gene", values=col).loc[enz, genes]       # noqa: E731
    A, B, Cc = pv("GeneArt_sites"), pv("CodonHarmony_mean_sites_per_row"), pv("CodonHarmony_pct_rows_with_site")
    hit = [e for e in enz if A.loc[e].any() or B.loc[e].any()]
    fig_h = 0.34 * max(len(hit), 2) + 3.4
    fig, axs = plt.subplots(1, 3, figsize=(15, fig_h))
    if not hit:
        for ax in axs:
            ax.axis("off")
        fig.text(0.5, 0.5, "None of the {} enzymes of the scorer has a site in any GeneArt or Codon Harmony sequence".format(len(enz)),
                 ha="center", fontsize=12)
    else:
        specs = [(A, "GeneArt: sites in the sequence", "{:.0f}", None),
                 (B, "Codon Harmony: mean sites per row", "{:.2f}", None),
                 (Cc, "Codon Harmony: % of rows with at least one site", "{:.0f}", (0, 100))]
        for ax, (M, ttl, f, lim) in zip(axs, specs):
            Mh = M.loc[hit].to_numpy(float)
            vmax = lim[1] if lim else max(float(Mh.max()), 1.0)
            ax.imshow(Mh, cmap="Blues", vmin=0, vmax=vmax, aspect="auto")
            for i in range(Mh.shape[0]):
                for j in range(Mh.shape[1]):
                    ax.text(j, i, f.format(Mh[i, j]), ha="center", va="center", fontsize=8,
                            color="white" if Mh[i, j] > 0.6 * vmax else "black")
            ax.set_xticks(range(len(genes)))
            ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=9)
            ax.set_yticks(range(len(hit)))
            ax.set_yticklabels(["{} ({})".format(e, cf.RESTRICTION_SITES[e]) for e in hit] if ax is axs[0] else [], fontsize=8.5)
            ax.set_title(ttl, fontsize=10, fontweight="bold", loc="left")
    none = [e for e in enz if e not in hit]
    fig.suptitle("Restriction sites per enzyme (both strands; enzymes of the scorer's list)", fontsize=12.5,
                 fontweight="bold", x=0.01, ha="left")
    fig.text(0.01, 0.01, textwrap.fill("Only enzymes with at least one site in any sequence are listed. No site in any sequence: {}.".format(
        ", ".join(none) if none else "none"), 200), fontsize=8, color="#333", va="bottom")
    fig.tight_layout(rect=[0, 0.05, 1, 0.94])
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------------------------------
# agreement with GeneArt (codon by codon)
# --------------------------------------------------------------------------------------------------
def match_matrix(ga_seq, seqs):
    """Rows = Codon Harmony sequences, columns = codons: 1 where the codon equals the GeneArt codon."""
    g = cf.strip_trailing_stop(ga_seq)
    n = len(g) // 3
    A = np.frombuffer(g.encode(), np.uint8)[:3 * n].reshape(n, 3)
    M = np.zeros((len(seqs), n), np.float32)
    ident = np.zeros(len(seqs))
    for r, s in enumerate(seqs):
        s = cf.strip_trailing_stop(s)
        m = min(n, len(s) // 3)
        B = np.frombuffer(s.encode(), np.uint8)[:3 * m].reshape(m, 3)
        M[r, :m] = (A[:m] == B).all(axis=1)
        ident[r] = 100.0 * (A[:m] == B).mean() if m else np.nan
    return M, ident


def plot_codon_agreement(genes, mats, path, dpi, win=9):
    fig, axs = plt.subplots(len(genes), 1, figsize=(13, 1.9 * len(genes) + 1.6), sharey=True)
    axs = np.atleast_1d(axs)
    for ax, g in zip(axs, genes):
        M = mats[g]
        raw = 100.0 * M.mean(axis=0)
        sm = pd.Series(raw).rolling(win, center=True, min_periods=1).mean().to_numpy()
        x = np.arange(1, len(raw) + 1)
        ax.fill_between(x, raw, color=C_DC, alpha=0.18, linewidth=0)
        ax.plot(x, sm, color=C_DC, lw=1.5)
        ax.axhline(100.0 * M.mean(), color="#444", ls="--", lw=1)
        ax.set_ylim(0, 100)
        ax.set_title("{}: on average {:.0f}% of codons are identical to GeneArt (dashed line); n = {:,} Codon Harmony rows".format(
            g, 100.0 * M.mean(), M.shape[0]), fontsize=9.5, loc="left", fontweight="bold")
        ax.set_ylabel("% rows\nidentical")
        ax.grid(alpha=0.3)
    axs[-1].set_xlabel("codon position")
    fig.suptitle("Where do the Codon Harmony rows use the same codon as GeneArt?  (light = per codon, line = {}-codon average)".format(win),
                 fontsize=12, fontweight="bold", x=0.01, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------------------------------
def pair_results(dna, ga, genes, metrics):
    scored = [m for m in metrics if DIRECTION[m] in ("higher", "lower", "band")]
    wins = {}
    rows = []
    for g in genes:
        idx = dna.index[dna["Name"] == g]
        cols = {}
        for m in scored:
            gv = float(ga.loc[ga["Name"] == g, m].iloc[0])
            cols[m] = row_winners(dna.loc[idx, m].to_numpy(float), gv, DIRECTION[m])
            rows.append({"metric": m, "gene": g, "pairs": len(idx), "CodonHarmony_wins": int((cols[m] == 1).sum()),
                         "GeneArt_wins": int((cols[m] == -1).sum()), "draws": int((cols[m] == 0).sum())})
        wins[g] = pd.DataFrame(cols, index=idx)
    return scored, wins, pd.DataFrame(rows)


def finite_metrics(names, *frames):
    ok, bad = [], []
    for m in names:
        v = np.concatenate([pd.to_numeric(f[m], errors="coerce").to_numpy(float) for f in frames])
        (ok if np.isfinite(v).all() else bad).append(m)
    return ok, bad


def main():
    ap = argparse.ArgumentParser(description="All plots: GeneArt vs every Codon Harmony row.")
    ap.add_argument("--scores", default=str(HERE / "outputs" / "scores" / "scores_codonharmony.xlsx"))
    ap.add_argument("--grid", default=None, help="grid output with the DNA (sheet Results); default: human_codonharmony_grid.xlsx in "
                    "../CodonHarmony/outputs (use the same file as score_codonharmony.py)")
    ap.add_argument("--trials", default=None, help="the ORIGINAL trials workbook (GeneArt sequences)")
    ap.add_argument("--out-dir", default=str(HERE / "outputs" / "plots_codonharmony"))
    ap.add_argument("--clusters", type=int, default=3, choices=[2, 3, 4, 5], help="K-means clusters per metric per gene")
    ap.add_argument("--workers", type=int, default=min(40, max(1, (os.cpu_count() or 2) - 1)))
    ap.add_argument("--pair-plots", type=int, default=10, help="per gene the N best-CAI distinct Codon Harmony sequences get "
                    "codon_compare overlay + scorecard + codon_features 9-panel figure (default 10; 0 = none; -1 = every row, "
                    "tens of thousands of files)")
    ap.add_argument("--no-profiles", action="store_true", help="only the plots that need scores_codonharmony.xlsx (no DNA, no scoring)")
    ap.add_argument("--no-host", action="store_true", help="no host-reference lines / plots")
    ap.add_argument("--max-rows-per-gene", type=int, default=0, help="TEST ONLY: random subset of rows per gene for every plot")
    ap.add_argument("--dpi", type=int, default=130)
    ap.add_argument("--stop-policy", choices=["strip", "asis"], default="strip",
                    help="terminal stop codon handling of the re-scored DNA; use the value used for score_codonharmony.py")
    a = ap.parse_args()
    k = a.clusters
    _norm = cf.strip_trailing_stop if a.stop_policy == "strip" else (lambda s: s)

    try:
        import seqfold  # noqa: F401
    except Exception:
        warnings.warn("seqfold is not installed: TIS_accessibility_5p, MFE_min_window and stable_structure_count have no "
                      "values and are left out; the MFE profile panel is missing (pip install seqfold)")

    # ---- score table --------------------------------------------------------------------------------------------------
    sc_dna = pd.read_excel(a.scores, sheet_name="CodonHarmony")
    ga = pd.read_excel(a.scores, sheet_name="GeneArt")
    genes = list(dict.fromkeys(ga["Name"]))
    if set(genes) != set(sc_dna["Name"]):
        sys.exit("Genes differ between sheets: GeneArt {} vs CodonHarmony {}".format(sorted(genes), sorted(set(sc_dna["Name"]))))
    if ga["Name"].duplicated().any():
        sys.exit("More than one GeneArt row per gene")
    in_table = [m for m, _ in cc.COMPARE_METRICS if m in sc_dna.columns and m in ga.columns]
    not_in = [m for m, _ in cc.COMPARE_METRICS if m not in in_table]
    metrics, empty = finite_metrics(in_table, sc_dna, ga)
    print("Metrics plotted: {} | not in the table: {} | skipped, empty or non-numeric column (seqfold missing?): {}".format(
        len(metrics), not_in or "none", empty or "none"))
    if a.max_rows_per_gene > 0:
        warnings.warn("--max-rows-per-gene: every plot of this run uses a random subset of the rows (TEST ONLY)")
        keep = []
        for g in genes:
            idx = sc_dna.index[sc_dna["Name"] == g].to_numpy()
            keep += list(np.random.RandomState(0).choice(idx, min(len(idx), a.max_rows_per_gene), replace=False))
        sc_dna = sc_dna.loc[sorted(keep)]
    dna = sc_dna.reset_index(drop=True)

    # ---- the DNA (needed for profiles / motifs / codon usage only) -----------------------------------------------------
    profiles = not a.no_profiles
    if profiles:
        grid_path = a.grid
        if grid_path is None:
            for n in ("human_codonharmony_grid.xlsx",):
                if (ROOT / "CodonHarmony" / "outputs" / n).exists():
                    grid_path = str(ROOT / "CodonHarmony" / "outputs" / n)
                    break
        if grid_path is None:
            sys.exit("grid workbook not found; pass --grid (or use --no-profiles)")
        res = pd.read_excel(grid_path, sheet_name="Results")
        res = res[res["Status"] == "OK"][["Name", "Variant", DNA]].copy()
        res[DNA] = res[DNA].fillna("").astype(str).map(cf.clean_seq).map(_norm)
        n0 = len(dna)
        dna = dna.merge(res.rename(columns={DNA: "_seq"}), on=["Name", "Variant"], how="left")
        if dna["_seq"].isna().any() or len(dna) != n0:
            sys.exit("{} rows of scores_codonharmony.xlsx have no (or an ambiguous) DNA in {} (different grid file?)".format(
                int(dna["_seq"].isna().sum()) + len(dna) - n0, grid_path))
        tr = sas.load_geneart(a.trials)
        ga = ga.merge(tr[["Name", "Protein", "GeneArt"]].rename(columns={"GeneArt": "_seq"}), on="Name", how="left")
        ga["_seq"] = ga["_seq"].map(lambda s: _norm(s) if isinstance(s, str) else s)
        if ga["_seq"].isna().any():
            sys.exit("GeneArt sequence missing for {}".format(list(ga[ga["_seq"].isna()]["Name"])))
        bad = [n for n, s, p in zip(ga["Name"], ga["_seq"], ga["Protein"]) if sas.translate(cf.strip_trailing_stop(s)) != p]
        if bad:
            sys.exit("GeneArt sequence does not translate to the input protein: {}".format(bad))

    # ---- human host reference ----------------------------------------------------------------------------------------
    host = hdin = None
    if not a.no_host:
        host = rc.load_host_counts()
        hdin = rc.load_host_dinuc()
        if host is None:
            warnings.warn("Human.txt not found/readable in host_files/Human: host lines and host plots skipped")
        else:
            for m, v, lab in (("GC_pct", host["gc"], "human host mean GC (Human.txt)"),
                              ("GC3_pct", host["gc3"], "human host GC3 (Human.txt)"),
                              ("GC1", host["gc1"], "human host GC1 (Human.txt)"),
                              ("GC2", host["gc2"], "human host GC2 (Human.txt)")):
                REF_LINES.setdefault(m, []).append((v, lab, "#444444"))
        REF_LINES.setdefault("ENC", []).append((rc.HOST_ENC, "human host ENc {:.1f} (CoCoPUTs)".format(rc.HOST_ENC), "#444444"))

    out = Path(a.out_dir)
    sub = ["distributions", "clusters", "summary", "extras"] + (["profiles", "geneart_features", "liabilities",
                                                                  "host_reference", "agreement"] if profiles else [])
    for s_ in sub:
        (out / s_).mkdir(parents=True, exist_ok=True)

    # ---- 1. one distribution + one cluster image per metric ------------------------------------------------------------
    st, cl, rowcl = analyse(dna, ga, metrics, genes, k)
    for m in metrics:                      # per-row cluster rank (for the colours)
        dna["_rank_" + m] = 0
        for g in genes:
            idx = dna.index[dna["Name"] == g]
            kk = int(st[(st.metric == m) & (st.gene == g)]["n_clusters_used"].iloc[0])
            dna.loc[idx, "_rank_" + m] = [NAMES[kk].index(c) for c in rowcl.loc[idx, m + "_cluster"]]
    for i, m in enumerate(metrics, 1):
        plot_distribution(m, dna, ga, genes, st, out / "distributions" / "{:02d}_{}.png".format(i, m), a.dpi)
        plot_clusters(m, dna, ga, genes, st, cl, k, out / "clusters" / "{:02d}_{}.png".format(i, m), a.dpi)
        print("  {:2d}/{} {}".format(i, len(metrics), m), flush=True)

    # ---- 2. summary plots from the score table ---------------------------------------------------------------------
    S = out / "summary"
    heatmap_cluster(st, genes, metrics, k, S / "heatmap_geneart_cluster.png", a.dpi)
    heatmap_beats(st, genes, metrics, S / "heatmap_geneart_beats_codonharmony.png", a.dpi)
    scored, wins, wc = pair_results(dna, ga, genes, metrics)
    plot_win_rate(wc, scored, S / "win_rate_per_metric.png", a.dpi)
    plot_win_heatmap(wc, genes, scored, S / "win_heatmap_per_gene.png", a.dpi)
    plot_tally(wins, genes, S / "tally_per_gene.png", a.dpi)
    if all(m in dna.columns for m in ("GC_pct", "CAI", "GC3_pct", "tAI")):
        plot_scatter(dna, ga, genes, S / "scatter_cai_gc_tai_gc3.png", a.dpi)
    if all(m in dna.columns for m in ("GC3_pct", "ENC")):
        plot_enc_gc3(dna, ga, genes, host, S / "enc_vs_gc3.png", a.dpi)

    tables = {}
    ex_dna = ex_ga = st_ex = rt = chk = None
    # ---- 3. plots that need the DNA ------------------------------------------------------------------------------------
    if profiles:
        data = score_all(list(dict.fromkeys(list(dna["_seq"]) + list(ga["_seq"]))), max(1, min(a.workers, 61)),
                         out / "profile_cache.pkl")
        chk = check_against_table(dna, ga, data, in_table)
        print("Check against scores_codonharmony.xlsx: {} metric(s) with differences".format(int((chk["differences_gt_1e-6"] > 0).sum())))
        if (chk["differences_gt_1e-6"] > 0).any():
            warnings.warn("scalar scores computed here differ from scores_codonharmony.xlsx for: {} (different reference files, a "
                          "different grid file, or seqfold present in only one of the runs?)".format(
                              list(chk[chk["differences_gt_1e-6"] > 0]["metric"])))
        rows = {g: list(dna.loc[dna["Name"] == g, "_seq"]) for g in genes}
        ga_res = {g: data[ga.loc[ga["Name"] == g, "_seq"].iloc[0]] for g in genes}

        for g in genes:
            plot_gene_profiles(g, ga_res[g], rows[g], data, out / "profiles" / "{}_profiles_GeneArt_vs_CodonHarmony.png".format(g), a.dpi)
            cf.plot_profiles(expand(ga_res[g], "{} GeneArt".format(g)), str(out / "geneart_features" / "{}_GeneArt_profiles.png".format(g)))
            print("  profiles", g, flush=True)

        L = out / "liabilities"
        plot_motif_map(genes, ga_res, rows, data, L / "motif_map_restriction_splice_homopolymer.png", a.dpi)
        rt = restriction_table(genes, ga_res, rows, data)
        plot_restriction_enzymes(rt, genes, L / "restriction_sites_per_enzyme.png", a.dpi)

        # per-codon agreement with GeneArt + identity per row
        mats = {}
        ex_dna = pd.DataFrame({"Name": dna["Name"], "Variant": dna["Variant"]})
        ex_ga = pd.DataFrame({"Name": ga["Name"], "Variant": "GeneArt"})
        ex_dna["pct_identity_nt"] = np.nan
        ex_dna["codons_differing"] = np.nan
        for g in genes:
            idx = dna.index[dna["Name"] == g]
            M, ident = match_matrix(ga.loc[ga["Name"] == g, "_seq"].iloc[0], rows[g])
            mats[g] = M
            ex_dna.loc[idx, "pct_identity_nt"] = ident
            ex_dna.loc[idx, "codons_differing"] = M.shape[1] - M.sum(axis=1)
        plot_codon_agreement(genes, mats, out / "agreement" / "codon_agreement_with_geneart.png", a.dpi)

        for df, src in ((ex_dna, dna), (ex_ga, ga)):
            e = [data[s]["extra"] for s in src["_seq"]]
            df["GC1"] = [x["gc"][0] for x in e]
            df["GC2"] = [x["gc"][1] for x in e]
            df["ENC_deviation"] = pd.to_numeric(pd.Series([data[s]["scores"].get("ENC_deviation") for s in src["_seq"]],
                                                          index=df.index), errors="coerce")
            if host:
                gc3 = np.array([x["gc"][2] for x in e])
                cod = np.vstack([x["cod"] for x in e])
                P = 1000.0 * cod / cod.sum(axis=1, keepdims=True)
                hv = np.array([host["per1000"][c] for c in rc.SENSE])
                df["codon_usage_L1_to_host"] = np.abs(P - hv).sum(axis=1) / 10.0
                df["GC_pos_L1_to_host"] = (np.abs(df["GC1"] - host["gc1"]) + np.abs(df["GC2"] - host["gc2"])
                                           + np.abs(gc3 - host["gc3"]))
            if not a.no_host:
                df["ENC_abs_diff_to_host"] = (pd.to_numeric(src["ENC"], errors="coerce") - rc.HOST_ENC).abs().to_numpy()
            if hdin:
                df["dinuc_total_L1_to_host"] = np.abs(np.vstack([x["dt"] for x in e]) - np.array([hdin[0][d] for d in rc.DINUC])).sum(axis=1)
                df["dinuc_junction_L1_to_host"] = np.abs(np.vstack([x["dj"] for x in e]) - np.array([hdin[1][d] for d in rc.DINUC])).sum(axis=1)

        H = out / "host_reference"
        if host:
            gene_mean = lambda g: np.mean([data[s]["extra"]["cod"] for s in rows[g]], axis=0)        # noqa: E731
            plot_host_codon(host, np.sum([ga_res[g]["extra"]["cod"] for g in genes], axis=0),
                            np.sum([gene_mean(g) for g in genes], axis=0), H / "codon_usage_per_1000_vs_host.png", a.dpi)
            plot_host_gc(host, dna, ga, ex_dna, ex_ga, H / "gc_by_codon_position_vs_host.png", a.dpi)
        if hdin:
            plot_host_dinuc(hdin, [np.vstack([ga_res[g]["extra"][kk_] for g in genes]) for kk_ in ("dt", "dj")],
                            [np.vstack([data[s]["extra"][kk_] for s in dna["_seq"]]) for kk_ in ("dt", "dj")],
                            H / "dinucleotides_vs_host.png", a.dpi)

        # distribution images of the extra columns (same plot as for the 22 metrics, host values drawn as lines)
        cand = [m for m in EXTRA_DIR if m in ex_dna.columns and m in ex_ga.columns]
        extras, _ = finite_metrics(cand, ex_dna, ex_ga)
        DIRECTION.update({m: EXTRA_DIR[m] for m in extras})
        st_ex, _, _ = analyse(ex_dna, ex_ga, extras, genes, 3)
        for i, m in enumerate(extras, 1):
            plot_distribution(m, ex_dna, ex_ga, genes, st_ex, out / "extras" / "X{:02d}_{}.png".format(i, m), a.dpi)
        print("  extras:", extras)

        # remaining metrics (Kozak, baseline-audit metrics, other scorer scalars): distribution + cluster image each
        cache = {}
        def mm(s):
            if s not in cache:
                cache[s] = more_metrics(s, data[s]["scores"])
            return cache[s]
        mdna = pd.concat([dna[["Name", "Variant"]], pd.DataFrame([mm(s) for s in dna["_seq"]], index=dna.index)], axis=1)
        mga = pd.concat([ga[["Name"]].assign(Variant="GeneArt"), pd.DataFrame([mm(s) for s in ga["_seq"]], index=ga.index)], axis=1)
        more, more_bad = finite_metrics(list(MORE_DIR), mdna, mga)
        print("  more metrics: {} | skipped (empty / not defined for these sequences): {}".format(more, more_bad or "none"))
        DIRECTION.update({m: MORE_DIR[m] for m in more})
        st_m, cl_m, rowcl_m = analyse(mdna, mga, more, genes, k)
        for m in more:
            mdna["_rank_" + m] = 0
            for g in genes:
                idx = mdna.index[mdna["Name"] == g]
                kk = int(st_m[(st_m.metric == m) & (st_m.gene == g)]["n_clusters_used"].iloc[0])
                mdna.loc[idx, "_rank_" + m] = [NAMES[kk].index(c) for c in rowcl_m.loc[idx, m + "_cluster"]]
        for d_ in ("distributions", "clusters"):
            (out / "more_metrics" / d_).mkdir(parents=True, exist_ok=True)
        for i, m in enumerate(more, 1):
            plot_distribution(m, mdna, mga, genes, st_m, out / "more_metrics" / "distributions" / "M{:02d}_{}.png".format(i, m), a.dpi)
            plot_clusters(m, mdna, mga, genes, st_m, cl_m, k, out / "more_metrics" / "clusters" / "M{:02d}_{}.png".format(i, m), a.dpi)

        # optional per-pair figures
        if a.pair_plots:
            cc._init_mpl()
            (out / "pairs").mkdir(exist_ok=True)
            (out / "pair_features").mkdir(exist_ok=True)
            sel = []
            for g in genes:
                d = dna[dna["Name"] == g]
                if a.pair_plots > 0:
                    d = d.sort_values("CAI", ascending=False, kind="stable").drop_duplicates("_seq").head(a.pair_plots)
                sel += [(g, r["Variant"].split("_")[0], r["_seq"], r["Variant"]) for _, r in d.iterrows()]
            for n, (g, vid, seq, variant) in enumerate(sel, 1):
                ra, rb = expand(ga_res[g], LA), expand(data[seq], LB)
                rws = pair_rows(ra["scores"], rb["scores"])
                tag = "{}_{}".format(g, vid)
                cc.overlay_profiles(ra, rb, LA, LB, str(out / "pairs" / (tag + "_overlay_profiles.png")))
                cc.scorecard(rws, LA, LB, str(out / "pairs" / (tag + "_scorecard.png")))
                cf.plot_profiles(dict(rb, name="{} {}".format(g, variant)), str(out / "pair_features" / (tag + "_profiles.png")))
                if n % 20 == 0:
                    print("  pair figures {}/{}".format(n, len(sel)), flush=True)

    # ---- 4. tables behind the plots ------------------------------------------------------------------------------------
    pub = st.drop(columns=["_code", "_bounds"])
    piv = lambda col: pub.pivot(index="gene", columns="metric", values=col).loc[genes, metrics]       # noqa: E731
    desc = [("GeneArt_vs_CodonHarmony", "one row per metric x gene: GeneArt value, Codon Harmony min/median/max, percentile, % of rows GeneArt beats, cluster of GeneArt"),
            ("Cluster_ranges", "one row per metric x gene x cluster: rows, % of gene rows, min, max, centre"),
            ("Win_counts", "per metric x gene: pairs where the Codon Harmony row wins / GeneArt wins / draw"),
            ("Pivot_percentile", "GeneArt percentile among the Codon Harmony rows (0 below all, 100 above all)"),
            ("Pivot_beats_pct", "% of Codon Harmony rows GeneArt is better than (blank for metrics without a rule)"),
            ("Pivot_cluster", "cluster name GeneArt falls into")]
    if profiles:
        desc += [("Extras_GeneArt_vs_CodonHarmony", "the same as GeneArt_vs_CodonHarmony for the extra columns (GC1, GC2, ENC_deviation, distances to the host)"),
                 ("Extras_CodonHarmony_rows", "per Codon Harmony row: GC1, GC2, ENC deviation, distances to the host, % identity and codons differing vs GeneArt"),
                 ("Extras_GeneArt", "the same extra columns for the GeneArt sequences"),
                 ("Restriction_per_enzyme", "per enzyme x gene: GeneArt sites, Codon Harmony mean sites per row and % of rows with a site"),
                 ("More_GeneArt_vs_CodonHarmony", "GeneArt_vs_CodonHarmony for the remaining metrics (Kozak, GC min/std 60 bp, CpG count, longest repeat, AT%, other scorer scalars)"),
                 ("More_CodonHarmony_rows", "per Codon Harmony row: the remaining metrics"), ("More_GeneArt", "the remaining metrics for the GeneArt sequences"),
                 ("Check_vs_scores_all", "scalar scores recomputed for the profile plots vs scores_codonharmony.xlsx (differences must be 0)")]
    with pd.ExcelWriter(out / "summary_tables.xlsx", engine="openpyxl") as w:
        pd.DataFrame(desc, columns=["sheet", "content"]).to_excel(w, sheet_name="Read_me", index=False)
        pub.to_excel(w, sheet_name="GeneArt_vs_CodonHarmony", index=False)
        cl.to_excel(w, sheet_name="Cluster_ranges", index=False)
        wc.to_excel(w, sheet_name="Win_counts", index=False)
        piv("GeneArt_percentile").to_excel(w, sheet_name="Pivot_percentile")
        piv("GeneArt_beats_pct").to_excel(w, sheet_name="Pivot_beats_pct")
        piv("GeneArt_cluster").to_excel(w, sheet_name="Pivot_cluster")
        if profiles:
            st_ex.drop(columns=["_code", "_bounds"]).to_excel(w, sheet_name="Extras_GeneArt_vs_CodonHarmony", index=False)
            ex_dna.to_excel(w, sheet_name="Extras_CodonHarmony_rows", index=False)
            ex_ga.to_excel(w, sheet_name="Extras_GeneArt", index=False)
            rt.to_excel(w, sheet_name="Restriction_per_enzyme", index=False)
            st_m.drop(columns=["_code", "_bounds"]).to_excel(w, sheet_name="More_GeneArt_vs_CodonHarmony", index=False)
            mdna.drop(columns=[c for c in mdna.columns if c.startswith("_rank_")]).to_excel(w, sheet_name="More_CodonHarmony_rows", index=False)
            mga.to_excel(w, sheet_name="More_GeneArt", index=False)
            chk.to_excel(w, sheet_name="Check_vs_scores_all", index=False)
        for ws in w.book.worksheets:
            ws.freeze_panes = "B2"
            for col in ws.columns:
                ws.column_dimensions[col[0].column_letter].width = min(60, max(10, len(str(col[0].value)) + 2))
    rowcl.to_csv(out / "row_clusters.csv", index=False)
    n_png = len(list(out.rglob("*.png")))
    print("Done. {} PNG files, summary_tables.xlsx and row_clusters.csv in {}".format(n_png, out))


if __name__ == "__main__":
    main()
