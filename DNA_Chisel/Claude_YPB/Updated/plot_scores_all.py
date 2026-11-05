# -*- coding: utf-8 -*-
"""
plot_scores_all.py  -  plots from scores_all.xlsx (no re-scoring)
=================================================================

Reads the score table written by score_all_sequences.py (sheets DNAChisel and GeneArt) and draws:

  1. distributions/<NN>_<metric>.png   ONE image per metric: every gene (rows of the image), every DNA Chisel
                                       row (grey dots) with box = median/IQR, and the GeneArt value as a red diamond.
  2. clusters/<NN>_<metric>.png        ONE image per metric: for each gene the DNA Chisel rows are split by
                                       K-means (sklearn, default 3 clusters = Low / Mid / High range of that metric
                                       in that gene); the GeneArt value is overlaid and the cluster it falls into
                                       is printed.
  3. summary/heatmap_geneart_cluster.png        gene x metric: which cluster GeneArt falls into (+ percentile)
     summary/heatmap_geneart_beats_dnachisel.png gene x scored metric: % of DNA Chisel rows GeneArt is better than
  4. score_summary.xlsx / row_clusters.csv   all numbers behind the plots (see README_plot_scores.md)

Nothing is re-computed from DNA: every value is read from scores_all.xlsx. The better/worse direction of each
metric and the GC band come from codon_compare.py (imported, not modified), so there is one definition.

Run from Sequence_scoring/ (pip install numpy pandas matplotlib openpyxl scikit-learn):
    python plot_scores_all.py
    python plot_scores_all.py --scores outputs/scores/scores_all.xlsx --out-dir outputs/score_plots --clusters 3
"""

from __future__ import annotations

import argparse
import sys
import textwrap
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap, to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from sklearn.cluster import KMeans

import codon_compare as cc            # given, unchanged: metric list, better-direction, GC band

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent if (HERE.parent / "host_files").exists() else HERE

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


def mlabel(m):
    return "{} {}".format(m, ARROW[DIRECTION[m]])


def fmt(v):
    return "{:.4g}".format(v)


def palette(k):
    cm = LinearSegmentedColormap.from_list("lmh", ["#0072B2", "#E69F00", "#D55E00"])
    return [cm(x) for x in np.linspace(0, 1, k)] if k > 1 else [cm(0.0)]


# --------------------------------------------------------------------------------------------------
# numbers
# --------------------------------------------------------------------------------------------------
def pct_below(x, v):
    """Percentile of v among x: % of x below v, ties count half (0 = below all, 100 = above all)."""
    return 100.0 * ((x < v - EPS).sum() + 0.5 * (np.abs(x - v) <= EPS).sum()) / len(x)


def beat_share(x, v, direction):
    """% of DNA Chisel rows that GeneArt is better than (ties count half). NaN for 'na' metrics."""
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
            if kk == 1:      # no spread in DNA Chisel: say where GeneArt is relative to that single value
                gname = "Same value" if abs(gv - x[0]) <= EPS else ("Above range" if gv > x[0] else "Below range")
            for i in range(kk):
                xs = x[lab == i]
                clus.append({"metric": m, "gene": g, "cluster": names[i], "rank": i, "n_clusters_used": kk,
                             "n_rows": len(xs), "pct_of_gene_rows": round(100 * len(xs) / len(x), 2),
                             "min": xs.min(), "max": xs.max(), "centre": cen[i],
                             "GeneArt_here": i == gr})
            stats.append({"metric": m, "gene": g, "better": DIRECTION[m], "GeneArt": gv,
                          "DNAChisel_n": len(x), "DNAChisel_min": x.min(), "DNAChisel_median": float(np.median(x)),
                          "DNAChisel_max": x.max(), "GeneArt_percentile": round(pct_below(x, gv), 2),
                          "GeneArt_beats_pct": round(beat_share(x, gv, DIRECTION[m]), 2),
                          "GeneArt_cluster": gname, "n_clusters_used": kk,
                          "GeneArt_outside_DNAChisel_range": bool(gv < x.min() - EPS or gv > x.max() + EPS),
                          "_code": heat_code(gr, kk, k), "_bounds": bounds})
    return pd.DataFrame(stats), pd.DataFrame(clus), rowcl


# --------------------------------------------------------------------------------------------------
# per-metric figures
# --------------------------------------------------------------------------------------------------
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
                          "Each grey dot = one DNA Chisel row of that gene (all Status-OK grid rows). Box = quartiles, thick bar = "
                          "median. Red diamond = GeneArt sequence of that gene. Same x-scale for every gene. Right column = GeneArt "
                          "percentile among that gene's DNA Chisel rows (0 = below all, 50 = median, 100 = above all).",
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
        star = "*" if sub.loc[g, "GeneArt_outside_DNAChisel_range"] else ""
        axt.text(0.05, i, "{:.0f}{}".format(sub.loc[g, "GeneArt_percentile"], star), va="center", fontsize=11,
                 fontweight="bold")
    ax.set_yticklabels(["{}\nn={:,}".format(g, int((dna["Name"] == g).sum())) for g in genes], fontsize=9)
    handles = [Line2D([], [], marker="o", ls="", color=C_DNA, alpha=0.6, markersize=5, label="DNA Chisel row"),
               Line2D([], [], color="#1f2d3a", lw=2.6, label="median (box = quartiles)"),
               Line2D([], [], marker="D", ls="", color=C_GA, markeredgecolor="white", markersize=9, label="GeneArt")]
    handles += ref_handles(metric)
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.45, 0.035), ncol=4, fontsize=8.5, frameon=False)
    fig.text(0.01, 0.008, "* GeneArt value lies outside the whole DNA Chisel range of that gene", fontsize=8, color="#444")
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
                          "K-means (k={}) on this metric, run separately inside each gene on its DNA Chisel rows. Dots are coloured by "
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
        star = " *" if r["GeneArt_outside_DNAChisel_range"] else ""
        axt.text(0.84, i, r["GeneArt_cluster"] + star, va="center", fontsize=11, fontweight="bold",
                 color="#555" if kk == 1 else gcol)
    ax.set_yticklabels(genes, fontsize=9.5)
    handles = [Line2D([], [], marker="o", ls="", color=cols[j], markersize=6, label=NAMES[k][j] + " cluster")
               for j in range(k)]
    handles += [Line2D([], [], color="#222", ls="--", label="cluster boundary"),
                Line2D([], [], marker="*", ls="", color="white", markeredgecolor="black", markersize=14, label="GeneArt")]
    handles += ref_handles(metric)
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.035), ncol=k + 2, fontsize=8.5, frameon=False)
    fig.text(0.01, 0.008, "* GeneArt value lies outside the whole DNA Chisel range of that gene.  Low/Mid/High describe the numeric "
             "range of the metric only, not good or bad: see the better/worse rule in the title.", fontsize=8, color="#444")
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


# --------------------------------------------------------------------------------------------------
# summary heatmaps
# --------------------------------------------------------------------------------------------------
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
    out = np.array([[bool(st[(st.metric == m) & (st.gene == g)]["GeneArt_outside_DNAChisel_range"].iloc[0])
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
                 "Colour = cluster (K-means inside each gene, on the DNA Chisel rows of that gene) of the metric's value range; "
                 "number = GeneArt percentile among those rows",
                 fontsize=11, fontweight="bold", loc="left", pad=48)
    handles = [Patch(color=colors[j], label="{} range".format(NAMES[k][j])) for j in range(k)]
    handles += [Patch(color=C_NOSPREAD, label="DNA Chisel has a single value (no spread)")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=k + 1, frameon=False, fontsize=9,
              title="Cell colour", title_fontsize=9)
    fig.text(0.01, 0.01, "Number in cell: percentile of GeneArt among the DNA Chisel rows of the same gene "
             "(0 = below all, 50 = median, 100 = above all; ties count half).  * = GeneArt is outside the DNA Chisel range of that gene.\n"
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
    ax.set_title("How often is GeneArt better than DNA Chisel?\n"
                 "Cell = % of the DNA Chisel rows of the same gene that the GeneArt sequence beats on that metric",
                 fontsize=11, fontweight="bold", loc="left", pad=48)
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.05, pad=0.30, aspect=40)
    cb.set_ticks([0, 25, 50, 75, 100])
    cb.set_ticklabels(["0\nGeneArt worse than\nall DNA Chisel rows", "25", "50\ntie / middle", "75",
                       "100\nGeneArt better than\nall DNA Chisel rows"], fontsize=8.5)
    note = ("Better/worse follows codon_compare: \u2191 higher wins, \u2193 lower wins, \u25CE GC% closer to {:g}% wins. Ties count "
            "half, so a metric where both sides have the same value in every row shows 50. Metrics with no better/worse rule "
            "(ENC, tAI_ramp, minmax_mean, codon_autocorr_lag1, GC3_pct) are not shown here (see the cluster heatmap)."
            ).format(BAND_MID)
    fig.text(0.01, 0.01, textwrap.fill(note, 200), fontsize=8, color="#333", va="bottom")
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# pair-wise win/draw/loss (each DNA Chisel row of a gene is one pair with the GeneArt sequence of that gene)
# --------------------------------------------------------------------------------------------------
C_DC = "#2B6CB0"          # DNA Chisel in the win/loss and scatter plots (GeneArt = C_GA)
GENE_COL = ["#1565C0", "#E65100", "#2E7D32", "#6A1B9A", "#C62828", "#00838F", "#AD1457", "#5D4037"]


def row_winners(x, v, direction):
    """+1 = DNA Chisel row wins, -1 = GeneArt wins, 0 = draw; same rule as codon_compare.winner_of."""
    if direction == "higher":
        d = x - v
    elif direction == "lower":
        d = v - x
    else:                                                    # band
        d = abs(v - BAND_MID) - np.abs(x - BAND_MID)
    return np.where(d > 1e-9, 1, np.where(d < -1e-9, -1, 0))


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
            rows.append({"metric": m, "gene": g, "pairs": len(idx), "DNAChisel_wins": int((cols[m] == 1).sum()),
                         "GeneArt_wins": int((cols[m] == -1).sum()), "draws": int((cols[m] == 0).sum())})
        wins[g] = pd.DataFrame(cols, index=idx)
    return scored, wins, pd.DataFrame(rows)


def plot_win_rate(wc, scored, path, dpi):
    tot = wc.groupby("metric")[["DNAChisel_wins", "draws", "GeneArt_wins"]].sum().loc[scored]
    pct = 100 * tot.div(tot.sum(axis=1), axis=0).sort_values("DNAChisel_wins")
    fig, ax = plt.subplots(figsize=(11, 0.42 * len(pct) + 2.4))
    y = np.arange(len(pct))
    ax.barh(y, pct["DNAChisel_wins"], color=C_DC, label="DNA Chisel row better")
    ax.barh(y, pct["draws"], left=pct["DNAChisel_wins"], color="#BDBDBD", label="draw")
    ax.barh(y, pct["GeneArt_wins"], left=pct["DNAChisel_wins"] + pct["draws"], color=C_GA, label="GeneArt better")
    ax.set_yticks(y)
    ax.set_yticklabels([mlabel(m) for m in pct.index], fontsize=9)
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of all pairs  (pair = one DNA Chisel row vs the GeneArt sequence of the same gene; n = {:,})"
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
    ax1.bar(x, mean["B"], color=C_DC, label="DNA Chisel better")
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
    ax2.set_xlabel("number of scored metrics (out of {}) on which the DNA Chisel row beats GeneArt".format(wins[genes[0]].shape[1]))
    ax2.grid(axis="x", alpha=0.3)
    ax2.set_title("Spread over DNA Chisel rows (one dot = one pair, bar = median)", fontsize=11, fontweight="bold", loc="left")
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
    h += [Line2D([], [], marker="o", ls="", color="#888", alpha=0.5, markersize=4, label="dot = DNA Chisel row"),
          Line2D([], [], marker="*", ls="", color="white", markeredgecolor="black", markersize=13,
                 label="star = GeneArt (gene colour)")]
    fig.legend(handles=h, loc="lower center", ncol=5, frameon=False, fontsize=9)
    fig.suptitle("Codon optimality against GC: every DNA Chisel row and the GeneArt sequence of each gene",
                 fontsize=12.5, fontweight="bold", x=0.01, ha="left")
    fig.tight_layout(rect=[0, 0.1, 1, 0.95])
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Plots from scores_all.xlsx: distributions, clusters, heatmaps.")
    ap.add_argument("--scores", default=str(HERE / "outputs" / "scores" / "scores_all.xlsx"))
    ap.add_argument("--out-dir", default=str(HERE / "outputs" / "score_plots"))
    ap.add_argument("--clusters", type=int, default=3, choices=[2, 3, 4, 5], help="K-means clusters per metric per gene")
    ap.add_argument("--dpi", type=int, default=130)
    a = ap.parse_args()

    dna = pd.read_excel(a.scores, sheet_name="DNAChisel")
    ga = pd.read_excel(a.scores, sheet_name="GeneArt")
    metrics = [m for m, _ in cc.COMPARE_METRICS if m in dna.columns and m in ga.columns]
    miss = [m for m, _ in cc.COMPARE_METRICS if m not in metrics]
    print("Metrics plotted: {} (not in the table: {})".format(len(metrics), miss or "none"))
    genes = list(dict.fromkeys(ga["Name"]))
    if set(genes) != set(dna["Name"]):
        sys.exit("Genes differ between sheets: GeneArt {} vs DNAChisel {}".format(sorted(genes), sorted(set(dna["Name"]))))
    if ga["Name"].duplicated().any():
        sys.exit("More than one GeneArt row per gene")
    dna = dna.reset_index(drop=True)
    k = a.clusters

    st, cl, rowcl = analyse(dna, ga, metrics, genes, k)
    # per-row cluster rank (for colouring) kept on a private copy
    for m in metrics:
        dna["_rank_" + m] = 0
        for g in genes:
            idx = dna.index[dna["Name"] == g]
            kk = int(st[(st.metric == m) & (st.gene == g)]["n_clusters_used"].iloc[0])
            dna.loc[idx, "_rank_" + m] = [NAMES[kk].index(c) for c in rowcl.loc[idx, m + "_cluster"]]

    out = Path(a.out_dir)
    for sub in ("distributions", "clusters", "summary"):
        (out / sub).mkdir(parents=True, exist_ok=True)
    for i, m in enumerate(metrics, 1):
        plot_distribution(m, dna, ga, genes, st, out / "distributions" / "{:02d}_{}.png".format(i, m), a.dpi)
        plot_clusters(m, dna, ga, genes, st, cl, k, out / "clusters" / "{:02d}_{}.png".format(i, m), a.dpi)
        print("  {:2d}/{} {}".format(i, len(metrics), m), flush=True)
    heatmap_cluster(st, genes, metrics, k, out / "summary" / "heatmap_geneart_cluster.png", a.dpi)
    heatmap_beats(st, genes, metrics, out / "summary" / "heatmap_geneart_beats_dnachisel.png", a.dpi)

    scored, wins, wc = pair_results(dna, ga, genes, metrics)
    plot_win_rate(wc, scored, out / "summary" / "win_rate_per_metric.png", a.dpi)
    plot_tally(wins, genes, out / "summary" / "tally_per_gene.png", a.dpi)
    plot_scatter(dna, ga, genes, out / "summary" / "scatter_cai_gc_tai_gc3.png", a.dpi)

    # tables behind the plots
    pub = st.drop(columns=["_code", "_bounds"])
    piv = lambda col: pub.pivot(index="gene", columns="metric", values=col).loc[genes, metrics]
    readme = pd.DataFrame({"sheet": ["GeneArt_vs_DNAChisel", "Cluster_ranges", "Pivot_percentile", "Pivot_beats_pct",
                                     "Pivot_cluster", "Win_counts"],
                           "content": ["one row per metric x gene: GeneArt value, DNA Chisel min/median/max, percentile, "
                                       "% of rows GeneArt beats, cluster of GeneArt",
                                       "one row per metric x gene x cluster: rows, % of gene rows, min, max, centre",
                                       "GeneArt percentile among DNA Chisel rows (0 below all, 100 above all)",
                                       "% of DNA Chisel rows GeneArt is better than (blank for metrics without a rule)",
                                       "cluster name GeneArt falls into",
                                       "per metric x gene: pairs where the DNA Chisel row wins / GeneArt wins / draw"]})
    with pd.ExcelWriter(out / "score_summary.xlsx", engine="openpyxl") as w:
        readme.to_excel(w, sheet_name="Read_me", index=False)
        pub.to_excel(w, sheet_name="GeneArt_vs_DNAChisel", index=False)
        cl.to_excel(w, sheet_name="Cluster_ranges", index=False)
        wc.to_excel(w, sheet_name="Win_counts", index=False)
        piv("GeneArt_percentile").to_excel(w, sheet_name="Pivot_percentile")
        piv("GeneArt_beats_pct").to_excel(w, sheet_name="Pivot_beats_pct")
        piv("GeneArt_cluster").to_excel(w, sheet_name="Pivot_cluster")
        for ws in w.book.worksheets:
            ws.freeze_panes = "B2"
            for col in ws.columns:
                ws.column_dimensions[col[0].column_letter].width = min(60, max(10, len(str(col[0].value)) + 2))
    rowcl.to_csv(out / "row_clusters.csv", index=False)
    print("Done. {} distribution + {} cluster plots + 2 heatmaps + 3 win/scatter plots in {}".format(len(metrics), len(metrics), out))


if __name__ == "__main__":
    main()
