# -*- coding: utf-8 -*-
"""
plot_profiles_host.py  -  positional profiles, motif map and host-reference plots: GeneArt vs ALL DNA Chisel rows
=================================================================================================================

Companion of plot_scores_all.py. That script needs only scores_all.xlsx. This one covers the plots that need the
DNA itself (positional profiles, motif positions, codon / dinucleotide usage), which scores_all.xlsx does not
contain. It therefore runs codon_features.compute_features once per distinct sequence (the same scorer, the same
human_codon.csv / human_tgcn.csv as score_all_sequences.py), keeps the profiles, and caches them
(profile_cache.pkl) so a second run does not score again. The scalar scores it gets back are compared with
scores_all.xlsx as a safeguard (sheet Check_vs_scores_all).

Imported, NOT modified: codon_features, codon_compare, score_all_sequences (GeneArt loader), run_comparison (host
reference loaders and per-sequence helpers), plot_scores_all (distribution plot, cluster analysis, colours).

Outputs (outputs/profile_plots/)
  profiles/<gene>_profiles_GeneArt_vs_DNAChisel.png   12 positional panels per gene: GeneArt line over the median and
                                                      5-95% / 25-75% bands of all DNA Chisel rows of that gene
  geneart_features/<gene>_GeneArt_profiles.png        codon_features 9-panel figure of each GeneArt sequence
  motif_map_restriction_splice_homopolymer.png        where the liabilities sit along each CDS (GeneArt vs DNA Chisel)
  host_reference/codon_usage_per_1000_vs_host.png     codon usage per 1000 codons: host vs GeneArt vs DNA Chisel
  host_reference/dinucleotides_vs_host.png            total and junction dinucleotides per 1000: host vs both
  extras/NN_<metric>.png                              distribution plots (as in plot_scores_all) for GC1, GC2 and the
                                                      distances to the host reference; host values drawn as lines
  pairs/ and pair_features/                           only with --pair-plots (per-pair overlay, scorecard, 9-panel)
  profile_host_summary.xlsx                           numbers behind the plots + Check_vs_scores_all

Run from Sequence_scoring/ (pip install numpy pandas matplotlib openpyxl scikit-learn seqfold tqdm):
    python plot_profiles_host.py
    python plot_profiles_host.py --workers 40 --pair-plots 10       # also overlay/scorecard/9-panel for the 10 best-CAI
                                                                    # distinct sequences of each gene
"""

from __future__ import annotations

import argparse
import os
import pickle
import sys
import warnings
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import codon_compare as cc            # given, unchanged
import codon_features as cf           # given, unchanged
import plot_scores_all as ps          # plots/analysis of the score table
import run_comparison as rc           # given: host loaders, codon/dinucleotide helpers
import score_all_sequences as sas     # given: GeneArt loader, reference file paths, translate

HERE = Path(__file__).resolve().parent
ROOT = sas.ROOT
DNA = sas.DNA
LA, LB = "GeneArt", "DNAChisel"
C_GA, C_DC = ps.C_GA, ps.C_DC
C_HOST = "#9E9E9E"

# (key, title, x label, path inside compute_features()["profiles"])
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


# --------------------------------------------------------------------------------------------------
# scoring of the sequences (once per distinct DNA, in worker processes) + cache
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
    """Scalar scores returned by the scorer here must equal the scores_all.xlsx values."""
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


# --------------------------------------------------------------------------------------------------
# plots
# --------------------------------------------------------------------------------------------------
def _aligned(arrs, ref_idx):
    out = np.full((len(arrs), len(ref_idx)), np.nan, np.float32)
    for i, (idx, val) in enumerate(arrs):
        if len(idx) == len(ref_idx) and np.array_equal(idx, ref_idx):
            out[i] = val
        else:
            _, a, b = np.intersect1d(ref_idx, idx, return_indices=True)
            out[i, a] = val[b]
    return out


def plot_gene_profiles(gene, ga_res, seqs, data, path, dpi):
    panels = [(k, t, x) for k, t, x, _ in PANELS] + [("acf", "Autocorrelation of per-codon tAI (ACF)", "lag"),
                                                    ("tm", "Fragment Tm (60-nt fragments, \u00b0C)", "fragment")]
    keep = [p for p in panels if ga_res["prof"][p[0]] is not None]
    ncol = 2
    nrow = int(np.ceil(len(keep) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(16, 2.7 * nrow + 1.2))
    axs = np.atleast_1d(axs).ravel()
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
            ax.bar(gi - 0.2, gv, 0.4, color=C_GA, label=LA)
            ax.bar(gi + 0.2, p50, 0.4, color=C_DC, yerr=[p50 - p25, p75 - p50], ecolor="#222", capsize=2)
        else:
            ax.fill_between(gi, p5, p95, color=C_DC, alpha=0.15, linewidth=0)
            ax.fill_between(gi, p25, p75, color=C_DC, alpha=0.30, linewidth=0)
            ax.plot(gi, p50, color=C_DC, lw=1.2)
            ax.plot(gi, gv, color=C_GA, lw=1.7)
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
         Line2D([], [], color=C_DC, lw=1.5, label="DNA Chisel median"),
         Patch(color=C_DC, alpha=0.30, label="DNA Chisel 25\u201375%"),
         Patch(color=C_DC, alpha=0.15, label="DNA Chisel 5\u201395%")]
    fig.legend(handles=h, loc="upper right", ncol=4, frameon=False, fontsize=9.5, bbox_to_anchor=(0.99, 0.995))
    fig.suptitle("{}: positional profiles, GeneArt vs all {:,} DNA Chisel rows (bars in ACF / Tm panels: median and "
                 "quartiles)".format(gene, len(seqs)), fontsize=12.5, fontweight="bold", x=0.01, ha="left", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


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
        ax.set_title("{}  ({:,} nt):  GeneArt has {} restriction / {} splice / {} homopolymer hits;  DNA Chisel rows with at "
                     "least one: {:.0f}% / {:.0f}% / {:.0f}%  (n = {:,})".format(
                         g, L, counts[0], counts[1], counts[2], *[100.0 * a / len(seqs) for a in anyhit], len(seqs)),
                     fontsize=9.5, loc="left", fontweight="bold", pad=14)
        ax.set_xlabel("nt position" if g == genes[-1] else "")
    fig.tight_layout(rect=[0, 0.06, 0.93, 0.96])
    cax = fig.add_axes([0.945, 0.2, 0.012, 0.6])
    cb = fig.colorbar(ScalarMappable(Normalize(0, 100), "Blues"), cax=cax)
    cb.set_label("% of the gene's DNA Chisel rows with a hit in that {}-nt bin".format(bin_nt), fontsize=8.5)
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
           label="DNA Chisel (mean over the rows of each gene, genes pooled)")
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
        ax.bar(x + 0.28, dc.mean(axis=0), 0.27, color=C_DC, label="DNA Chisel (mean over all rows; whiskers = quartiles)")
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
def main():
    ap = argparse.ArgumentParser(description="Positional profiles, motif map, host-reference plots.")
    ap.add_argument("--scores", default=str(HERE / "outputs" / "scores" / "scores_all.xlsx"))
    ap.add_argument("--grid", default=None, help="grid output with the DNA (sheet Results); default: human_grid_1.xlsx, "
                    "else human_grid.xlsx in ../DNAChisel/outputs")
    ap.add_argument("--trials", default=None, help="the ORIGINAL trials workbook (GeneArt sequences)")
    ap.add_argument("--out-dir", default=str(HERE / "outputs" / "profile_plots"))
    ap.add_argument("--workers", type=int, default=min(40, max(1, (os.cpu_count() or 2) - 1)))
    ap.add_argument("--pair-plots", type=int, default=0, help="0 = none; N>0 = per gene the N best-CAI distinct DNA Chisel "
                    "sequences get overlay + scorecard + 9-panel figure; -1 = every row (tens of thousands of files)")
    ap.add_argument("--max-rows-per-gene", type=int, default=0, help="TEST ONLY: random subset of rows per gene")
    ap.add_argument("--no-host", action="store_true", help="skip the host-reference plots")
    ap.add_argument("--dpi", type=int, default=130)
    a = ap.parse_args()

    try:
        import seqfold  # noqa: F401
    except Exception:
        warnings.warn("seqfold is not installed: the MFE profile panel will be missing (pip install seqfold)")

    sc_dna = pd.read_excel(a.scores, sheet_name="DNAChisel")
    sc_ga = pd.read_excel(a.scores, sheet_name="GeneArt")
    metrics = [m for m, _ in cc.COMPARE_METRICS if m in sc_dna.columns and m in sc_ga.columns]
    genes = list(dict.fromkeys(sc_ga["Name"]))

    grid_path = a.grid
    if grid_path is None:
        for n in ("human_grid_1.xlsx", "human_grid.xlsx"):
            if (ROOT / "DNAChisel" / "outputs" / n).exists():
                grid_path = str(ROOT / "DNAChisel" / "outputs" / n)
                break
    if grid_path is None:
        sys.exit("grid workbook not found; pass --grid")
    res = pd.read_excel(grid_path, sheet_name="Results")
    res = res[res["Status"] == "OK"][["Name", "Variant", DNA]].copy()
    res[DNA] = res[DNA].fillna("").astype(str).map(cf.clean_seq)
    dna = sc_dna.merge(res.rename(columns={DNA: "_seq"}), on=["Name", "Variant"], how="left")
    if dna["_seq"].isna().any() or len(dna) != len(sc_dna):
        sys.exit("{} rows of scores_all.xlsx have no DNA in {} (different grid file?)".format(
            int(dna["_seq"].isna().sum()), grid_path))
    tr = sas.load_geneart(a.trials)
    ga = sc_ga.merge(tr[["Name", "Protein", "GeneArt"]].rename(columns={"GeneArt": "_seq"}), on="Name", how="left")
    if ga["_seq"].isna().any():
        sys.exit("GeneArt sequence missing for {}".format(list(ga[ga["_seq"].isna()]["Name"])))
    bad = [n for n, s, p in zip(ga["Name"], ga["_seq"], ga["Protein"]) if sas.translate(cf.strip_trailing_stop(s)) != p]
    if bad:
        sys.exit("GeneArt sequence does not translate to the input protein: {}".format(bad))

    if a.max_rows_per_gene > 0:
        warnings.warn("--max-rows-per-gene: every plot of this run uses a random subset of the rows (TEST ONLY)")
        dna = dna.groupby("Name", group_keys=False).apply(
            lambda d: d.sample(min(len(d), a.max_rows_per_gene), random_state=0)).reset_index(drop=True)
    dna = dna.reset_index(drop=True)

    out = Path(a.out_dir)
    for sub in ("profiles", "geneart_features", "host_reference", "extras"):
        (out / sub).mkdir(parents=True, exist_ok=True)
    data = score_all(list(dict.fromkeys(list(dna["_seq"]) + list(ga["_seq"]))), max(1, min(a.workers, 61)),
                     out / "profile_cache.pkl")

    chk = check_against_table(dna, ga, data, metrics)
    print("Check against scores_all.xlsx: {} metric(s) with differences".format(int((chk["differences_gt_1e-6"] > 0).sum())))
    if (chk["differences_gt_1e-6"] > 0).any():
        warnings.warn("scalar scores computed here differ from scores_all.xlsx for: {} (different reference files or "
                      "seqfold missing in one of the runs?)".format(list(chk[chk["differences_gt_1e-6"] > 0]["metric"])))

    rows = {g: list(dna.loc[dna["Name"] == g, "_seq"]) for g in genes}
    ga_res = {g: data[ga.loc[ga["Name"] == g, "_seq"].iloc[0]] for g in genes}

    # ---- positional profiles per gene + GeneArt 9-panel figure ------------------------------------------------------
    for g in genes:
        plot_gene_profiles(g, ga_res[g], rows[g], data, out / "profiles" / "{}_profiles_GeneArt_vs_DNAChisel.png".format(g), a.dpi)
        cf.plot_profiles(expand(ga_res[g], "{} GeneArt".format(g)), str(out / "geneart_features" / "{}_GeneArt_profiles.png".format(g)))
        print("  profiles", g, flush=True)
    plot_motif_map(genes, ga_res, rows, data, out / "motif_map_restriction_splice_homopolymer.png", a.dpi)

    # ---- host reference ----------------------------------------------------------------------------------------------
    host = None if a.no_host else rc.load_host_counts()
    hdin = None if a.no_host else rc.load_host_dinuc()
    if not a.no_host and host is None:
        warnings.warn("Human.txt not found/readable in host_files/Human: host plots and host distances skipped")
    ex_dna = pd.DataFrame({"Name": dna["Name"], "Variant": dna["Variant"]})
    ex_ga = pd.DataFrame({"Name": ga["Name"], "Variant": "GeneArt"})
    for df, src in ((ex_dna, dna), (ex_ga, ga)):
        e = [data[s]["extra"] for s in src["_seq"]]
        df["GC1"] = [x["gc"][0] for x in e]
        df["GC2"] = [x["gc"][1] for x in e]
        gc3 = np.array([x["gc"][2] for x in e])
        cod = np.vstack([x["cod"] for x in e])
        if host:
            P = 1000.0 * cod / cod.sum(axis=1, keepdims=True)
            hv = np.array([host["per1000"][c] for c in rc.SENSE])
            df["codon_usage_L1_to_host"] = np.abs(P - hv).sum(axis=1) / 10.0
            df["GC_pos_L1_to_host"] = np.abs(df["GC1"] - host["gc1"]) + np.abs(df["GC2"] - host["gc2"]) + np.abs(gc3 - host["gc3"])
            df["ENC_abs_diff_to_host"] = (src["ENC"] - rc.HOST_ENC).abs()
        if hdin:
            df["dinuc_total_L1_to_host"] = np.abs(np.vstack([x["dt"] for x in e]) - np.array([hdin[0][d] for d in rc.DINUC])).sum(axis=1)
            df["dinuc_junction_L1_to_host"] = np.abs(np.vstack([x["dj"] for x in e]) - np.array([hdin[1][d] for d in rc.DINUC])).sum(axis=1)
    if host:
        gene_mean = lambda g: np.mean([data[s]["extra"]["cod"] for s in rows[g]], axis=0)        # noqa: E731
        plot_host_codon(host, np.sum([ga_res[g]["extra"]["cod"] for g in genes], axis=0),
                        np.sum([gene_mean(g) for g in genes], axis=0), out / "host_reference" / "codon_usage_per_1000_vs_host.png", a.dpi)
        for m, v, lab in (("GC_pct", host["gc"], "human host mean GC (Human.txt)"), ("GC3_pct", host["gc3"], "human host GC3 (Human.txt)"),
                          ("GC1", host["gc1"], "human host GC1 (Human.txt)"), ("GC2", host["gc2"], "human host GC2 (Human.txt)"),
                          ("ENC", rc.HOST_ENC, "human host ENc {:.1f} (CoCoPUTs)".format(rc.HOST_ENC))):
            ps.REF_LINES.setdefault(m, []).append((v, lab, "#444444"))
    if hdin:
        plot_host_dinuc(hdin, [np.vstack([ga_res[g]["extra"][k] for g in genes]) for k in ("dt", "dj")],
                        [np.vstack([data[s]["extra"][k] for s in dna["_seq"]]) for k in ("dt", "dj")],
                        out / "host_reference" / "dinucleotides_vs_host.png", a.dpi)

    # ---- distribution plots of the extra columns (same plot as plot_scores_all) ----------------------------------------
    extras = [m for m in EXTRA_DIR if m in ex_dna.columns]
    ps.DIRECTION.update({m: EXTRA_DIR[m] for m in extras})
    st, _, _ = ps.analyse(ex_dna, ex_ga, extras, genes, 3)
    if "ENC" not in ps.DIRECTION:
        pass
    # ENC itself (with the host line) is already in plot_scores_all's distribution folder; the host line is added here
    # only for the ENC / GC_pct / GC3_pct images of THIS run: re-draw them so the host reference is visible.
    redo = [m for m in ("GC_pct", "GC3_pct", "ENC") if m in sc_dna.columns]
    st_main, _, _ = ps.analyse(sc_dna, sc_ga, redo, genes, 3)
    (out / "extras").mkdir(exist_ok=True)
    for i, m in enumerate(extras, 1):
        ps.plot_distribution(m, ex_dna, ex_ga, genes, st, out / "extras" / "X{:02d}_{}.png".format(i, m), a.dpi)
    for m in redo:
        ps.plot_distribution(m, sc_dna, sc_ga, genes, st_main, out / "extras" / "host_line_{}.png".format(m), a.dpi)

    # ---- optional per-pair figures --------------------------------------------------------------------------------------
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

    # ---- tables --------------------------------------------------------------------------------------------------------
    with pd.ExcelWriter(out / "profile_host_summary.xlsx", engine="openpyxl") as w:
        pd.DataFrame({"sheet": ["Check_vs_scores_all", "Extras_GeneArt_vs_DNAChisel", "DNAChisel_extras", "GeneArt_extras"],
                      "content": ["scalar scores recomputed by this script vs scores_all.xlsx (differences must be 0)",
                                  "per extra column x gene: GeneArt value, DNA Chisel min/median/max, percentile, cluster",
                                  "per DNA Chisel row: GC1, GC2 and distances to the host reference",
                                  "the same for the GeneArt sequences"]}).to_excel(w, sheet_name="Read_me", index=False)
        chk.to_excel(w, sheet_name="Check_vs_scores_all", index=False)
        st.drop(columns=["_code", "_bounds"]).to_excel(w, sheet_name="Extras_GeneArt_vs_DNAChisel", index=False)
        ex_dna.to_excel(w, sheet_name="DNAChisel_extras", index=False)
        ex_ga.to_excel(w, sheet_name="GeneArt_extras", index=False)
        for ws in w.book.worksheets:
            ws.freeze_panes = "A2"
            for col in ws.columns:
                ws.column_dimensions[col[0].column_letter].width = min(40, max(10, len(str(col[0].value)) + 2))
    print("Done. Outputs in {}".format(out))


if __name__ == "__main__":
    main()
