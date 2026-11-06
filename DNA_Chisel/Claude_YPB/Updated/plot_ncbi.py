# -*- coding: utf-8 -*-
"""
plot_ncbi.py  -  plots of the NCBI reference CDS of each gene (NCBI only)
=========================================================================

Only the NCBI sequences are used: no GeneArt, no DNA Chisel. Same figure types, layout and colours as plot_all.py
(imported from it, not modified), so the images look the same as the ones already made.

Reads
  scores_ncbi*.xlsx   sheet NCBI: the metrics of the NCBI sequence of each gene (same columns as scores_all.xlsx)
  trials workbook     column "NCBI sequence" (DNA for the positional profiles and motif positions)

Outputs (outputs/plots_ncbi/)
  metrics/NN_<metric>.png            one image per score in scores_ncbi, all genes in the image
  profiles/<gene>_profiles_NCBI.png  12 positional panels per gene (GC, MFE, %MinMax, rolling CAI / tAI / rare codons,
                                     GC3, windowed ENC, CpG, hairpin, ACF, fragment Tm)
  liabilities/motif_map_NCBI.png     restriction / splice / homopolymer hits along each CDS
  host_reference/                    codon usage, GC by codon position, dinucleotides: NCBI vs the human host
  summary/enc_vs_gc3.png             ENC against GC3 with Wright's expected curve
  ncbi_tables.xlsx                   the numbers behind the plots + check against scores_ncbi

Run from Sequence_scoring/   (next to plot_all.py and its helper scripts)
    python plot_ncbi.py
    python plot_ncbi.py --ncbi-scores outputs/scores/scores_ncbi_1.xlsx --trials "../Codon_optimization_trials_27102026.xlsx"
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

try:
    import plot_all as pa                   # given, unchanged
except ImportError:                         # the file was saved under another name
    import plot_all_1 as pa                 # type: ignore

cf, cc, rc, sas = pa.cf, pa.cc, pa.rc, pa.sas
HERE = Path(__file__).resolve().parent
ROOT = pa.ROOT
C_NC = "#2E7D32"                            # NCBI colour


def plot_metric(m, genes, nc, path, dpi):
    sub = ("Green square = the NCBI reference CDS of that gene (value printed above it), one row per gene, same x-scale for "
           "every gene. Dashed lines / shaded band = reference values for this metric (audit thresholds, human GC band, human host).")
    fig, ax, axt = pa._frame(m, genes, sub, [4.2, 0.05], [])
    for i, g in enumerate(genes):
        v = float(nc.loc[nc["Name"] == g, m].iloc[0])
        ax.scatter([v], [i], marker="s", s=80, color=C_NC, edgecolor="white", linewidths=1, zorder=6)
        ax.annotate(pa.fmt(v), (v, i), xytext=(0, 12), textcoords="offset points", ha="center", fontsize=8, color=C_NC,
                    fontweight="bold", zorder=7)
    ax.set_yticklabels(genes, fontsize=9.5)
    ax.margins(x=0.10)
    h = [Line2D([], [], marker="s", ls="", color=C_NC, markeredgecolor="white", markersize=9, label="NCBI")] + pa.ref_handles(m)
    fig.legend(handles=h, loc="lower center", bbox_to_anchor=(0.45, 0.035), ncol=4, fontsize=8.5, frameon=False)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def panel_list():
    return [(k, t, x) for k, t, x, _ in pa.PANELS] + [("acf", "Autocorrelation of per-codon tAI (ACF)", "lag"),
                                                      ("tm", "Fragment Tm (60-nt fragments, \u00b0C)", "fragment")]


def plot_profiles(gene, res, path, dpi):
    keep = [p for p in panel_list() if res["prof"][p[0]] is not None]
    ncol = 2
    nrow = int(np.ceil(len(keep) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(16, 2.7 * nrow + 1.8))
    axs = np.atleast_1d(axs).ravel()
    refs = set()
    for ax, (key, title, xl) in zip(axs, keep):
        gi, gv = res["prof"][key]
        if key in ("acf", "tm"):
            ax.bar(gi, gv, 0.6, color=C_NC)
        else:
            ax.plot(gi, gv, color=C_NC, lw=1.7)
        if key == "gc":
            ax.axhspan(pa.BAND[0], pa.BAND[1], color="#2E7D32", alpha=0.07, zorder=0)
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
    h = [Line2D([], [], color=C_NC, lw=2, label="NCBI")]
    if "band" in refs:
        h.append(Patch(color="#2E7D32", alpha=0.2, label="human GC band {:g}\u2013{:g}%".format(*pa.BAND)))
    if "cai" in refs:
        h.append(Line2D([], [], color="#6A1B9A", ls="--", label="audit: CAI > 0.85"))
    if "rare" in refs:
        h.append(Line2D([], [], color="#6A1B9A", ls="--", label="audit: rare codons > 5% problematic"))
    top = 0.935
    fig.suptitle("{}: positional profiles, NCBI sequence".format(gene), fontsize=13, fontweight="bold", x=0.01, ha="left", y=0.995)
    sub = "Green = NCBI reference CDS. Same panels as the positional-profile figure of plot_all.py; axes scaled to this sequence."
    fig.text(0.01, 0.962, sub, fontsize=8.5, color="#444", va="center")
    fig.legend(handles=h, loc="upper left", ncol=5, frameon=False, fontsize=8.5, bbox_to_anchor=(0.005, top + 0.003))
    fig.tight_layout(rect=[0, 0, 1, top - 0.012])
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_motif_map(genes, res, path, dpi, bin_nt=30):
    cats = [("restriction", "Restriction site"), ("splice", "Splice motif"), ("homopolymer", "Homopolymer \u22656")]
    fig, axs = plt.subplots(len(genes), 1, figsize=(14, 2.0 * len(genes) + 1.6))
    axs = np.atleast_1d(axs)
    for ax, g in zip(axs, genes):
        L = res[g]["length"]
        nb = int(np.ceil(L / bin_nt))
        ax.set_xlim(0, nb * bin_nt)
        ax.set_ylim(2.5, -0.5)
        counts = []
        for ci, (k, _) in enumerate(cats):
            hits = res[g]["pos"][k]
            counts.append(len(hits))
            if hits:
                ax.scatter([h["pos"] for h in hits], [ci] * len(hits), marker="D", s=46, color=C_NC, edgecolor="white",
                           linewidths=0.8, zorder=5)
            if k == "restriction":
                for h in hits:
                    ax.annotate(h["name"], (h["pos"], ci), xytext=(0, 8), textcoords="offset points", rotation=60,
                                fontsize=6.5, color=C_NC, ha="left", va="bottom")
            ax.axhline(ci, color="#ddd", lw=0.8, zorder=0)
        ax.set_yticks(range(3))
        ax.set_yticklabels([c[1] for c in cats], fontsize=8.5)
        ax.set_title("{}  ({:,} nt):  NCBI has {} restriction / {} splice / {} homopolymer hits".format(g, res[g]["length"], *counts),
                     fontsize=9.5, loc="left", fontweight="bold", pad=14)
        ax.set_xlabel("nt position" if g == genes[-1] else "")
    fig.tight_layout(rect=[0, 0.06, 0.93, 0.96])
    fig.legend(handles=[Line2D([], [], marker="D", ls="", color=C_NC, markeredgecolor="white", markersize=8,
                               label="NCBI hit (restriction: enzyme name; either strand, start position)")],
               loc="lower center", frameon=False, fontsize=9)
    fig.suptitle("Where the sequence liabilities sit along each CDS (NCBI sequences)", fontsize=13, fontweight="bold", x=0.01, ha="left")
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_host_codon(host, vec, path, dpi):
    order = [c for c in [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"] if c in rc.SENSE]
    order.sort(key=lambda c: (cf.GENETIC_CODE[c], c))
    pos = {c: i for i, c in enumerate(rc.SENSE)}
    v = 1000.0 * vec / vec.sum()
    x = np.arange(len(order))
    fig, ax = plt.subplots(figsize=(19, 5.4))
    ax.bar(x - 0.28, [host["per1000"][c] for c in order], 0.27, color=pa.C_HOST, label="Human host (CoCoPUTs, Human.txt)")
    ax.bar(x, [v[pos[c]] for c in order], 0.27, color=C_NC, label="NCBI (8 genes pooled)")
    ax.set_xticks(x)
    ax.set_xticklabels(["{}\n{}".format(cf.GENETIC_CODE[c], c) for c in order], fontsize=6.5)
    ax.set_ylabel("per 1000 sense codons")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.3)
    ax.set_title("Codon usage per 1000 codons, grouped by amino acid", fontsize=12, fontweight="bold", loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_host_gc(host, nc, ex, path, dpi):
    names = ["GC overall", "GC1", "GC2", "GC3"]
    vals = [nc["GC_pct"].to_numpy(float), ex["GC1"].to_numpy(float), ex["GC2"].to_numpy(float), nc["GC3_pct"].to_numpy(float)]
    hv = [host["gc"], host["gc1"], host["gc2"], host["gc3"]]
    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    x = np.arange(4)
    ax.bar(x - 0.28, hv, 0.27, color=pa.C_HOST, label="Human host (Human.txt)")
    ax.bar(x, [v.mean() for v in vals], 0.27, color=C_NC, label="NCBI (mean of the genes; black dots = each gene)")
    for i, v in enumerate(vals):
        ax.scatter(np.full(len(v), i), v, s=14, color="black", zorder=4)
    ax.axhspan(pa.BAND[0], pa.BAND[1], xmin=0, xmax=0.25, color="#2E7D32", alpha=0.08, zorder=0)
    ax.text(-0.45, pa.BAND[1] + 0.6, "human GC band {:g}\u2013{:g}% (overall GC only)".format(*pa.BAND), fontsize=8, color="#2E7D32")
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("GC %")
    ax.set_ylim(0, max([85.0] + [float(v.max()) + 6 for v in vals]))
    ax.legend(frameon=False, fontsize=9, loc="upper left", bbox_to_anchor=(0, -0.08), ncol=1)
    ax.grid(axis="y", alpha=0.3)
    ax.set_title("GC by codon position against the human host", fontsize=12, fontweight="bold", loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_host_dinuc(hd, arr, path, dpi):
    fig, axs = plt.subplots(2, 1, figsize=(13, 8.5))
    for ax, k, ttl in zip(axs, (0, 1), ("Total dinucleotides", "Junction dinucleotides (3rd base of a codon + 1st base of the next)")):
        x = np.arange(16)
        a = arr[k]
        ax.bar(x - 0.28, [hd[k][d] for d in rc.DINUC], 0.27, color=pa.C_HOST, label="Human host")
        ax.bar(x, a.mean(axis=0), 0.27, color=C_NC, label="NCBI (mean of the genes; dots = each gene)")
        ax.scatter(np.repeat(x, len(a)), a.T.ravel(), s=9, color="black", zorder=4)
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


def plot_enc_gc3(nc, genes, host, path, dpi):
    fig, ax = plt.subplots(figsize=(10, 6.4))
    s = np.linspace(0.02, 0.98, 200)
    ax.plot(100 * s, [cf.enc_expected_from_gc3(v) for v in s], color="#333", lw=1.6, zorder=1)
    for gi, g in enumerate(genes):
        r = nc[nc["Name"] == g].iloc[0]
        ax.scatter([r["GC3_pct"]], [r["ENC"]], marker="*", s=240, color=pa.GENE_COL[gi % 8], edgecolor="black", linewidths=1.2, zorder=5)
    if host:
        ax.scatter([host["gc3"]], [rc.HOST_ENC], marker="X", s=150, color="black", zorder=6)
        ax.annotate("human host (GC3 {:.0f}%, ENc {:.1f})".format(host["gc3"], rc.HOST_ENC), (host["gc3"], rc.HOST_ENC),
                    xytext=(8, 8), textcoords="offset points", fontsize=9)
    ax.set_xlabel("GC3 (%)")
    ax.set_ylabel("ENC (Nc): 20 = strongest codon bias, 61 = none")
    ax.set_ylim(15, 63)
    ax.grid(alpha=0.3)
    h = [Line2D([], [], color="#333", lw=1.6, label="Nc expected from GC3 alone (Wright 1990): points below = more bias than GC3 explains")]
    h += [Line2D([], [], marker="*", ls="", color=pa.GENE_COL[i % 8], markeredgecolor="black", markersize=13, label="NCBI " + g)
          for i, g in enumerate(genes)]
    if host:
        h.append(Line2D([], [], marker="X", ls="", color="black", markersize=9, label="human host"))
    fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, fontsize=8.5)
    fig.suptitle("Codon bias (ENC) against GC3: the NCBI sequence of each gene", fontsize=12, fontweight="bold", x=0.01, ha="left")
    fig.tight_layout(rect=[0, 0.13, 1, 0.96])
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def trials_table(path):
    t = pd.read_excel(path, header=2).iloc[:, 1:6]
    t.columns = ["Name", "Protein", "Host", "GeneArt", "NCBI"]
    t = t.dropna(subset=["Name", "Protein"])
    t = t[t["Host"].astype(str).str.lower().str.contains("homo")].reset_index(drop=True)
    t["Protein"] = t["Protein"].astype(str).str.strip().str.upper()
    t["GeneArt"] = t["GeneArt"].map(lambda s: cf.clean_seq(s) if isinstance(s, str) else "")
    t["NCBI"] = t["NCBI"].map(lambda s: cf.clean_seq(s) if isinstance(s, str) else "")
    return t


def find_trials(arg):
    if arg:
        return Path(arg)
    found = sorted(ROOT.glob("Codon*ptimization*trials*.xlsx"), key=lambda p: p.stat().st_mtime)
    if not found:
        sys.exit("trials workbook not found in {}; pass --trials".format(ROOT))
    if len(found) > 1:
        print("several trials workbooks found; using the newest: {}".format(found[-1].name))
    return found[-1]


def main():
    ap = argparse.ArgumentParser(description="Plots of the NCBI sequences only.")
    ap.add_argument("--ncbi-scores", default=None, help="scores_ncbi*.xlsx (sheet NCBI); default: outputs/scores/scores_ncbi*.xlsx")
    ap.add_argument("--trials", default=None, help="trials workbook with the 'NCBI sequence' column (newest match by default)")
    ap.add_argument("--out-dir", default=str(HERE / "outputs" / "plots_ncbi"))
    ap.add_argument("--no-profiles", action="store_true", help="only the metric images (no DNA needed)")
    ap.add_argument("--no-host", action="store_true", help="no human-host reference lines and no host_reference/ plots")
    ap.add_argument("--dpi", type=int, default=130)
    a = ap.parse_args()

    p = a.ncbi_scores
    if p is None:
        cand = sorted((HERE / "outputs" / "scores").glob("scores_ncbi*.xlsx")) + sorted(ROOT.glob("scores_ncbi*.xlsx"))
        if not cand:
            sys.exit("scores_ncbi*.xlsx not found; pass --ncbi-scores")
        p = str(cand[0])
    nc = pd.read_excel(p, sheet_name="NCBI")
    genes = list(dict.fromkeys(nc["Name"]))
    if nc["Name"].duplicated().any():
        sys.exit("More than one NCBI row per gene")
    in_table = [m for m, _ in cc.COMPARE_METRICS if m in nc.columns and m != "kozak_score"]
    ok = [m for m in in_table if np.isfinite(pd.to_numeric(nc[m], errors="coerce").to_numpy(float)).all()]
    print("Metrics: {} plotted | no data (seqfold missing?): {}".format(len(ok), [m for m in in_table if m not in ok] or "none"))

    host = hdin = None
    if not a.no_host:
        host, hdin = rc.load_host_counts(), rc.load_host_dinuc()
        if host is not None:
            for m, v, lab in (("GC_pct", host["gc"], "human host mean GC (Human.txt)"), ("GC3_pct", host["gc3"], "human host GC3 (Human.txt)")):
                pa.REF_LINES.setdefault(m, []).append((v, lab, "#444444"))
        pa.REF_LINES.setdefault("ENC", []).append((rc.HOST_ENC, "human host ENc {:.1f} (CoCoPUTs)".format(rc.HOST_ENC), "#444444"))

    out = Path(a.out_dir)
    for d in ["metrics", "summary"] + ([] if a.no_profiles else ["profiles", "liabilities", "host_reference"]):
        (out / d).mkdir(parents=True, exist_ok=True)

    for i, m in enumerate(in_table, 1):
        if m in ok:
            plot_metric(m, genes, nc, out / "metrics" / "{:02d}_{}.png".format(i, m), a.dpi)
            print("  {:2d}/{} {}".format(i, len(in_table), m), flush=True)
    tables = {"NCBI_scores": nc}
    if host and all(m in nc.columns for m in ("GC3_pct", "ENC")):
        plot_enc_gc3(nc, genes, host, out / "summary" / "enc_vs_gc3.png", a.dpi)

    if not a.no_profiles:
        tr = trials_table(find_trials(a.trials))
        miss = [g for g in genes if g not in set(tr["Name"]) or not tr.loc[tr["Name"] == g, "NCBI"].iloc[0]]
        if miss:
            sys.exit("no NCBI sequence in the trials workbook for {}".format(miss))
        raw = {g: tr.loc[tr["Name"] == g, "NCBI"].iloc[0] for g in genes}
        for g in genes:
            if sas.translate(cf.strip_trailing_stop(raw[g])) != tr.loc[tr["Name"] == g, "Protein"].iloc[0]:
                warnings.warn("NCBI sequence of {} does not translate to the input protein".format(g))
        # scores_ncbi was made from one of two forms of the sequence; use the one that reproduces it
        best = None
        for label, seqs in (("as in the workbook (stop codon kept)", dict(raw)),
                            ("trailing stop removed", {g: cf.strip_trailing_stop(s) for g, s in raw.items()})):
            data = pa.score_all(list(dict.fromkeys(seqs.values())), 1, out / "ncbi_cache.pkl")
            frame = nc[["Name"] + in_table].copy()
            frame["_seq"] = frame["Name"].map(seqs)
            chk = pa.check_against_table(frame, frame.iloc[0:0], data, in_table)
            nbad = int((chk["differences_gt_1e-6"] > 0).sum())
            print("  NCBI sequence {}: {} metric(s) differ from the score table".format(label, nbad))
            if best is None or nbad < best[0]:
                best = (nbad, label, seqs, data, chk)
            if nbad == 0:
                break
        nbad, label, seqs, data, chk = best
        if nbad:
            warnings.warn("scores computed here still differ from scores_ncbi for {} metric(s); check human_codon.csv, "
                          "human_tgcn.csv and seqfold".format(nbad))
        print("Profiles use the NCBI sequence {}.".format(label))
        res = {g: data[seqs[g]] for g in genes}
        for g in genes:
            plot_profiles(g, res[g], out / "profiles" / "{}_profiles_NCBI.png".format(g), a.dpi)
            print("  profiles", g, flush=True)
        plot_motif_map(genes, res, out / "liabilities" / "motif_map_NCBI.png", a.dpi)
        tables["Check_vs_scores_ncbi"] = chk
        ex = pd.DataFrame({"Name": genes, "GC1": [res[g]["extra"]["gc"][0] for g in genes], "GC2": [res[g]["extra"]["gc"][1] for g in genes]})
        tables["GC_by_position"] = ex
        if host:
            plot_host_codon(host, np.sum([res[g]["extra"]["cod"] for g in genes], axis=0),
                            out / "host_reference" / "codon_usage_per_1000_vs_host.png", a.dpi)
            plot_host_gc(host, nc.set_index("Name").loc[genes].reset_index(), ex, out / "host_reference" / "gc_by_codon_position_vs_host.png", a.dpi)
        if hdin:
            plot_host_dinuc(hdin, [np.vstack([res[g]["extra"][k] for g in genes]) for k in ("dt", "dj")],
                            out / "host_reference" / "dinucleotides_vs_host.png", a.dpi)
        tables["Restriction_per_enzyme"] = pd.DataFrame([{"enzyme": e, "site": cf.RESTRICTION_SITES[e], "gene": g,
                                                          "NCBI_sites": int(res[g]["scores"].get("RS_" + e, 0) or 0)}
                                                         for e in cf.RESTRICTION_SITES for g in genes])
        tables["Motif_hits"] = pd.DataFrame([{"gene": g, "type": k, **h} for g in genes for k in ("restriction", "splice", "homopolymer")
                                            for h in res[g]["pos"][k]])

    with pd.ExcelWriter(out / "ncbi_tables.xlsx", engine="openpyxl") as w:
        for name, df in tables.items():
            df.to_excel(w, sheet_name=name[:31], index=False)
        for ws in w.book.worksheets:
            ws.freeze_panes = "B2"
            for col in ws.columns:
                ws.column_dimensions[col[0].column_letter].width = min(50, max(10, len(str(col[0].value)) + 2))
    print("Done. {} PNG files and ncbi_tables.xlsx in {}".format(len(list(out.rglob("*.png"))), out))


if __name__ == "__main__":
    main()
