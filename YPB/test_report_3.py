"""
PowerPoint report for drug_synergy_2.py
=======================================
Reads the OUTPUT of drug_synergy_2.py (Drug_combination_results.xlsx) and the
INPUT module genes (Drug_target.xlsx) and writes Drug_Combination_Report.pptx
with one slide per drug pair (+ one summary slide).

Every change compared with the old test_report_2.py is marked  # [CHANGED n]
and listed in the table at the bottom of this file.
"""
import os
import re

import networkx as nx
import pandas as pd
import plotly.graph_objects as go

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.util import Inches, Pt

# ======================================================
# SETTINGS
# ======================================================
NETWORK_FILE = "STRING_PPI_filtered.txt"
RESULTS_FILE = "Drug_combination_results.xlsx"   # written by drug_synergy_2.py
INPUT_FILE = "Drug_target.xlsx"                  # same input as drug_synergy_2.py
OUT_PPTX = "Drug_Combination_Report.pptx"
IMG_DIR = "report_images"
ALPHA = 0.05                                     # same alpha as drug_synergy_2.py
SUMMARY_SLIDE = True

# ======================================================
# LOAD NETWORK  (same loading as drug_synergy_2.py main())
# ======================================================
string_network = pd.read_csv(NETWORK_FILE, sep="\t")
# [CHANGED 1] edge_attr=True (as in drug_synergy_2.py) instead of only 'Weight'
G = nx.from_pandas_edgelist(string_network, source="source", target="target", edge_attr=True)

# ======================================================
# LOAD INPUTS
# ======================================================
# [CHANGED 2] results and module genes now come from two different files.
# Before, both sheets were read from Drug_combination_results.xlsx, but that file
# now holds only the results sheet; Module_genes lives in Drug_target.xlsx.
results = pd.read_excel(RESULTS_FILE, sheet_name="Drug_combination_results")
modules = pd.read_excel(INPUT_FILE, sheet_name="Module_genes")

os.makedirs(IMG_DIR, exist_ok=True)

# ======================================================
# HELPERS
# ======================================================
def split_genes(cell):
    """'MAP2K1, MAP2K2' -> ['MAP2K1', 'MAP2K2']  (same cleaning as drug_synergy_2.py, item 8)"""
    if pd.isna(cell):
        return []
    return [g.strip() for g in str(cell).split(",") if g.strip()]


def to_float(x):
    """Result cells can be 'ERROR' / 'NA' for failed pairs -> NaN."""
    return pd.to_numeric(x, errors="coerce")


def fmt(x, nd=3):
    return "n/a" if pd.isna(x) else f"{x:.{nd}f}"


def fmt_p(p):
    # p can never be below 1/(1+1000) = 0.000999 (drug_synergy_2.py, item 7)
    if pd.isna(p):
        return "n/a"
    return "<0.001" if p < 0.001 else f"{p:.3f}"


def safe_name(text):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text)).strip("_")


def nearest_path(graph, source, targets):
    """Shortest path from `source` to the closest node in `targets` (None if unreachable)."""
    # [CHANGED 12] the source itself counts: a target that is also a disease gene has
    # distance 0 (no connector path). Before, the source was removed first, so a 1-hop
    # path to a neighbour was drawn although the proximity distance is 0.
    if source in targets:
        return [source]
    targets = set(targets)
    best = None
    # BFS layer by layer until the first target is met
    seen = {source}
    frontier = {source: [source]}
    while frontier and best is None:
        nxt = {}
        for node, path in frontier.items():
            for nb in graph.neighbors(node):
                if nb in seen:
                    continue
                seen.add(nb)
                nxt[nb] = path + [nb]
                if nb in targets and best is None:
                    best = path + [nb]
        frontier = nxt
    return best


# ======================================================
# FIGURE  [CHANGED 3]
# ======================================================
# Old code: each drug target + ALL its STRING neighbours was drawn as a "cluster".
# drug_synergy_2.py (item 3) no longer expands targets, so the picture must not
# either, otherwise it would show something the numbers were not calculated on.
# New picture shows exactly the sets used in the analysis:
#   red    = disease module genes present in the network
#   blue   = Drug A targets      orange = Drug B targets   (large, labelled)
#   grey   = genes on the shortest path from each target to its nearest disease
#            gene (this is what the proximity z-score measures)
COL_DISEASE = "rgb(220,50,47)"
COL_A = "rgb(31,119,180)"
COL_B = "rgb(255,127,14)"
COL_PATH = "rgb(150,150,150)"


def build_figure(disease, drug_a, drug_b, targets_a, targets_b, disease_genes):
    ta = [g for g in targets_a if g in G]
    tb = [g for g in targets_b if g in G]
    dis = [g for g in disease_genes if g in G]
    dis_set = set(dis)

    # connector nodes: shortest path target -> nearest disease gene
    path_nodes = set()
    for t in ta + tb:
        p = nearest_path(G, t, dis_set)
        if p:
            path_nodes.update(p)

    nodes = set(dis) | set(ta) | set(tb) | path_nodes
    sub = G.subgraph(nodes)
    pos = nx.spring_layout(sub, seed=42, k=0.5)

    edge_x, edge_y = [], []
    for u, v in sub.edges():
        edge_x += [pos[u][0], pos[v][0], None]
        edge_y += [pos[u][1], pos[v][1], None]

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=edge_x, y=edge_y, mode="lines",
                             line=dict(width=0.4, color="#CCCCCC"),
                             hoverinfo="skip", showlegend=False))

    def layer(node_list, color, size, name, labels=False, border=None):
        node_list = [n for n in node_list if n in pos]
        if not node_list:
            return
        fig.add_trace(go.Scatter(
            x=[pos[n][0] for n in node_list],
            y=[pos[n][1] for n in node_list],
            mode="markers+text" if labels else "markers",
            text=node_list if labels else None,
            textposition="top center",
            textfont=dict(size=18, color="black"),
            marker=dict(size=size, color=color, opacity=0.9,
                        line=dict(width=border[0], color=border[1]) if border else dict(width=0)),
            name=name, showlegend=True))

    targets_all = set(ta) | set(tb)
    layer([n for n in path_nodes if n not in dis_set and n not in targets_all],
          COL_PATH, 8, f"Path to disease (n={len(path_nodes - dis_set - targets_all)})")
    layer([n for n in dis if n not in targets_all],
          COL_DISEASE, 11, f"Disease module (n={len(dis)})")
    # targets drawn last (on top). A red outline = target is also a disease gene.
    for lst, col, name in ((ta, COL_A, f"{drug_a} targets"), (tb, COL_B, f"{drug_b} targets")):
        layer(lst, col, 24, f"{name} (n={len(lst)})", labels=True, border=(3, "rgb(220,50,47)"))

    fig.update_layout(
        title=dict(text=f"{disease}: {drug_a} + {drug_b}", font=dict(size=24)),
        legend=dict(font=dict(size=18)),
        width=1300, height=760,
        plot_bgcolor="white",
        margin=dict(l=20, r=20, t=70, b=20),
        xaxis=dict(visible=False), yaxis=dict(visible=False))
    return fig, len(dis), ta, tb


# ======================================================
# PPT HELPERS
# ======================================================
prs = Presentation()
prs.slide_width = Inches(13.33)
prs.slide_height = Inches(7.5)
# [CHANGED 4] layout 6 = truly blank. Layout 5 is "Title Only" and left an empty
# title placeholder on every slide ("Click to add title") behind the text box.
BLANK = prs.slide_layouts[6]

GREEN = RGBColor(180, 238, 180)
AMBER = RGBColor(255, 224, 150)
GREY = RGBColor(235, 235, 235)
CLASS_FILL = {                                   # [CHANGED 5]
    "Complementary Exposure": RGBColor(198, 239, 206),
    "Overlapping Exposure": RGBColor(255, 235, 156),
}


def add_text(slide, text, x, y, w, h, size, bold=False, align=None, color=None):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(size)
    p.font.bold = bold
    if color:
        p.font.color.rgb = color
    if align:
        p.alignment = align
    return tb


def set_cell(cell, text, size=14, bold=False, fill=None, align=PP_ALIGN.CENTER):
    cell.text = text
    p = cell.text_frame.paragraphs[0]
    p.font.size = Pt(size)
    p.font.bold = bold
    p.alignment = align
    cell.vertical_anchor = MSO_ANCHOR.MIDDLE
    if fill is not None:
        cell.fill.solid()
        cell.fill.fore_color.rgb = fill


def close_fill(z, p):
    """green = close (z<0, p<alpha); amber = z<0 but weak p (borderline, as the log warned)"""
    if pd.isna(z) or z >= 0:
        return None
    return GREEN if (not pd.isna(p) and p < ALPHA) else AMBER


# ======================================================
# SUMMARY SLIDE  [CHANGED 6]
# ======================================================
def add_summary_slide(df):
    slide = prs.slides.add_slide(BLANK)
    add_text(slide, "Drug combination summary", 0.5, 0.2, 12, 0.6, 28, bold=True)
    cols = ["Disease", "Drug A", "Drug B", "Z (A)", "Z (B)", "Sab", "Classification"]
    tbl = slide.shapes.add_table(len(df) + 1, len(cols), Inches(0.5), Inches(1.1),
                                 Inches(12.3), Inches(0.5 * (len(df) + 1))).table
    widths = [1.6, 2.0, 2.0, 1.3, 1.3, 1.1, 3.0]
    for i, w in enumerate(widths):
        tbl.columns[i].width = Inches(w)
    for i, c in enumerate(cols):
        set_cell(tbl.cell(0, i), c, 16, bold=True)
    for r, (_, row) in enumerate(df.iterrows(), start=1):
        za, zb = to_float(row["Z_score_drugA"]), to_float(row["Z_score_drugB"])
        pa, pb = to_float(row["Pvalue_drugA"]), to_float(row["Pvalue_drugB"])
        sab = to_float(row["Sab_score"])
        cls = str(row["Classification"])
        vals = [row["Disease"], row["Drug_A"], row["Drug_B"], fmt(za, 2), fmt(zb, 2), fmt(sab, 2), cls]
        fills = [None, None, None, close_fill(za, pa), close_fill(zb, pb), None, CLASS_FILL.get(cls)]
        for c, (v, f) in enumerate(zip(vals, fills)):
            set_cell(tbl.cell(r, c), str(v), 14, fill=f)
    add_text(slide,
             "Z: drug-to-disease proximity (z < 0 = close; green p < 0.05, amber p >= 0.05 borderline).  "
             "Sab: separation of the two drugs' targets (< 0 overlapping, >= 0 separated).  "
             "Targets only (no neighbour expansion).",
             0.5, 6.4, 12.3, 0.8, 12, color=RGBColor(90, 90, 90))


# ======================================================
# MAIN LOOP
# ======================================================
# [CHANGED 7] loop over the result rows (one per pair) and look the disease genes
# up in Module_genes. Before it looped over modules and filtered 'combination'.
if SUMMARY_SLIDE and len(results):
    add_summary_slide(results)

for _, row in results.iterrows():
    disease = row["Disease"]
    # [CHANGED 8] drug names: old column 'Drugs' no longer exists -> Drug_A / Drug_B
    drug_a = str(row.get("Drug_A", "") or "Drug A")
    drug_b = str(row.get("Drug_B", "") or "Drug B")
    if drug_a in ("nan", ""):
        drug_a = "Drug A"
    if drug_b in ("nan", ""):
        drug_b = "Drug B"

    if disease not in modules.columns:
        print(f"WARNING: disease '{disease}' not found in Module_genes sheet -> slide skipped")
        continue
    # [CHANGED 9] strip names exactly like drug_synergy_2.py (item 8)
    disease_genes = [str(g).strip() for g in modules[disease].dropna().to_list()]

    targets_a = split_genes(row["Drug_target1"])
    targets_b = split_genes(row["Drug_target2"])

    za, zb = to_float(row["Z_score_drugA"]), to_float(row["Z_score_drugB"])
    pa, pb = to_float(row["Pvalue_drugA"]), to_float(row["Pvalue_drugB"])
    sab = to_float(row["Sab_score"])
    classification = str(row["Classification"])
    failed = classification.startswith("Error")

    slide = prs.slides.add_slide(BLANK)
    add_text(slide, f"{disease}: {drug_a} + {drug_b}", 0.5, 0.2, 12.3, 0.6, 28, bold=True)

    # --- Table [CHANGED 10]: drug names, targets, z with p, Sab ---
    headers = ["Drug A", "Targets A", "Drug B", "Targets B", "Z (A) [p]", "Z (B) [p]", "Sab"]
    values = [drug_a, ", ".join(targets_a), drug_b, ", ".join(targets_b),
              f"{fmt(za)}  [{fmt_p(pa)}]", f"{fmt(zb)}  [{fmt_p(pb)}]", fmt(sab)]
    fills = [None, None, None, None, close_fill(za, pa), close_fill(zb, pb), None]
    tbl = slide.shapes.add_table(2, 7, Inches(0.5), Inches(1.0), Inches(12.3), Inches(1.1)).table
    for i, w in enumerate([1.6, 1.9, 1.6, 1.9, 2.0, 2.0, 1.3]):
        tbl.columns[i].width = Inches(w)
    for i, h in enumerate(headers):
        set_cell(tbl.cell(0, i), h, 15, bold=True)
    for i, (v, f) in enumerate(zip(values, fills)):
        set_cell(tbl.cell(1, i), v, 15, fill=f)

    # --- Graph [CHANGED 3] ---
    if failed:
        add_text(slide, "Analysis failed for this pair - see log of drug_synergy_2.py",
                 0.5, 3.5, 12.3, 1, 20, align=PP_ALIGN.CENTER)
    else:
        fig, n_dis, ta, tb = build_figure(disease, drug_a, drug_b, targets_a, targets_b, disease_genes)
        img = os.path.join(IMG_DIR, safe_name(f"{disease}_{drug_a}_{drug_b}") + ".png")
        fig.write_image(img, scale=3)
        # [CHANGED 11] picture sized by HEIGHT so it never runs into the classification
        # box (old: width 7.5in on a 1400x900 figure = 4.8in tall from y=2.3 -> 7.1in,
        # overlapping the box at y=6.7in)
        pic_h = 4.2
        pic_w = pic_h * 1300 / 760
        slide.shapes.add_picture(img, Inches((13.33 - pic_w) / 2), Inches(2.3), height=Inches(pic_h))

        missing = [g for g in targets_a + targets_b if g not in G]
        if missing:
            add_text(slide, f"Not in STRING network: {', '.join(missing)}", 0.5, 6.5, 4, 0.4, 11,
                     color=RGBColor(180, 0, 0))

    # --- Classification (bottom centre), coloured by class [CHANGED 5] ---
    box = slide.shapes.add_shape(1, Inches(3.4), Inches(6.7), Inches(6.5), Inches(0.6))
    box.fill.solid()
    box.fill.fore_color.rgb = CLASS_FILL.get(classification, GREY)
    box.line.fill.background()
    tf = box.text_frame
    tf.paragraphs[0].text = f"Classification: {classification}"
    tf.paragraphs[0].font.size = Pt(20)
    tf.paragraphs[0].font.bold = True
    tf.paragraphs[0].font.color.rgb = RGBColor(0, 0, 0)
    tf.paragraphs[0].alignment = PP_ALIGN.CENTER
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE

# ======================================================
# SAVE
# ======================================================
prs.save(OUT_PPTX)
print(f"PPT generated -> {OUT_PPTX}")

# ----------------------------------------------------------------------------
# WHAT CHANGED vs test_report_2.py
# 1  edge_attr=True, same network loading as drug_synergy_2.py.
# 2  Module genes now read from Drug_target.xlsx (results file has only 1 sheet).
# 3  Graph: no neighbour expansion (drug_synergy_2.py item 3). Shows disease
#    module + drug targets + shortest path target -> nearest disease gene.
#    Old colour-averaging of overlapping clusters removed (no clusters now).
# 4  Blank layout (6) instead of Title-Only (5): no empty title placeholder.
# 5  Classification box coloured by class; old code highlighted Sab using
#    P_value_separation, which no longer exists (item 4 removed it).
# 6  New summary slide with all pairs (switch off: SUMMARY_SLIDE = False).
# 7  Loop over result rows; module genes looked up by Disease.
# 8  'Drugs' column -> Drug_A / Drug_B (item 11).
# 9  Gene names stripped like drug_synergy_2.py (item 8).
# 10 Table shows names, targets, z AND p (Pvalue_drugA/B, item 7); z cells green
#    (z<0, p<0.05) or amber (z<0, p>=0.05 = the "borderline" warning in the log).
# 12 Connector path: a target that is itself a disease gene gets no path (distance 0).
# 11 Picture size fixed so it does not overlap the classification box; failed
#    pairs ("Error: ...") and genes missing from STRING are shown, not crashed on.
# ----------------------------------------------------------------------------
