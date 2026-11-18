import pandas as pd
import numpy as np
import networkx as nx
import re
import plotly.graph_objects as go

from itertools import combinations
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.enum.text import PP_ALIGN
from pptx.dml.color import RGBColor

# ======================================================
# LOAD NETWORK
# ======================================================
string_network = pd.read_table("STRING_PPI_filtered.txt")
G = nx.from_pandas_edgelist(string_network, 'source', 'target', edge_attr='Weight')

# ======================================================
# LOAD INPUTS
# ======================================================
combination = pd.read_excel(
    'Drug_combination_results.xlsx',
    sheet_name='Drug_combination_results'
)
modules = pd.read_excel(
    'Drug_combination_results.xlsx',
    sheet_name='Module_genes'
)

# ======================================================
# COLORS (UNCHANGED)
# ======================================================
def get_distinct_colors(n):
    return [
        'rgb(0, 0, 255)', 'rgb(0, 128, 0)', 'rgb(0, 255, 255)',
        'rgb(128, 0, 128)', 'rgb(255, 165, 0)', 'rgb(0, 255, 0)',
        'rgb(75, 0, 130)', 'rgb(255, 192, 203)',
    ][:n]

# ======================================================
# TARGET PARSER (UNCHANGED)
# ======================================================
def extract_targets(row):
    targets = []
    for col in ["Drug_target1", "Drug_target2", "Drug_target3"]:
        if col in row and pd.notna(row[col]):
            targets.append([x.strip() for x in str(row[col]).split(",") if x.strip()])
    return targets

# ======================================================
# PPT SETUP
# ======================================================
prs = Presentation()
prs.slide_width = Inches(13.33)
prs.slide_height = Inches(7.5)

layout = prs.slide_layouts[5]  # blank slide

# ======================================================
# MAIN LOOP
# ======================================================
for module in modules:

    protein = modules[module].dropna().to_list()
    subset = combination[combination["Disease"] == module]

    for _, row in subset.iterrows():

        seed_genes = extract_targets(row)
        seed_colors = get_distinct_colors(len(seed_genes))

        # --- Build clusters ---
        clusters = {}
        for i, genes in enumerate(seed_genes):
            nodes = []
            for g in genes:
                if g in G:
                    nodes += [g] + list(G.neighbors(g))
            clusters[f'C{i+1}'] = {
                'name': ",".join(genes),
                'nodes': list(set(nodes)),
                'color': seed_colors[i]
            }

        disease_nodes = [g for g in protein if g in G]
        clusters["Disease"] = {
            'name': 'Disease Module',
            'nodes': disease_nodes,
            'color': 'rgb(255,0,0)'
        }

        # --- Subgraph ---
        all_nodes = set()
        for c in clusters.values():
            all_nodes.update(c['nodes'])

        subgraph = G.subgraph(all_nodes)
        pos = nx.spring_layout(subgraph, seed=42, k=0.5)

        # --- Edges ---
        edge_x, edge_y = [], []
        for u, v in subgraph.edges():
            x0, y0 = pos[u]
            x1, y1 = pos[v]
            edge_x += [x0, x1, None]
            edge_y += [y0, y1, None]

        # --- Nodes ---
        node_x, node_y, node_color, node_size = [], [], [], []
        cluster_rgb = {
            k: tuple(map(int, re.findall(r'\d+', v['color'])))
            for k, v in clusters.items()
        }

        for node in subgraph.nodes():
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)

            membership = [c for c in clusters if node in clusters[c]['nodes']]
            if len(membership) > 1:
                r = sum(cluster_rgb[c][0] for c in membership) // len(membership)
                g = sum(cluster_rgb[c][1] for c in membership) // len(membership)
                b = sum(cluster_rgb[c][2] for c in membership) // len(membership)
                node_color.append(f"rgb({r},{g},{b})")
            else:
                node_color.append(clusters[membership[0]]['color'])

            node_size.append(20 if node in disease_nodes else 12)

        # --- Plot ---
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=edge_x, y=edge_y, mode='lines',
                                 line=dict(width=0.5, color='#AAAAAA'),
                                 showlegend=False))
        fig.add_trace(go.Scatter(x=node_x, y=node_y, mode='markers',
                                 marker=dict(size=node_size,
                                             color=node_color,
                                             opacity=0.85),
                                 showlegend=False))

        # --- Legend (larger font) ---
        for c in clusters.values():
            fig.add_trace(go.Scatter(
                x=[None], y=[None],
                mode='markers',
                marker=dict(size=14, color=c['color']),
                name=f"{c['name']} (n={len(c['nodes'])})",
                showlegend=True
            ))

        fig.update_layout(
            title=f"{module} - {row['Drugs']}",
            legend=dict(font=dict(size=16)),
            width=1400,
            height=900,
            xaxis=dict(visible=False),
            yaxis=dict(visible=False)
        )

        img_name = f"{module}_{row['Drugs']}.png".replace(" ", "_")
        fig.write_image(img_name, scale=3)  # HIGH QUALITY

        # ======================================================
        # PPT SLIDE
        # ======================================================
        slide = prs.slides.add_slide(layout)

        # --- Title ---
        title = slide.shapes.add_textbox(Inches(0.5), Inches(0.2), Inches(12), Inches(0.6))
        p = title.text_frame.paragraphs[0]
        p.text = f"{module} - {row['Drugs']}"
        p.font.size = Pt(28)
        p.font.bold = True

        # --- Table ---
        table = slide.shapes.add_table(
            2, 5, Inches(0.5), Inches(1.0), Inches(12), Inches(1.2)
        ).table

        headers = ["Drug Target 1", "Drug Target 2", "Z (Drug A)", "Z (Drug B)", "Sab"]
        values = [
            row["Drug_target1"], row["Drug_target2"],
            f"{row['Z_score_drugA']:.3f}",
            f"{row['Z_score_drugB']:.3f}",
            f"{row['Sab_score']:.3f}"
        ]

        for i, h in enumerate(headers):
            cell = table.cell(0, i)
            cell.text = h
            cell.text_frame.paragraphs[0].font.bold = True
            cell.text_frame.paragraphs[0].font.size = Pt(16)

        for i, v in enumerate(values):
            cell = table.cell(1, i)
            cell.text = v
            cell.text_frame.paragraphs[0].font.size = Pt(16)

            if i == 4 and row["P_value_separation"] < 0.05:
                cell.fill.solid()
                cell.fill.fore_color.rgb = RGBColor(180, 238, 180)

        # --- Graph ---
        slide.shapes.add_picture(
            img_name, Inches(3.0), Inches(2.3), width=Inches(7.5)
        )

        # --- Classification (bottom center) ---
        box = slide.shapes.add_textbox(
            Inches(3), Inches(6.7), Inches(7), Inches(0.6)
        )
        p = box.text_frame.paragraphs[0]
        p.text = f"Classification: {row['Classification']}"
        p.font.size = Pt(20)
        p.font.bold = True
        p.alignment = PP_ALIGN.CENTER

# ======================================================
# SAVE
# ======================================================
prs.save("Drug_Combination_Report.pptx")
print("PPT generated ? Drug_Combination_Report.pptx")
