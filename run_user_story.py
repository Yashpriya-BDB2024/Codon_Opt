"""
run_user_story.py  --  drives the existing Network_proximity scripts for the user story:

    Dabrafenib (BRAF) + Trametinib (MAP2K1,MAP2K2) in Melanoma
    Cetuximab (EGFR)  + Vemurafenib (BRAF)         in CRC

Put this file in the Network_proximity folder (next to drug_synergy.py) together with:
    Melanoma_Genes.xlsx, CRC_Genes.xlsx, STRING_PPI_filtered.txt
(STRING_PPI_filtered.txt is created by make_string_ppi.py)
It does NOT modify drug_synergy.py or precompute_shortest_paths.py.

Usage:
    python run_user_story.py inputs   # step 1: build Drug_target.xlsx only
    python run_user_story.py run      # steps 1-3: inputs + cache + proximity analysis (slow)
    python run_user_story.py excel    # step 4: CSV -> formatted Excel (after a run)
    python run_user_story.py all      # everything
"""
import os
import subprocess
import sys

import pandas as pd
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

# ======================= CONFIG =======================
BASE = os.path.dirname(os.path.abspath(__file__))

DISEASE_FILES = {                       # disease name -> gene list (first column)
    "Melanoma": "Melanoma_Genes.xlsx",
    "CRC": "CRC_Genes.xlsx",
}

# (disease, drug A, drug A target genes, drug B, drug B target genes)
# Multiple genes: comma, NO spaces (drug_synergy.py does not strip them).
COMBINATIONS = [
    ("Melanoma", "Dabrafenib", "BRAF", "Trametinib", "MAP2K1,MAP2K2"),
    ("CRC", "Cetuximab", "EGFR", "Vemurafenib", "BRAF"),
]

# Gene symbols renamed by HGNC that STRING still lists under the old name.
# (disease-list symbol -> symbol used in STRING).  Verify on genenames.org.
SYMBOL_FIXES = {"POU2AF2": "C11orf53", "POU2AF3": "COLCA2"}

STRING_FILE = "STRING_PPI_filtered.txt"
INFO_FILE = "9606.protein.info.v12.5.txt.gz"   # optional: only used to explain excluded genes
CACHE_FILE = "shortest_paths_cache.pkl"
INPUT_XLSX = "Drug_target.xlsx"
RESULT_CSV = "Drug_combination_results.csv"    # written by drug_synergy.py
OUTPUT_XLSX = "Drug_combination_results.xlsx"

# Documented on the Run_info sheet -- EDIT to match what you actually used / lead confirmed.
RUN_NOTES = {
    "STRING version": "v12.5 (Homo sapiens, 9606)",
    "Network type": "Full network, AB pairs only",
    "Confidence cutoff": "combined_score >= 700 (high confidence)",
    "Edge weight": "TO BE CONFIRMED (see make_string_ppi.py WEIGHT_MODE)",
    "Randomisations": "1000 per test, no fixed random seed",
    "Scripts": "drug_synergy.py / precompute_shortest_paths.py (unmodified)",
}
# ======================================================

MEANING = {
    "Complementary Exposure": "Both drugs act near the disease module, in different network regions",
    "Overlapping Exposure": "Both drugs act near the disease module, in overlapping network regions",
    "Indirect Exposure": "Only one drug acts near the disease module; the two drugs' targets overlap",
    "Single Exposure": "Only one drug acts near the disease module; the two drugs' targets are separate",
    "Independent Action": "Neither drug acts near the disease module; targets are separate",
    "Non-exposure": "Neither drug acts near the disease module; targets overlap",
}


def path(name):
    return os.path.join(BASE, name)


# ---------------------------------------------------------------- step 1
def load_gene_list(file):
    genes = pd.read_excel(path(file)).iloc[:, 0].dropna().astype(str).str.strip()
    return list(dict.fromkeys(genes))   # unique, keep order


def build_input_workbook():
    modules = {}
    for disease, file in DISEASE_FILES.items():
        genes = [SYMBOL_FIXES.get(g, g) for g in load_gene_list(file)]
        changed = [f"{g}->{SYMBOL_FIXES[g]}" for g in load_gene_list(file) if g in SYMBOL_FIXES]
        if changed:
            print(f"[{disease}] renamed to STRING symbols: {', '.join(changed)}")
        modules[disease] = pd.Series(genes)
    module_df = pd.DataFrame(modules)   # shorter columns are padded with blanks

    combo_df = pd.DataFrame(
        [{"Disease": d, "Drug_target1": ta, "Drug_target2": tb}
         for d, _, ta, _, tb in COMBINATIONS])
    for col in ("Drug_target1", "Drug_target2"):
        assert not combo_df[col].str.contains(", ").any(), "No spaces after commas in targets"

    with pd.ExcelWriter(path(INPUT_XLSX)) as xw:
        module_df.to_excel(xw, sheet_name="Module_genes", index=False)
        combo_df.to_excel(xw, sheet_name="Drug_combination", index=False)
    print(f"Built {INPUT_XLSX}: " +
          ", ".join(f"{k}={v.shape[0]} genes" for k, v in modules.items()))


# ---------------------------------------------------------------- step 2-3
def run_script(script):
    print(f"\n>>> running {script} ...")
    r = subprocess.run([sys.executable, script], cwd=BASE)
    if r.returncode != 0:
        sys.exit(f"{script} failed (exit code {r.returncode}). See messages above.")


def run_analysis():
    if not os.path.exists(path(STRING_FILE)):
        sys.exit(f"{STRING_FILE} not found. Run make_string_ppi.py first.")
    if not os.path.exists(path(CACHE_FILE)):
        print(f"{CACHE_FILE} not found -> building it (slow, memory heavy)")
        run_script("precompute_shortest_paths.py")
    else:
        print(f"Using existing {CACHE_FILE}")
    run_script("drug_synergy.py")


# ---------------------------------------------------------------- step 4
def excluded_genes_table():
    """Disease genes that could not be used, with the reason."""
    net = pd.read_csv(path(STRING_FILE), sep="\t")
    nodes = set(net["source"].astype(str)) | set(net["target"].astype(str))
    in_info = None
    if os.path.exists(path(INFO_FILE)):
        in_info = set(pd.read_csv(path(INFO_FILE), sep="\t")["preferred_name"].astype(str))
    rows = []
    for disease, file in DISEASE_FILES.items():
        for orig in load_gene_list(file):
            g = SYMBOL_FIXES.get(orig, orig)
            if g in nodes:
                continue
            if in_info is None:
                reason = "Not in network at the chosen cutoff"
            elif g in in_info:
                reason = "In STRING but has no interactions at the chosen cutoff (isolated)"
            else:
                reason = "Symbol not found in STRING protein info file"
            rows.append({"Disease": disease, "Gene (as in disease list)": orig,
                         "Symbol used": g, "Reason excluded": reason})
    return pd.DataFrame(rows, columns=["Disease", "Gene (as in disease list)",
                                       "Symbol used", "Reason excluded"])


def style_sheet(ws, widths, number_formats=None):
    head_fill = PatternFill("solid", start_color="1F3864")
    for c in ws[1]:
        c.font = Font(name="Arial", bold=True, color="FFFFFF", size=10)
        c.fill = head_fill
        c.alignment = Alignment(wrap_text=True, vertical="center")
    for row in ws.iter_rows(min_row=2):
        for c in row:
            c.font = Font(name="Arial", size=10)
            c.alignment = Alignment(wrap_text=True, vertical="top")
    for i, w in enumerate(widths, 1):
        ws.column_dimensions[get_column_letter(i)].width = w
    for col_letter, fmt in (number_formats or {}).items():
        for c in ws[col_letter][1:]:
            c.number_format = fmt
    ws.freeze_panes = "A2"


def build_excel():
    if not os.path.exists(path(RESULT_CSV)):
        sys.exit(f"{RESULT_CSV} not found - run the analysis first.")
    res = pd.read_csv(path(RESULT_CSV))
    for col in ("Drug_target1", "Drug_target2"):
        res[col] = res[col].astype(str).str.strip()

    names = {(d, ta, tb): (na, nb) for d, na, ta, nb, tb in COMBINATIONS}
    res["Drug_A"] = [names.get((r.Disease, r.Drug_target1, r.Drug_target2), ("", ""))[0]
                     for r in res.itertuples()]
    res["Drug_B"] = [names.get((r.Disease, r.Drug_target1, r.Drug_target2), ("", ""))[1]
                     for r in res.itertuples()]

    def meaning(cls):
        for key, text in MEANING.items():
            if str(cls).replace("Significant ", "").startswith(key):
                return text
        return ""
    res["Meaning"] = res["Classification"].map(meaning)

    ordered = ["Disease", "Drug_A", "Drug_target1", "Drug_B", "Drug_target2",
               "Z_score_drugA", "Z_score_drugB", "Sab_score", "P_value_separation",
               "Classification", "Meaning"]
    res = res[ordered].rename(columns={
        "Drug_A": "Drug A", "Drug_target1": "Drug A target(s)",
        "Drug_B": "Drug B", "Drug_target2": "Drug B target(s)",
        "Z_score_drugA": "Z-score drug A", "Z_score_drugB": "Z-score drug B",
        "Sab_score": "Separation (s_AB)", "P_value_separation": "p-value (separation)"})

    legend = pd.DataFrame({
        "Term": ["Z-score (drug)", "Separation (s_AB)", "p-value (separation)"],
        "Explanation": [
            "Negative = the drug's targets are closer to the disease genes than random (below 0 counts as 'near')",
            "Negative = the two drugs' targets overlap in the network; positive = they are separate",
            "Empirical p-value from 1000 random target sets; < 0.05 = significant"]})
    notes = pd.DataFrame({"Setting": list(RUN_NOTES), "Value": list(RUN_NOTES.values())})

    with pd.ExcelWriter(path(OUTPUT_XLSX), engine="openpyxl") as xw:
        res.to_excel(xw, sheet_name="Results", index=False)
        excluded_genes_table().to_excel(xw, sheet_name="Excluded_genes", index=False)
        notes.to_excel(xw, sheet_name="Run_info", index=False)
        legend.to_excel(xw, sheet_name="Run_info", index=False, startrow=len(notes) + 3)
        style_sheet(xw.sheets["Results"], [12, 13, 16, 13, 16, 12, 12, 13, 14, 38, 60],
                    {"F": "0.00", "G": "0.00", "H": "0.000", "I": "0.000"})
        style_sheet(xw.sheets["Excluded_genes"], [12, 24, 16, 60])
        style_sheet(xw.sheets["Run_info"], [26, 90])
    print(f"Saved {OUTPUT_XLSX}")


# ---------------------------------------------------------------- main
if __name__ == "__main__":
    step = sys.argv[1] if len(sys.argv) > 1 else "all"
    if step in ("inputs", "run", "all"):
        build_input_workbook()
    if step in ("run", "all"):
        run_analysis()
    if step in ("excel", "all"):
        build_excel()
    if step not in ("inputs", "run", "excel", "all"):
        sys.exit(__doc__)
