"""
Check that disease genes and drug targets exist in STRING_PPI_filtered.txt.

Run from the folder containing the files:
    python check_genes_in_string.py

Edit the CONFIG block below if your file names differ.
"""
import pandas as pd

# ---------------- CONFIG ----------------
STRING_FILE = "STRING_PPI_filtered.txt"

# disease name -> Excel file with the gene list in the first column
DISEASE_FILES = {
    "Melanoma": "Melanoma Genes.xlsx",
    "CRC": "CRC_Genes.xlsx",
}

# disease -> (drug A targets, drug B targets)  -- gene symbols, no spaces
COMBINATIONS = {
    "Melanoma": [("Dabrafenib", ["BRAF"]), ("Trametinib", ["MAP2K1", "MAP2K2"])],
    "CRC": [("Cetuximab", ["EGFR"]), ("Vemurafenib", ["BRAF"])],
}
OUTPUT_CSV = "gene_coverage_check.csv"
# ----------------------------------------


def main():
    net = pd.read_csv(STRING_FILE, sep="\t")
    nodes = set(net["source"].astype(str)) | set(net["target"].astype(str))
    # number of direct neighbours per gene (used to flag poorly connected genes)
    degree = pd.concat([net["source"], net["target"]]).astype(str).value_counts()
    print(f"Network: {len(nodes)} genes, {len(net)} edges\n")

    rows = []
    for disease, path in DISEASE_FILES.items():
        genes = pd.read_excel(path).iloc[:, 0].dropna().astype(str).str.strip()
        genes = list(dict.fromkeys(genes))  # unique, keep order
        missing = [g for g in genes if g not in nodes]
        print(f"[{disease}] disease genes: {len(genes)} | in network: "
              f"{len(genes) - len(missing)} | missing: {len(missing)}")
        if missing:
            print("   missing:", ", ".join(missing))
        for g in genes:
            rows.append({"Disease": disease, "Type": "Disease gene", "Name": "",
                         "Gene": g, "In_network": g in nodes,
                         "Degree": int(degree.get(g, 0))})

        for drug, targets in COMBINATIONS[disease]:
            for t in targets:
                ok = t in nodes
                print(f"   drug target {drug} -> {t}: "
                      f"{'OK (degree ' + str(int(degree[t])) + ')' if ok else 'MISSING'}")
                rows.append({"Disease": disease, "Type": "Drug target", "Name": drug,
                             "Gene": t, "In_network": ok,
                             "Degree": int(degree.get(t, 0))})
        print()

    out = pd.DataFrame(rows)
    out.to_csv(OUTPUT_CSV, index=False)
    print(f"Full report saved to {OUTPUT_CSV}")
    if not out[(out.Type == "Drug target")]["In_network"].all():
        print("WARNING: at least one drug target is missing from the network - "
              "the proximity run will fail or be unreliable for that drug.")


if __name__ == "__main__":
    main()