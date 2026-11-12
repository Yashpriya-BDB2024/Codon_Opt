"""
Build STRING_PPI_filtered.txt (source, target, Weight) for drug_synergy.py

Inputs  (download from string-db.org -> Download, Homo sapiens 9606):
  9606.protein.links.v12.5.min700.onlyAB.tsv.gz   (full network, >=0.700, AB only, TSV)
  9606.protein.info.v12.5.txt.gz                  (protein ID -> gene symbol)
Output:
  STRING_PPI_filtered.txt   tab-separated: source, target, Weight

Weight is used by precompute_shortest_paths.py as an edge DISTANCE
(shortest path = smallest total weight), so STRONGER interactions must get SMALLER values.

WEIGHT_MODE options (confirm which one your lead used):
  "hops"  -> Weight = 1.0 for every edge. Distance = number of interaction steps.
             Unweighted shortest paths, the classic network-proximity setup.
  "score" -> Weight = 1 - combined_score/1000  (score 900 -> 0.1, score 700 -> 0.3)
  "log"   -> Weight = -ln(combined_score/1000) (score 900 -> 0.105, score 700 -> 0.357)
"""
import pandas as pd

# ---------------- CONFIG ----------------
LINKS_FILE = "9606.protein.links.v12.5.min700.onlyAB.tsv.gz"
INFO_FILE = "9606.protein.info.v12.5.txt.gz"
OUTPUT_FILE = "STRING_PPI_filtered.txt"
MIN_SCORE = 700            # already applied on the website; kept as a safety net
WEIGHT_MODE = "hops"       # "hops", "score" or "log"
# ----------------------------------------

import numpy as np


def make_weight(scores: pd.Series) -> pd.Series:
    if WEIGHT_MODE == "hops":
        return pd.Series(1.0, index=scores.index)
    if WEIGHT_MODE == "score":
        return 1.0 - scores / 1000.0
    if WEIGHT_MODE == "log":
        return -np.log(scores / 1000.0)
    raise ValueError(f"Unknown WEIGHT_MODE: {WEIGHT_MODE}")


def main():
    links = pd.read_csv(LINKS_FILE, sep="\t")
    info = pd.read_csv(INFO_FILE, sep="\t")
    print(f"Links read: {len(links)} | proteins in info file: {len(info)}")

    links = links[links["combined_score"] >= MIN_SCORE].copy()

    id2name = dict(zip(info["#string_protein_id"], info["preferred_name"]))
    links["source"] = links["protein1"].map(id2name)
    links["target"] = links["protein2"].map(id2name)

    unmapped = links["source"].isna() | links["target"].isna()
    print(f"Edges with an unmapped protein ID (dropped): {int(unmapped.sum())}")
    links = links[~unmapped]

    # drop self-loops and any remaining duplicate A-B / B-A pairs
    links = links[links["source"] != links["target"]]
    key = links[["source", "target"]].apply(lambda r: tuple(sorted(r)), axis=1)
    links = links.loc[~key.duplicated()]

    links["Weight"] = make_weight(links["combined_score"]).astype(float)
    assert links["Weight"].notna().all(), "Weight contains blanks"
    assert pd.api.types.is_float_dtype(links["Weight"]), "Weight must be numeric"

    out = links[["source", "target", "Weight"]]
    out.to_csv(OUTPUT_FILE, sep="\t", index=False)

    genes = set(out["source"]) | set(out["target"])
    print(f"Saved {OUTPUT_FILE}: {len(out)} edges, {len(genes)} genes, "
          f"WEIGHT_MODE={WEIGHT_MODE}")
    print(out.head())


if __name__ == "__main__":
    main()
