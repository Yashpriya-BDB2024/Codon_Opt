import pandas as pd

links = pd.read_csv("9606.protein.links.v12.5.min700.onlyAB.tsv.gz", sep="\t")
info  = pd.read_csv("9606.protein.info.v12.5.txt.gz", sep="\t")

id2name = dict(zip(info["#string_protein_id"], info["preferred_name"]))
links["source"] = links["protein1"].map(id2name)
links["target"] = links["protein2"].map(id2name)
links = links.dropna(subset=["source", "target"])

links["Weight"] = ...   # still to confirm with your lead

links[["source", "target", "Weight"]].to_csv("STRING_PPI_filtered.txt", sep="\t", index=False)