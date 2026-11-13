import networkx as nx
import pandas as pd
import pickle
from typing import Dict
from tqdm import tqdm


def precompute_shortest_paths(
    tsv_file: str,
    output_file: str
) -> Dict:

    print(f"Loading network from {tsv_file}...")

    df = pd.read_csv(tsv_file, sep="\t")

    network = nx.Graph()

    for _, row in df.iterrows():

        network.add_edge(
            str(row["source"]),
            str(row["target"])
        )

    print(
        f"Network loaded: "
        f"{network.number_of_nodes()} nodes, "
        f"{network.number_of_edges()} edges"
    )

    print("Calculating shortest paths...")

    shortest_paths = {}

    for source in tqdm(
        network.nodes(),
        desc="Shortest paths progress"
    ):

        shortest_paths[source] = (
            nx.single_source_shortest_path_length(
                network,
                source
            )
        )

    with open(output_file, "wb") as f:
        pickle.dump(shortest_paths, f)

    print(
        f"Shortest paths computed and saved to "
        f"{output_file}"
    )

    return shortest_paths


if __name__ == "__main__":

    precompute_shortest_paths(
        tsv_file="STRING_PPI_filtered.txt",
        output_file="shortest_paths_cache.pkl"
    )
