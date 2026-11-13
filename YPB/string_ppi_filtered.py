import pandas as pd


def create_string_ppi_filtered(
    protein_file,
    string_file,
    output_file="STRING_PPI_filtered.txt"
):
    """
    Create STRING_PPI_filtered.txt from STRING v12.5 files.

    Output columns:
        source
        target
        Weight

    Weight = combined_score / 1000
    """

    print("Loading protein info...")
    protein = pd.read_csv(
        protein_file,
        sep="\t"
    )[["#string_protein_id", "preferred_name"]]

    print(f"Protein mappings loaded: {len(protein):,}")

    print("Loading STRING interactions...")
    string = pd.read_csv(
        string_file,
        sep="\t"
    )

    original_interactions = len(string)

    print(f"STRING interactions loaded: {original_interactions:,}")

    print("Converting STRING scores to 0-1 scale...")
    string["combined_score"] = (
        string["combined_score"] / 1000.0
    )

    print("Mapping protein1 -> source...")
    string = pd.merge(
        string,
        protein,
        how="left",
        left_on="protein1",
        right_on="#string_protein_id"
    )

    string = string.rename(
        columns={"preferred_name": "source"}
    ).drop(
        columns=["#string_protein_id"]
    )

    print("Mapping protein2 -> target...")
    string = pd.merge(
        string,
        protein,
        how="left",
        left_on="protein2",
        right_on="#string_protein_id"
    )

    string = string.rename(
        columns={"preferred_name": "target"}
    ).drop(
        columns=["#string_protein_id"]
    )

    # --------------------------------------------------
    # Validation checks
    # --------------------------------------------------

    print("\n========== VALIDATION ==========")

    missing_source = string["source"].isna().sum()
    missing_target = string["target"].isna().sum()

    print(f"Missing source mappings : {missing_source:,}")
    print(f"Missing target mappings : {missing_target:,}")

    if missing_source > 0 or missing_target > 0:
        print("WARNING: Some STRING IDs could not be mapped!")

    print("\nExample mappings:")

    cols_to_show = [
        "protein1",
        "source",
        "protein2",
        "target"
    ]

    print(
        string[cols_to_show]
        .head(10)
        .to_string(index=False)
    )

    # --------------------------------------------------
    # Create final file
    # --------------------------------------------------

    ppi = string[
        ["source", "target", "combined_score"]
    ].copy()

    ppi.rename(
        columns={"combined_score": "Weight"},
        inplace=True
    )

    # Duplicate check
    duplicate_pairs = ppi.duplicated(
        subset=["source", "target"]
    ).sum()

    print(f"\nDuplicate gene pairs: {duplicate_pairs:,}")

    # Weight checks
    print(
        f"Weight range: "
        f"{ppi['Weight'].min():.3f} "
        f"to "
        f"{ppi['Weight'].max():.3f}"
    )

    if ppi["Weight"].min() < 0.9:
        print(
            "WARNING: Found interactions with Weight < 0.9"
        )

    final_rows = len(ppi)

    print(f"Input interactions : {original_interactions:,}")
    print(f"Output interactions: {final_rows:,}")

    if original_interactions == final_rows:
        print("PASS: No interactions lost during processing.")
    else:
        print("WARNING: Interaction count changed!")

    print("\nFirst 10 rows of final file:")
    print(ppi.head(10).to_string(index=False))

    ppi.to_csv(
        output_file,
        sep="\t",
        index=False
    )

    print("\n========== COMPLETE ==========")
    print(f"Created: {output_file}")
    print(f"Interactions: {len(ppi):,}")
    print(
        f"Unique genes: "
        f"{len(set(ppi['source']).union(set(ppi['target']))):,}"
    )


if __name__ == "__main__":

    create_string_ppi_filtered(
        protein_file="9606.protein.info.v12.5.txt.gz",
        string_file="9606.protein.links.detailed.v12.5.min900.onlyAB.tsv.gz",
        output_file="STRING_PPI_filtered.txt"
    )






################################################# OUTPUT #######################################################
# Protein mappings loaded: 19,699
# STRING interactions: 125,344
# Missing source mappings : 0
# Missing target mappings : 0
# Duplicate gene pairs: 0
# Weight range: 0.900 to 0.999
# PASS: No interactions lost during processing.
# Unique genes: 14,301



