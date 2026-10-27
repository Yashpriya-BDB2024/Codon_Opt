from pathlib import Path
import pandas as pd
from python_codon_tables import get_codons_table
from reverse_translation import frequency_weighted_sampling


def main():
    """Generate DNA sequences for Homo sapiens proteins in the Excel file."""
    base_dir = Path(__file__).resolve().parent
    input_file = base_dir / "codon_optimization_trials.xlsx"
    output_file = base_dir / "human_codon_sequences.xlsx"

    host_col = "Codon optimization- Host"
    protein_col = "Amino acid sequence"
    output_col = "Generated DNA Sequence"

    # Read the input workbook.
    df = pd.read_excel(input_file)

    for col in (host_col, protein_col):
        if col not in df.columns:
            raise ValueError(f"Missing required Excel column: {col}")

    # Keep only Homo sapiens records.
    hosts = (
        df[host_col]
        .fillna("")
        .astype(str)
        .str.lower()
        .str.replace(r"\s+", "", regex=True)
    )

    human_df = df.loc[hosts == "homosapiens"].copy()

    if human_df.empty:
        raise ValueError("No Homo sapiens records found.")

    # Load human codon usage data.
    human_table = get_codons_table("h_sapiens")

    # Generate a DNA sequence for each protein.
    dna_sequences = []

    for row_idx, protein_seq in human_df[protein_col].items():
        if pd.isna(protein_seq):
            raise ValueError(
                f"Missing protein sequence at Excel row {row_idx + 2}."
            )

        try:
            dna_seq = frequency_weighted_sampling(
                protein_seq=str(protein_seq),
                codon_usage_table=human_table,
                seed=42,
            )
        except (TypeError, ValueError, RuntimeError) as error:
            raise ValueError(
                f"Could not process Excel row {row_idx + 2}: {error}"
            ) from error

        dna_sequences.append(dna_seq)

        if "Name" in human_df.columns:
            print(f"Processed {human_df.at[row_idx, 'Name']}")

    human_df[output_col] = dna_sequences

    # Save the human records without modifying the original workbook.
    human_df.to_excel(output_file, index=False)

    print(f"Processed {len(human_df)} human protein records.")
    print(f"Results saved to: {output_file}")


if __name__ == "__main__":
    main()
