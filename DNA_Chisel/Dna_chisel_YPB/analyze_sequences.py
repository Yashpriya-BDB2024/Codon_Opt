import argparse
from pathlib import Path
import pandas as pd
from Bio.Seq import Seq


def analyze_protein_sequence(sequence):
    """Validate a protein sequence and report its basic properties."""

    protein = (
        ""
        if pd.isna(sequence)
        else "".join(str(sequence).split()).upper()
    )

    # A single terminal '*' is accepted as stop-symbol notation.
    protein_for_check = (
        protein[:-1] if protein.endswith("*") else protein
    )

    standard_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
    invalid_residues = sorted(
        set(protein_for_check) - standard_amino_acids
    )

    valid = bool(protein_for_check) and not invalid_residues

    amino_acid_counts = {
        aa: protein_for_check.count(aa)
        for aa in sorted(standard_amino_acids)
        if protein_for_check.count(aa) > 0
    }

    return {
        "Protein length": len(protein_for_check),
        "Protein invalid residues": ", ".join(invalid_residues),
        "Protein valid": valid,
        "Protein starts with methionine": (
            protein_for_check.startswith("M") if valid else None
        ),
        "Protein amino acid counts": str(amino_acid_counts),
    }


def analyze_dna_sequence(dna_sequence, protein_sequence):
    """Validate DNA and compare its translation with the intended protein."""

    dna = (
        ""
        if pd.isna(dna_sequence)
        else "".join(str(dna_sequence).split()).upper()
    )
    protein = (
        ""
        if pd.isna(protein_sequence)
        else "".join(str(protein_sequence).split()).upper()
    )

    # Permit a terminal '*' in the protein as stop-symbol notation.
    protein_for_check = protein[:-1] if protein.endswith("*") else protein

    valid_bases = set("ACGT")
    invalid_bases = sorted(set(dna) - valid_bases)
    dna_valid = bool(dna) and not invalid_bases
    divisible_by_three = bool(dna) and len(dna) % 3 == 0

    base_counts = {base: dna.count(base) for base in "ACGT"}
    gc_percent = (
        round((base_counts["G"] + base_counts["C"]) / len(dna) * 100, 2)
        if dna_valid
        else None
    )

    longest_run = 0
    current_run = 0
    previous_base = None

    for base in dna:
        if base not in valid_bases:
            current_run = 0
            previous_base = None
            continue

        if base == previous_base:
            current_run += 1
        else:
            current_run = 1
            previous_base = base

        longest_run = max(longest_run, current_run)

    protein_valid = (
        bool(protein_for_check)
        and not (set(protein_for_check) - set("ACDEFGHIKLMNPQRSTVWY"))
    )

    translation_matches = None

    if dna_valid and divisible_by_three and protein_valid:
        translated = str(Seq(dna).translate(table=1))

        # A final DNA stop codon is optional in this comparison.
        if translated.endswith("*"):
            translated = translated[:-1]

        translation_matches = translated == protein_for_check

    return {
        "DNA length": len(dna),
        "DNA length divisible by 3": divisible_by_three,
        "DNA invalid bases": ", ".join(invalid_bases),
        "DNA valid": dna_valid,
        "A count": base_counts["A"],
        "C count": base_counts["C"],
        "G count": base_counts["G"],
        "T count": base_counts["T"],
        "DNA GC percentage": gc_percent,
        "Longest homopolymer run": longest_run,
        "DNA starts with ATG": dna.startswith("ATG") if dna_valid else None,
        "DNA ends with stop codon": (
            dna[-3:] in {"TAA", "TAG", "TGA"}
            if dna_valid and len(dna) >= 3
            else None
        ),
        "Translation matches protein": translation_matches,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Analyse protein and generated DNA sequences in Excel."
    )
    parser.add_argument("--input", required=True, help="Input Excel file")
    parser.add_argument("--output", required=True, help="Output Excel report")
    args = parser.parse_args()

    input_path = Path(args.input).resolve()
    output_path = Path(args.output).resolve()

    if not input_path.is_file():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    if input_path == output_path:
        raise ValueError("Input and output files must be different.")

    protein_column = "Amino acid sequence"
    dna_column = "Generated DNA Sequence"
    host_column = "Codon optimization- Host"

    data = pd.read_excel(input_path)

    required_columns = {protein_column, dna_column}
    missing_columns = required_columns - set(data.columns)

    if missing_columns:
        raise ValueError(
            f"Missing required columns: {sorted(missing_columns)}. "
            f"Available columns: {data.columns.tolist()}"
        )

    if host_column in data.columns:
        data = data[
            data[host_column].astype("string").str.strip().str.casefold()
            == "homosapiens"
        ].copy()

    if data.empty:
        raise ValueError("No records available for analysis.")

    protein_results = [
        analyze_protein_sequence(sequence)
        for sequence in data[protein_column]
    ]

    dna_results = [
        analyze_dna_sequence(dna, protein)
        for dna, protein in zip(
            data[dna_column], data[protein_column]
        )
    ]

    report = data.copy()

    protein_report = pd.DataFrame(protein_results, index=data.index)
    dna_report = pd.DataFrame(dna_results, index=data.index)

    report = pd.concat([report, protein_report, dna_report], axis=1)

    report["Overall basic status"] = [
        (
            "Missing protein or DNA"
            if pd.isna(protein) or not str(protein).strip()
            or pd.isna(dna) or not str(dna).strip()
            else "Invalid protein"
            if not protein_result["Protein valid"]
            else "Invalid DNA"
            if not dna_result["DNA valid"]
            else "DNA length not divisible by 3"
            if not dna_result["DNA length divisible by 3"]
            else "Translation mismatch"
            if dna_result["Translation matches protein"] is False
            else "Checks passed"
            if dna_result["Translation matches protein"] is True
            else "Translation not verified"
        )
        for protein, dna, protein_result, dna_result in zip(
            data[protein_column],
            data[dna_column],
            protein_results,
            dna_results,
        )
    ]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    report.to_excel(output_path, index=False)

    print(f"Records analysed: {len(report)}")
    print("Overall status counts:")
    print(report["Overall basic status"].value_counts(dropna=False).to_string())
    print(f"Report saved to: {output_path}")


if __name__ == "__main__":
    main()
