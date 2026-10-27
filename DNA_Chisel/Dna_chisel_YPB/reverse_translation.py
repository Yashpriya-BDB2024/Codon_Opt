import random
import math
from Bio.Data.CodonTable import standard_dna_table


def frequency_weighted_sampling(protein_seq, codon_usage_table, seed=42):
    """Generate a DNA sequence using host-specific codon frequencies.
    Each amino acid is replaced by a synonymous codon selected according to its relative usage frequency.

    Args:
        protein_seq (str): Protein sequence using standard one-letter amino acid codes.
        codon_usage_table (dict): Codon frequencies grouped by amino acid.
        seed (int): Random seed for reproducible results.

    Returns:
        str: Generated DNA coding sequence without a stop codon.

    Raises:
        TypeError: If the input types are invalid.
        ValueError: If the protein sequence or codon table is invalid.
        RuntimeError: If the generated DNA fails translation validation.
    """
    if not isinstance(protein_seq, str) or not protein_seq.strip():
        raise ValueError("Protein sequence cannot be empty.")

    protein_seq = "".join(protein_seq.split()).upper()

    valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
    invalid_amino_acids = set(protein_seq) - valid_amino_acids

    if invalid_amino_acids:
        raise ValueError(
            f"Invalid amino acid symbols: {sorted(invalid_amino_acids)}"
        )

    if not isinstance(codon_usage_table, dict):
        raise TypeError("Codon usage table must be a dictionary.")

    if isinstance(seed, bool) or not isinstance(seed, int):
        raise TypeError("Seed must be an integer.")

    rng = random.Random(seed)

    codon_to_aa = {
        codon: aa
        for codon, aa in standard_dna_table.forward_table.items()
    }

    dna_codons = []

    for position, aa in enumerate(protein_seq, start=1):
        if aa not in codon_usage_table:
            raise ValueError(
                f"Missing codon frequencies for {aa} at position {position}."
            )

        codon_freq = codon_usage_table[aa]

        if not isinstance(codon_freq, dict) or not codon_freq:
            raise ValueError(f"Invalid codon data for amino acid {aa}.")

        codons = []
        weights = []

        for codon, freq in codon_freq.items():
            if not isinstance(codon, str):
                raise ValueError(f"Invalid codon: {codon!r}")

            codon = codon.upper()

            if (
                len(codon) != 3
                or any(base not in "ACGT" for base in codon)
            ):
                raise ValueError(f"Invalid DNA codon: {codon!r}")

            if codon_to_aa.get(codon) != aa:
                raise ValueError(
                    f"Codon {codon} does not encode amino acid {aa}."
                )

            if (
                isinstance(freq, bool)
                or not isinstance(freq, (int, float))
                or not math.isfinite(freq)
                or freq < 0
            ):
                raise ValueError(
                    f"Invalid frequency for codon {codon}: {freq!r}"
                )

            if freq > 0:
                codons.append(codon)
                weights.append(float(freq))

        if not codons:
            raise ValueError(f"No usable codons for amino acid {aa}.")

        dna_codons.append(
            rng.choices(codons, weights=weights, k=1)[0]
        )

    dna_seq = "".join(dna_codons)

    translated_seq = "".join(
        codon_to_aa[dna_seq[i:i + 3]]
        for i in range(0, len(dna_seq), 3)
    )

    if translated_seq != protein_seq:
        raise RuntimeError("Generated DNA does not match the input protein.")

    return dna_seq
