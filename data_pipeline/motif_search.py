"""
GeneWeaver Custom DNA Motif Search Engine

Features:
- Search custom DNA motifs in multiple FASTA sequences.
- Detect overlapping motif occurrences.
- Report positions using 1-based indexing.
- Count motif occurrences.
- Summarize search results.
- Export results to CSV.
"""

import csv

VALID_DNA_BASES = set("ACGT")


def validate_motif(motif):
    """Validate and normalize a DNA motif."""

    if not isinstance(motif, str):
        raise TypeError("Motif must be a string.")

    motif = motif.strip().upper()

    if not motif:
        raise ValueError("Motif cannot be empty.")

    invalid_bases = set(motif) - VALID_DNA_BASES

    if invalid_bases:
        raise ValueError(f"Invalid DNA characters: {sorted(invalid_bases)}")

    return motif


def find_motif_positions(sequence, motif):
    """Find all motif positions, including overlapping matches."""

    motif = validate_motif(motif)
    sequence = str(sequence).strip().upper()

    positions = []
    start = 0

    while True:
        index = sequence.find(motif, start)

        if index == -1:
            break

        positions.append(index + 1)
        start = index + 1

    return positions


def search_motif(records, motif):
    """Search a custom motif across multiple Biopython records."""

    motif = validate_motif(motif)
    results = []

    for record in records:
        sequence = str(record.seq).upper()
        positions = find_motif_positions(sequence, motif)

        results.append(
            {
                "sequence_id": record.id,
                "sequence_length": len(sequence),
                "motif": motif,
                "occurrence_count": len(positions),
                "positions": positions,
            }
        )

    return results


def summarize_motif_search(results):
    """Generate a summary of motif search results."""

    total_occurrences = sum(result["occurrence_count"] for result in results)

    sequences_with_matches = sum(result["occurrence_count"] > 0 for result in results)

    return {
        "total_sequences": len(results),
        "sequences_with_matches": sequences_with_matches,
        "sequences_without_matches": (len(results) - sequences_with_matches),
        "total_occurrences": total_occurrences,
    }


def export_motif_results(results, output_file):
    """Export motif search results to a CSV file."""

    with open(output_file, "w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=[
                "sequence_id",
                "sequence_length",
                "motif",
                "occurrence_count",
                "positions",
            ],
        )

        writer.writeheader()

        for result in results:
            row = result.copy()
            row["positions"] = ", ".join(
                str(position) for position in result["positions"]
            )
            writer.writerow(row)

    return output_file
