import json
from pathlib import Path

# Supported DNA bases (standard IUPAC DNA alphabet)
VALID_DNA_BASES = set("ACGTRYSWKMBDHVN")


def validate_comparison_sequence(sequence, name="Sequence"):
    """Validate and normalize a DNA sequence."""
    if not isinstance(sequence, str):
        raise TypeError(f"{name} must be a string.")

    sequence = sequence.strip().upper()

    if not sequence:
        raise ValueError(f"{name} cannot be empty.")

    invalid_bases = sorted(set(sequence) - VALID_DNA_BASES)

    if invalid_bases:
        raise ValueError(
            f"{name} contains invalid DNA bases: {', '.join(invalid_bases)}"
        )

    return sequence


def align_sequences(
    sequence1,
    sequence2,
    match_score=2,
    mismatch_score=-1,
    gap_score=-2,
):
    """
    Perform global DNA sequence alignment using the
    Needleman-Wunsch algorithm.

    Tie-breaking order: diagonal, then up, then left.
    """
    sequence1 = validate_comparison_sequence(sequence1, "Sequence 1")
    sequence2 = validate_comparison_sequence(sequence2, "Sequence 2")

    rows = len(sequence1) + 1
    cols = len(sequence2) + 1

    scores = [[0] * cols for _ in range(rows)]

    # Initialize gap penalties
    for i in range(1, rows):
        scores[i][0] = i * gap_score

    for j in range(1, cols):
        scores[0][j] = j * gap_score

    # Fill the scoring matrix
    for i in range(1, rows):
        for j in range(1, cols):
            diagonal = scores[i - 1][j - 1] + (
                match_score if sequence1[i - 1] == sequence2[j - 1] else mismatch_score
            )

            up = scores[i - 1][j] + gap_score
            left = scores[i][j - 1] + gap_score

            scores[i][j] = max(diagonal, up, left)

    # Trace back to construct the alignment
    aligned1 = []
    aligned2 = []

    i = len(sequence1)
    j = len(sequence2)

    while i > 0 or j > 0:
        if i > 0 and j > 0:
            score = (
                match_score
                if (sequence1[i - 1] == sequence2[j - 1])
                else mismatch_score
            )

            if scores[i][j] == scores[i - 1][j - 1] + score:
                aligned1.append(sequence1[i - 1])
                aligned2.append(sequence2[j - 1])
                i -= 1
                j -= 1
                continue

        if i > 0 and scores[i][j] == scores[i - 1][j] + gap_score:
            aligned1.append(sequence1[i - 1])
            aligned2.append("-")
            i -= 1
        else:
            aligned1.append("-")
            aligned2.append(sequence2[j - 1])
            j -= 1

    aligned1.reverse()
    aligned2.reverse()

    return (
        "".join(aligned1),
        "".join(aligned2),
        scores[-1][-1],
    )


def compare_sequences(sequence1, sequence2):
    """Align two DNA sequences and calculate comparison statistics."""
    sequence1 = validate_comparison_sequence(sequence1, "Sequence 1")
    sequence2 = validate_comparison_sequence(sequence2, "Sequence 2")

    aligned1, aligned2, alignment_score = align_sequences(sequence1, sequence2)

    matches = 0
    mismatches = 0
    gaps = 0

    mismatch_positions = []
    gap_positions = []

    for position, (base1, base2) in enumerate(zip(aligned1, aligned2), start=1):
        if base1 == "-" or base2 == "-":
            gaps += 1
            gap_positions.append(position)
        elif base1 == base2:
            matches += 1
        else:
            mismatches += 1
            mismatch_positions.append(position)

    alignment_length = len(aligned1)

    identity_percentage = matches / alignment_length * 100 if alignment_length else 0.0

    return {
        "sequence_1": sequence1,
        "sequence_2": sequence2,
        "aligned_sequence_1": aligned1,
        "aligned_sequence_2": aligned2,
        "alignment_score": alignment_score,
        "alignment_length": alignment_length,
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
        "identity_percentage": round(identity_percentage, 2),
        "mismatch_positions": mismatch_positions,
        "gap_positions": gap_positions,
    }


def summarize_comparison(result):
    """Create a readable summary of the comparison."""
    return {
        "alignment_length": result["alignment_length"],
        "matches": result["matches"],
        "mismatches": result["mismatches"],
        "gaps": result["gaps"],
        "identity_percentage": result["identity_percentage"],
        "alignment_score": result["alignment_score"],
    }


def export_comparison_report(result, output_file):
    """Export the comparison result as a JSON report."""
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as file:
        json.dump(result, file, indent=4)

    return str(output_path)
