"""
GeneWeaver DNA Sequence Quality Control Module

Checks DNA sequence validity, GC content, ambiguous bases,
sequence length and overall quality status.
"""

from data_pipeline.fasta_processor import (
    VALID_DNA_BASES,
    find_invalid_bases,
    calculate_sequence_statistics,
)


def analyze_sequence_quality(
    sequence,
    min_length=1,
    max_length=None,
    min_gc=20.0,
    max_gc=80.0,
    max_ambiguous_percent=10.0,
):
    """
    Perform quality checks on one DNA sequence.

    Returns a dictionary containing quality metrics,
    warnings, errors and final quality status.
    """

    sequence = str(sequence).strip().upper()

    warnings = []
    errors = []

    # 1. Empty sequence check
    if not sequence:
        return {
            "status": "FAIL",
            "length": 0,
            "GC_percentage": 0.0,
            "ambiguous_bases": 0,
            "ambiguous_percentage": 0.0,
            "invalid_bases": [],
            "warnings": [],
            "errors": ["Sequence is empty."],
        }

    length = len(sequence)

    # 2. Invalid nucleotide check
    invalid_bases = find_invalid_bases(sequence)

    if invalid_bases:
        errors.append(f"Invalid DNA bases detected: {invalid_bases}")

    # 3. Sequence length check
    if length < min_length:
        warnings.append(f"Sequence is shorter than minimum length ({min_length}).")

    if max_length is not None and length > max_length:
        warnings.append(f"Sequence exceeds maximum length ({max_length}).")

    # 4. Calculate statistics only for valid DNA
    if invalid_bases:
        gc_percentage = None
        ambiguous_bases = sum(1 for base in sequence if base not in "ATGC")
        ambiguous_percentage = round((ambiguous_bases / length) * 100, 2)

    else:
        statistics = calculate_sequence_statistics(sequence)

        gc_percentage = statistics["GC_percentage"]
        ambiguous_bases = statistics["ambiguous_bases"]

        ambiguous_percentage = round((ambiguous_bases / length) * 100, 2)

        # 5. GC content check
        if gc_percentage < min_gc:
            warnings.append(f"GC content is below {min_gc}%.")

        elif gc_percentage > max_gc:
            warnings.append(f"GC content is above {max_gc}%.")

        # 6. Ambiguous nucleotide check
        if ambiguous_percentage > max_ambiguous_percent:
            warnings.append(
                "Ambiguous nucleotide percentage exceeds " f"{max_ambiguous_percent}%."
            )

    # 7. Determine final status
    if errors:
        status = "FAIL"
    elif warnings:
        status = "WARNING"
    else:
        status = "PASS"

    return {
        "status": status,
        "length": length,
        "GC_percentage": gc_percentage,
        "ambiguous_bases": ambiguous_bases,
        "ambiguous_percentage": ambiguous_percentage,
        "invalid_bases": invalid_bases,
        "warnings": warnings,
        "errors": errors,
    }


def run_quality_control(
    records,
    min_length=1,
    max_length=None,
    min_gc=20.0,
    max_gc=80.0,
    max_ambiguous_percent=10.0,
):
    """
    Run quality control on multiple BioPython FASTA records.
    """

    if not records:
        raise ValueError("No sequences provided for quality control.")

    results = []

    for record in records:
        sequence = str(record.seq)

        quality = analyze_sequence_quality(
            sequence,
            min_length=min_length,
            max_length=max_length,
            min_gc=min_gc,
            max_gc=max_gc,
            max_ambiguous_percent=max_ambiguous_percent,
        )

        results.append(
            {
                "sequence_id": record.id,
                **quality,
            }
        )

    # 8. Overall quality summary
    summary = {
        "total_sequences": len(results),
        "passed": sum(1 for item in results if item["status"] == "PASS"),
        "warnings": sum(1 for item in results if item["status"] == "WARNING"),
        "failed": sum(1 for item in results if item["status"] == "FAIL"),
    }

    return {
        "project": "GeneWeaver",
        "report_type": "DNA Sequence Quality Control",
        "summary": summary,
        "results": results,
    }
