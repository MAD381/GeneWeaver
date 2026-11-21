"""Generate consolidated JSON and CSV reports for GeneWeaver."""

import csv
import json
from datetime import datetime, timezone
from pathlib import Path

from data_pipeline.fasta_processor import (
    calculate_sequence_statistics,
    generate_reverse_complement,
)

from data_pipeline.quality_control import run_quality_control

from ai_engine.dna_intelligence import generate_dna_insights


def build_analysis_report(
    records,
    motif_results=None,
    comparison_result=None,
):
    """
    Build a consolidated report for uploaded DNA sequences.

    Includes:
    - DNA statistics
    - AI-assisted insights
    - Reverse complements
    - DNA quality control
    - Optional motif search results
    - Optional pairwise comparison results
    """

    if not records:
        raise ValueError("At least one DNA sequence is required.")

    # Run quality control on all sequences.
    qc_report = run_quality_control(records)

    sequence_reports = []

    for index, record in enumerate(records):
        sequence = str(record.seq).strip().upper()

        statistics = calculate_sequence_statistics(sequence)
        analysis = generate_dna_insights(sequence)
        reverse_complement = generate_reverse_complement(sequence)

        # Get QC results for this sequence.
        quality = qc_report["results"][index]

        sequence_reports.append(
            {
                "sequence_id": str(record.id),
                "description": str(record.description),
                "sequence_length": len(sequence),
                "statistics": statistics,
                "patterns": analysis["patterns"],
                "insights": analysis["insights"],
                "reverse_complement": str(reverse_complement),
                "quality_control": quality,
            }
        )

    report = {
        "report_title": "GeneWeaver Consolidated DNA Analysis Report",
        "report_version": "1.1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "project": "GeneWeaver",
        "sequence_count": len(sequence_reports),
        "sequences": sequence_reports,
        "quality_control": {
            "summary": qc_report["summary"],
            "results": qc_report["results"],
        },
        "motif_search_results": motif_results or [],
        "pairwise_comparison": comparison_result,
        "methodology_notes": [
            "DNA statistics are calculated using GeneWeaver's existing "
            "sequence statistics function.",
            "GC and AT percentages use total sequence length as the denominator.",
            "Ambiguous bases are counted separately.",
            "AI-assisted insights are explainable, rule-based observations.",
            "Pattern observations do not establish gene function or biological significance.",
            "Quality control checks sequence validity, length, GC content "
            "and ambiguous nucleotide percentage.",
            "FASTA sequence-level quality indicators are not Phred sequencing "
            "quality scores.",
            "Pairwise identity, when included, follows GeneWeaver's existing "
            "alignment calculation.",
        ],
    }

    return report


def export_json_report(report, output_file):
    """Save a consolidated report as formatted JSON."""

    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as file:
        json.dump(report, file, indent=4, ensure_ascii=False)

    return output_path


def export_csv_report(report, output_file):
    """
    Save sequence-level analysis and quality-control results as CSV.

    Nested motif, comparison, pattern and QC details are stored
    as JSON strings in their respective CSV columns.
    """

    if not isinstance(report, dict) or "sequences" not in report:
        raise ValueError("A valid GeneWeaver analysis report is required.")

    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "sequence_id",
        "description",
        "sequence_length",
        "A",
        "T",
        "G",
        "C",
        "GC_percentage",
        "AT_percentage",
        "ambiguous_bases",
        "patterns",
        "insights",
        "reverse_complement",
        "quality_status",
        "ambiguous_percentage",
        "invalid_bases",
        "quality_warnings",
        "quality_errors",
        "motif_search_results",
        "pairwise_comparison",
    ]

    with output_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as file:

        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()

        for sequence in report["sequences"]:
            statistics = sequence["statistics"]
            quality = sequence.get("quality_control", {})

            writer.writerow(
                {
                    "sequence_id": sequence["sequence_id"],
                    "description": sequence["description"],
                    "sequence_length": sequence["sequence_length"],
                    "A": statistics["A"],
                    "T": statistics["T"],
                    "G": statistics["G"],
                    "C": statistics["C"],
                    "GC_percentage": statistics["GC_percentage"],
                    "AT_percentage": statistics["AT_percentage"],
                    "ambiguous_bases": statistics["ambiguous_bases"],
                    "patterns": json.dumps(
                        sequence["patterns"],
                        ensure_ascii=False,
                    ),
                    "insights": json.dumps(
                        sequence["insights"],
                        ensure_ascii=False,
                    ),
                    "reverse_complement": sequence["reverse_complement"],
                    "quality_status": quality.get("status", ""),
                    "ambiguous_percentage": quality.get("ambiguous_percentage", ""),
                    "invalid_bases": json.dumps(
                        quality.get("invalid_bases", []),
                        ensure_ascii=False,
                    ),
                    "quality_warnings": json.dumps(
                        quality.get("warnings", []),
                        ensure_ascii=False,
                    ),
                    "quality_errors": json.dumps(
                        quality.get("errors", []),
                        ensure_ascii=False,
                    ),
                    "motif_search_results": json.dumps(
                        report.get("motif_search_results", []),
                        ensure_ascii=False,
                    ),
                    "pairwise_comparison": json.dumps(
                        report.get("pairwise_comparison"),
                        ensure_ascii=False,
                    ),
                }
            )

    return output_path
