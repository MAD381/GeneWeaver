import csv
import hashlib
import io
import json
import os
import tempfile

import streamlit as st

from data_pipeline.fasta_processor import (
    parse_fasta,
    calculate_sequence_statistics,
    generate_reverse_complement,
)

from data_pipeline.motif_search import (
    search_motif,
    summarize_motif_search,
)

from data_pipeline.pairwise_comparison import compare_sequences

from ai_engine.dna_intelligence import generate_dna_insights
from dashboard.gc_analytics import display_gc_analytics
from data_pipeline.analysis_report import (
    build_analysis_report,
    export_json_report,
    export_csv_report,
)
from data_pipeline.quality_control import run_quality_control

st.set_page_config(
    page_title="GeneWeaver DNA Intelligence",
    page_icon="🧬",
    layout="wide",
)

st.title("🧬 GeneWeaver DNA Intelligence Dashboard")

st.write(
    "An AI-assisted bioinformatics interface for FASTA processing, "
    "DNA statistics, pattern detection, and sequence exploration."
)

st.info(
    "This dashboard uses GeneWeaver's existing DNA analysis engine. "
    "Pattern detections are research observations, not clinical or "
    "experimental conclusions."
)

uploaded_file = st.file_uploader(
    "Upload a FASTA file",
    type=["fasta", "fa", "fna", "txt"],
)

if uploaded_file is None:
    st.warning("Upload a FASTA file to begin DNA analysis.")
    st.stop()

# Keep optional analysis results associated with the uploaded FASTA file.
upload_bytes = uploaded_file.getvalue()
upload_fingerprint = hashlib.sha256(upload_bytes).hexdigest()

if st.session_state.get("geneweaver_upload_fingerprint") != upload_fingerprint:
    st.session_state["geneweaver_upload_fingerprint"] = upload_fingerprint
    st.session_state["geneweaver_motif_results"] = None
    st.session_state["geneweaver_pairwise_comparison"] = None

temp_path = None

try:
    with tempfile.NamedTemporaryFile(
        delete=False,
        suffix=".fasta",
    ) as temp_file:
        temp_file.write(upload_bytes)
        temp_path = temp_file.name

    records = parse_fasta(temp_path)

    if not records:
        st.error("No DNA sequences were found in this FASTA file.")
        st.stop()

    st.success(f"Successfully loaded {len(records)} DNA sequence(s).")

    selected_index = st.selectbox(
        "Select a DNA sequence",
        options=list(range(len(records))),
        format_func=lambda i: (f"{records[i].id} — {len(records[i].seq)} bases"),
    )

    selected_record = records[selected_index]
    sequence = str(selected_record.seq).upper()

    # GC Content Analytics
    st.divider()
    display_gc_analytics(records)
    st.divider()

    # Custom DNA Motif Search
    st.header("🔬 Custom DNA Motif Search")

    st.write(
        "Enter a DNA motif to search across all sequences "
        "in your uploaded FASTA file."
    )

    with st.form("motif_search_form"):
        motif_input = st.text_input(
            "Enter DNA motif (example: ATG)",
            placeholder="Enter A, T, G, C sequence",
        )

        search_submitted = st.form_submit_button("Search DNA Motif")

    if search_submitted:
        st.session_state["geneweaver_motif_results"] = None

        try:
            motif_results = search_motif(records, motif_input)
            st.session_state["geneweaver_motif_results"] = motif_results
            summary = summarize_motif_search(motif_results)

            st.success("Motif search completed successfully.")

            metric_columns = st.columns(3)

            metric_columns[0].metric(
                "Total Sequences",
                summary["total_sequences"],
            )

            metric_columns[1].metric(
                "Sequences with Matches",
                summary["sequences_with_matches"],
            )

            metric_columns[2].metric(
                "Total Occurrences",
                summary["total_occurrences"],
            )

            st.subheader("Motif Search Results")

            display_results = []

            for result in motif_results:
                display_results.append(
                    {
                        "Sequence ID": result["sequence_id"],
                        "Sequence Length": result["sequence_length"],
                        "Motif": result["motif"],
                        "Occurrence Count": result["occurrence_count"],
                        "Positions (1-based)": ", ".join(
                            str(position) for position in result["positions"]
                        ),
                    }
                )

            st.dataframe(
                display_results,
                use_container_width=True,
            )

            csv_buffer = io.StringIO()

            writer = csv.DictWriter(
                csv_buffer,
                fieldnames=[
                    "sequence_id",
                    "sequence_length",
                    "motif",
                    "occurrence_count",
                    "positions",
                ],
            )

            writer.writeheader()

            for result in motif_results:
                writer.writerow(
                    {
                        "sequence_id": result["sequence_id"],
                        "sequence_length": result["sequence_length"],
                        "motif": result["motif"],
                        "occurrence_count": result["occurrence_count"],
                        "positions": ", ".join(
                            str(position) for position in result["positions"]
                        ),
                    }
                )

            st.download_button(
                label="Download Motif Search Results (CSV)",
                data=csv_buffer.getvalue(),
                file_name="geneweaver_motif_results.csv",
                mime="text/csv",
            )

            st.caption(
                "Positions use 1-based indexing. Overlapping motif "
                "occurrences are included. Only A, T, G and C motifs "
                "are accepted by this search engine."
            )

        except (ValueError, TypeError) as error:
            st.error(f"Motif search error: {error}")

    st.divider()

    # Pairwise DNA Sequence Comparison
    st.header("🧬 Pairwise DNA Sequence Comparison")

    st.write(
        "Select two sequences from your FASTA file to align them "
        "and examine their matching and differing bases."
    )

    if len(records) < 2:
        st.warning(
            "Upload a FASTA file containing at least two sequences "
            "to use pairwise comparison."
        )

    else:
        sequence_options = list(range(len(records)))

        with st.form("pairwise_comparison_form"):
            comparison_col1, comparison_col2 = st.columns(2)

            with comparison_col1:
                first_index = st.selectbox(
                    "First DNA sequence",
                    options=sequence_options,
                    format_func=lambda i: (
                        f"{records[i].id} — {len(records[i].seq)} bases"
                    ),
                    key="comparison_sequence_1",
                )

            with comparison_col2:
                second_index = st.selectbox(
                    "Second DNA sequence",
                    options=sequence_options,
                    index=1,
                    format_func=lambda i: (
                        f"{records[i].id} — {len(records[i].seq)} bases"
                    ),
                    key="comparison_sequence_2",
                )

            comparison_submitted = st.form_submit_button("Compare DNA Sequences")

        if comparison_submitted:
            st.session_state["geneweaver_pairwise_comparison"] = None

            try:
                first_record = records[first_index]
                second_record = records[second_index]

                comparison = compare_sequences(
                    str(first_record.seq),
                    str(second_record.seq),
                )

                st.success("DNA sequence comparison completed.")

                st.subheader("Comparison Summary")

                metric_columns = st.columns(4)

                metric_columns[0].metric(
                    "Sequence Identity",
                    f"{comparison['identity_percentage']}%",
                )

                metric_columns[1].metric(
                    "Matching Bases",
                    comparison["matches"],
                )

                metric_columns[2].metric(
                    "Mismatches",
                    comparison["mismatches"],
                )

                metric_columns[3].metric(
                    "Gaps",
                    comparison["gaps"],
                )

                st.subheader("Global Sequence Alignment")

                st.text(f"Sequence 1 ({first_record.id}):")
                st.code(comparison["aligned_sequence_1"])

                st.text(f"Sequence 2 ({second_record.id}):")
                st.code(comparison["aligned_sequence_2"])

                st.write(f"**Alignment score:** {comparison['alignment_score']}")

                st.write(
                    f"**Alignment length:** "
                    f"{comparison['alignment_length']} columns"
                )

                st.write(
                    "**Mismatch positions (1-based alignment columns):** "
                    + (
                        ", ".join(
                            str(position)
                            for position in comparison["mismatch_positions"]
                        )
                        if comparison["mismatch_positions"]
                        else "None"
                    )
                )

                st.write(
                    "**Gap positions (1-based alignment columns):** "
                    + (
                        ", ".join(
                            str(position) for position in comparison["gap_positions"]
                        )
                        if comparison["gap_positions"]
                        else "None"
                    )
                )

                st.subheader("Detailed Comparison")

                comparison_table = [
                    {
                        "Measurement": "Matches",
                        "Value": comparison["matches"],
                    },
                    {
                        "Measurement": "Mismatches",
                        "Value": comparison["mismatches"],
                    },
                    {
                        "Measurement": "Gaps",
                        "Value": comparison["gaps"],
                    },
                    {
                        "Measurement": "Identity Percentage",
                        "Value": (f"{comparison['identity_percentage']}%"),
                    },
                    {
                        "Measurement": "Alignment Score",
                        "Value": comparison["alignment_score"],
                    },
                ]

                st.dataframe(
                    comparison_table,
                    use_container_width=True,
                    hide_index=True,
                )

                comparison_report = {
                    "sequence_1_id": first_record.id,
                    "sequence_2_id": second_record.id,
                    **comparison,
                }

                st.session_state["geneweaver_pairwise_comparison"] = comparison_report

                st.download_button(
                    label="Download Pairwise Comparison Report (JSON)",
                    data=json.dumps(
                        comparison_report,
                        indent=4,
                    ),
                    file_name="geneweaver_pairwise_comparison.json",
                    mime="application/json",
                )

                st.caption(
                    "Global alignment uses Needleman-Wunsch scoring. "
                    "Identity is calculated as matching bases divided "
                    "by total alignment columns, including gaps. "
                    "Mismatch and gap positions refer to alignment "
                    "columns, not original sequence coordinates. "
                    "Ambiguous IUPAC bases are accepted, but are counted "
                    "as matches only when their symbols are identical."
                )

            except (ValueError, TypeError) as error:
                st.error(f"Sequence comparison error: {error}")

    st.divider()

    # Selected Sequence Overview
    st.header("Sequence Overview")

    col1, col2, col3 = st.columns(3)

    col1.metric("Sequence ID", selected_record.id)
    col2.metric("Sequence Length", f"{len(sequence):,} bases")
    col3.metric("Total Sequences", len(records))

    # DNA Statistics
    st.subheader("DNA Statistics")

    statistics = calculate_sequence_statistics(sequence)

    stat_columns = st.columns(4)

    stat_columns[0].metric(
        "GC Content",
        f"{statistics['GC_percentage']}%",
    )

    stat_columns[1].metric(
        "AT Content",
        f"{statistics['AT_percentage']}%",
    )

    stat_columns[2].metric("A Bases", statistics["A"])
    stat_columns[3].metric("T Bases", statistics["T"])

    st.write("Complete statistics")
    st.json(statistics)

    # DNA Pattern Detection
    st.subheader("DNA Pattern Detection")

    analysis = generate_dna_insights(sequence)
    patterns = analysis["patterns"]

    pattern_columns = st.columns(4)

    for column, (pattern_name, count) in zip(
        pattern_columns,
        patterns.items(),
    ):
        column.metric(
            pattern_name.replace("_", " "),
            count,
        )

    # AI-Assisted Insights
    st.subheader("AI-Assisted DNA Insights")

    for insight in analysis["insights"]:
        st.write(f"• {insight}")

    st.caption(
        "These insights are generated using explainable, "
        "rule-based analysis. They do not establish gene function "
        "or biological significance."
    )

    # Reverse Complement
    st.subheader("Reverse Complement")

    reverse_complement = generate_reverse_complement(sequence)

    st.text_area(
        "Reverse-complement sequence",
        value=str(reverse_complement),
        height=150,
    )

    # Original Sequence
    st.subheader("Original DNA Sequence")

    st.text_area(
        "Selected sequence",
        value=sequence,
        height=150,
    )

    # JSON Report
    report = {
        "sequence_id": selected_record.id,
        "sequence_length": len(sequence),
        "statistics": statistics,
        "patterns": patterns,
        "insights": analysis["insights"],
        "reverse_complement": str(reverse_complement),
    }

    st.download_button(
        label="Download Selected Sequence Report (JSON)",
        data=json.dumps(report, indent=4),
        file_name=f"{selected_record.id}_geneweaver_report.json",
        mime="application/json",
    )
    # DNA Quality Control Dashboard
    st.divider()
    st.header("🧬 DNA Quality Control")

    st.write(
        "Check sequence length, GC content, ambiguous bases, and invalid bases "
        "for every sequence in the uploaded FASTA file."
    )

    qc_report = run_quality_control(records)
    qc_summary = qc_report["summary"]
    qc_results = qc_report["results"]

    # Quality Control Summary
    st.subheader("Quality Control Summary")

    qc_col1, qc_col2, qc_col3, qc_col4 = st.columns(4)

    qc_col1.metric(
        "Total Sequences",
        qc_summary.get("total_sequences", len(qc_results)),
    )

    qc_col2.metric(
        "Passed",
        qc_summary.get("passed", 0),
    )

    qc_col3.metric(
        "Warnings",
        qc_summary.get("warnings", 0),
    )

    qc_col4.metric(
        "Failed",
        qc_summary.get("failed", 0),
    )

    # Individual Sequence Quality Results
    st.subheader("Individual Sequence Quality Results")

    qc_table = []

    for index, result in enumerate(qc_results):
        qc_table.append(
            {
                "Sequence ID": result.get(
                    "sequence_id",
                    (
                        records[index].id
                        if index < len(records)
                        else f"Sequence {index + 1}"
                    ),
                ),
                "Status": result.get("status", "Unknown"),
                "Length": result.get("length", 0),
                "GC Percentage": result.get("GC_percentage"),
                "Ambiguous Percentage": result.get("ambiguous_percentage"),
                "Invalid Bases": str(result.get("invalid_bases", [])),
                "Warnings": "; ".join(str(item) for item in result.get("warnings", [])),
                "Errors": "; ".join(str(item) for item in result.get("errors", [])),
            }
        )

    st.dataframe(
        qc_table,
        use_container_width=True,
        hide_index=True,
    )

    # Detailed Quality Control Information
    st.subheader("Detailed Quality Control Results")

    for index, result in enumerate(qc_results):
        sequence_id = result.get(
            "sequence_id",
            records[index].id if index < len(records) else f"Sequence {index + 1}",
        )

        status = str(result.get("status", "Unknown")).upper()

        with st.expander(f"{sequence_id} — {status}"):
            st.write(f"**Sequence length:** {result.get('length', 0)}")

            st.write(f"**GC percentage:** {result.get('GC_percentage')}")

            st.write(f"**Ambiguous bases:** " f"{result.get('ambiguous_bases', 0)}")

            st.write(
                f"**Ambiguous percentage:** " f"{result.get('ambiguous_percentage')}"
            )

            st.write(f"**Invalid bases:** " f"{result.get('invalid_bases', [])}")

            warnings = result.get("warnings", [])
            errors = result.get("errors", [])

            if errors:
                st.error("\n".join(str(item) for item in errors))

            if warnings:
                st.warning("\n".join(str(item) for item in warnings))

            if not warnings and not errors:
                st.success("No quality warnings or errors detected.")

    st.caption(
        "These are computational sequence-quality checks, not clinical or "
        "experimental validation. FASTA files do not necessarily contain "
        "Phred quality scores."
    )

    # Consolidated Analysis Report
    st.divider()
    st.header("📄 Consolidated Analysis Report")

    st.write(
        "Download a combined report containing sequence statistics, "
        "AI-assisted insights, reverse complements, and any motif-search "
        "or pairwise-comparison results completed for this FASTA file."
    )

    consolidated_report = build_analysis_report(
        records,
        motif_results=st.session_state.get("geneweaver_motif_results"),
        comparison_result=st.session_state.get("geneweaver_pairwise_comparison"),
    )

    with tempfile.TemporaryDirectory() as report_dir:
        json_output_path = os.path.join(
            report_dir,
            "geneweaver_consolidated_report.json",
        )

        csv_output_path = os.path.join(
            report_dir,
            "geneweaver_consolidated_report.csv",
        )

        export_json_report(
            consolidated_report,
            json_output_path,
        )

        export_csv_report(
            consolidated_report,
            csv_output_path,
        )

        with open(json_output_path, "rb") as json_file:
            consolidated_json_data = json_file.read()

        with open(csv_output_path, "rb") as csv_file:
            consolidated_csv_data = csv_file.read()

    report_col1, report_col2 = st.columns(2)

    with report_col1:
        st.download_button(
            label="Download Complete Analysis Report (JSON)",
            data=consolidated_json_data,
            file_name="geneweaver_consolidated_report.json",
            mime="application/json",
            use_container_width=True,
        )

    with report_col2:
        st.download_button(
            label="Download Complete Analysis Report (CSV)",
            data=consolidated_csv_data,
            file_name="geneweaver_consolidated_report.csv",
            mime="text/csv",
            use_container_width=True,
        )


except Exception as error:
    st.error(f"Analysis failed: {error}")

finally:
    if temp_path and os.path.exists(temp_path):
        os.remove(temp_path)
