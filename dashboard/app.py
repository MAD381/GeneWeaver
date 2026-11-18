import csv
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

temp_path = None

try:
    with tempfile.NamedTemporaryFile(
        delete=False,
        suffix=".fasta",
    ) as temp_file:
        temp_file.write(uploaded_file.getvalue())
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
        try:
            motif_results = search_motif(records, motif_input)
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

                st.write(f"**Alignment score:** " f"{comparison['alignment_score']}")

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

except Exception as error:
    st.error(f"Analysis failed: {error}")

finally:
    if temp_path and os.path.exists(temp_path):
        os.remove(temp_path)
