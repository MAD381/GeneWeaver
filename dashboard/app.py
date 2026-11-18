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

    # GC Content Analytics for all uploaded sequences
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

            # Prepare downloadable CSV report
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
