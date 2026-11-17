import json
import os
import tempfile

import streamlit as st

from data_pipeline.fasta_processor import (
    parse_fasta,
    calculate_sequence_statistics,
    generate_reverse_complement,
)
from ai_engine.dna_intelligence import generate_dna_insights

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

    st.header("Sequence Overview")

    col1, col2, col3 = st.columns(3)

    col1.metric("Sequence ID", selected_record.id)
    col2.metric("Sequence Length", f"{len(sequence):,} bases")
    col3.metric("Total Sequences", len(records))

    st.subheader("DNA Statistics")

    statistics = calculate_sequence_statistics(sequence)

    stat_columns = st.columns(4)

    stat_columns[0].metric("GC Content", f"{statistics['GC_percentage']}%")
    stat_columns[1].metric("AT Content", f"{statistics['AT_percentage']}%")
    stat_columns[2].metric("A Bases", statistics["A"])
    stat_columns[3].metric("T Bases", statistics["T"])

    st.write("Complete statistics")
    st.json(statistics)

    st.subheader("DNA Pattern Detection")

    analysis = generate_dna_insights(sequence)
    patterns = analysis["patterns"]

    pattern_columns = st.columns(4)

    for column, (pattern_name, count) in zip(
        pattern_columns,
        patterns.items(),
    ):
        column.metric(pattern_name.replace("_", " "), count)

    st.subheader("AI-Assisted DNA Insights")

    for insight in analysis["insights"]:
        st.write(f"• {insight}")

    st.caption(
        "These insights are generated using explainable, "
        "rule-based analysis. They do not establish gene function "
        "or biological significance."
    )

    st.subheader("Reverse Complement")

    reverse_complement = generate_reverse_complement(sequence)

    st.text_area(
        "Reverse-complement sequence",
        value=str(reverse_complement),
        height=150,
    )

    st.subheader("Original DNA Sequence")

    st.text_area(
        "Selected sequence",
        value=sequence,
        height=150,
    )

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
