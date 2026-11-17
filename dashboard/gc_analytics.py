import pandas as pd
import plotly.express as px
import streamlit as st


def calculate_gc_analytics(records):
    """
    Calculate GC content and nucleotide statistics
    for every FASTA sequence.
    """

    results = []

    for record in records:
        sequence = str(record.seq).upper()

        valid_bases = [base for base in sequence if base in "ATGC"]

        length = len(sequence)
        valid_length = len(valid_bases)

        gc_count = valid_bases.count("G") + valid_bases.count("C")
        at_count = valid_bases.count("A") + valid_bases.count("T")

        gc_percentage = (gc_count / valid_length) * 100 if valid_length else 0

        at_percentage = (at_count / valid_length) * 100 if valid_length else 0

        results.append(
            {
                "Sequence ID": record.id,
                "Sequence Length": length,
                "Valid Bases": valid_length,
                "GC Count": gc_count,
                "AT Count": at_count,
                "GC Percentage": round(gc_percentage, 2),
                "AT Percentage": round(at_percentage, 2),
                "Ambiguous Bases": length - valid_length,
            }
        )

    return pd.DataFrame(results)


def display_gc_analytics(records):

    st.header("📊 Interactive GC Content Analytics")

    st.write(
        "Explore GC composition, compare DNA sequences "
        "and identify differences in nucleotide distribution."
    )

    if not records:
        st.warning("No DNA sequences available.")
        return

    df = calculate_gc_analytics(records)

    st.subheader("Sequence Summary")

    col1, col2, col3 = st.columns(3)

    col1.metric("Total Sequences", len(df))
    col2.metric("Average GC Content", f"{df['GC Percentage'].mean():.2f}%")
    col3.metric("Highest GC Content", f"{df['GC Percentage'].max():.2f}%")

    st.subheader("GC Content Comparison")

    figure = px.bar(
        df,
        x="Sequence ID",
        y="GC Percentage",
        color="GC Percentage",
        hover_data=[
            "Sequence Length",
            "Valid Bases",
            "GC Count",
            "Ambiguous Bases",
        ],
        title="GC Content Across DNA Sequences",
        labels={
            "GC Percentage": "GC Content (%)",
            "Sequence ID": "DNA Sequence",
        },
        color_continuous_scale="Viridis",
    )

    figure.update_layout(
        yaxis_range=[0, 100],
        xaxis_tickangle=-30,
    )

    st.plotly_chart(figure, use_container_width=True)

    st.subheader("Nucleotide Composition")

    composition = df[["Sequence ID", "GC Count", "AT Count"]].melt(
        id_vars="Sequence ID",
        var_name="Nucleotide Group",
        value_name="Count",
    )

    composition_chart = px.bar(
        composition,
        x="Sequence ID",
        y="Count",
        color="Nucleotide Group",
        barmode="group",
        title="AT and GC Base Distribution",
        color_discrete_sequence=["#636EFA", "#EF553B"],
    )

    st.plotly_chart(
        composition_chart,
        use_container_width=True,
    )

    st.subheader("Detailed GC Analysis")

    st.dataframe(
        df,
        use_container_width=True,
        hide_index=True,
    )

    csv_data = df.to_csv(index=False).encode("utf-8")

    st.download_button(
        label="Download GC Analytics CSV",
        data=csv_data,
        file_name="geneweaver_gc_analytics.csv",
        mime="text/csv",
    )

    highest = df.loc[df["GC Percentage"].idxmax()]
    lowest = df.loc[df["GC Percentage"].idxmin()]

    st.subheader("Automated Observations")

    st.write(
        f"• Highest GC content: {highest['Sequence ID']} "
        f"({highest['GC Percentage']}%)."
    )

    st.write(
        f"• Lowest GC content: {lowest['Sequence ID']} "
        f"({lowest['GC Percentage']}%)."
    )

    st.caption(
        "GC percentages are calculated using valid A, T, G and C "
        "bases only. Ambiguous bases are excluded from the "
        "percentage denominator."
    )
