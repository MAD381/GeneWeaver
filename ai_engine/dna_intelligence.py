import csv
import json
from pathlib import Path

from data_pipeline.fasta_processor import (
    calculate_sequence_statistics,
    generate_reverse_complement,
)

# ==========================================================
# GENEWEAVER - HYBRID DNA INTELLIGENCE ENGINE
# Rule-based bioinformatics + optional Machine Learning
# ==========================================================

FEATURE_COLUMNS = [
    "length",
    "A",
    "T",
    "G",
    "C",
    "GC_percentage",
    "AT_percentage",
    "ambiguous_bases",
]


def count_motif(sequence, motif):
    """Count motif occurrences, including overlapping matches."""
    sequence = sequence.upper()
    motif = motif.upper()

    if not motif:
        return 0

    return sum(
        sequence.startswith(motif, index)
        for index in range(len(sequence) - len(motif) + 1)
    )


def detect_dna_patterns(sequence):
    """Detect selected DNA patterns. These are observations, not diagnoses."""
    sequence = sequence.upper()

    patterns = {
        "ATG_start_codon": count_motif(sequence, "ATG"),
        "TATAAA_motif": count_motif(sequence, "TATAAA"),
        "CpG_sites": count_motif(sequence, "CG"),
        "NGG_PAM_like_sites": 0,
    }

    # Count canonical N-G-G patterns (A/G/C/T followed by GG).
    for index in range(len(sequence) - 2):
        if sequence[index] in "ACGT" and sequence[index + 1 : index + 3] == "GG":
            patterns["NGG_PAM_like_sites"] += 1

    return patterns


def generate_dna_insights(sequence):
    """Generate explainable observations from DNA sequence statistics."""
    statistics = calculate_sequence_statistics(sequence)
    patterns = detect_dna_patterns(sequence)

    insights = []

    gc_content = statistics.get("GC_percentage", 0)
    ambiguous = statistics.get("ambiguous_bases", 0)
    length = statistics.get("length", len(sequence))

    if gc_content >= 60:
        insights.append("High GC content detected.")
    elif gc_content <= 40:
        insights.append("Low GC content detected.")
    else:
        insights.append("GC content is in the moderate range.")

    if ambiguous:
        insights.append(f"The sequence contains {ambiguous} ambiguous DNA bases.")
    else:
        insights.append("No ambiguous DNA bases detected.")

    if length < 100:
        insights.append("This is a relatively short DNA sequence.")
    else:
        insights.append("The sequence contains at least 100 bases.")

    if patterns["NGG_PAM_like_sites"]:
        insights.append(
            "NGG-like patterns detected; this does not establish "
            "CRISPR target suitability."
        )

    return {
        "statistics": statistics,
        "patterns": patterns,
        "insights": insights,
    }


def extract_ml_features(sequence):
    """Convert a DNA sequence into numeric features for machine learning."""
    stats = calculate_sequence_statistics(sequence)

    return {
        "length": stats.get("length", len(sequence)),
        "A": stats.get("A", 0),
        "T": stats.get("T", 0),
        "G": stats.get("G", 0),
        "C": stats.get("C", 0),
        "GC_percentage": stats.get("GC_percentage", 0),
        "AT_percentage": stats.get("AT_percentage", 0),
        "ambiguous_bases": stats.get("ambiguous_bases", 0),
    }


def run_ml_classification(records, training_csv):
    """
    Train a Random Forest using a user-provided labelled CSV.
    Predictions are only as reliable as the training data.
    """
    try:
        from sklearn.ensemble import RandomForestClassifier
    except ImportError as error:
        raise RuntimeError(
            "scikit-learn is required for ML. Install it using: "
            "python -m pip install scikit-learn"
        ) from error

    training_path = Path(training_csv)

    if not training_path.is_file():
        raise FileNotFoundError(f"Training CSV not found: {training_path}")

    features = []
    labels = []

    with training_path.open("r", newline="", encoding="utf-8-sig") as file:
        reader = csv.DictReader(file)

        if not reader.fieldnames:
            raise ValueError("Training CSV is empty or has no header.")

        required_columns = FEATURE_COLUMNS + ["label"]
        missing = [
            column for column in required_columns if column not in reader.fieldnames
        ]

        if missing:
            raise ValueError("Training CSV is missing columns: " + ", ".join(missing))

        for row_number, row in enumerate(reader, start=2):
            try:
                values = [float(row[column]) for column in FEATURE_COLUMNS]
            except (TypeError, ValueError):
                raise ValueError(f"Invalid numeric feature in CSV row {row_number}.")

            if any(value < 0 for value in values):
                raise ValueError(f"Negative feature value in CSV row {row_number}.")

            label = (row.get("label") or "").strip()

            if not label:
                raise ValueError(f"Missing label in CSV row {row_number}.")

            features.append(values)
            labels.append(label)

    if len(features) < 4:
        raise ValueError("Please provide at least four labelled training examples.")

    unique_labels = set(labels)

    if len(unique_labels) < 2:
        raise ValueError("Training data must contain at least two different labels.")

    if any(labels.count(label) < 2 for label in unique_labels):
        raise ValueError("Each label must have at least two training examples.")

    model = RandomForestClassifier(
        n_estimators=100,
        random_state=42,
        class_weight="balanced",
    )

    model.fit(features, labels)

    predictions = []

    for record in records:
        sequence = str(record.seq)
        row = extract_ml_features(sequence)
        values = [[row[column] for column in FEATURE_COLUMNS]]

        prediction = model.predict(values)[0]
        probabilities = model.predict_proba(values)[0]
        confidence = float(max(probabilities))

        predictions.append(
            {
                "sequence_id": record.id,
                "predicted_label": str(prediction),
                "model_probability": round(confidence, 4),
            }
        )

    return {
        "model": "Random Forest",
        "training_examples": len(features),
        "classes": sorted(unique_labels),
        "predictions": predictions,
        "warning": (
            "Model probabilities are not validated accuracy scores. "
            "Use representative, scientifically labelled data and "
            "independent validation before relying on predictions."
        ),
    }


def run_ai_analysis(records, output_dir, training_csv=None):
    """Run DNA intelligence and save a JSON report."""
    results = []

    for record in records:
        sequence = str(record.seq).upper()

        analysis = generate_dna_insights(sequence)

        results.append(
            {
                "sequence_id": record.id,
                "sequence_length": len(sequence),
                "reverse_complement": generate_reverse_complement(sequence),
                "analysis": analysis,
            }
        )

    report = {
        "project": "GeneWeaver",
        "module": "Hybrid DNA Intelligence",
        "analysis_type": "Rule-based bioinformatics",
        "sequences_analyzed": len(results),
        "results": results,
        "machine_learning": {
            "status": "Not run",
            "reason": "No labelled training CSV was supplied.",
        },
    }

    if training_csv:
        report["machine_learning"] = {
            "status": "Completed",
            **run_ml_classification(records, training_csv),
        }

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    report_file = output_path / "geneweaver_ai_report.json"

    with report_file.open("w", encoding="utf-8") as file:
        json.dump(report, file, indent=4)

    print("\nGENEWEAVER AI INTELLIGENCE")
    print("-------------------------")
    print(f"Sequences analyzed : {len(results)}")
    print(f"ML status          : {report['machine_learning']['status']}")
    print(f"JSON report        : {report_file}")

    return report
