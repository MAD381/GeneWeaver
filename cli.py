
import argparse
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from data_pipeline.fasta_processor import (
    parse_fasta,
    chunk_sequence,
)


# ==========================================================
# COMMAND-LINE INPUT VALIDATION
# ==========================================================

def positive_integer(value):
    """Accept only positive integers."""
    try:
        number = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Enter a valid integer.")

    if number <= 0:
        raise argparse.ArgumentTypeError("Value must be greater than zero.")

    return number


def non_negative_integer(value):
    """Accept zero or positive integers."""
    try:
        number = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Enter a valid integer.")

    if number < 0:
        raise argparse.ArgumentTypeError("Value cannot be negative.")

    return number


# ==========================================================
# WEEK 2: FASTA PROCESSING PIPELINE
# ==========================================================

def run_pipeline(
    input_file,
    output_dir,
    chunk_size,
    overlap,
    ai_analysis=False,
    training_csv=None,
):
    """Process FASTA records and optionally run DNA intelligence."""

    input_path = Path(input_file)
    output_path = Path(output_dir)

    if not input_path.is_file():
        raise FileNotFoundError(
            f"Input FASTA file not found: {input_path}"
        )

    if chunk_size <= 0:
        raise ValueError("Chunk size must be greater than zero.")

    if overlap < 0 or overlap >= chunk_size:
        raise ValueError(
            "Overlap must be non-negative and smaller than chunk size."
        )

    output_path.mkdir(parents=True, exist_ok=True)

    # Use the existing Week 1 FASTA parser.
    records = parse_fasta(str(input_path))

    if not records:
        raise ValueError("No FASTA sequences found.")

    processed_records = []
    total_chunks = 0

    print("\n====================================")
    print("       GENEWEAVER - WEEK 1 & 2")
    print("====================================")
    print("FASTA Processing Pipeline")
    print("------------------------------------")
    print("Processing sequences...\n")

    # Generate overlapping chunks.
    for record in tqdm(records, desc="Processing", unit="sequence"):
        record_id = record.id
        sequence = str(record.seq)

        chunks = chunk_sequence(
            sequence,
            chunk_size,
            overlap,
        )

        for index, (start, end, chunk) in enumerate(chunks, start=1):
            processed_record = SeqRecord(
                Seq(chunk),
                id=f"{record_id}_chunk_{index}",
                description=f"start={start} end={end}",
            )

            processed_records.append(processed_record)
            total_chunks += 1

    # Save the processed FASTA file.
    output_file = output_path / "processed_chunks.fasta"

    SeqIO.write(
        processed_records,
        str(output_file),
        "fasta",
    )

    print("\nFASTA processing completed!")
    print("------------------------------------")
    print(f"Sequences processed : {len(records)}")
    print(f"Chunks generated    : {total_chunks}")
    print(f"Output file         : {output_file}")
    print("Status              : SUCCESS")

    # ======================================================
    # HYBRID AI INTEGRATION
    # ======================================================

    if ai_analysis or training_csv:

        print("\nStarting Hybrid DNA Intelligence...")

        from ai_engine.dna_intelligence import run_ai_analysis

        run_ai_analysis(
            records=records,
            output_dir=output_path,
            training_csv=training_csv,
        )

    return {
        "sequences_processed": len(records),
        "chunks_generated": total_chunks,
        "output_file": str(output_file),
    }


# ==========================================================
# COMMAND-LINE INTERFACE
# ==========================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "GeneWeaver: FASTA Processing and "
            "Hybrid DNA Intelligence Pipeline"
        )
    )

    # Existing Week 2 arguments.
    parser.add_argument(
        "--input",
        default="data/sample.fasta",
        help="Input FASTA file path",
    )

    parser.add_argument(
        "--output-dir",
        default="outputs",
        help="Output directory",
    )

    parser.add_argument(
        "--chunk-size",
        type=positive_integer,
        default=1000,
        help="Number of DNA bases in each chunk",
    )

    parser.add_argument(
        "--overlap",
        type=non_negative_integer,
        default=100,
        help="Overlapping bases between chunks",
    )

    # New Hybrid AI arguments.
    parser.add_argument(
        "--ai-analysis",
        action="store_true",
        help="Enable DNA intelligence analysis",
    )

    parser.add_argument(
        "--training-csv",
        default=None,
        help="Optional labelled CSV file for Random Forest training",
    )

    args = parser.parse_args()

    try:
        run_pipeline(
            input_file=args.input,
            output_dir=args.output_dir,
            chunk_size=args.chunk_size,
            overlap=args.overlap,
            ai_analysis=args.ai_analysis,
            training_csv=args.training_csv,
        )

    except (ValueError, OSError, RuntimeError) as error:
        print(
            f"\nERROR: {error}",
            file=sys.stderr,
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
