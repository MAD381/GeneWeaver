import argparse
import sys
from pathlib import Path
from tqdm import tqdm

from data_pipeline.fasta_processor import process_fasta

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
    """Run FASTA processing and optional DNA intelligence."""

    input_path = Path(input_file)
    output_path = Path(output_dir)

    if not input_path.is_file():
        raise FileNotFoundError(f"Input FASTA file not found: {input_path}")

    if overlap >= chunk_size:
        raise ValueError("Overlap must be smaller than chunk size.")

    print("\n====================================")
    print("       GENEWEAVER - WEEK 1 & 2")
    print("====================================")
    print("FASTA Processing Pipeline")
    print("------------------------------------")

    progress_bar = None

    def update_progress(done, total, chunks):
        nonlocal progress_bar

        if progress_bar is None:
            progress_bar = tqdm(
                total=total,
                desc="Processing",
                unit="sequence",
            )

        progress_bar.update(1)
        progress_bar.set_postfix(chunks=chunks)

    try:
        result = process_fasta(
            input_file=input_path,
            output_dir=output_path,
            chunk_size=chunk_size,
            overlap=overlap,
            progress_callback=update_progress,
        )
    finally:
        if progress_bar is not None:
            progress_bar.close()

    print("\nFASTA processing completed!")
    print("------------------------------------")
    print(f"Sequences processed : {result['total_sequences']}")
    print(f"Chunks generated    : {result['total_chunks']}")
    print(f"Processing time     : " f"{result['processing_time_seconds']} seconds")
    print(f"Output file         : {result['output_file']}")
    print("Status              : SUCCESS")

    # ======================================================
    # HYBRID AI INTEGRATION
    # ======================================================

    if ai_analysis or training_csv:

        print("\nStarting Hybrid DNA Intelligence...")

        from ai_engine.dna_intelligence import run_ai_analysis

        records = None

        # AI analysis uses the original FASTA records.
        from data_pipeline.fasta_processor import parse_fasta

        records = parse_fasta(input_path)

        run_ai_analysis(
            records=records,
            output_dir=output_path,
            training_csv=training_csv,
        )

    return result


# ==========================================================
# COMMAND-LINE INTERFACE
# ==========================================================


def main():

    parser = argparse.ArgumentParser(
        description=(
            "GeneWeaver: FASTA Processing and " "Hybrid DNA Intelligence Pipeline"
        )
    )

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
