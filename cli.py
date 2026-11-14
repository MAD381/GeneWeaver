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

# --------------------------------------------------
# 1. COMMAND-LINE INPUT VALIDATION
# --------------------------------------------------


def positive_integer(value):
    """Accept only positive integers."""

    try:
        number = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Value must be an integer.")

    if number <= 0:
        raise argparse.ArgumentTypeError("Value must be greater than zero.")

    return number


def non_negative_integer(value):
    """Accept zero or positive integers."""

    try:
        number = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Value must be zero or greater.")

    if number < 0:
        raise argparse.ArgumentTypeError("Value cannot be negative.")

    return number


# --------------------------------------------------
# 2. MAIN FASTA PROCESSING PIPELINE
# --------------------------------------------------


def run_pipeline(input_file, output_dir, chunk_size, overlap):
    """Process FASTA sequences with validation and progress reporting."""

    input_path = Path(input_file)
    output_path = Path(output_dir)

    # Validate input file
    if not input_path.exists():
        raise FileNotFoundError(f"Input FASTA file not found: {input_path}")

    if not input_path.is_file():
        raise ValueError(f"Input path is not a file: {input_path}")

    # Validate chunk settings
    if not isinstance(chunk_size, int) or chunk_size <= 0:
        raise ValueError("Chunk size must be a positive integer.")

    if not isinstance(overlap, int) or overlap < 0:
        raise ValueError("Overlap must be a non-negative integer.")

    if overlap >= chunk_size:
        raise ValueError("Overlap must be smaller than chunk size.")

    # Validate output directory
    if output_path.exists() and not output_path.is_dir():
        raise ValueError(f"Output path exists but is not a directory: {output_path}")

    # Create output directory if necessary
    output_path.mkdir(parents=True, exist_ok=True)

    # Parse FASTA using the existing Week 1 function
    records = parse_fasta(str(input_path))

    if not records:
        raise ValueError("No FASTA sequences found in the input file.")

    processed_records = []
    total_chunks = 0

    print("\nGENEWEAVER - WEEK 2")
    print("-------------------")
    print("FASTA Processing Pipeline")
    print("------------------------")
    print(f"Input file          : {input_path}")
    print(f"Chunk size          : {chunk_size}")
    print(f"Overlap             : {overlap}")
    print(f"Sequences found     : {len(records)}")
    print("\nProcessing sequences...\n")

    # Process sequences with progress bar
    for record in tqdm(records, desc="Processing", unit="sequence", total=len(records)):

        record_id = record.id
        sequence = str(record.seq)

        # Generate chunks using Week 1 function
        chunks = chunk_sequence(sequence, chunk_size=chunk_size, overlap=overlap)

        # Convert chunks into FASTA records
        for index, (start, end, chunk) in enumerate(chunks, start=1):

            processed_record = SeqRecord(
                Seq(chunk),
                id=f"{record_id}_chunk_{index}",
                description=f"start={start} end={end}",
            )

            processed_records.append(processed_record)
            total_chunks += 1

    # Save processed FASTA records
    output_file = output_path / "processed_chunks.fasta"

    SeqIO.write(processed_records, str(output_file), "fasta")

    # Final processing summary
    summary = {
        "total_sequences": len(records),
        "total_chunks": total_chunks,
        "output_file": str(output_file),
        "status": "SUCCESS",
    }

    print("\nProcessing completed successfully!")
    print("-----------------------------------")
    print(f"Sequences processed : {summary['total_sequences']}")
    print(f"Chunks generated    : {summary['total_chunks']}")
    print(f"Output file         : {summary['output_file']}")
    print(f"Status              : {summary['status']}")

    return summary


# --------------------------------------------------
# 3. COMMAND-LINE INTERFACE
# --------------------------------------------------


def main():

    parser = argparse.ArgumentParser(
        description=(
            "GeneWeaver - FASTA Processing Pipeline "
            "with sequence chunking and progress reporting"
        ),
        epilog=(
            "Example: python cli.py --input data/sample.fasta "
            "--output-dir outputs --chunk-size 1000 --overlap 100"
        ),
    )

    parser.add_argument(
        "--input", default="data/sample.fasta", help="Input FASTA file path"
    )

    parser.add_argument(
        "--output-dir",
        default="outputs",
        help="Directory where processed FASTA will be saved",
    )

    parser.add_argument(
        "--chunk-size",
        type=positive_integer,
        default=1000,
        help="Number of bases in each chunk (must be positive)",
    )

    parser.add_argument(
        "--overlap",
        type=non_negative_integer,
        default=100,
        help="Overlapping bases between chunks (default: 100)",
    )

    args = parser.parse_args()

    try:

        run_pipeline(
            input_file=args.input,
            output_dir=args.output_dir,
            chunk_size=args.chunk_size,
            overlap=args.overlap,
        )

    except (ValueError, OSError) as error:

        print(f"\nERROR: {error}", file=sys.stderr)

        sys.exit(1)


# --------------------------------------------------
# 4. PROGRAM ENTRY POINT
# --------------------------------------------------

if __name__ == "__main__":
    main()
