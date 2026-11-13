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


def run_pipeline(input_file, output_dir, chunk_size, overlap):
    """Process FASTA sequences with a progress bar."""

    input_path = Path(input_file)
    output_path = Path(output_dir)

    # Validate input file
    if not input_path.exists():
        raise FileNotFoundError(f"Input FASTA file not found: {input_path}")

    # Validate chunk settings
    if chunk_size <= 0:
        raise ValueError("Chunk size must be greater than zero.")

    if overlap < 0 or overlap >= chunk_size:
        raise ValueError("Overlap must be non-negative and smaller than chunk size.")

    # Create output directory
    output_path.mkdir(parents=True, exist_ok=True)

    # Parse FASTA records using Week 1 function
    records = parse_fasta(str(input_path))

    if not records:
        raise ValueError("No FASTA sequences found.")

    processed_records = []
    total_chunks = 0

    print("\nGENEWEAVER - WEEK 2")
    print("-------------------")
    print("FASTA Processing Pipeline")
    print("------------------------")
    print("Processing sequences...\n")

    # Process records with progress bar
    for record in tqdm(records, desc="Processing", unit="sequence"):

        record_id = record.id
        sequence = str(record.seq)

        # Generate chunks using Week 1 function
        chunks = chunk_sequence(sequence, chunk_size, overlap)

        # Each chunk contains start, end and DNA sequence
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

    # Display final summary
    print("\nProcessing completed successfully!")
    print("-----------------------------------")
    print(f"Sequences processed : {len(records)}")
    print(f"Chunks generated    : {total_chunks}")
    print(f"Output file         : {output_file}")
    print("Status              : SUCCESS")


def main():

    parser = argparse.ArgumentParser(description="GeneWeaver FASTA Processing Pipeline")

    parser.add_argument(
        "--input", default="data/sample.fasta", help="Input FASTA file path"
    )

    parser.add_argument("--output-dir", default="outputs", help="Output directory")

    parser.add_argument(
        "--chunk-size", type=int, default=1000, help="Number of bases in each chunk"
    )

    parser.add_argument(
        "--overlap", type=int, default=100, help="Overlapping bases between chunks"
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


if __name__ == "__main__":
    main()
