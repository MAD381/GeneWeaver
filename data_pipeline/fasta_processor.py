from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

VALID_DNA_BASES = set("ACGTRYSWKMBDHVN")


def validate_sequence(sequence):
    sequence = str(sequence).strip().upper()
    return bool(sequence) and all(base in VALID_DNA_BASES for base in sequence)


def parse_fasta(file_path):
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    records = []

    for record in SeqIO.parse(file_path, "fasta"):
        sequence = str(record.seq).upper()

        if not validate_sequence(sequence):
            raise ValueError(f"Invalid DNA sequence: {record.id}")

        record.seq = Seq(sequence)
        records.append(record)

    if not records:
        raise ValueError("No FASTA sequences found.")

    return records


def chunk_sequence(sequence, chunk_size=1000, overlap=100):
    sequence = str(sequence).strip().upper()

    if not validate_sequence(sequence):
        raise ValueError("Invalid DNA sequence.")

    if chunk_size <= 0:
        raise ValueError("Chunk size must be positive.")

    if overlap < 0 or overlap >= chunk_size:
        raise ValueError("Overlap must be smaller than chunk size.")

    chunks = []
    step = chunk_size - overlap

    for start in range(0, len(sequence), step):
        end = min(start + chunk_size, len(sequence))
        chunks.append((start, end, sequence[start:end]))

        if end == len(sequence):
            break

    return chunks


def process_fasta(input_file, output_dir, chunk_size=1000, overlap=100):
    records = parse_fasta(input_file)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / "processed_chunks.fasta"
    processed_records = []

    for record in records:
        chunks = chunk_sequence(record.seq, chunk_size, overlap)

        for index, (start, end, chunk) in enumerate(chunks, start=1):
            processed_records.append(
                SeqRecord(
                    Seq(chunk),
                    id=f"{record.id}_chunk_{index}",
                    description=f"start={start} end={end}",
                )
            )

    SeqIO.write(processed_records, output_file, "fasta")

    print("\nGENEWEAVER - WEEK 1")
    print("-------------------")
    print(f"Sequences processed: {len(records)}")
    print(f"Chunks generated: {len(processed_records)}")
    print(f"Output file: {output_file}")
    print("Status: SUCCESS")

    return processed_records


if __name__ == "__main__":
    process_fasta("data/sample.fasta", "outputs", chunk_size=8, overlap=2)
