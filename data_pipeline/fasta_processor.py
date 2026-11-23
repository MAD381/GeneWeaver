from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import hashlib
import time

VALID_DNA_BASES = set("ACGTRYSWKMBDHVN")


# --------------------------------------------------
# 1. DNA VALIDATION
# --------------------------------------------------


def find_invalid_bases(sequence):
    """Return invalid DNA characters and their 1-based positions."""

    sequence = str(sequence).strip().upper()
    invalid_bases = []

    for position, base in enumerate(sequence, start=1):
        if base not in VALID_DNA_BASES:
            invalid_bases.append((base, position))

    return invalid_bases


def validate_sequence(sequence):
    """Check whether a DNA sequence is valid and non-empty."""

    sequence = str(sequence).strip().upper()

    if not sequence:
        return False

    return len(find_invalid_bases(sequence)) == 0


# --------------------------------------------------
# 2. DNA SEQUENCE STATISTICS
# --------------------------------------------------


def calculate_sequence_statistics(sequence):
    """Calculate nucleotide counts and sequence composition."""

    sequence = str(sequence).strip().upper()

    if not validate_sequence(sequence):
        raise ValueError("Invalid or empty DNA sequence.")

    length = len(sequence)

    counts = {
        "A": sequence.count("A"),
        "T": sequence.count("T"),
        "G": sequence.count("G"),
        "C": sequence.count("C"),
    }

    gc_percentage = ((counts["G"] + counts["C"]) / length) * 100
    at_percentage = ((counts["A"] + counts["T"]) / length) * 100

    ambiguous_bases = sum(1 for base in sequence if base not in "ATGC")

    return {
        "length": length,
        "A": counts["A"],
        "T": counts["T"],
        "G": counts["G"],
        "C": counts["C"],
        "GC_percentage": round(gc_percentage, 2),
        "AT_percentage": round(at_percentage, 2),
        "ambiguous_bases": ambiguous_bases,
    }


# --------------------------------------------------
# 3. SEQUENCE FINGERPRINT
# --------------------------------------------------


def sequence_fingerprint(sequence):
    """Generate a SHA-256 fingerprint for a DNA sequence."""

    sequence = str(sequence).strip().upper()

    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


# --------------------------------------------------
# 4. REVERSE COMPLEMENT GENERATION
# --------------------------------------------------


def generate_reverse_complement(sequence):
    """Generate the reverse complement of a valid DNA sequence."""

    sequence = str(sequence).strip().upper()

    if not validate_sequence(sequence):
        raise ValueError(
            "Cannot generate reverse complement for invalid or empty DNA sequence."
        )

    return str(Seq(sequence).reverse_complement())


# --------------------------------------------------
# 5. DUPLICATE DETECTION
# --------------------------------------------------


def detect_duplicates(records):
    """Detect duplicate FASTA IDs and identical sequences."""

    seen_ids = set()
    duplicate_ids = []
    sequence_groups = {}

    for record in records:
        record_id = record.id
        sequence = str(record.seq).upper()

        if record_id in seen_ids:
            duplicate_ids.append(record_id)

        seen_ids.add(record_id)

        sequence_groups.setdefault(sequence, []).append(record_id)

    duplicate_sequences = [ids for ids in sequence_groups.values() if len(ids) > 1]

    return {
        "duplicate_ids": duplicate_ids,
        "duplicate_sequences": duplicate_sequences,
    }


# --------------------------------------------------
# 6. FASTA PARSING
# --------------------------------------------------


def parse_fasta(file_path):
    """Parse FASTA file and validate every sequence."""

    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {file_path}")

    with open(file_path, "r", encoding="utf-8") as fasta_file:
        records = list(SeqIO.parse(fasta_file, "fasta"))

    if not records:
        raise ValueError("No FASTA sequences found.")

    for record in records:
        sequence = str(record.seq).upper()

        if not validate_sequence(sequence):
            invalid = find_invalid_bases(sequence)

            raise ValueError(f"Invalid DNA sequence in record '{record.id}': {invalid}")

    return records


# --------------------------------------------------
# 7. SEQUENCE CHUNKING
# --------------------------------------------------


def chunk_sequence(sequence, chunk_size=1000, overlap=100):
    """Split a sequence into overlapping chunks."""

    sequence = str(sequence).strip().upper()

    if not validate_sequence(sequence):
        raise ValueError("Cannot chunk invalid or empty DNA sequence.")

    if not isinstance(chunk_size, int) or chunk_size <= 0:
        raise ValueError("chunk_size must be a positive integer.")

    if not isinstance(overlap, int) or overlap < 0:
        raise ValueError("overlap must be a non-negative integer.")

    if overlap >= chunk_size:
        raise ValueError("overlap must be smaller than chunk_size.")

    chunks = []
    step = chunk_size - overlap

    for start in range(0, len(sequence), step):
        end = min(start + chunk_size, len(sequence))
        chunk = sequence[start:end]

        chunks.append((start, end, chunk))

        if end == len(sequence):
            break

    return chunks


# --------------------------------------------------
# 8. CHUNK ANALYSIS
# --------------------------------------------------


def analyze_chunks(chunks):
    """Calculate statistics for every sequence chunk."""

    results = []

    for start, end, sequence in chunks:
        statistics = calculate_sequence_statistics(sequence)

        results.append(
            {
                "start": start,
                "end": end,
                "length": statistics["length"],
                "GC_percentage": statistics["GC_percentage"],
                "ambiguous_bases": statistics["ambiguous_bases"],
            }
        )

    return results


# --------------------------------------------------
# 9. COMPLETE FASTA PROCESSING
# --------------------------------------------------


def process_fasta(
    input_file,
    output_dir,
    chunk_size=1000,
    overlap=100,
    progress_callback=None,
):
    """Process FASTA sequences and save chunks efficiently."""

    start_time = time.perf_counter()

    records = parse_fasta(input_file)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / "processed_chunks.fasta"

    total_chunks = 0
    total_sequences = len(records)

    with open(output_file, "w", encoding="utf-8") as fasta_output:
        for sequence_index, record in enumerate(records, start=1):
            sequence = str(record.seq).upper()

            chunks = chunk_sequence(
                sequence,
                chunk_size=chunk_size,
                overlap=overlap,
            )

            for chunk_index, (start, end, chunk) in enumerate(chunks, start=1):
                chunk_record = SeqRecord(
                    Seq(chunk),
                    id=f"{record.id}_chunk_{chunk_index}",
                    description=f"start={start} end={end}",
                )

                SeqIO.write(chunk_record, fasta_output, "fasta")
                total_chunks += 1

            if progress_callback is not None:
                progress_callback(
                    sequence_index,
                    total_sequences,
                    total_chunks,
                )

    processing_time = time.perf_counter() - start_time

    return {
        "total_sequences": total_sequences,
        "total_chunks": total_chunks,
        "output_file": str(output_file),
        "processing_time_seconds": round(processing_time, 4),
    }
