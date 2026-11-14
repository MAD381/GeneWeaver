from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import hashlib

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
            invalid_bases.append((position, base))

    return invalid_bases


def validate_sequence(sequence):
    """Check whether a non-empty DNA sequence is valid."""
    sequence = str(sequence).strip().upper()
    return bool(sequence) and not find_invalid_bases(sequence)


# --------------------------------------------------
# 2. DNA SEQUENCE STATISTICS
# --------------------------------------------------


def calculate_sequence_statistics(sequence):
    """Calculate nucleotide counts, GC%, AT% and sequence length."""
    sequence = str(sequence).strip().upper()

    if not validate_sequence(sequence):
        raise ValueError(
            "Cannot calculate statistics for an invalid or empty sequence."
        )

    counts = {
        "A": sequence.count("A"),
        "T": sequence.count("T"),
        "G": sequence.count("G"),
        "C": sequence.count("C"),
    }

    length = len(sequence)
    gc_count = counts["G"] + counts["C"]
    at_count = counts["A"] + counts["T"]
    ambiguous_count = length - sum(counts.values())

    return {
        "length": length,
        "A": counts["A"],
        "T": counts["T"],
        "G": counts["G"],
        "C": counts["C"],
        "GC_percentage": round((gc_count / length) * 100, 2),
        "AT_percentage": round((at_count / length) * 100, 2),
        "ambiguous_bases": ambiguous_count,
    }


# --------------------------------------------------
# 3. SEQUENCE FINGERPRINTING
# --------------------------------------------------


def sequence_fingerprint(sequence):
    """Generate a SHA-256 fingerprint for a DNA sequence."""
    sequence = str(sequence).strip().upper()

    if not validate_sequence(sequence):
        raise ValueError("Cannot generate fingerprint for an invalid sequence.")

    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


# --------------------------------------------------
# 4. DUPLICATE DETECTION
# --------------------------------------------------


def detect_duplicates(records):
    """Find duplicate FASTA identifiers and identical sequences."""

    id_counts = {}
    sequence_ids = {}

    for record in records:
        id_counts[record.id] = id_counts.get(record.id, 0) + 1

        sequence = str(record.seq).upper()
        sequence_ids.setdefault(sequence, []).append(record.id)

    duplicate_ids = [record_id for record_id, count in id_counts.items() if count > 1]

    duplicate_sequences = [ids for ids in sequence_ids.values() if len(ids) > 1]

    return {
        "duplicate_ids": duplicate_ids,
        "duplicate_sequences": duplicate_sequences,
    }


# --------------------------------------------------
# 5. FASTA FILE PARSING
# --------------------------------------------------


def parse_fasta(file_path):
    """Read FASTA records and validate every DNA sequence."""

    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    records = []

    with file_path.open("r", encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence = str(record.seq).upper()

            if not sequence:
                raise ValueError(f"Empty DNA sequence: {record.id}")

            invalid_bases = find_invalid_bases(sequence)

            if invalid_bases:
                details = ", ".join(
                    f"'{base}' at position {position}"
                    for position, base in invalid_bases
                )
                raise ValueError(f"Invalid DNA sequence in {record.id}: {details}")

            record.seq = Seq(sequence)
            records.append(record)

    if not records:
        raise ValueError("No FASTA sequences found.")

    return records


# --------------------------------------------------
# 6. OVERLAPPING DNA CHUNKING
# --------------------------------------------------


def chunk_sequence(sequence, chunk_size=1000, overlap=100):
    """Split a DNA sequence into overlapping chunks."""

    sequence = str(sequence).strip().upper()

    if not validate_sequence(sequence):
        raise ValueError("Invalid or empty DNA sequence.")

    if not isinstance(chunk_size, int) or isinstance(chunk_size, bool):
        raise ValueError("Chunk size must be an integer.")

    if not isinstance(overlap, int) or isinstance(overlap, bool):
        raise ValueError("Overlap must be an integer.")

    if chunk_size <= 0:
        raise ValueError("Chunk size must be positive.")

    if overlap < 0 or overlap >= chunk_size:
        raise ValueError("Overlap must be non-negative and smaller than chunk size.")

    chunks = []
    step = chunk_size - overlap

    for start in range(0, len(sequence), step):
        end = min(start + chunk_size, len(sequence))

        chunks.append((start, end, sequence[start:end]))

        if end == len(sequence):
            break

    return chunks


# --------------------------------------------------
# 7. CHUNK ANALYSIS
# --------------------------------------------------


def analyze_chunks(chunks):
    """Calculate length and GC content for each DNA chunk."""

    results = []

    for index, (start, end, chunk) in enumerate(chunks, start=1):
        statistics = calculate_sequence_statistics(chunk)

        results.append(
            {
                "chunk_number": index,
                "start": start,
                "end": end,
                "length": statistics["length"],
                "GC_percentage": statistics["GC_percentage"],
            }
        )

    return results


# --------------------------------------------------
# 8. COMPLETE FASTA PROCESSING
# --------------------------------------------------


def process_fasta(input_file, output_dir, chunk_size=1000, overlap=100):
    """Process FASTA records, analyze sequences and save DNA chunks."""

    records = parse_fasta(input_file)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / "processed_chunks.fasta"

    processed_records = []
    duplicate_report = detect_duplicates(records)

    total_bases = 0
    total_gc = 0
    total_at = 0
    total_ambiguous = 0

    print("\nGENEWEAVER - WEEK 1 ADVANCED PROCESSING")
    print("---------------------------------------")

    for record in records:
        sequence = str(record.seq)

        statistics = calculate_sequence_statistics(sequence)
        fingerprint = sequence_fingerprint(sequence)

        total_bases += statistics["length"]
        total_gc += statistics["G"] + statistics["C"]
        total_at += statistics["A"] + statistics["T"]
        total_ambiguous += statistics["ambiguous_bases"]

        chunks = chunk_sequence(sequence, chunk_size, overlap)
        chunk_analysis = analyze_chunks(chunks)

        print(f"\nSequence ID: {record.id}")
        print(f"Length: {statistics['length']}")
        print(f"GC content: {statistics['GC_percentage']}%")
        print(f"AT content: {statistics['AT_percentage']}%")
        print(f"Ambiguous bases: {statistics['ambiguous_bases']}")
        print(f"SHA-256: {fingerprint}")
        print(f"Chunks generated: {len(chunks)}")

        for index, (start, end, chunk) in enumerate(chunks, start=1):
            gc_content = chunk_analysis[index - 1]["GC_percentage"]

            processed_records.append(
                SeqRecord(
                    Seq(chunk),
                    id=f"{record.id}_chunk_{index}",
                    description=(f"start={start} end={end} GC={gc_content}%"),
                )
            )

    SeqIO.write(processed_records, output_file, "fasta")

    overall_gc = round((total_gc / total_bases) * 100, 2) if total_bases else 0

    overall_at = round((total_at / total_bases) * 100, 2) if total_bases else 0

    print("\nPROCESSING SUMMARY")
    print("------------------")
    print(f"Sequences processed: {len(records)}")
    print(f"Total bases analyzed: {total_bases}")
    print(f"Overall GC content: {overall_gc}%")
    print(f"Overall AT content: {overall_at}%")
    print(f"Ambiguous bases: {total_ambiguous}")
    print(f"Chunks generated: {len(processed_records)}")

    print("\nDUPLICATE DETECTION")
    print("-------------------")
    print(f"Duplicate IDs: {duplicate_report['duplicate_ids']}")
    print(f"Identical sequences: {duplicate_report['duplicate_sequences']}")

    print(f"\nOutput file: {output_file}")
    print("Status: SUCCESS")

    return processed_records


# --------------------------------------------------
# 9. RUN THE PROCESSOR
# --------------------------------------------------

if __name__ == "__main__":
    process_fasta("data/sample.fasta", "outputs", chunk_size=8, overlap=2)
