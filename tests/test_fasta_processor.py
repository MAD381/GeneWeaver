import unittest
import tempfile
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from data_pipeline.fasta_processor import (
    validate_sequence,
    chunk_sequence,
    calculate_sequence_statistics,
    sequence_fingerprint,
    detect_duplicates,
    analyze_chunks,
    generate_reverse_complement,
    process_fasta,
)


class TestFastaProcessor(unittest.TestCase):

    # Original Week 1 tests

    def test_valid_dna(self):
        self.assertTrue(validate_sequence("ATGC"))

    def test_invalid_dna(self):
        self.assertFalse(validate_sequence("ATGX"))

    def test_empty_sequence(self):
        self.assertFalse(validate_sequence(""))

    def test_lowercase_sequence(self):
        self.assertTrue(validate_sequence("atgc"))

    def test_sequence_chunking(self):
        chunks = chunk_sequence("ATGCGTACGT", chunk_size=6, overlap=2)

        self.assertEqual(len(chunks), 2)
        self.assertEqual(chunks[0][2], "ATGCGT")
        self.assertEqual(chunks[1][2], "GTACGT")

    def test_invalid_chunk_size(self):
        with self.assertRaises(ValueError):
            chunk_sequence("ATGC", chunk_size=0)

    def test_invalid_overlap(self):
        with self.assertRaises(ValueError):
            chunk_sequence("ATGC", chunk_size=4, overlap=4)

    # Sequence statistics tests

    def test_sequence_statistics(self):
        result = calculate_sequence_statistics("ATGCCGNN")

        self.assertEqual(result["length"], 8)
        self.assertEqual(result["A"], 1)
        self.assertEqual(result["T"], 1)
        self.assertEqual(result["G"], 2)
        self.assertEqual(result["C"], 2)
        self.assertEqual(result["GC_percentage"], 50.0)
        self.assertEqual(result["AT_percentage"], 25.0)
        self.assertEqual(result["ambiguous_bases"], 2)

    def test_sequence_fingerprint(self):
        first = sequence_fingerprint("ATGC")
        second = sequence_fingerprint("ATGC")

        self.assertEqual(first, second)
        self.assertEqual(len(first), 64)

    def test_duplicate_detection(self):
        records = [
            SeqRecord(Seq("ATGC"), id="sequence_1"),
            SeqRecord(Seq("ATGC"), id="sequence_2"),
            SeqRecord(Seq("GGCC"), id="sequence_1"),
        ]

        result = detect_duplicates(records)

        self.assertEqual(result["duplicate_ids"], ["sequence_1"])
        self.assertEqual(
            result["duplicate_sequences"],
            [["sequence_1", "sequence_2"]],
        )

    def test_chunk_analysis(self):
        chunks = [
            (0, 4, "ATGC"),
            (2, 6, "GCGC"),
        ]

        result = analyze_chunks(chunks)

        self.assertEqual(len(result), 2)
        self.assertEqual(result[0]["length"], 4)
        self.assertEqual(result[0]["GC_percentage"], 50.0)
        self.assertEqual(result[1]["GC_percentage"], 100.0)

    # Reverse complement tests

    def test_reverse_complement(self):
        result = generate_reverse_complement("ATGC")
        self.assertEqual(result, "GCAT")

    def test_reverse_complement_lowercase(self):
        result = generate_reverse_complement("atgc")
        self.assertEqual(result, "GCAT")

    def test_reverse_complement_invalid_dna(self):
        with self.assertRaises(ValueError):
            generate_reverse_complement("ATGX")

    # FASTA processing progress callback test

    def test_process_fasta_progress_callback(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_file = Path(temp_dir) / "sample.fasta"
            output_dir = Path(temp_dir) / "output"

            input_file.write_text(
                ">seq1\nATGCGT\n>seq2\nGGCCTA\n",
                encoding="utf-8",
            )

            progress_updates = []

            def callback(done, total, chunks):
                progress_updates.append((done, total, chunks))

            result = process_fasta(
                input_file,
                output_dir,
                chunk_size=4,
                overlap=1,
                progress_callback=callback,
            )

            self.assertEqual(result["total_sequences"], 2)
            self.assertEqual(len(progress_updates), 2)
            self.assertEqual(progress_updates[-1], (2, 2, 4))
            self.assertTrue(Path(result["output_file"]).exists())

    # FASTA processing time test

    def test_process_fasta_processing_time(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_file = Path(temp_dir) / "sample.fasta"
            output_dir = Path(temp_dir) / "output"

            input_file.write_text(
                ">seq1\nATGCGT\n",
                encoding="utf-8",
            )

            result = process_fasta(input_file, output_dir)

            self.assertIn("processing_time_seconds", result)
            self.assertGreaterEqual(
                result["processing_time_seconds"],
                0,
            )


if __name__ == "__main__":
    unittest.main()
