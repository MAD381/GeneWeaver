import unittest

from data_pipeline.fasta_processor import (
    validate_sequence,
    chunk_sequence,
)


class TestFastaProcessor(unittest.TestCase):

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


if __name__ == "__main__":
    unittest.main()

