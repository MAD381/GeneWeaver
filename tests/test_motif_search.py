import unittest
import tempfile
import os

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from data_pipeline.motif_search import (
    validate_motif,
    find_motif_positions,
    search_motif,
    summarize_motif_search,
    export_motif_results,
)


class TestMotifSearch(unittest.TestCase):

    def test_valid_motif(self):
        self.assertEqual(validate_motif("atg"), "ATG")

    def test_empty_motif(self):
        with self.assertRaises(ValueError):
            validate_motif("")

    def test_invalid_motif(self):
        with self.assertRaises(ValueError):
            validate_motif("ATGX")

    def test_overlapping_matches(self):
        self.assertEqual(find_motif_positions("AAAAA", "AAA"), [1, 2, 3])

    def test_lowercase_sequence(self):
        self.assertEqual(find_motif_positions("atgatg", "ATG"), [1, 4])

    def test_no_match(self):
        self.assertEqual(find_motif_positions("CCCC", "ATG"), [])

    def test_multiple_sequences(self):
        records = [
            SeqRecord(Seq("ATGATG"), id="Gene_1"),
            SeqRecord(Seq("CCCATG"), id="Gene_2"),
        ]

        results = search_motif(records, "ATG")

        self.assertEqual(results[0]["occurrence_count"], 2)
        self.assertEqual(results[1]["occurrence_count"], 1)

    def test_summary(self):
        records = [
            SeqRecord(Seq("ATGATG"), id="Gene_1"),
            SeqRecord(Seq("CCCC"), id="Gene_2"),
        ]

        results = search_motif(records, "ATG")
        summary = summarize_motif_search(results)

        self.assertEqual(summary["total_sequences"], 2)
        self.assertEqual(summary["sequences_with_matches"], 1)
        self.assertEqual(summary["sequences_without_matches"], 1)
        self.assertEqual(summary["total_occurrences"], 2)

    def test_csv_export(self):
        records = [SeqRecord(Seq("ATGATG"), id="Gene_1")]

        results = search_motif(records, "ATG")

        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = os.path.join(temp_dir, "motif_results.csv")

            export_motif_results(results, output_file)

            self.assertTrue(os.path.exists(output_file))

            with open(output_file, encoding="utf-8") as file:
                content = file.read()

            self.assertIn("Gene_1", content)
            self.assertIn("ATG", content)


if __name__ == "__main__":
    unittest.main()
