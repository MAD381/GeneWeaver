import json
import tempfile
import unittest
from pathlib import Path

from data_pipeline.pairwise_comparison import (
    validate_comparison_sequence,
    align_sequences,
    compare_sequences,
    summarize_comparison,
    export_comparison_report,
)


class TestPairwiseComparison(unittest.TestCase):

    def test_valid_sequence(self):
        self.assertEqual(validate_comparison_sequence("atgc"), "ATGC")

    def test_empty_sequence(self):
        with self.assertRaises(ValueError):
            validate_comparison_sequence("")

    def test_invalid_sequence(self):
        with self.assertRaises(ValueError):
            validate_comparison_sequence("ATGX")

    def test_non_string_sequence(self):
        with self.assertRaises(TypeError):
            validate_comparison_sequence(123)

    def test_identical_sequences(self):
        result = compare_sequences("ATGC", "ATGC")

        self.assertEqual(result["matches"], 4)
        self.assertEqual(result["mismatches"], 0)
        self.assertEqual(result["gaps"], 0)
        self.assertEqual(result["identity_percentage"], 100.0)

    def test_mismatch_detection(self):
        result = compare_sequences("ATGCGT", "ATGAGT")

        self.assertEqual(result["matches"], 5)
        self.assertEqual(result["mismatches"], 1)
        self.assertEqual(result["mismatch_positions"], [4])

    def test_different_sequence_lengths(self):
        result = compare_sequences("ATGCGT", "ATGGT")

        self.assertEqual(result["gaps"], 1)
        self.assertEqual(result["gap_positions"], [4])
        self.assertEqual(
            len(result["aligned_sequence_1"]), len(result["aligned_sequence_2"])
        )

    def test_lowercase_sequences(self):
        result = compare_sequences("atgc", "ATGC")

        self.assertEqual(result["identity_percentage"], 100.0)

    def test_alignment_score(self):
        aligned1, aligned2, score = align_sequences("ATGC", "ATGC")

        self.assertEqual(aligned1, "ATGC")
        self.assertEqual(aligned2, "ATGC")
        self.assertEqual(score, 8)

    def test_summary(self):
        result = compare_sequences("ATGC", "ATGA")
        summary = summarize_comparison(result)

        self.assertEqual(summary["matches"], 3)
        self.assertEqual(summary["mismatches"], 1)
        self.assertEqual(summary["identity_percentage"], 75.0)

    def test_json_report_export(self):
        result = compare_sequences("ATGC", "ATGA")

        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = Path(temp_dir) / "comparison.json"

            saved_path = export_comparison_report(result, output_file)

            self.assertTrue(Path(saved_path).exists())

            with open(saved_path, "r", encoding="utf-8") as file:
                saved_data = json.load(file)

            self.assertEqual(saved_data["identity_percentage"], 75.0)


if __name__ == "__main__":
    unittest.main()
