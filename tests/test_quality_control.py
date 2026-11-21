import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from data_pipeline.quality_control import (
    analyze_sequence_quality,
    run_quality_control,
)


class TestQualityControl(unittest.TestCase):

    def test_valid_sequence_passes(self):
        result = analyze_sequence_quality("ATGCGCAT")

        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["length"], 8)

    def test_invalid_bases_fail(self):
        result = analyze_sequence_quality("ATGCXYZ")

        self.assertEqual(result["status"], "FAIL")
        self.assertTrue(len(result["invalid_bases"]) > 0)

    def test_high_ambiguous_content_warning(self):
        result = analyze_sequence_quality("ATGCNNNNNN")

        self.assertEqual(result["status"], "WARNING")
        self.assertGreater(result["ambiguous_percentage"], 10)

    def test_short_sequence_warning(self):
        result = analyze_sequence_quality(
            "ATGC",
            min_length=10,
        )

        self.assertEqual(result["status"], "WARNING")

    def test_gc_content_calculation(self):
        result = analyze_sequence_quality("GGCCATAT")

        self.assertEqual(result["GC_percentage"], 50.0)

    def test_empty_sequence_fails(self):
        result = analyze_sequence_quality("")

        self.assertEqual(result["status"], "FAIL")

    def test_multiple_sequence_report(self):
        records = [
            SeqRecord(Seq("ATGCATGC"), id="Gene_A"),
            SeqRecord(Seq("GGCCGGCC"), id="Gene_B"),
        ]

        report = run_quality_control(records)

        self.assertEqual(
            report["summary"]["total_sequences"],
            2,
        )

        self.assertEqual(len(report["results"]), 2)

    def test_empty_records_raise_error(self):
        with self.assertRaises(ValueError):
            run_quality_control([])


if __name__ == "__main__":
    unittest.main()
