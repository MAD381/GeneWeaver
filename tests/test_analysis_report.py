import csv
import json
import tempfile
import unittest
from pathlib import Path

from data_pipeline.fasta_processor import parse_fasta
from data_pipeline.analysis_report import (
    build_analysis_report,
    export_json_report,
    export_csv_report,
)


class TestAnalysisReport(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.records = parse_fasta("data/sample.fasta")

    def test_build_analysis_report(self):
        report = build_analysis_report(self.records)

        self.assertEqual(
            report["sequence_count"],
            len(self.records),
        )
        self.assertEqual(
            len(report["sequences"]),
            len(self.records),
        )
        self.assertEqual(
            report["project"],
            "GeneWeaver",
        )

    def test_report_contains_sequence_statistics(self):
        report = build_analysis_report(self.records)

        for sequence in report["sequences"]:
            self.assertIn("sequence_id", sequence)
            self.assertIn("sequence_length", sequence)
            self.assertIn("statistics", sequence)
            self.assertIn("GC_percentage", sequence["statistics"])
            self.assertIn("reverse_complement", sequence)

    def test_json_export(self):
        report = build_analysis_report(self.records)

        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = Path(temp_dir) / "report.json"

            export_json_report(report, output_file)

            with output_file.open(encoding="utf-8") as file:
                saved_report = json.load(file)

            self.assertEqual(
                saved_report["sequence_count"],
                len(self.records),
            )

    def test_csv_export(self):
        report = build_analysis_report(self.records)

        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = Path(temp_dir) / "report.csv"

            export_csv_report(report, output_file)

            with output_file.open(
                encoding="utf-8",
                newline="",
            ) as file:
                rows = list(csv.DictReader(file))

            self.assertEqual(len(rows), len(self.records))
            self.assertIn("GC_percentage", rows[0])
            self.assertIn("sequence_id", rows[0])

    def test_empty_records_raise_error(self):
        with self.assertRaises(ValueError):
            build_analysis_report([])


if __name__ == "__main__":
    unittest.main()
