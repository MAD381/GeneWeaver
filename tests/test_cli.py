import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import cli


class TestGeneWeaverCLI(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.test_folder = Path(self.temp_dir.name)

        self.input_file = self.test_folder / "sample.fasta"
        self.output_dir = self.test_folder / "output"

        self.input_file.write_text(
            ">sequence_1\n"
            "ATGCGTACGTAGCTAGCTAG\n"
            ">sequence_2\n"
            "GCTAGCTAGCATCGATCGTA\n"
        )

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_positive_integer(self):
        self.assertEqual(cli.positive_integer("100"), 100)

    def test_positive_integer_rejects_zero(self):
        with self.assertRaises(Exception):
            cli.positive_integer("0")

    def test_non_negative_integer(self):
        self.assertEqual(cli.non_negative_integer("0"), 0)

    def test_non_negative_integer_rejects_negative(self):
        with self.assertRaises(Exception):
            cli.non_negative_integer("-1")

    def test_run_pipeline_success(self):
        result = cli.run_pipeline(
            input_file=self.input_file,
            output_dir=self.output_dir,
            chunk_size=10,
            overlap=2,
        )

        self.assertEqual(result["total_sequences"], 2)
        self.assertGreater(result["total_chunks"], 0)
        self.assertGreaterEqual(result["processing_time_seconds"], 0)

        output_file = self.output_dir / "processed_chunks.fasta"
        self.assertTrue(output_file.exists())

    def test_missing_input_file(self):
        missing_file = self.test_folder / "missing.fasta"

        with self.assertRaises(FileNotFoundError):
            cli.run_pipeline(
                input_file=missing_file,
                output_dir=self.output_dir,
                chunk_size=10,
                overlap=2,
            )

    def test_invalid_overlap(self):
        with self.assertRaises(ValueError):
            cli.run_pipeline(
                input_file=self.input_file,
                output_dir=self.output_dir,
                chunk_size=10,
                overlap=10,
            )

    def test_progress_callback_is_used(self):
        from data_pipeline.fasta_processor import process_fasta

        with patch(
            "cli.process_fasta",
            wraps=process_fasta,
        ) as mocked_process:

            cli.run_pipeline(
                input_file=self.input_file,
                output_dir=self.output_dir,
                chunk_size=10,
                overlap=2,
            )

            self.assertTrue(mocked_process.called)

            callback = mocked_process.call_args.kwargs["progress_callback"]

            self.assertTrue(callable(callback))


if __name__ == "__main__":
    unittest.main()
