"""End-to-end command-line regression tests."""

from __future__ import annotations

import shutil
import unittest
import json
from pathlib import Path

from treetip_placer.cli import main


class CliTests(unittest.TestCase):
    """Verify that a MLST profile run produces its documented artifacts."""

    def test_profile_placement_writes_outputs(self) -> None:
        fixtures = Path("tests/fixtures")
        output = Path.cwd() / ".treetip-test-output"
        if output.exists():
            shutil.rmtree(output)
        try:
            exit_code = main(
                [
                    "place",
                    "--tree",
                    str(fixtures / "tree.nwk"),
                    "--reference-profiles",
                    str(fixtures / "reference_profiles.tsv"),
                    "--query-profiles",
                    str(fixtures / "query_profiles.tsv"),
                    "--output",
                    str(output),
                    "--preset",
                    "thorough",
                    "--threads",
                    "2",
                    "--unresolved-weight",
                    "0",
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertTrue((output / "trees" / "placed_tree.nwk").is_file())
            self.assertTrue((output / "placements" / "placements.tsv").is_file())
            self.assertTrue((output / "placements" / "placements.jplace").is_file())
            self.assertTrue((output / "qc" / "query_qc.tsv").is_file())
            self.assertIn("query_A", (output / "trees" / "placed_tree.nwk").read_text(encoding="utf-8"))
            self.assertIn("query_B", (output / "trees" / "placed_tree.nwk").read_text(encoding="utf-8"))
        finally:
            if output.exists():
                shutil.rmtree(output)

    def test_alignment_placement_writes_outputs(self) -> None:
        fixtures = Path("tests/fixtures")
        output = Path.cwd() / ".treetip-alignment-test-output"
        if output.exists():
            shutil.rmtree(output)
        try:
            exit_code = main(
                [
                    "place",
                    "--tree",
                    str(fixtures / "tree.nwk"),
                    "--reference-alignment",
                    str(fixtures / "reference_alignment_full.fasta"),
                    "--query-alignment",
                    str(fixtures / "query_alignment_full.fasta"),
                    "--output",
                    str(output),
                    "--preset",
                    "thorough",
                    "--threads",
                    "2",
                    "--unresolved-weight",
                    "0",
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertIn("query_A", (output / "trees" / "placed_tree.nwk").read_text(encoding="utf-8"))
        finally:
            if output.exists():
                shutil.rmtree(output)

    def test_mlst_leave_one_out_writes_validation_outputs(self) -> None:
        fixtures = Path("tests/fixtures")
        output = Path.cwd() / ".treetip-validation-test-output"
        if output.exists():
            shutil.rmtree(output)
        try:
            exit_code = main(
                [
                    "validate-mlst",
                    "--tree",
                    str(fixtures / "validation_tree.nwk"),
                    "--reference-profiles",
                    str(fixtures / "validation_reference_profiles.tsv"),
                    "--output",
                    str(output),
                    "--max-holdouts",
                    "2",
                    "--preset",
                    "thorough",
                    "--threads",
                    "2",
                    "--unresolved-weight",
                    "0",
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertTrue((output / "leave_one_out.tsv").is_file())
            self.assertTrue((output / "leave_one_out.json").is_file())
            self.assertTrue((output / "validation_summary.json").is_file())
            summary = json.loads((output / "validation_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["holdouts_requested"], 2)
        finally:
            if output.exists():
                shutil.rmtree(output)


if __name__ == "__main__":
    unittest.main()
