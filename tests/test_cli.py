"""End-to-end command-line regression tests."""

from __future__ import annotations

import shutil
import unittest
from pathlib import Path

from treetip_placer.cli import main, select_clade_calibration_controls
from treetip_placer.models import AlignedSequence


class CliTests(unittest.TestCase):
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
                    "--min-comparable-sites",
                    "2",
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertIn("query_A", (output / "trees" / "placed_tree.nwk").read_text(encoding="utf-8"))
        finally:
            if output.exists():
                shutil.rmtree(output)

    def test_clade_calibration_writes_auditable_outputs(self) -> None:
        root = Path.cwd() / ".treetip-clade-calibration-test"
        if root.exists():
            shutil.rmtree(root)
        root.mkdir()
        try:
            tree = root / "tree.nwk"
            alignment = root / "reference.fasta"
            clades = root / "tip_clades.tsv"
            output = root / "output"
            tree.write_text("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);\n", encoding="utf-8")
            alignment.write_text(
                ">A\nAAAAAA\n>B\nAAAAAT\n>C\nCCCCCC\n>D\nCCCCCG\n", encoding="utf-8"
            )
            clades.write_text(
                "reference_tip_id\tclade\nA\tC-I\nB\tC-I\nC\tC-II\nD\tC-II\n",
                encoding="utf-8",
            )
            exit_code = main(
                [
                    "calibrate-clade-assignment",
                    "--tree", str(tree),
                    "--reference-alignment", str(alignment),
                    "--tip-clades", str(clades),
                    "--output", str(output),
                    "--preset", "thorough",
                    "--min-comparable-sites", "2",
                    "--unresolved-weight", "0",
                ]
            )
            self.assertEqual(exit_code, 0)
            records = (output / "clade_calibration_records.tsv").read_text(encoding="utf-8")
            self.assertIn("truth_clade", records)
            self.assertIn("relative_fit_weight", records)
            self.assertTrue((output / "clade_confusion.tsv").is_file())
            self.assertTrue((output / "clade_calibration_exclusions.tsv").is_file())
            self.assertTrue((output / "clade_calibration_summary.json").is_file())
            self.assertTrue((output / "report.html").is_file())
            report = (output / "report.html").read_text(encoding="utf-8")
            self.assertIn("Clade-assignment calibration", report)
        finally:
            if root.exists():
                shutil.rmtree(root)

    def test_clade_controls_exclude_singleton_labelled_clades(self) -> None:
        controls = select_clade_calibration_controls(
            (
                AlignedSequence("A", ("A", "A")),
                AlignedSequence("B", ("A", "T")),
                AlignedSequence("C", ("C", "C")),
            ),
            {"A": "C-I", "B": "C-I", "C": "C-VI"},
            10,
        )

        self.assertEqual(controls, ("A", "B"))



if __name__ == "__main__":
    unittest.main()
