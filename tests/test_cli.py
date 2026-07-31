"""End-to-end command-line regression tests."""

from __future__ import annotations

import shutil
import unittest
from pathlib import Path

from treetip_placer.cli import main


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



if __name__ == "__main__":
    unittest.main()
