"""End-to-end command-line regression tests."""

from __future__ import annotations

import shutil
import unittest
from pathlib import Path

from datetime import datetime, timezone

from treetip_placer.cli import (
    ORGANISM_REFERENCE_PACKAGES,
    default_classify_output_directory,
    main,
    resolve_mafft_executable,
    resolve_organism_reference_package,
    select_clade_calibration_controls,
)
from treetip_placer.models import AlignedSequence


class CliTests(unittest.TestCase):
    def test_package_version_matches_pyproject(self) -> None:
        """Pin __init__.__version__ to pyproject.toml so a version bump cannot update one and not the other.

        __version__ is stamped into --version and into the provenance block of
        every calibration report and summary JSON; a mismatch silently mislabels
        every artifact a released version produces.
        """

        import tomllib

        from treetip_placer import __version__

        pyproject = tomllib.loads(
            (Path(__file__).resolve().parent.parent / "pyproject.toml").read_text(encoding="utf-8")
        )

        self.assertEqual(__version__, pyproject["project"]["version"])

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

    def test_resolve_organism_reference_package_finds_the_bundled_package(self) -> None:
        self.assertIn("cdifficile-mlst", ORGANISM_REFERENCE_PACKAGES)

        resolved = resolve_organism_reference_package("cdifficile-mlst")

        self.assertTrue(resolved.is_file())
        self.assertEqual(resolved.name, "paired_reference.json")

    def test_resolve_organism_reference_package_rejects_a_missing_file(self) -> None:
        # ORGANISM_REFERENCE_PACKAGES maps a key to a directory name under
        # reference-packages/; if that directory doesn't exist (a broken
        # checkout, or a future registry entry added before its data lands),
        # the error must name what was expected rather than raising a bare
        # FileNotFoundError from deep inside argparse dispatch.
        ORGANISM_REFERENCE_PACKAGES["nonexistent-test-organism"] = "does-not-exist-on-disk"
        try:
            with self.assertRaisesRegex(ValueError, "does-not-exist-on-disk"):
                resolve_organism_reference_package("nonexistent-test-organism")
        finally:
            del ORGANISM_REFERENCE_PACKAGES["nonexistent-test-organism"]

    def test_resolve_mafft_executable_prefers_the_explicit_path(self) -> None:
        self.assertEqual(resolve_mafft_executable("/custom/path/to/mafft"), "/custom/path/to/mafft")

    def test_resolve_mafft_executable_raises_a_clear_error_when_not_found(self) -> None:
        with self.assertRaisesRegex(ValueError, "MAFFT"):
            resolve_mafft_executable(None, path_lookup=lambda _name: None)

    def test_resolve_mafft_executable_falls_back_to_path_lookup(self) -> None:
        self.assertEqual(
            resolve_mafft_executable(None, path_lookup=lambda _name: "/usr/local/bin/mafft"),
            "/usr/local/bin/mafft",
        )

    def test_default_classify_output_directory_is_stem_based_and_stable(self) -> None:
        fixed_now = datetime(2026, 8, 6, 12, 30, 45, tzinfo=timezone.utc)

        name = default_classify_output_directory("my_isolates.fasta", now=fixed_now)

        self.assertEqual(name, "my_isolates_classify_20260806T123045Z")


if __name__ == "__main__":
    unittest.main()
