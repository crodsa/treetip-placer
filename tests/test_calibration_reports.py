"""Tests for the shared calibration HTML report."""

from __future__ import annotations

import shutil
import unittest
from dataclasses import dataclass
from pathlib import Path

from treetip_placer.calibration_reports import write_calibration_report
from treetip_placer.cli import build_panel_gate_summary, write_calibration_records


SUMMARY = {
    "version": "0.7.0",
    "method": "fixed_tree_leave_one_unique_sequence_out_placement_accuracy",
    "warning": "Internal topology diagnostics only: this is not phylogenetic support.",
    "tree_path": "/private/tree.nwk",
    "tree_sha256": "abc123",
    "summary": {"controls": 2, "outcome_correct": 1, "outcome_far": 1},
}

COLUMNS = ("source_sample_id", "outcome", "edge_hops_from_truth")

RECORDS = (
    {"source_sample_id": "safe_sample", "outcome": "CORRECT", "edge_hops_from_truth": 0},
    {"source_sample_id": "sample<unsafe>", "outcome": "FAR", "edge_hops_from_truth": 7},
)

ARTIFACTS = (("placement_accuracy_records.tsv", "Per-control records (TSV)"),)


class CalibrationReportTests(unittest.TestCase):
    def setUp(self) -> None:
        self.root = Path.cwd() / ".treetip-calibration-report-test"
        if self.root.exists():
            shutil.rmtree(self.root)
        self.root.mkdir()

    def tearDown(self) -> None:
        if self.root.exists():
            shutil.rmtree(self.root)

    def render(self, **overrides) -> str:
        parameters = {
            "output_directory": self.root,
            "title": "Placement-accuracy calibration",
            "summary": SUMMARY,
            "records": RECORDS,
            "record_columns": COLUMNS,
            "outcome_column": "outcome",
            "problem_outcomes": frozenset({"FAR"}),
            "artifacts": ARTIFACTS,
        }
        parameters.update(overrides)
        write_calibration_report(**parameters)
        return (self.root / "report.html").read_text(encoding="utf-8")

    def test_escapes_sample_identifiers(self) -> None:
        document = self.render()

        self.assertIn("sample&lt;unsafe&gt;", document)
        self.assertNotIn("sample<unsafe>", document)

    def test_marks_only_problem_rows(self) -> None:
        document = self.render()

        self.assertEqual(document.count('class="problem"'), 1)

    def test_renders_the_scope_warning_and_provenance(self) -> None:
        document = self.render()

        self.assertIn("not phylogenetic support", document)
        self.assertIn("abc123", document)
        self.assertIn("0.7.0", document)

    def test_threshold_sweep_is_omitted_when_empty(self) -> None:
        document = self.render()

        self.assertNotIn("Threshold sweep", document)

    def test_threshold_sweep_highlights_the_default_row(self) -> None:
        document = self.render(
            summary={**SUMMARY, "default_threshold": 0.01},
            threshold_sweep=(
                {"threshold": 0.0, "controls": 2, "called": 2, "called_accuracy": 0.5},
                {"threshold": 0.01, "controls": 2, "called": 1, "called_accuracy": 1.0},
            ),
        )

        self.assertIn("Threshold sweep", document)
        self.assertEqual(document.count('class="default-threshold"'), 1)

    def test_links_every_artifact(self) -> None:
        document = self.render()

        self.assertIn('href="placement_accuracy_records.tsv"', document)

    def test_real_panel_gate_summary_renders_outcome_counts(self) -> None:
        """Regression guard: drives the real production summary-building path.

        calibrate-panel-gate has no CLI test (it needs MAFFT), so build_panel_gate_summary
        is called directly here — not re-implemented — so that deleting its "summary" key
        makes this test fail rather than a hand-written fixture masking the regression.
        """

        @dataclass
        class FakeDecision:
            sample_id: str
            assignment: str
            decision_status: str
            classical_nearest_distance: float | None
            cryptic_nearest_distance: float | None
            distance_margin: float | None

        decisions = (
            FakeDecision("classical::a", "CLASSICAL", "PASS", 0.0, 0.2, 0.2),
            FakeDecision("classical::b", "CLASSICAL", "PASS", 0.0, 0.15, 0.15),
            FakeDecision("cryptic::c", "CLASSICAL", "PASS", 0.05, 0.2, -0.15),
        )
        truth_by_id = {"classical::a": "CLASSICAL", "classical::b": "CLASSICAL", "cryptic::c": "CRYPTIC"}
        source_id_by_id = {"classical::a": "a", "classical::b": "b", "cryptic::c": "c"}

        conflict_count, calibration_rows, fields = write_calibration_records(
            self.root, decisions, truth_by_id, source_id_by_id
        )
        summary = build_panel_gate_summary(
            classical_controls_count=2,
            cryptic_controls_count=1,
            default_threshold=0.01,
            conflict_count=conflict_count,
            sweep=({"threshold": 0.01, "controls": 3, "called": 3, "called_accuracy": 1.0},),
            calibration_rows=calibration_rows,
            classical_reference_alignment_sha256="def456",
            cryptic_reference_alignment_sha256="ghi789",
        )

        write_calibration_report(
            output_directory=self.root,
            title="Panel-gate calibration",
            summary=summary,
            records=calibration_rows,
            record_columns=fields,
            outcome_column="outcome",
            problem_outcomes=frozenset({"INCORRECT_CALL"}),
            artifacts=(("calibration_records.tsv", "Per-control records (TSV)"),),
            threshold_sweep=summary["threshold_sweep"],
        )
        document = (self.root / "report.html").read_text(encoding="utf-8")

        self.assertNotIn("No outcome counts were recorded.", document)
        self.assertIn("<dt>classical_controls</dt>", document)
        self.assertIn("<dt>cryptic_controls</dt>", document)
        self.assertIn("<dt>flagged_panel_label_conflicts</dt>", document)
        self.assertIn("<title>Panel-gate calibration</title>", document)


if __name__ == "__main__":
    unittest.main()
