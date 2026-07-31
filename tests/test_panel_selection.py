"""Unit tests for conservative classical/cryptic panel selection."""

from __future__ import annotations

import unittest

from treetip_placer.models import PlacementResult, QueryQc
from treetip_placer.cli import calibration_queries, panel_threshold_sweep, select_unique_sequence_controls
from treetip_placer.models import AlignedSequence
from treetip_placer.panel_selection import choose_reference_panel


def result(sample_id: str, distance: float | None, *, qc_status: str = "PASS") -> PlacementResult:
    """Build the small result shape needed to test panel selection."""

    qc = QueryQc(
        sample_id=sample_id,
        status=qc_status,
        coverage=1.0,
        observed_sites=100,
        comparable_sites=100,
        nearest_reference="reference",
        nearest_distance=distance,
    )
    return PlacementResult(sample_id=sample_id, qc=qc, final_status="PASS")


class PanelSelectionTests(unittest.TestCase):
    def test_selects_clearly_closer_panel(self) -> None:
        decision = choose_reference_panel(
            [result("query", 0.01)], [result("query", 0.10)], 0.02
        )[0]

        self.assertEqual(decision.assignment, "CLASSICAL")
        self.assertEqual(decision.selected_panel, "classical")
        self.assertAlmostEqual(decision.distance_margin or 0.0, 0.09)

    def test_marks_small_distance_margin_ambiguous(self) -> None:
        decision = choose_reference_panel(
            [result("query", 0.01)], [result("query", 0.025)], 0.02
        )[0]

        self.assertEqual(decision.assignment, "AMBIGUOUS")
        self.assertEqual(decision.decision_status, "UNRESOLVED")
        self.assertIsNone(decision.selected_panel)

    def test_zero_distance_is_not_treated_as_missing(self) -> None:
        decision = choose_reference_panel(
            [result("query", 0.045)], [result("query", 0.0)], 0.02
        )[0]

        self.assertEqual(decision.assignment, "CRYPTIC")
        self.assertEqual(decision.selected_panel, "cryptic")

    def test_rejects_query_when_both_panels_fail_qc(self) -> None:
        decision = choose_reference_panel(
            [result("query", 0.01, qc_status="FAIL")],
            [result("query", 0.02, qc_status="FAIL")],
            0.02,
        )[0]

        self.assertEqual(decision.assignment, "UNRESOLVED")
        self.assertEqual(decision.decision_status, "FAIL")
        self.assertIn("NO_PANEL_PASSED_QC", decision.findings)

    def test_unique_control_selection_excludes_duplicate_sequences(self) -> None:
        controls = select_unique_sequence_controls(
            (
                AlignedSequence("B", tuple("AAAA")),
                AlignedSequence("A", tuple("AAAA")),
                AlignedSequence("C", tuple("AAAT")),
            ),
            10,
        )

        self.assertEqual([control.sample_id for control in controls], ["A", "C"])
        self.assertEqual(
            [control.sample_id for control in calibration_queries("classical", controls)],
            ["classical::A", "classical::C"],
        )

    def test_threshold_sweep_includes_requested_default(self) -> None:
        rows = panel_threshold_sweep(
            [result("query", 0.01)],
            [result("query", 0.10)],
            {"query": "CLASSICAL"},
            0.006,
        )

        self.assertTrue(any(row["threshold"] == 0.006 for row in rows))


if __name__ == "__main__":
    unittest.main()
