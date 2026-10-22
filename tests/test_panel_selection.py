"""Unit tests for conservative classical/cryptic panel selection."""

from __future__ import annotations

import unittest

from treetip_placer.models import PlacementResult, QcFinding, QueryQc
from treetip_placer.cli import calibration_queries, panel_threshold_sweep, select_unique_sequence_controls
from treetip_placer.models import AlignedSequence
from treetip_placer.panel_selection import choose_reference_panel


def result(
    sample_id: str,
    distance: float | None,
    *,
    qc_status: str = "PASS",
    final_status: str = "PASS",
    qc_findings: tuple[QcFinding, ...] = (),
) -> PlacementResult:
    """Build the small result shape needed to test panel selection."""

    qc = QueryQc(
        sample_id=sample_id,
        status=qc_status,
        coverage=1.0,
        observed_sites=100,
        comparable_sites=100,
        nearest_reference="reference",
        nearest_distance=distance,
        findings=list(qc_findings),
    )
    return PlacementResult(sample_id=sample_id, qc=qc, final_status=final_status)


class PanelSelectionTests(unittest.TestCase):
    def test_selects_clearly_closer_panel(self) -> None:
        decision = choose_reference_panel([result("query", 0.01)], [result("query", 0.10)], 0.02)[0]

        self.assertEqual(decision.assignment, "CLASSICAL")
        self.assertEqual(decision.selected_panel, "classical")
        self.assertAlmostEqual(decision.distance_margin or 0.0, 0.09)

    def test_marks_small_distance_margin_ambiguous(self) -> None:
        decision = choose_reference_panel([result("query", 0.01)], [result("query", 0.025)], 0.02)[0]

        self.assertEqual(decision.assignment, "AMBIGUOUS")
        self.assertEqual(decision.decision_status, "UNRESOLVED")
        self.assertIsNone(decision.selected_panel)

    def test_zero_distance_is_not_treated_as_missing(self) -> None:
        decision = choose_reference_panel([result("query", 0.045)], [result("query", 0.0)], 0.02)[0]

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

    def test_rejects_panel_whose_placement_failed_post_qc(self) -> None:
        decision = choose_reference_panel(
            [result("query", 0.01, qc_status="PASS", final_status="FAIL")],
            [result("query", 0.08, qc_status="PASS", final_status="PASS")],
            0.02,
        )[0]

        self.assertEqual(decision.assignment, "CRYPTIC")
        self.assertEqual(decision.selected_panel, "cryptic")
        self.assertEqual(decision.decision_status, "WARN")
        self.assertIn("ONLY_ONE_PANEL_PASSED_QC", decision.findings)

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

    def test_warn_decision_carries_the_selected_placement_s_own_findings(self) -> None:
        # A "Decision: WARN" / "Selected placement: WARN" with an empty
        # findings column is a black box: the reader can see something was
        # flagged but not what. The chosen placement's own QC findings (which
        # already carry a human-readable message, not just a code) must be
        # copied onto the decision so the report can show the real reason.
        warn_finding = QcFinding(
            "DIFFUSE_DISTANCE_WEIGHT", "warning", "Best relative-fit weight is 0.013, below threshold."
        )
        decision = choose_reference_panel(
            [result("query", 0.01, final_status="WARN", qc_findings=(warn_finding,))],
            [result("query", 0.10)],
            0.02,
        )[0]

        self.assertEqual(decision.decision_status, "WARN")
        self.assertEqual(decision.selected_placement_status, "WARN")
        self.assertIn("DIFFUSE_DISTANCE_WEIGHT", decision.findings)

    def test_pass_decision_has_no_findings_to_report(self) -> None:
        # Regression guard: a clean PASS must not gain spurious findings from
        # this change -- there is nothing to explain.
        decision = choose_reference_panel([result("query", 0.01, final_status="PASS")], [result("query", 0.10)], 0.02)[
            0
        ]

        self.assertEqual(decision.decision_status, "PASS")
        self.assertEqual(decision.findings, ())

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
