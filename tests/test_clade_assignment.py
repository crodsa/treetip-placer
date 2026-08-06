"""Tests for fixed-tree discrete clade assignment."""

from __future__ import annotations

import shutil
import unittest
from pathlib import Path

from treetip_placer.clade_assignment import (
    assign_phylogenetic_clades,
    load_tip_clades,
    write_clade_assignment_outputs,
)
from treetip_placer.models import EdgeFit, PlacementResult, QueryQc
from treetip_placer.panel_selection import PanelDecision
from treetip_placer.tree import Tree


class CladeAssignmentTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tree = Tree.from_newick("((A:0.1,B:0.1):0.1,C:0.2);")
        tip_a = self.tree.tips[self.tree.tip_index["A"]]
        self.result = PlacementResult(
            sample_id="query",
            qc=QueryQc("query", "PASS", 1.0, 10, 10, "A", 0.01),
            final_status="PASS",
            best_fit=EdgeFit(
                edge_id=self.tree.edge_by_child[tip_a.node_id].edge_id,
                parent_id=tip_a.parent.node_id if tip_a.parent else -1,
                child_id=tip_a.node_id,
                distal_length=0.05,
                pendant_length=0.01,
                objective=0.0,
                rmse=0.0,
            ),
        )
        self.decision = PanelDecision(
            sample_id="query",
            assignment="CRYPTIC",
            decision_status="PASS",
            selected_panel="cryptic",
            selected_placement_status="PASS",
            classical_nearest_reference="classical_reference",
            classical_nearest_distance=0.2,
            cryptic_nearest_reference="A",
            cryptic_nearest_distance=0.01,
            distance_margin=0.19,
            minimum_distance_separation=0.01,
            findings=(),
        )

    def test_assigns_single_clade_below_selected_edge(self) -> None:
        assignments = assign_phylogenetic_clades(
            [self.decision],
            {"cryptic": (self.tree, [self.result])},
            {"cryptic": {"A": "C-II", "B": "C-II", "C": "C-III"}},
        )

        assignment = assignments[0]
        self.assertEqual(assignment.phylogenetic_clade, "C-II")
        self.assertEqual(assignment.clade_assignment_status, "PASS")
        self.assertEqual(assignment.assignment_method, "selected_edge_descendant_tip_labels")
        self.assertEqual(assignment.supporting_reference_tips, ("A",))

    def test_marks_missing_metadata_as_unresolved(self) -> None:
        assignment = assign_phylogenetic_clades(
            [self.decision],
            {"cryptic": (self.tree, [self.result])},
            {"cryptic": {}},
        )[0]

        self.assertEqual(assignment.clade_assignment_status, "UNRESOLVED")
        self.assertIn("NO_TIP_CLADE_METADATA_FOR_SELECTED_PANEL", assignment.findings)

    def test_unresolved_placement_cannot_emit_a_clade_call(self) -> None:
        unresolved_result = PlacementResult(
            sample_id="query",
            qc=self.result.qc,
            final_status="UNRESOLVED",
            best_fit=self.result.best_fit,
        )
        assignment = assign_phylogenetic_clades(
            [self.decision],
            {"cryptic": (self.tree, [unresolved_result])},
            {"cryptic": {"A": "C-II", "B": "C-II", "C": "C-III"}},
        )[0]

        self.assertIsNone(assignment.phylogenetic_clade)
        self.assertEqual(assignment.clade_assignment_status, "UNRESOLVED")
        self.assertEqual(assignment.findings, ("SELECTED_PLACEMENT_UNRESOLVED",))

    def test_nearest_reference_arbitrates_when_edge_consensus_disagrees(self) -> None:
        # Selected edge is C's own edge: its descendant subtree is just {C},
        # so consensus trivially reads C's label ("C-III"). But the query's
        # single nearest reference by real genetic distance is A ("C-II") --
        # a real disagreement between tree-edge topology and raw distance,
        # the exact failure calibrate-clade-assignment hit on the 2026-08-06
        # panel: a selected edge's descendant subtree can span tips that are
        # topologically nearby but genetically distant on a weakly resolved
        # tree, while the nearest-reference distance is unaffected by that.
        tip_c = self.tree.tips[self.tree.tip_index["C"]]
        conflicting_result = PlacementResult(
            sample_id="query",
            qc=QueryQc("query", "PASS", 1.0, 10, 10, "A", 0.01),
            final_status="PASS",
            best_fit=EdgeFit(
                edge_id=self.tree.edge_by_child[tip_c.node_id].edge_id,
                parent_id=tip_c.parent.node_id if tip_c.parent else -1,
                child_id=tip_c.node_id,
                distal_length=0.1,
                pendant_length=0.01,
                objective=0.0,
                rmse=0.0,
            ),
        )

        assignment = assign_phylogenetic_clades(
            [self.decision],
            {"cryptic": (self.tree, [conflicting_result])},
            {"cryptic": {"A": "C-II", "B": "C-II", "C": "C-III"}},
        )[0]

        self.assertEqual(assignment.phylogenetic_clade, "C-II")
        self.assertEqual(assignment.clade_assignment_status, "WARN")
        self.assertEqual(
            assignment.assignment_method,
            "nearest_labelled_reference_edge_consensus_conflict",
        )
        self.assertIn(
            "SELECTED_EDGE_CONSENSUS_DISAGREES_WITH_NEAREST_REFERENCE", assignment.findings
        )
        self.assertEqual(assignment.supporting_reference_tips, ("A",))

    def test_edge_consensus_agreeing_with_nearest_reference_still_passes(self) -> None:
        # Same scenario as test_assigns_single_clade_below_selected_edge, made
        # explicit as a regression guard for the new arbitration branch: when
        # the edge consensus and the nearest reference agree, PASS must still
        # come from the ordinary edge-consensus path, not the conflict branch.
        assignment = assign_phylogenetic_clades(
            [self.decision],
            {"cryptic": (self.tree, [self.result])},
            {"cryptic": {"A": "C-II", "B": "C-II", "C": "C-III"}},
        )[0]

        self.assertEqual(assignment.phylogenetic_clade, "C-II")
        self.assertEqual(assignment.clade_assignment_status, "PASS")
        self.assertEqual(assignment.assignment_method, "selected_edge_descendant_tip_labels")

    def test_unlabelled_nearest_reference_cannot_arbitrate(self) -> None:
        # The query's nearest reference by distance is C, but C carries no
        # label in this TSV. With nothing to cross-check against, the edge
        # consensus (A's label) must still stand as before -- the arbitration
        # only overrides when it has a labelled opinion to override with.
        tip_c = self.tree.tips[self.tree.tip_index["C"]]
        result_with_unlabelled_nearest = PlacementResult(
            sample_id="query",
            qc=QueryQc("query", "PASS", 1.0, 10, 10, tip_c.name, 0.02),
            final_status="PASS",
            best_fit=self.result.best_fit,
        )

        assignment = assign_phylogenetic_clades(
            [self.decision],
            {"cryptic": (self.tree, [result_with_unlabelled_nearest])},
            {"cryptic": {"A": "C-II", "B": "C-II"}},
        )[0]

        self.assertEqual(assignment.phylogenetic_clade, "C-II")
        self.assertEqual(assignment.clade_assignment_status, "PASS")
        self.assertEqual(assignment.assignment_method, "selected_edge_descendant_tip_labels")

    def test_mixed_edge_uses_nearest_labelled_reference_as_warn_fallback(self) -> None:
        mixed_child = self.tree.root.children[0]
        mixed_result = PlacementResult(
            sample_id="query",
            qc=self.result.qc,
            final_status="PASS",
            best_fit=EdgeFit(
                edge_id=self.tree.edge_by_child[mixed_child.node_id].edge_id,
                parent_id=mixed_child.parent.node_id if mixed_child.parent else -1,
                child_id=mixed_child.node_id,
                distal_length=0.05,
                pendant_length=0.01,
                objective=0.0,
                rmse=0.0,
            ),
        )
        assignment = assign_phylogenetic_clades(
            [self.decision],
            {"cryptic": (self.tree, [mixed_result])},
            {"cryptic": {"A": "C-II", "B": "C-III", "C": "C-III"}},
        )[0]

        self.assertEqual(assignment.phylogenetic_clade, "C-II")
        self.assertEqual(assignment.clade_assignment_status, "WARN")
        self.assertEqual(
            assignment.assignment_method,
            "nearest_labelled_reference_mixed_edge_fallback",
        )

    def test_loader_and_output_reject_unknown_tip_ids(self) -> None:
        root = Path.cwd() / ".treetip-clade-assignment-test"
        if root.exists():
            shutil.rmtree(root)
        root.mkdir()
        try:
            metadata = root / "tip_clades.tsv"
            metadata.write_text("reference_tip_id\tclade\nA\tC-II\n", encoding="utf-8")
            self.assertEqual(load_tip_clades(metadata, self.tree, "cryptic"), {"A": "C-II"})
            metadata.write_text("reference_tip_id\tclade\nUNKNOWN\tC-II\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unknown reference tip"):
                load_tip_clades(metadata, self.tree, "cryptic")

            output_assignment = assign_phylogenetic_clades(
                [self.decision],
                {"cryptic": (self.tree, [self.result])},
                {"cryptic": {"A": "C-II"}},
            )
            write_clade_assignment_outputs(root, output_assignment, {"method": "test"})
            contents = (root / "clade_assignment.tsv").read_text(encoding="utf-8")
            self.assertIn("phylogenetic_clade", contents)
            self.assertIn("C-II", contents)
        finally:
            if root.exists():
                shutil.rmtree(root)


if __name__ == "__main__":
    unittest.main()
