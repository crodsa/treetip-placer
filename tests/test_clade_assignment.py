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
