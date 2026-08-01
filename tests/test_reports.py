"""Tests for the human-readable automatic-placement report."""

from __future__ import annotations

import shutil
import unittest
from pathlib import Path

from treetip_placer.panel_selection import PanelDecision
from treetip_placer.clade_assignment import CladeAssignment
from treetip_placer.reports import (
    render_clade_context,
    render_tree_figure,
    render_tree_zoom,
    write_auto_report,
)
from treetip_placer.tree import Tree


class ReportTests(unittest.TestCase):
    def test_report_escapes_sample_identifiers_and_links_artifacts(self) -> None:
        root = Path.cwd() / ".treetip-report-test"
        if root.exists():
            shutil.rmtree(root)
        root.mkdir()
        try:
            trees = root / "selected_trees"
            trees.mkdir()
            (trees / "classical.selected_placed_tree.nwk").write_text(
                "(classical_reference:0.1);\n", encoding="utf-8"
            )
            (trees / "cryptic.selected_placed_tree.nwk").write_text(
                "(cryptic_reference:0.1,sample<unsafe>:0.1);\n", encoding="utf-8"
            )
            decision = PanelDecision(
                sample_id="sample<unsafe>",
                assignment="CRYPTIC",
                decision_status="PASS",
                selected_panel="cryptic",
                selected_placement_status="PASS",
                classical_nearest_reference="classical",
                classical_nearest_distance=0.1,
                cryptic_nearest_reference="cryptic",
                cryptic_nearest_distance=0.01,
                distance_margin=0.09,
                minimum_distance_separation=0.01,
                findings=(),
            )
            write_auto_report(
                root,
                [decision],
                {"minimum_distance_separation": 0.01, "reference_package_path": "package.json"},
                [
                    CladeAssignment(
                        sample_id="sample<unsafe>",
                        selected_panel="cryptic",
                        phylogenetic_clade="C-II",
                        clade_assignment_status="PASS",
                        assignment_method="test",
                        nearest_reference="cryptic_reference",
                        supporting_reference_tips=("cryptic_reference",),
                        findings=(),
                    )
                ],
            )
            contents = (root / "report.html").read_text(encoding="utf-8")

            self.assertIn("sample&lt;unsafe&gt;", contents)
            self.assertIn("panel_assignment.tsv", contents)
            self.assertIn("cryptic.selected_placed_tree.nwk", contents)
            self.assertIn("Cryptic selected tree", contents)
            self.assertNotIn("Classical selected tree", contents)
            self.assertIn("sample&lt;unsafe&gt; placement zoom", contents)
            self.assertEqual(contents.count('class="dendrogram"'), 2)
            self.assertIn("C-II", contents)
            self.assertIn("clade_assignment.tsv", contents)
        finally:
            if root.exists():
                shutil.rmtree(root)

    def test_dendrogram_highlights_only_selected_queries(self) -> None:
        figure = render_tree_figure(
            "cryptic",
            Tree.from_newick("(reference:0.1,query:0.2);"),
            {"query"},
        )

        self.assertIn("<svg", figure)
        self.assertIn("query-tip", figure)
        self.assertIn(">query</text>", figure)
        self.assertIn(">reference</text>", figure)

    def test_placement_zoom_uses_the_selected_query(self) -> None:
        figure = render_tree_zoom(
            "cryptic",
            Tree.from_newick("(reference:0.1,query:0.2);"),
            "query",
        )

        self.assertIn("placement zoom", figure)
        self.assertIn(">query</text>", figure)

    def test_clade_context_collapses_non_target_pure_subtrees(self) -> None:
        figure = render_clade_context(
            "cryptic",
            Tree.from_newick("((A:0.1,B:0.1):0.1,(query:0.1,(C:0.1,D:0.1):0.1):0.1);"),
            "query",
            "C-II",
            {"A": "C-I", "B": "C-I", "C": "C-III", "D": "C-III"},
        )

        self.assertIn("C-I (2 tips)", figure)
        self.assertIn("C-III (2 tips)", figure)
        self.assertIn(">query</text>", figure)
        self.assertIn("target C-II expanded", figure)

    def test_report_allows_unresolved_clade_assignments(self) -> None:
        root = Path.cwd() / ".treetip-report-unresolved-clade-test"
        if root.exists():
            shutil.rmtree(root)
        root.mkdir()
        try:
            decision = PanelDecision(
                sample_id="query",
                assignment="CRYPTIC",
                decision_status="PASS",
                selected_panel="cryptic",
                selected_placement_status="PASS",
                classical_nearest_reference="classical",
                classical_nearest_distance=0.2,
                cryptic_nearest_reference="cryptic",
                cryptic_nearest_distance=0.01,
                distance_margin=0.19,
                minimum_distance_separation=0.01,
                findings=(),
            )
            unresolved = CladeAssignment(
                sample_id="query",
                selected_panel="cryptic",
                phylogenetic_clade=None,
                clade_assignment_status="UNRESOLVED",
                assignment_method=None,
                nearest_reference="cryptic",
                supporting_reference_tips=(),
                findings=("NO_TIP_CLADE_METADATA_FOR_SELECTED_PANEL",),
            )
            write_auto_report(root, [decision], {"minimum_distance_separation": 0.01}, [unresolved])
            self.assertIn("UNRESOLVED", (root / "report.html").read_text(encoding="utf-8"))
        finally:
            if root.exists():
                shutil.rmtree(root)


if __name__ == "__main__":
    unittest.main()
