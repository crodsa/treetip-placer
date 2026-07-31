"""Tests for the human-readable automatic-placement report."""

from __future__ import annotations

import shutil
import unittest
from pathlib import Path

from treetip_placer.panel_selection import PanelDecision
from treetip_placer.reports import render_tree_figure, render_tree_zoom, write_auto_report
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
            )
            contents = (root / "report.html").read_text(encoding="utf-8")

            self.assertIn("sample&lt;unsafe&gt;", contents)
            self.assertIn("panel_assignment.tsv", contents)
            self.assertIn("cryptic.selected_placed_tree.nwk", contents)
            self.assertIn("Cryptic selected tree", contents)
            self.assertNotIn("Classical selected tree", contents)
            self.assertIn("sample&lt;unsafe&gt; placement zoom", contents)
            self.assertEqual(contents.count('class="dendrogram"'), 2)
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
        self.assertNotIn(">reference</text>", figure)

    def test_placement_zoom_uses_the_selected_query(self) -> None:
        figure = render_tree_zoom(
            "cryptic",
            Tree.from_newick("(reference:0.1,query:0.2);"),
            "query",
        )

        self.assertIn("placement zoom", figure)
        self.assertIn(">query</text>", figure)


if __name__ == "__main__":
    unittest.main()
