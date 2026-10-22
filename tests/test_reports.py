"""Tests for the human-readable automatic-placement report."""

from __future__ import annotations

import shutil
import unittest
from pathlib import Path

from treetip_placer.panel_selection import PanelDecision
from treetip_placer.clade_assignment import CladeAssignment
from treetip_placer.reports import (
    render_clade_context,
    render_decision_row,
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
            (trees / "classical.selected_placed_tree.nwk").write_text("(classical_reference:0.1);\n", encoding="utf-8")
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
                {"cryptic": {"cryptic_reference": "C-II"}},
            )
            contents = (root / "report.html").read_text(encoding="utf-8")

            self.assertIn("sample&lt;unsafe&gt;", contents)
            self.assertIn("panel_assignment.tsv", contents)
            self.assertIn("cryptic.selected_placed_tree.nwk", contents)
            self.assertIn("global clade context", contents)
            self.assertNotIn("Classical panel", contents)
            self.assertIn("sample&lt;unsafe&gt; placement branch", contents)
            self.assertEqual(contents.count('class="dendrogram"'), 2)
            self.assertIn("C-II", contents)
            self.assertIn("clade_assignment.tsv", contents)
        finally:
            if root.exists():
                shutil.rmtree(root)

    def test_decision_row_merges_panel_and_clade_findings(self) -> None:
        # A reader looking at one row of the report table should see every
        # reason a WARN was raised, not just whichever layer (panel decision
        # vs. clade assignment) happened to have one -- the two findings
        # sources must not silently overwrite or hide each other.
        decision = PanelDecision(
            sample_id="sample1",
            assignment="CRYPTIC",
            decision_status="WARN",
            selected_panel="cryptic",
            selected_placement_status="WARN",
            classical_nearest_reference="classical",
            classical_nearest_distance=0.1,
            cryptic_nearest_reference="cryptic",
            cryptic_nearest_distance=0.01,
            distance_margin=0.09,
            minimum_distance_separation=0.01,
            findings=("DIFFUSE_DISTANCE_WEIGHT",),
        )
        clade_assignment = CladeAssignment(
            sample_id="sample1",
            selected_panel="cryptic",
            phylogenetic_clade="C-II",
            clade_assignment_status="WARN",
            assignment_method="nearest_labelled_reference_edge_consensus_conflict",
            nearest_reference="cryptic_reference",
            supporting_reference_tips=("cryptic_reference",),
            findings=("SELECTED_EDGE_CONSENSUS_DISAGREES_WITH_NEAREST_REFERENCE",),
        )

        row = render_decision_row(decision, clade_assignment)

        self.assertIn("DIFFUSE_DISTANCE_WEIGHT", row)
        self.assertIn("SELECTED_EDGE_CONSENSUS_DISAGREES_WITH_NEAREST_REFERENCE", row)

    def test_decision_row_with_no_clade_assignment_still_shows_decision_findings(self) -> None:
        decision = PanelDecision(
            sample_id="sample1",
            assignment="AMBIGUOUS",
            decision_status="UNRESOLVED",
            selected_panel=None,
            selected_placement_status=None,
            classical_nearest_reference="classical",
            classical_nearest_distance=0.1,
            cryptic_nearest_reference="cryptic",
            cryptic_nearest_distance=0.095,
            distance_margin=0.005,
            minimum_distance_separation=0.01,
            findings=("PANEL_DISTANCE_MARGIN_TOO_SMALL",),
        )

        row = render_decision_row(decision, None)

        self.assertIn("PANEL_DISTANCE_MARGIN_TOO_SMALL", row)

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

        self.assertIn("placement branch", figure)
        self.assertIn(">query</text>", figure)
        self.assertIn(">reference</text>", figure)
        self.assertIn("placement-branch", figure)

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
        self.assertIn("clade C-II", figure)
        self.assertIn("global clade context", figure)
        self.assertNotIn(">A [", figure)
        self.assertNotIn(">B [", figure)
        self.assertNotIn(">C [", figure)
        self.assertNotIn(">D [", figure)

    def test_clade_context_persists_without_metadata_and_omits_reference_ids(self) -> None:
        figure = render_clade_context(
            "cryptic",
            Tree.from_newick("((A:0.1,query:0.1):0.1,B:0.1);"),
            "query",
            None,
            {},
        )

        self.assertIn("Unlabelled (1 tips)", figure)
        self.assertIn(">query</text>", figure)
        self.assertNotIn(">A</text>", figure)
        self.assertNotIn(">B</text>", figure)

    def test_placement_branch_collapses_excess_topological_neighbours(self) -> None:
        figure = render_tree_zoom(
            "cryptic",
            Tree.from_newick(
                "((sister:0.1,query:0.1):0.1,((far_1:0.1,far_2:0.1):0.1,(far_3:0.1,far_4:0.1):0.1):0.1);"
            ),
            "query",
            max_reference_tips=2,
        )

        self.assertIn(">sister</text>", figure)
        self.assertIn("… (", figure)
        self.assertIn("placement-branch", figure)

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

    def test_multi_query_batch_isolates_other_queries_from_reference_counts_and_zoom(self) -> None:
        tree = Tree.from_newick("((ref_1:0.1,ref_2:0.1):0.1,((ref_3:0.1,query_A:0.01):0.05,query_B:0.02):0.1);")
        tip_clades = {"ref_1": "C-I", "ref_2": "C-I", "ref_3": "C-II"}

        context_a = render_clade_context(
            "cryptic",
            tree,
            "query_A",
            "C-II",
            tip_clades,
            query_ids={"query_A", "query_B"},
        )
        zoom_a = render_tree_zoom(
            "cryptic",
            tree,
            "query_A",
            query_ids={"query_A", "query_B"},
        )

        self.assertIn("C-I (2 tips)", context_a)
        self.assertIn("C-II (1 tips)", context_a)
        self.assertNotIn("unlabelled", context_a.lower())
        self.assertNotIn("query_B", context_a)
        self.assertIn(">query_A</text>", zoom_a)
        self.assertIn(">ref_3</text>", zoom_a)
        self.assertNotIn("query_B", zoom_a)


if __name__ == "__main__":
    unittest.main()
