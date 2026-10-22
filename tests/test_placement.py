"""Regression tests for the fixed-tree distance placement MVP."""

from __future__ import annotations

import unittest
from pathlib import Path

from treetip_placer.alignment import read_fasta_alignment, require_matching_alignment_widths
from treetip_placer.models import AlignedSequence, PlacementConfig
from treetip_placer.outputs import render_placed_tree
from treetip_placer.placement import place_query
from treetip_placer.tree import Tree


class PlacementTests(unittest.TestCase):
    """Exercise placement math, QC, and interoperable tree rendering."""

    def setUp(self) -> None:
        self.tree = Tree.from_newick("(A:0.5,B:0,C:0.5);")
        self.references = (
            AlignedSequence("A", tuple("AAAACCCC")),
            AlignedSequence("B", tuple("AAAAGGGG")),
            AlignedSequence("C", tuple("CCCCGGGG")),
        )
        self.config = PlacementConfig(
            preset="thorough",
            min_coverage=0.5,
            min_comparable_sites=2,
            unresolved_weight=0.0,
        )

    def test_exact_sequence_joins_its_reference_terminal_edge(self) -> None:
        result = place_query(
            AlignedSequence("query_A", self.references[0].characters), self.references, self.tree, self.config
        )

        self.assertEqual(result.final_status, "PASS")
        self.assertIsNotNone(result.best_fit)
        assert result.best_fit is not None
        self.assertEqual(self.tree.nodes[result.best_fit.child_id].name, "A")
        self.assertAlmostEqual(result.best_fit.distal_length, 0.5)
        self.assertAlmostEqual(result.best_fit.pendant_length, 0.0)
        self.assertAlmostEqual(result.best_fit.objective, 0.0)

    def test_missing_sequence_fails_preplacement_qc(self) -> None:
        result = place_query(AlignedSequence("empty", (None,) * 8), self.references, self.tree, self.config)

        self.assertEqual(result.final_status, "FAIL")
        self.assertIsNone(result.best_fit)
        self.assertTrue(any(finding.code == "LOW_COVERAGE" for finding in result.qc.findings))

    def test_diffuse_weight_without_bootstrap_is_a_warning(self) -> None:
        configuration = PlacementConfig(
            preset="thorough",
            min_coverage=0.5,
            min_comparable_sites=2,
            unresolved_weight=1.0,
        )
        result = place_query(
            AlignedSequence("intermediate", tuple("AAAACCGG")),
            self.references,
            self.tree,
            configuration,
        )

        self.assertEqual(result.final_status, "WARN")
        self.assertTrue(any(finding.code == "DIFFUSE_DISTANCE_WEIGHT" for finding in result.qc.findings))

    def test_rendered_tree_is_parseable_and_contains_query(self) -> None:
        result = place_query(
            AlignedSequence("query_A", self.references[0].characters), self.references, self.tree, self.config
        )
        rendered = render_placed_tree(self.tree, [result])

        reparsed = Tree.from_newick(rendered)
        self.assertIn("query_A", reparsed.tip_index)
        self.assertIn("A", reparsed.tip_index)

    def test_fasta_adapter_requires_identical_column_counts(self) -> None:
        fixtures = Path("tests/fixtures")
        reference = read_fasta_alignment(fixtures / "reference_alignment.fasta")
        query = read_fasta_alignment(fixtures / "query_alignment.fasta")
        require_matching_alignment_widths(reference, query)
        self.assertEqual(reference.site_count, 4)
        self.assertEqual(query.sequences[0].characters, ("A", "C", "G", "T"))

    def test_pruning_records_the_original_tip_attachment(self) -> None:
        tree = Tree.from_newick("((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);")

        pruned, attachment = tree.without_tip("A")

        self.assertNotIn("A", pruned.tip_index)
        self.assertEqual(set(pruned.tip_index), {"B", "C", "D"})
        self.assertEqual(attachment.kind, "edge")
        assert attachment.edge_child_id is not None
        expected_edge = pruned.edge_by_child[attachment.edge_child_id]
        self.assertEqual(pruned.edge_distance_to_attachment(expected_edge, attachment), 0)

    def test_patristic_distance_is_zero_at_the_exact_edge_attachment(self) -> None:
        # Pruning A from ((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3) collapses the
        # (A,B) ancestor: B now hangs directly off the root at 0.4, and the
        # recorded attachment sits 0.3 from the root along that edge -- the
        # same position the (A,B) ancestor used to occupy.
        tree = Tree.from_newick("((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);")
        pruned, attachment = tree.without_tip("A")
        edge = pruned.edge_by_child[attachment.edge_child_id]

        distance = pruned.patristic_distance_to_attachment(edge, attachment.distal_length, attachment)

        self.assertAlmostEqual(distance, 0.0)

    def test_patristic_distance_on_the_same_edge_is_the_distal_length_difference(self) -> None:
        tree = Tree.from_newick("((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);")
        pruned, attachment = tree.without_tip("A")
        edge = pruned.edge_by_child[attachment.edge_child_id]

        distance = pruned.patristic_distance_to_attachment(edge, 0.0, attachment)

        self.assertAlmostEqual(distance, 0.3)

    def test_patristic_distance_routes_through_the_shared_root_for_a_distant_edge(self) -> None:
        # True attachment: 0.3 from the root along the root-B edge.
        # Placed point: 0.05 from Y along the Y-C edge, where Y=(C,D) hangs
        # 0.3 off the root. Path: attachment -> root (0.3) -> Y (0.3) -> point (0.05).
        tree = Tree.from_newick("((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);")
        pruned, attachment = tree.without_tip("A")
        c_edge = pruned.edge_by_child[pruned.tips[pruned.tip_index["C"]].node_id]

        distance = pruned.patristic_distance_to_attachment(c_edge, 0.05, attachment)

        self.assertAlmostEqual(distance, 0.65)

    def test_patristic_distance_to_a_node_attachment_from_its_own_edge(self) -> None:
        # A polytomy leaves the pruned tip's parent as a surviving node (X)
        # rather than collapsing it, exercising the "node" attachment kind.
        tree = Tree.from_newick("((A:0.1,B:0.1,E:0.1):0.3,(C:0.1,D:0.1):0.3);")
        pruned, attachment = tree.without_tip("A")
        self.assertEqual(attachment.kind, "node")
        b_edge = pruned.edge_by_child[pruned.tips[pruned.tip_index["B"]].node_id]

        at_the_node = pruned.patristic_distance_to_attachment(b_edge, 0.0, attachment)
        along_the_branch = pruned.patristic_distance_to_attachment(b_edge, 0.06, attachment)

        self.assertAlmostEqual(at_the_node, 0.0)
        self.assertAlmostEqual(along_the_branch, 0.06)

    def test_patristic_distance_to_a_node_attachment_from_a_distant_edge(self) -> None:
        # X (the attachment node) hangs 0.3 off the root; Y=(C,D) hangs 0.3
        # off the root on the other side. Path: point -> Y (0.1) -> root (0.3)
        # -> X (0.3).
        tree = Tree.from_newick("((A:0.1,B:0.1,E:0.1):0.3,(C:0.1,D:0.1):0.3);")
        pruned, attachment = tree.without_tip("A")
        c_edge = pruned.edge_by_child[pruned.tips[pruned.tip_index["C"]].node_id]

        distance = pruned.patristic_distance_to_attachment(c_edge, 0.1, attachment)

        self.assertAlmostEqual(distance, 0.7)

    def test_sequence_qc_prefers_reference_meeting_min_comparable_sites(self) -> None:
        from treetip_placer.sequence import evaluate_sequence_qc

        query = AlignedSequence("q", ("A", "C", "G", "T"))
        sparse_ref = AlignedSequence("sparse_zero_dist", ("A", None, None, None))
        complete_ref = AlignedSequence("complete_small_dist", ("A", "C", "G", "A"))
        qc = evaluate_sequence_qc(
            query,
            (sparse_ref, complete_ref),
            min_coverage=0.5,
            warn_coverage=0.9,
            min_comparable_sites=3,
            outlier_distance=0.5,
        )

        self.assertEqual(qc.status, "PASS")
        self.assertEqual(qc.nearest_reference, "complete_small_dist")
        self.assertEqual(qc.comparable_sites, 4)


if __name__ == "__main__":
    unittest.main()
