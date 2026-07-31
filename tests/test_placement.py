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
        result = place_query(AlignedSequence("query_A", self.references[0].characters), self.references, self.tree, self.config)

        self.assertEqual(result.final_status, "PASS")
        self.assertIsNotNone(result.best_fit)
        assert result.best_fit is not None
        self.assertEqual(self.tree.nodes[result.best_fit.child_id].name, "A")
        self.assertAlmostEqual(result.best_fit.distal_length, 0.5)
        self.assertAlmostEqual(result.best_fit.pendant_length, 0.0)
        self.assertAlmostEqual(result.best_fit.objective, 0.0)

    def test_missing_sequence_fails_preplacement_qc(self) -> None:
        result = place_query(
            AlignedSequence("empty", (None,) * 8), self.references, self.tree, self.config
        )

        self.assertEqual(result.final_status, "FAIL")
        self.assertIsNone(result.best_fit)
        self.assertTrue(any(finding.code == "LOW_COVERAGE" for finding in result.qc.findings))

    def test_rendered_tree_is_parseable_and_contains_query(self) -> None:
        result = place_query(AlignedSequence("query_A", self.references[0].characters), self.references, self.tree, self.config)
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


if __name__ == "__main__":
    unittest.main()
