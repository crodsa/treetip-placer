"""Tests for leave-one-out placement-accuracy grading."""

from __future__ import annotations

import unittest

from treetip_placer.models import AlignedSequence, PlacementConfig, PlacementResult
from treetip_placer.placement import place_query
from treetip_placer.placement_accuracy import (
    grade_placement_accuracy,
    hop_distribution,
    summarize_placement_accuracy,
)
from treetip_placer.tree import Tree


CONCORDANT_SEQUENCES = {
    "A": "AAAAAA",
    "B": "AAAAAT",
    "C": "CCCCCC",
    "D": "CCCCCG",
}

# 'A' carries a C-clade sequence while the topology keeps it beside 'B',
# so its leave-one-out placement lands two edges from its true attachment.
DISCORDANT_SEQUENCES = {
    "A": "CCCCCC",
    "B": "AAAAAA",
    "C": "CCCCCC",
    "D": "CCCCCG",
}


def grade_holdout(newick: str, sequences: dict[str, object], held_out: str, near_hop_threshold: int):
    """Run one leave-one-out control and grade it.

    A sequence may be given as a string of characters or, to model absent data,
    as a tuple already containing None. Reading FASTA maps '-' and 'N' to None
    (`alignment.SEQUENCE_MISSING`), but an AlignedSequence built directly in a
    test keeps whatever characters it is handed — so a literal '-' string here
    would be scored as real data and grade WARN rather than failing QC.
    """

    tree = Tree.from_newick(newick)
    by_id = {
        name: AlignedSequence(name, tuple(text) if isinstance(text, str) else text)
        for name, text in sequences.items()
    }
    pruned, attachment = tree.without_tip(held_out)
    references = tuple(by_id[tip.name or ""] for tip in pruned.tips)
    config = PlacementConfig(
        preset="thorough", min_coverage=0.5, min_comparable_sites=2, unresolved_weight=0.0
    )
    result = place_query(by_id[held_out], references, pruned, config)
    return grade_placement_accuracy(result, pruned, attachment, near_hop_threshold), attachment


class PlacementAccuracyTests(unittest.TestCase):
    def test_recovered_attachment_grades_correct_with_zero_hops(self) -> None:
        grade, attachment = grade_holdout(
            "((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);", CONCORDANT_SEQUENCES, "A", 2
        )

        self.assertEqual(grade.outcome, "CORRECT")
        self.assertEqual(grade.edge_hops_from_truth, 0)
        self.assertEqual(attachment.kind, "edge")

    def test_node_kind_attachment_is_graded(self) -> None:
        # A polytomy leaves a node with two surviving children after pruning,
        # exercising the 'node' branch of edge_distance_to_attachment.
        grade, attachment = grade_holdout(
            "((A:0.1,B:0.1,E:0.1):0.3,(C:0.1,D:0.1):0.3);",
            {**CONCORDANT_SEQUENCES, "E": "AAAATT"},
            "A",
            2,
        )

        self.assertEqual(attachment.kind, "node")
        self.assertEqual(grade.outcome, "CORRECT")
        self.assertEqual(grade.edge_hops_from_truth, 0)

    def test_displaced_placement_grades_near_within_the_threshold(self) -> None:
        grade, _ = grade_holdout(
            "((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);", DISCORDANT_SEQUENCES, "A", 2
        )

        self.assertEqual(grade.edge_hops_from_truth, 2)
        self.assertEqual(grade.outcome, "NEAR")

    def test_the_same_placement_grades_far_under_a_stricter_threshold(self) -> None:
        grade, _ = grade_holdout(
            "((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);", DISCORDANT_SEQUENCES, "A", 1
        )

        self.assertEqual(grade.edge_hops_from_truth, 2)
        self.assertEqual(grade.outcome, "FAR")

    def test_failed_qc_grades_unresolved_without_hops(self) -> None:
        grade, _ = grade_holdout(
            "((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);",
            {**CONCORDANT_SEQUENCES, "A": (None,) * 6},
            "A",
            2,
        )

        self.assertEqual(grade.outcome, "UNRESOLVED")
        self.assertIsNone(grade.edge_hops_from_truth)

    def test_unresolved_status_is_not_graded_even_with_a_best_fit(self) -> None:
        """The regression guarding the status-before-best_fit ordering.

        place_query populates best_fit even when it returns UNRESOLVED, so an
        implementation testing only `best_fit is None` passes every other test
        here and fails exactly this one.
        """

        tree = Tree.from_newick("((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);")
        pruned, attachment = tree.without_tip("A")
        placed = place_query(
            AlignedSequence("A", tuple("AAAAAA")),
            tuple(AlignedSequence(tip.name or "", tuple(CONCORDANT_SEQUENCES[tip.name or ""])) for tip in pruned.tips),
            pruned,
            PlacementConfig(preset="thorough", min_coverage=0.5, min_comparable_sites=2, unresolved_weight=0.0),
        )
        self.assertIsNotNone(placed.best_fit)
        unresolved = PlacementResult(
            sample_id=placed.sample_id,
            qc=placed.qc,
            final_status="UNRESOLVED",
            best_fit=placed.best_fit,
            bootstrap_support=0.4,
        )

        grade = grade_placement_accuracy(unresolved, pruned, attachment, 2)

        self.assertEqual(grade.outcome, "UNRESOLVED")
        self.assertIsNone(grade.edge_hops_from_truth)

    def test_summaries_count_outcomes_and_hops(self) -> None:
        records = [
            {"outcome": "CORRECT", "edge_hops_from_truth": 0},
            {"outcome": "CORRECT", "edge_hops_from_truth": 0},
            {"outcome": "NEAR", "edge_hops_from_truth": 2},
            {"outcome": "UNRESOLVED", "edge_hops_from_truth": None},
        ]

        self.assertEqual(
            summarize_placement_accuracy(records),
            {"controls": 4, "outcome_correct": 2, "outcome_near": 1, "outcome_unresolved": 1},
        )
        self.assertEqual(hop_distribution(records), {"0": 2, "2": 1, "unresolved": 1})


class PlacementAccuracyCliTests(unittest.TestCase):
    def test_command_writes_records_and_summary(self) -> None:
        import json
        import shutil
        from pathlib import Path

        from treetip_placer.cli import main

        root = Path.cwd() / ".treetip-placement-accuracy-test"
        if root.exists():
            shutil.rmtree(root)
        root.mkdir()
        try:
            tree = root / "tree.nwk"
            alignment = root / "reference.fasta"
            output = root / "output"
            tree.write_text("((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);\n", encoding="utf-8")
            alignment.write_text(
                ">A\nAAAAAA\n>B\nAAAAAT\n>C\nCCCCCC\n>D\nCCCCCG\n", encoding="utf-8"
            )

            exit_code = main(
                [
                    "calibrate-placement-accuracy",
                    "--tree", str(tree),
                    "--reference-alignment", str(alignment),
                    "--output", str(output),
                    "--preset", "thorough",
                    "--min-comparable-sites", "2",
                    "--unresolved-weight", "0",
                ]
            )

            self.assertEqual(exit_code, 0)
            records = (output / "placement_accuracy_records.tsv").read_text(encoding="utf-8")
            self.assertIn("edge_hops_from_truth", records)
            self.assertIn("true_attachment_kind", records)
            summary = json.loads((output / "placement_accuracy_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["summary"]["controls"], 4)
            self.assertEqual(summary["summary"]["outcome_correct"], 4)
            self.assertEqual(summary["near_hop_threshold"], 2)
            self.assertEqual(sum(summary["hop_distribution"].values()), 4)
            self.assertIn("not phylogenetic support", summary["warning"])
        finally:
            if root.exists():
                shutil.rmtree(root)


if __name__ == "__main__":
    unittest.main()
