"""Grade leave-one-out placements against their recorded true attachment.

The grade is an unweighted topology diagnostic on a fixed backbone. It is not
phylogenetic support, not a branch-length error, and not evidence of external
generalization.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

from .models import PlacementResult
from .tree import PrunedTipAttachment, Tree


ACCURACY_RECORD_COLUMNS: tuple[str, ...] = (
    "source_sample_id",
    "outcome",
    "edge_hops_from_truth",
    "patristic_distance_from_truth",
    "true_attachment_kind",
    "placement_status",
    "placement_edge_id",
    "coverage",
    "comparable_sites",
    "nearest_reference_distance",
    "relative_fit_weight",
    "best_objective",
    "best_rmse",
    "bootstrap_support",
)

# A placement with either status must never receive an accuracy grade, matching
# the 0.6.3 policy that an unresolved placement cannot emit a discrete call.
UNGRADABLE_STATUSES = frozenset({"FAIL", "UNRESOLVED"})


@dataclass(frozen=True)
class AccuracyGrade:
    """One control's discrete accuracy outcome, its edge-hop distance, and
    its branch-length-weighted (patristic) distance from the true attachment.

    The hop count and the patristic distance can disagree on a poorly
    resolved tree: a placement many hops away along a chain of near-zero-
    length branches can still sit at essentially zero patristic distance
    from the truth. Neither is phylogenetic support.
    """

    outcome: str
    edge_hops_from_truth: int | None
    patristic_distance_from_truth: float | None


def grade_placement_accuracy(
    result: PlacementResult,
    pruned_tree: Tree,
    attachment: PrunedTipAttachment,
    near_hop_threshold: int,
) -> AccuracyGrade:
    """Compare a placed edge against the held-out tip's true attachment.

    Status is checked before best_fit: place_query populates best_fit even when
    it returns UNRESOLVED, so grading on a missing best_fit alone would score
    placements the project treats as having no defensible discrete call.
    """

    if result.final_status in UNGRADABLE_STATUSES or result.best_fit is None:
        return AccuracyGrade(outcome="UNRESOLVED", edge_hops_from_truth=None, patristic_distance_from_truth=None)
    selected = pruned_tree.edge_by_child.get(result.best_fit.child_id)
    if selected is None:
        return AccuracyGrade(outcome="UNRESOLVED", edge_hops_from_truth=None, patristic_distance_from_truth=None)
    hops = pruned_tree.edge_distance_to_attachment(selected, attachment)
    patristic_distance = pruned_tree.patristic_distance_to_attachment(
        selected, result.best_fit.distal_length, attachment
    )
    if hops == 0:
        outcome = "CORRECT"
    elif hops <= near_hop_threshold:
        outcome = "NEAR"
    else:
        outcome = "FAR"
    return AccuracyGrade(
        outcome=outcome, edge_hops_from_truth=hops, patristic_distance_from_truth=patristic_distance
    )


def build_accuracy_record(
    result: PlacementResult,
    grade: AccuracyGrade,
    attachment: PrunedTipAttachment,
    source_sample_id: str,
) -> dict[str, object]:
    """Build one auditable record row for a graded control."""

    best = result.best_fit
    return {
        "source_sample_id": source_sample_id,
        "outcome": grade.outcome,
        "edge_hops_from_truth": grade.edge_hops_from_truth,
        "patristic_distance_from_truth": grade.patristic_distance_from_truth,
        "true_attachment_kind": attachment.kind,
        "placement_status": result.final_status,
        "placement_edge_id": best.edge_id if best else None,
        "coverage": result.qc.coverage,
        "comparable_sites": result.qc.comparable_sites,
        "nearest_reference_distance": result.qc.nearest_distance,
        "relative_fit_weight": best.relative_fit_weight if best else None,
        "best_objective": best.objective if best else None,
        "best_rmse": best.rmse if best else None,
        "bootstrap_support": result.bootstrap_support,
    }


def summarize_placement_accuracy(
    records: Sequence[Mapping[str, object]],
) -> dict[str, int]:
    """Count controls per outcome using stable machine-readable keys."""

    counts: dict[str, int] = {"controls": len(records)}
    for record in records:
        key = f"outcome_{str(record['outcome']).lower()}"
        counts[key] = counts.get(key, 0) + 1
    return counts


def hop_distribution(records: Sequence[Mapping[str, object]]) -> dict[str, int]:
    """Count controls per observed hop distance.

    Retaining the raw distribution lets a near-hop threshold be re-evaluated
    from a completed run without repeating the placements.
    """

    counts: dict[str, int] = {}
    for record in records:
        hops = record["edge_hops_from_truth"]
        key = "unresolved" if hops is None else str(hops)
        counts[key] = counts.get(key, 0) + 1
    # Sort numerically, not lexically, so 10 hops follows 2 rather than preceding it.
    return dict(
        sorted(
            counts.items(),
            key=lambda item: (item[0] == "unresolved", int(item[0]) if item[0] != "unresolved" else 0),
        )
    )


def summarize_patristic_distance(records: Sequence[Mapping[str, object]]) -> dict[str, float | int | None]:
    """Summarize the branch-length-weighted distance over graded controls.

    Complements hop_distribution: on a poorly resolved tree the hop count can
    look bad while the underlying genetic distance recovered stays small, and
    this summary is what makes that visible without re-running the calibration.
    """

    distances = sorted(
        float(record["patristic_distance_from_truth"])
        for record in records
        if record["patristic_distance_from_truth"] is not None
    )
    if not distances:
        return {"graded_controls": 0, "mean": None, "median": None, "max": None}
    midpoint = len(distances) // 2
    if len(distances) % 2 == 0:
        median = (distances[midpoint - 1] + distances[midpoint]) / 2
    else:
        median = distances[midpoint]
    return {
        "graded_controls": len(distances),
        "mean": sum(distances) / len(distances),
        "median": median,
        "max": distances[-1],
    }
