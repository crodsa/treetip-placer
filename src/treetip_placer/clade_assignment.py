"""Conservative clade assignment from labelled tips on immutable reference trees."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping, Sequence

from .models import PlacementResult
from .panel_selection import PanelDecision
from .tree import Tree


@dataclass(frozen=True)
class CladeAssignment:
    """A discrete clade label inferred from a fixed-tree placement."""

    sample_id: str
    selected_panel: str | None
    phylogenetic_clade: str | None
    clade_assignment_status: str
    assignment_method: str | None
    nearest_reference: str | None
    supporting_reference_tips: tuple[str, ...]
    findings: tuple[str, ...]


def load_tip_clades(path: str | Path | None, tree: Tree, panel_name: str) -> dict[str, str]:
    """Load a strict TSV mapping of reference tip IDs to discrete clade labels.

    The expected columns are ``reference_tip_id`` and ``clade``.  Every declared
    reference ID must be present in the supplied fixed tree, preventing a silent
    mismatch between annotation data and the placement backbone.
    """

    if path is None:
        return {}
    metadata_path = Path(path)
    with metadata_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"reference_tip_id", "clade"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(
                f"The {panel_name} tip-clade TSV must contain reference_tip_id and clade columns"
            )
        labels: dict[str, str] = {}
        for row_number, row in enumerate(reader, start=2):
            tip_id = (row["reference_tip_id"] or "").strip()
            clade = (row["clade"] or "").strip()
            if not tip_id or not clade:
                raise ValueError(
                    f"The {panel_name} tip-clade TSV has an empty reference_tip_id or clade at row {row_number}"
                )
            if tip_id not in tree.tip_index:
                raise ValueError(
                    f"The {panel_name} tip-clade TSV declares unknown reference tip '{tip_id}'"
                )
            previous = labels.get(tip_id)
            if previous is not None and previous != clade:
                raise ValueError(
                    f"The {panel_name} tip-clade TSV assigns conflicting clades to reference tip '{tip_id}'"
                )
            labels[tip_id] = clade
    return labels


def assign_phylogenetic_clades(
    decisions: Sequence[PanelDecision],
    panels: Mapping[str, tuple[Tree, Sequence[PlacementResult]]],
    tip_clades: Mapping[str, Mapping[str, str]],
) -> list[CladeAssignment]:
    """Assign a clade only when a selected fixed-tree placement supports one label.

    A single clade among all labelled reference tips below the selected edge is a
    direct fixed-tree clade call. If the selected edge is mixed or unlabelled, a
    nearest-labelled-reference fallback is retained as ``WARN`` rather than being
    silently promoted to a direct edge-derived call. The function never invents a
    new topology or recomputes bootstrap support.
    """

    result_by_panel = {
        panel_name: {result.sample_id: result for result in results}
        for panel_name, (_, results) in panels.items()
    }
    assignments: list[CladeAssignment] = []
    for decision in decisions:
        if decision.selected_panel is None:
            assignments.append(
                CladeAssignment(
                    sample_id=decision.sample_id,
                    selected_panel=None,
                    phylogenetic_clade=None,
                    clade_assignment_status="UNRESOLVED",
                    assignment_method=None,
                    nearest_reference=None,
                    supporting_reference_tips=(),
                    findings=("NO_PANEL_SELECTED",),
                )
            )
            continue

        panel_name = decision.selected_panel
        tree, _ = panels[panel_name]
        labels = tip_clades.get(panel_name, {})
        result = result_by_panel[panel_name][decision.sample_id]
        nearest_reference = result.qc.nearest_reference
        if not labels:
            assignments.append(
                CladeAssignment(
                    sample_id=decision.sample_id,
                    selected_panel=panel_name,
                    phylogenetic_clade=None,
                    clade_assignment_status="UNRESOLVED",
                    assignment_method=None,
                    nearest_reference=nearest_reference,
                    supporting_reference_tips=(),
                    findings=("NO_TIP_CLADE_METADATA_FOR_SELECTED_PANEL",),
                )
            )
            continue
        if result.best_fit is None or result.final_status == "FAIL":
            assignments.append(
                CladeAssignment(
                    sample_id=decision.sample_id,
                    selected_panel=panel_name,
                    phylogenetic_clade=None,
                    clade_assignment_status="UNRESOLVED",
                    assignment_method=None,
                    nearest_reference=nearest_reference,
                    supporting_reference_tips=(),
                    findings=("NO_USABLE_SELECTED_PLACEMENT",),
                )
            )
            continue
        if result.final_status == "UNRESOLVED":
            assignments.append(
                CladeAssignment(
                    sample_id=decision.sample_id,
                    selected_panel=panel_name,
                    phylogenetic_clade=None,
                    clade_assignment_status="UNRESOLVED",
                    assignment_method=None,
                    nearest_reference=nearest_reference,
                    supporting_reference_tips=(),
                    findings=("SELECTED_PLACEMENT_UNRESOLVED",),
                )
            )
            continue

        edge = tree.edge_by_child[result.best_fit.child_id]
        subtree_tips = [
            tree.tips[index].name or ""
            for index in sorted(tree.edge_subtree_tip_indices(edge))
        ]
        labelled_tips = tuple(tip_id for tip_id in subtree_tips if tip_id in labels)
        subtree_clades = {labels[tip_id] for tip_id in labelled_tips}
        if len(subtree_clades) == 1:
            status = "PASS" if result.final_status == "PASS" else "WARN"
            findings = () if status == "PASS" else ("SELECTED_PLACEMENT_NOT_PASS",)
            assignments.append(
                CladeAssignment(
                    sample_id=decision.sample_id,
                    selected_panel=panel_name,
                    phylogenetic_clade=next(iter(subtree_clades)),
                    clade_assignment_status=status,
                    assignment_method="selected_edge_descendant_tip_labels",
                    nearest_reference=nearest_reference,
                    supporting_reference_tips=labelled_tips[:20],
                    findings=findings,
                )
            )
            continue
        if len(subtree_clades) > 1:
            if nearest_reference in labels:
                assignments.append(
                    CladeAssignment(
                        sample_id=decision.sample_id,
                        selected_panel=panel_name,
                        phylogenetic_clade=labels[nearest_reference],
                        clade_assignment_status="WARN",
                        assignment_method="nearest_labelled_reference_mixed_edge_fallback",
                        nearest_reference=nearest_reference,
                        supporting_reference_tips=(nearest_reference,),
                        findings=(
                            "MIXED_CLADE_LABELS_ON_SELECTED_EDGE",
                            "NEAREST_LABELLED_REFERENCE_FALLBACK",
                        ),
                    )
                )
                continue
            assignments.append(
                CladeAssignment(
                    sample_id=decision.sample_id,
                    selected_panel=panel_name,
                    phylogenetic_clade=None,
                    clade_assignment_status="UNRESOLVED",
                    assignment_method="selected_edge_descendant_tip_labels",
                    nearest_reference=nearest_reference,
                    supporting_reference_tips=labelled_tips[:20],
                    findings=("MIXED_CLADE_LABELS_ON_SELECTED_EDGE",),
                )
            )
            continue
        if nearest_reference in labels:
            assignments.append(
                CladeAssignment(
                    sample_id=decision.sample_id,
                    selected_panel=panel_name,
                    phylogenetic_clade=labels[nearest_reference],
                    clade_assignment_status="WARN",
                    assignment_method="nearest_labelled_reference_fallback",
                    nearest_reference=nearest_reference,
                    supporting_reference_tips=(nearest_reference,),
                    findings=("SELECTED_EDGE_HAS_NO_LABELLED_DESCENDANT_TIPS",),
                )
            )
            continue
        assignments.append(
            CladeAssignment(
                sample_id=decision.sample_id,
                selected_panel=panel_name,
                phylogenetic_clade=None,
                clade_assignment_status="UNRESOLVED",
                assignment_method=None,
                nearest_reference=nearest_reference,
                supporting_reference_tips=(),
                findings=("NO_LABELLED_REFERENCE_FOR_SELECTED_PLACEMENT",),
            )
        )
    return assignments


def write_clade_assignment_outputs(
    output_directory: str | Path,
    assignments: Sequence[CladeAssignment],
    metadata: Mapping[str, object],
) -> None:
    """Write authoritative clade-assignment TSV and JSON sidecars."""

    root = Path(output_directory)
    fields = [
        "sample_id",
        "selected_panel",
        "phylogenetic_clade",
        "clade_assignment_status",
        "assignment_method",
        "nearest_reference",
        "supporting_reference_tips",
        "findings",
    ]
    with (root / "clade_assignment.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for assignment in assignments:
            row = asdict(assignment)
            row["supporting_reference_tips"] = " | ".join(assignment.supporting_reference_tips)
            row["findings"] = " | ".join(assignment.findings)
            writer.writerow(row)
    (root / "clade_assignment.json").write_text(
        json.dumps({"metadata": dict(metadata), "assignments": [asdict(item) for item in assignments]}, indent=2, sort_keys=True)
        + "\n",
        encoding="utf-8",
    )
