"""Interoperable run outputs for trees, tables, and machine-readable provenance."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path
from typing import Mapping, Sequence

from .models import EdgeFit, PlacementResult, RunMetadata
from .tree import Node, Tree, format_length, quote_label


def write_run_outputs(
    output_directory: str | Path,
    tree: Tree,
    results: Sequence[PlacementResult],
    metadata: RunMetadata,
    events: Sequence[Mapping[str, object]],
) -> None:
    """Write all documented run artifacts to an existing empty directory."""

    root = Path(output_directory)
    trees = root / "trees"
    placements = root / "placements"
    qc = root / "qc"
    for directory in (trees, placements, qc):
        directory.mkdir(parents=True, exist_ok=True)

    placed_newick = render_placed_tree(tree, results)
    (trees / "placed_tree.nwk").write_text(placed_newick, encoding="utf-8")
    (trees / "reference_tree.normalized.nwk").write_text(tree.newick(), encoding="utf-8")
    (trees / "placed_tree.nexus").write_text(
        f"#NEXUS\n\nBEGIN TREES;\n  TREE placed_tree = {placed_newick.strip()}\nEND;\n",
        encoding="utf-8",
    )

    write_placements_tsv(placements / "placements.tsv", results)
    write_alternatives_tsv(placements / "alternative_edges.tsv", results)
    write_qc_tsv(qc / "query_qc.tsv", results)
    write_json(placements / "placements.json", [result_to_dict(result) for result in results])
    write_json(placements / "placements.jplace", jplace_document(tree, results))
    write_json(root / "run_manifest.json", asdict(metadata))
    with (root / "events.jsonl").open("w", encoding="utf-8", newline="\n") as handle:
        for event in events:
            handle.write(json.dumps(event, sort_keys=True) + "\n")


def render_placed_tree(tree: Tree, results: Sequence[PlacementResult]) -> str:
    """Render independent fixed-backbone placements as a single Newick tree."""

    attachments: dict[int, list[tuple[str, EdgeFit]]] = {}
    for result in results:
        if result.best_fit is not None and result.final_status != "FAIL":
            attachments.setdefault(result.best_fit.child_id, []).append((result.sample_id, result.best_fit))

    def render_node(node: Node) -> str:
        label = quote_label(node.name) if node.name else ""
        if node.is_leaf:
            return label
        rendered_children = [render_edge(child) for child in node.children]
        return "(" + ",".join(rendered_children) + ")" + label

    def render_edge(child: Node) -> str:
        length = child.length or 0.0
        downstream = length
        subtree = render_node(child)
        by_position: dict[float, list[tuple[str, EdgeFit]]] = {}
        for sample_id, fit in attachments.get(child.node_id, []):
            position = min(max(fit.distal_length, 0.0), length)
            by_position.setdefault(round(position, 12), []).append((sample_id, fit))
        for rounded_position in sorted(by_position, reverse=True):
            position = rounded_position
            leaves = [
                f"{quote_label(sample_id)}:{format_length(fit.pendant_length)}"
                for sample_id, fit in sorted(by_position[rounded_position])
            ]
            subtree = "(" + ",".join([f"{subtree}:{format_length(downstream - position)}", *leaves]) + ")"
            downstream = position
        return f"{subtree}:{format_length(downstream)}"

    return render_node(tree.root) + ";\n"


def write_placements_tsv(path: Path, results: Sequence[PlacementResult]) -> None:
    """Write one authoritative summary row per query."""

    fields = [
        "sample_id",
        "final_status",
        "qc_status",
        "coverage",
        "observed_sites",
        "comparable_sites",
        "nearest_reference",
        "nearest_distance",
        "edge_id",
        "parent_node_id",
        "child_node_id",
        "distal_length",
        "pendant_length",
        "distance_objective",
        "rmse",
        "relative_fit_weight",
        "bootstrap_edge_support",
        "search_scope",
        "findings",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for result in results:
            fit = result.best_fit
            writer.writerow(
                {
                    "sample_id": result.sample_id,
                    "final_status": result.final_status,
                    "qc_status": result.qc.status,
                    "coverage": f"{result.qc.coverage:.6f}",
                    "observed_sites": result.qc.observed_sites,
                    "comparable_sites": result.qc.comparable_sites,
                    "nearest_reference": result.qc.nearest_reference or "",
                    "nearest_distance": value_or_empty(result.qc.nearest_distance),
                    "edge_id": fit.edge_id if fit else "",
                    "parent_node_id": fit.parent_id if fit else "",
                    "child_node_id": fit.child_id if fit else "",
                    "distal_length": value_or_empty(fit.distal_length if fit else None),
                    "pendant_length": value_or_empty(fit.pendant_length if fit else None),
                    "distance_objective": value_or_empty(fit.objective if fit else None),
                    "rmse": value_or_empty(fit.rmse if fit else None),
                    "relative_fit_weight": value_or_empty(fit.relative_fit_weight if fit else None),
                    "bootstrap_edge_support": value_or_empty(result.bootstrap_support),
                    "search_scope": result.search_scope or "",
                    "findings": " | ".join(f"{finding.code}: {finding.message}" for finding in result.qc.findings),
                }
            )


def write_alternatives_tsv(path: Path, results: Sequence[PlacementResult]) -> None:
    """Write ranked alternative edge fits for uncertainty inspection."""

    fields = [
        "sample_id",
        "rank",
        "edge_id",
        "parent_node_id",
        "child_node_id",
        "distal_length",
        "pendant_length",
        "distance_objective",
        "rmse",
        "relative_fit_weight",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for result in results:
            for rank, fit in enumerate(result.alternatives, start=1):
                writer.writerow(
                    {
                        "sample_id": result.sample_id,
                        "rank": rank,
                        "edge_id": fit.edge_id,
                        "parent_node_id": fit.parent_id,
                        "child_node_id": fit.child_id,
                        "distal_length": value_or_empty(fit.distal_length),
                        "pendant_length": value_or_empty(fit.pendant_length),
                        "distance_objective": value_or_empty(fit.objective),
                        "rmse": value_or_empty(fit.rmse),
                        "relative_fit_weight": value_or_empty(fit.relative_fit_weight),
                    }
                )


def write_qc_tsv(path: Path, results: Sequence[PlacementResult]) -> None:
    """Write one row per QC finding; queries without findings receive an empty row."""

    fields = ["sample_id", "preplacement_status", "code", "severity", "message"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for result in results:
            findings = result.qc.findings or [None]
            for finding in findings:
                writer.writerow(
                    {
                        "sample_id": result.sample_id,
                        "preplacement_status": result.qc.status,
                        "code": finding.code if finding else "",
                        "severity": finding.severity if finding else "",
                        "message": finding.message if finding else "",
                    }
                )


def jplace_document(tree: Tree, results: Sequence[PlacementResult]) -> dict[str, object]:
    """Return an explicitly distance-based, experimental JPlace-shaped document."""

    placements: list[dict[str, object]] = []
    for result in results:
        if not result.alternatives:
            continue
        placements.append(
            {
                "p": [
                    [
                        fit.edge_id,
                        fit.objective,
                        fit.relative_fit_weight,
                        fit.distal_length,
                        fit.pendant_length,
                    ]
                    for fit in result.alternatives
                ],
                "n": [result.sample_id],
            }
        )
    return {
        "tree": render_jplace_tree(tree),
        "fields": [
            "edge_num",
            "distance_objective",
            "relative_fit_weight",
            "distal_length",
            "pendant_length",
        ],
        "placements": placements,
        "metadata": {
            "placement_method": "weighted_least_squares_categorical_distance",
            "warning": "This is not a maximum-likelihood placement result.",
        },
        "version": 3,
    }


def render_jplace_tree(tree: Tree) -> str:
    """Render a reference tree with stable edge numbers in JPlace notation."""

    def render_node(node: Node) -> str:
        label = quote_label(node.name) if node.name else ""
        if node.is_leaf:
            return label
        children = []
        for child in node.children:
            edge = tree.edge_by_child[child.node_id]
            children.append(f"{render_node(child)}:{format_length(edge.length)}{{{edge.edge_id}}}")
        return "(" + ",".join(children) + ")" + label

    return render_node(tree.root) + ";"


def result_to_dict(result: PlacementResult) -> dict[str, object]:
    """Serialize result dataclasses without leaking Python implementation details."""

    return {
        "sample_id": result.sample_id,
        "final_status": result.final_status,
        "search_scope": result.search_scope,
        "note": result.note,
        "qc": {
            "status": result.qc.status,
            "coverage": result.qc.coverage,
            "observed_sites": result.qc.observed_sites,
            "comparable_sites": result.qc.comparable_sites,
            "nearest_reference": result.qc.nearest_reference,
            "nearest_distance": result.qc.nearest_distance,
            "findings": [asdict(finding) for finding in result.qc.findings],
        },
        "best_fit": asdict(result.best_fit) if result.best_fit else None,
        "alternatives": [asdict(fit) for fit in result.alternatives],
        "bootstrap_support": result.bootstrap_support,
    }


def write_json(path: Path, contents: object) -> None:
    """Write stable, UTF-8 JSON for downstream tooling."""

    path.write_text(json.dumps(contents, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def value_or_empty(value: float | None) -> str:
    """Format optional numerical values for TSV without locale changes."""

    return "" if value is None else format(value, ".12g")
