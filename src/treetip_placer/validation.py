"""Leave-one-out validation for categorical MLST fixed-backbone placement."""

from __future__ import annotations

import csv
import json
import random
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, Mapping, Sequence

from .models import PlacementConfig, PlacementResult, Profile
from .placement import place_query
from .tree import PrunedTipAttachment, Tree


@dataclass(frozen=True)
class ValidationRecord:
    """One held-out profile placement and its known backbone attachment."""

    held_out_id: str
    placement_status: str
    qc_status: str
    true_attachment: PrunedTipAttachment
    selected_edge_id: int | None
    selected_parent_node_id: int | None
    selected_child_node_id: int | None
    topological_distance_edges: int | None
    recovered_true_attachment: bool | None
    nearest_reference: str | None
    nearest_distance: float | None
    relative_fit_weight: float | None
    bootstrap_edge_support: float | None
    findings: tuple[str, ...]


_WORKER_TREE_NEWICK: str | None = None
_WORKER_PROFILES: tuple[Profile, ...] = ()
_WORKER_CONFIG: PlacementConfig | None = None
_WORKER_MAX_CACHE_MEGABYTES = 256


def initialize_validation_worker(
    tree_newick: str,
    profiles_in_tip_order: tuple[Profile, ...],
    config: PlacementConfig,
    max_cache_megabytes: int,
) -> None:
    """Store immutable validation inputs once per process worker."""

    global _WORKER_TREE_NEWICK, _WORKER_PROFILES, _WORKER_CONFIG, _WORKER_MAX_CACHE_MEGABYTES
    _WORKER_TREE_NEWICK = tree_newick
    _WORKER_PROFILES = profiles_in_tip_order
    _WORKER_CONFIG = config
    _WORKER_MAX_CACHE_MEGABYTES = max_cache_megabytes


def validate_holdout_worker(held_out_id: str) -> ValidationRecord:
    """Run one validation task in an initialized process worker."""

    if _WORKER_TREE_NEWICK is None or _WORKER_CONFIG is None:
        raise RuntimeError("Validation worker was not initialized")
    tree = Tree.from_newick(_WORKER_TREE_NEWICK, _WORKER_MAX_CACHE_MEGABYTES)
    return validate_holdout(held_out_id, tree, _WORKER_PROFILES, _WORKER_CONFIG)


def select_holdout_ids(
    profiles_in_tip_order: Sequence[Profile], maximum: int, seed: int
) -> tuple[str, ...]:
    """Choose a deterministic, tree-order-preserving subset of profiles."""

    if maximum == 0 or maximum >= len(profiles_in_tip_order):
        return tuple(profile.sample_id for profile in profiles_in_tip_order)
    chosen = {
        profile.sample_id
        for profile in random.Random(seed).sample(list(profiles_in_tip_order), maximum)
    }
    return tuple(profile.sample_id for profile in profiles_in_tip_order if profile.sample_id in chosen)


def validate_holdout(
    held_out_id: str,
    full_tree: Tree,
    profiles_in_tip_order: Sequence[Profile],
    config: PlacementConfig,
) -> ValidationRecord:
    """Prune one known sample, place it, and compare to its original attachment."""

    profiles_by_id = {profile.sample_id: profile for profile in profiles_in_tip_order}
    query = profiles_by_id.get(held_out_id)
    if query is None:
        raise ValueError(f"Held-out profile '{held_out_id}' is absent from the reference panel")
    pruned_tree, true_attachment = full_tree.without_tip(held_out_id)
    pruned_references = tuple(
        profiles_by_id[tip.name or ""] for tip in pruned_tree.tips
    )
    result = place_query(query, pruned_references, pruned_tree, config)
    return validation_record(held_out_id, result, pruned_tree, true_attachment)


def validation_record(
    held_out_id: str,
    result: PlacementResult,
    pruned_tree: Tree,
    true_attachment: PrunedTipAttachment,
) -> ValidationRecord:
    """Convert an internal placement result into a stable validation record."""

    fit = result.best_fit
    topological_distance = None
    recovered = None
    if fit is not None:
        selected_edge = pruned_tree.edge_by_child[fit.child_id]
        topological_distance = pruned_tree.edge_distance_to_attachment(
            selected_edge, true_attachment
        )
        recovered = topological_distance == 0
    return ValidationRecord(
        held_out_id=held_out_id,
        placement_status=result.final_status,
        qc_status=result.qc.status,
        true_attachment=true_attachment,
        selected_edge_id=fit.edge_id if fit else None,
        selected_parent_node_id=fit.parent_id if fit else None,
        selected_child_node_id=fit.child_id if fit else None,
        topological_distance_edges=topological_distance,
        recovered_true_attachment=recovered,
        nearest_reference=result.qc.nearest_reference,
        nearest_distance=result.qc.nearest_distance,
        relative_fit_weight=fit.relative_fit_weight if fit else None,
        bootstrap_edge_support=result.bootstrap_support,
        findings=tuple(f"{finding.code}: {finding.message}" for finding in result.qc.findings),
    )


def run_serial_validation(
    held_out_ids: Sequence[str],
    tree: Tree,
    profiles_in_tip_order: Sequence[Profile],
    config: PlacementConfig,
    progress: Callable[[int, int], None],
) -> list[ValidationRecord]:
    """Validate in deterministic order without process workers."""

    records: list[ValidationRecord] = []
    for completed, held_out_id in enumerate(held_out_ids, start=1):
        records.append(validate_holdout(held_out_id, tree, profiles_in_tip_order, config))
        progress(completed, len(held_out_ids))
    return records


def run_parallel_validation(
    held_out_ids: Sequence[str],
    tree_newick: str,
    profiles_in_tip_order: tuple[Profile, ...],
    config: PlacementConfig,
    max_cache_megabytes: int,
    workers: int,
    progress: Callable[[int, int], None],
) -> list[ValidationRecord]:
    """Validate independent held-outs in parallel while preserving output order."""

    records: list[ValidationRecord | None] = [None] * len(held_out_ids)
    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=initialize_validation_worker,
        initargs=(tree_newick, profiles_in_tip_order, config, max_cache_megabytes),
    ) as executor:
        futures = {
            executor.submit(validate_holdout_worker, held_out_id): index
            for index, held_out_id in enumerate(held_out_ids)
        }
        completed = 0
        for future in as_completed(futures):
            records[futures[future]] = future.result()
            completed += 1
            progress(completed, len(held_out_ids))
    return [record for record in records if record is not None]


def summarize_validation(records: Sequence[ValidationRecord]) -> dict[str, object]:
    """Calculate transparent recovery and placement-status summaries."""

    status_counts: dict[str, int] = {}
    for record in records:
        status_counts[record.placement_status] = status_counts.get(record.placement_status, 0) + 1
    placed = [record for record in records if record.topological_distance_edges is not None]
    recovered = [record for record in placed if record.recovered_true_attachment]
    distances = [record.topological_distance_edges for record in placed]
    return {
        "holdouts_requested": len(records),
        "placements_with_selected_edge": len(placed),
        "recovered_true_attachment": len(recovered),
        "recovery_rate": len(recovered) / len(placed) if placed else None,
        "mean_topological_distance_edges": sum(distances) / len(distances) if distances else None,
        "placement_status_counts": dict(sorted(status_counts.items())),
        "metric_definition": (
            "A recovery has graph-edge distance zero: the selected edge touches a true "
            "node attachment or exactly matches a true attachment within a pruned edge."
        ),
    }


def write_validation_outputs(
    output_directory: str | Path,
    records: Sequence[ValidationRecord],
    summary: Mapping[str, object],
    manifest: Mapping[str, object],
    events: Sequence[Mapping[str, object]],
) -> None:
    """Write portable leave-one-out tables, JSON, and structured progress events."""

    root = Path(output_directory)
    write_validation_tsv(root / "leave_one_out.tsv", records)
    write_json(root / "leave_one_out.json", [validation_record_to_dict(record) for record in records])
    write_json(root / "validation_summary.json", dict(summary))
    write_json(root / "validation_manifest.json", dict(manifest))
    with (root / "events.jsonl").open("w", encoding="utf-8", newline="\n") as handle:
        for event in events:
            handle.write(json.dumps(event, sort_keys=True) + "\n")


def write_validation_tsv(path: Path, records: Sequence[ValidationRecord]) -> None:
    """Write one directly importable row per held-out reference profile."""

    fields = [
        "held_out_id",
        "placement_status",
        "qc_status",
        "true_attachment_kind",
        "true_attachment_node_id",
        "true_attachment_edge_parent_id",
        "true_attachment_edge_child_id",
        "true_attachment_distal_length",
        "selected_edge_id",
        "selected_parent_node_id",
        "selected_child_node_id",
        "topological_distance_edges",
        "recovered_true_attachment",
        "nearest_reference",
        "nearest_distance",
        "relative_fit_weight",
        "bootstrap_edge_support",
        "findings",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record in records:
            attachment = record.true_attachment
            writer.writerow(
                {
                    "held_out_id": record.held_out_id,
                    "placement_status": record.placement_status,
                    "qc_status": record.qc_status,
                    "true_attachment_kind": attachment.kind,
                    "true_attachment_node_id": value_or_empty(attachment.node_id),
                    "true_attachment_edge_parent_id": value_or_empty(attachment.edge_parent_id),
                    "true_attachment_edge_child_id": value_or_empty(attachment.edge_child_id),
                    "true_attachment_distal_length": value_or_empty(attachment.distal_length),
                    "selected_edge_id": value_or_empty(record.selected_edge_id),
                    "selected_parent_node_id": value_or_empty(record.selected_parent_node_id),
                    "selected_child_node_id": value_or_empty(record.selected_child_node_id),
                    "topological_distance_edges": value_or_empty(record.topological_distance_edges),
                    "recovered_true_attachment": value_or_empty(record.recovered_true_attachment),
                    "nearest_reference": record.nearest_reference or "",
                    "nearest_distance": value_or_empty(record.nearest_distance),
                    "relative_fit_weight": value_or_empty(record.relative_fit_weight),
                    "bootstrap_edge_support": value_or_empty(record.bootstrap_edge_support),
                    "findings": " | ".join(record.findings),
                }
            )


def validation_record_to_dict(record: ValidationRecord) -> dict[str, object]:
    """Return a JSON-compatible representation without Python model details."""

    return {
        **asdict(record),
        "findings": list(record.findings),
    }


def write_json(path: Path, contents: object) -> None:
    """Write stable, UTF-8 JSON for downstream processing."""

    path.write_text(json.dumps(contents, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def value_or_empty(value: object | None) -> str:
    """Format optional scalar table values without locale changes."""

    if value is None:
        return ""
    if isinstance(value, float):
        return format(value, ".12g")
    return str(value).lower() if isinstance(value, bool) else str(value)
