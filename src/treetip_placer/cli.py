"""Command-line interface for the TreeTip Placer MVP."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

from . import __version__
from .alignment import (
    read_fasta_alignment,
    read_fasta_records,
    require_matching_alignment_widths,
    write_fasta_records,
)
from .calibration_reports import write_calibration_report
from .clade_assignment import (
    assign_phylogenetic_clades,
    load_tip_clades,
    write_clade_assignment_outputs,
)
from .models import AlignedSequence, PlacementConfig, PlacementResult, RunMetadata
from .outputs import render_placed_tree, write_run_outputs
from .panel_selection import PanelDecision, choose_reference_panel, write_panel_selection_outputs
from .placement import place_query
from .placement_accuracy import (
    ACCURACY_RECORD_COLUMNS,
    build_accuracy_record,
    grade_placement_accuracy,
    hop_distribution,
    summarize_placement_accuracy,
)
from .reference_package import (
    PanelReference,
    PairedReferencePackage,
    create_paired_reference_package,
    load_paired_reference_package,
    validate_paired_reference_package,
)
from .reports import write_auto_report
from .sequence import evaluate_sequence_qc
from .tree import Tree


_WORKER_TREE: Tree | None = None
_WORKER_REFERENCES: tuple[AlignedSequence, ...] = ()
_WORKER_CONFIG: PlacementConfig | None = None


@dataclass(frozen=True)
class PanelRun:
    """Placement results and immutable tree for one named reference panel."""

    name: str
    tree: Tree
    results: tuple[PlacementResult, ...]


def worker_initialize(
    newick: str,
    references: tuple[AlignedSequence, ...],
    config: PlacementConfig,
    max_cache_megabytes: int,
) -> None:
    """Create isolated worker state so caches do not need to be pickled per query."""

    global _WORKER_TREE, _WORKER_REFERENCES, _WORKER_CONFIG
    _WORKER_TREE = Tree.from_newick(newick, max_cache_megabytes)
    _WORKER_REFERENCES = references
    _WORKER_CONFIG = config


def worker_place(query: AlignedSequence):
    """Place a single query in a process worker."""

    if _WORKER_TREE is None or _WORKER_CONFIG is None:
        raise RuntimeError("Placement worker was not initialized")
    return place_query(
        query,
        _WORKER_REFERENCES,
        _WORKER_TREE,
        _WORKER_CONFIG,
    )


class EventLogger:
    """Human-readable stderr progress plus structured run events."""

    def __init__(self, verbosity: int) -> None:
        self.verbosity = verbosity
        self.started_at = time.perf_counter()
        self.events: list[dict[str, object]] = []

    def event(self, level: int, stage: str, message: str, **data: object) -> None:
        elapsed_seconds = time.perf_counter() - self.started_at
        event = {
            "elapsed_seconds": round(elapsed_seconds, 6),
            "level": level,
            "stage": stage,
            "message": message,
            **data,
        }
        self.events.append(event)
        if self.verbosity >= level:
            suffix = "" if not data else " " + json.dumps(data, sort_keys=True)
            print(f"[{elapsed_seconds:8.2f}s] {stage}: {message}{suffix}", file=sys.stderr)


def build_parser() -> argparse.ArgumentParser:
    """Create the English-only command-line schema."""

    parser = argparse.ArgumentParser(
        prog="treetip",
        description="Quality-controlled fixed-tree placement for immutable aligned sequences.",
    )
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    place = subparsers.add_parser(
        "place", help="Place aligned query sequences on a fixed reference tree."
    )
    place.add_argument("--tree", required=True, help="Reference branch-length Newick tree.")
    place.add_argument(
        "--reference-alignment", required=True, help="Reference fixed-column FASTA alignment."
    )
    place.add_argument("--query-alignment", required=True, help="Query fixed-column FASTA alignment.")
    place.add_argument("--output", required=True, help="New, empty output directory.")
    add_placement_controls(place)

    auto = subparsers.add_parser(
        "place-auto",
        help="Map raw query FASTA to classical and cryptic panels, then select conservatively.",
    )
    auto.add_argument("--reference-package", help="Schema-v1 paired classical/cryptic reference package JSON.")
    auto.add_argument("--classical-tree", help="Classical reference branch-length Newick tree.")
    auto.add_argument("--classical-reference-alignment", help="Classical reference fixed-column FASTA alignment.")
    auto.add_argument("--cryptic-tree", help="Cryptic reference branch-length Newick tree.")
    auto.add_argument("--cryptic-reference-alignment", help="Cryptic reference fixed-column FASTA alignment.")
    auto.add_argument(
        "--classical-tip-clades",
        help="Optional TSV with reference_tip_id and clade columns for classical fixed-tree clade calls.",
    )
    auto.add_argument(
        "--cryptic-tip-clades",
        help="Optional TSV with reference_tip_id and clade columns for cryptic fixed-tree clade calls.",
    )
    auto.add_argument("--query-fasta", required=True, help="Raw query FASTA; mapped separately to both panels.")
    auto.add_argument("--mafft", required=True, help="MAFFT executable or wrapper for --add --keeplength mapping.")
    auto.add_argument("--mafft-threads", type=positive_integer, default=1, help="MAFFT threads per mapping.")
    auto.add_argument(
        "--min-panel-distance-separation",
        type=nonnegative_float,
        default=None,
        help="Minimum nearest-reference distance advantage required for a panel assignment.",
    )
    auto.add_argument("--output", required=True, help="New, empty output directory.")
    add_placement_controls(auto)

    calibrate = subparsers.add_parser(
        "calibrate-panel-gate",
        help="Calibrate classical/cryptic distance-gate behavior with internal leave-one-out controls.",
    )
    calibrate.add_argument(
        "--classical-reference-alignment", required=True, help="Classical reference fixed-column FASTA alignment."
    )
    calibrate.add_argument(
        "--cryptic-reference-alignment", required=True, help="Cryptic reference fixed-column FASTA alignment."
    )

    clade_calibration = subparsers.add_parser(
        "calibrate-clade-assignment",
        help="Validate fixed-tree discrete clade calls with internal leave-one-out controls.",
    )
    clade_calibration.add_argument("--tree", required=True, help="Reference branch-length Newick tree.")
    clade_calibration.add_argument(
        "--reference-alignment", required=True, help="Reference fixed-column FASTA alignment."
    )
    clade_calibration.add_argument(
        "--tip-clades", required=True, help="Validated TSV with reference_tip_id and clade columns."
    )
    clade_calibration.add_argument(
        "--max-controls-per-clade",
        type=positive_integer,
        default=20,
        help="Maximum deterministic unique-sequence controls selected per clade.",
    )
    clade_calibration.add_argument("--output", required=True, help="New, empty output directory.")
    add_placement_controls(clade_calibration)

    accuracy = subparsers.add_parser(
        "calibrate-placement-accuracy",
        help="Measure whether leave-one-out placement recovers a held-out reference tip's position.",
    )
    accuracy.add_argument("--tree", required=True, help="Reference branch-length Newick tree.")
    accuracy.add_argument(
        "--reference-alignment", required=True, help="Reference fixed-column FASTA alignment."
    )
    accuracy.add_argument(
        "--max-controls",
        type=positive_integer,
        default=50,
        help="Maximum deterministic unique-sequence controls evaluated.",
    )
    accuracy.add_argument(
        "--near-hop-threshold",
        type=positive_integer,
        default=2,
        help="Largest edge-hop distance still reported as NEAR; provisional and uncalibrated.",
    )
    accuracy.add_argument("--output", required=True, help="New, empty output directory.")
    add_placement_controls(accuracy)

    calibrate.add_argument("--mafft", required=True, help="MAFFT executable or wrapper for cross-panel mapping.")
    calibrate.add_argument("--mafft-threads", type=positive_integer, default=1, help="MAFFT threads per mapping.")
    calibrate.add_argument(
        "--max-controls-per-panel",
        type=positive_integer,
        default=50,
        help="Maximum unique-sequence controls selected deterministically from each panel.",
    )
    calibrate.add_argument(
        "--min-panel-distance-separation",
        type=nonnegative_float,
        default=0.01,
        help="Distance margin evaluated as the default decision threshold.",
    )
    calibrate.add_argument("--output", required=True, help="New, empty output directory.")
    calibrate.add_argument("--min-coverage", type=probability, default=0.70)
    calibrate.add_argument("--warn-coverage", type=probability, default=0.90)
    calibrate.add_argument("--min-comparable-sites", type=positive_integer, default=100)
    calibrate.add_argument("--outlier-distance", type=probability, default=0.50)
    calibrate.add_argument(
        "-v", "--verbose", action="count", default=0, help="Increase stderr progress detail; repeat for more detail."
    )

    package = subparsers.add_parser(
        "create-reference-package",
        help="Create a checksummed portable paired classical/cryptic reference manifest.",
    )
    package.add_argument("--classical-tree", required=True, help="Classical reference branch-length Newick tree.")
    package.add_argument("--classical-reference-alignment", required=True, help="Classical reference fixed-column FASTA alignment.")
    package.add_argument("--cryptic-tree", required=True, help="Cryptic reference branch-length Newick tree.")
    package.add_argument("--cryptic-reference-alignment", required=True, help="Cryptic reference fixed-column FASTA alignment.")
    package.add_argument(
        "--classical-tip-clades",
        help="Optional validated classical reference_tip_id/clade TSV to checksum in the package.",
    )
    package.add_argument(
        "--cryptic-tip-clades",
        help="Optional validated cryptic reference_tip_id/clade TSV to checksum in the package.",
    )
    package.add_argument(
        "--classical-expected-clade-count",
        type=positive_integer,
        help="Require a complete classical tip-clade partition with this many distinct labels.",
    )
    package.add_argument(
        "--cryptic-expected-clade-count",
        type=positive_integer,
        help="Require a complete cryptic tip-clade partition with this many distinct labels.",
    )
    package.add_argument("--output", required=True, help="New paired reference package JSON path.")
    package.add_argument("--name", help="Optional human-readable package name; defaults to the output filename stem.")
    package.add_argument(
        "--min-panel-distance-separation",
        type=nonnegative_float,
        default=0.01,
        help="Calibrated minimum nearest-reference distance advantage for panel assignment.",
    )
    package.add_argument(
        "-v", "--verbose", action="count", default=0, help="Increase stderr progress detail; repeat for more detail."
    )

    validate_package = subparsers.add_parser(
        "validate-reference-package",
        help="Validate paired-reference checksums and tree/alignment identity contracts.",
    )
    validate_package.add_argument("--reference-package", required=True, help="Schema-v1 paired reference package JSON.")
    return parser


def add_placement_controls(parser: argparse.ArgumentParser) -> None:
    """Add the shared placement and QC controls to a command parser."""

    parser.add_argument(
        "--preset",
        choices=("fast", "balanced", "thorough"),
        default="balanced",
        help="Search scope: fast local, balanced local, or exhaustive thorough search.",
    )
    parser.add_argument(
        "--threads",
        default="auto",
        help="Parallel query workers; use 'auto' or a positive integer.",
    )
    parser.add_argument(
        "--bootstrap",
        type=nonnegative_integer,
        default=0,
        help="Alignment-site bootstrap replicates per successfully placed query.",
    )
    parser.add_argument("--seed", type=int, default=1, help="Bootstrap random seed.")
    parser.add_argument("--candidate-count", type=positive_integer, default=16)
    parser.add_argument("--candidate-hops", type=positive_integer, default=4)
    parser.add_argument("--min-coverage", type=probability, default=0.70)
    parser.add_argument("--warn-coverage", type=probability, default=0.90)
    parser.add_argument("--min-comparable-sites", type=positive_integer, default=100)
    parser.add_argument("--outlier-distance", type=probability, default=0.50)
    parser.add_argument("--unresolved-weight", type=probability, default=0.60)
    parser.add_argument(
        "--force",
        action="store_true",
        help="Attempt placement after a pre-placement QC failure; failures remain annotated.",
    )
    parser.add_argument(
        "--max-cache-megabytes",
        type=positive_integer,
        default=256,
        help="Per-process limit for cached node-to-tip distance rows.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase stderr progress detail; repeat for more detail.",
    )


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command line and convert expected errors to a clear exit status."""

    parser = build_parser()
    arguments = parser.parse_args(argv)
    try:
        if arguments.command == "place":
            return run_place(arguments)
        if arguments.command == "place-auto":
            return run_auto_place(arguments)
        if arguments.command == "calibrate-panel-gate":
            return run_calibrate_panel_gate(arguments)
        if arguments.command == "calibrate-clade-assignment":
            return run_calibrate_clade_assignment(arguments)
        if arguments.command == "calibrate-placement-accuracy":
            return run_calibrate_placement_accuracy(arguments)
        if arguments.command == "create-reference-package":
            return run_create_reference_package(arguments)
        if arguments.command == "validate-reference-package":
            return run_validate_reference_package(arguments)
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("ERROR: interrupted by user", file=sys.stderr)
        return 130
    return 2


def run_place(arguments: argparse.Namespace) -> int:
    """Execute one complete placement run."""

    logger = EventLogger(arguments.verbose)
    logger.event(1, "start", "Starting TreeTip Placer", version=__version__)
    output = prepare_output_directory(arguments.output)
    tree = Tree.from_path(arguments.tree, arguments.max_cache_megabytes)
    if len(tree.tips) < 2:
        raise ValueError("Reference tree must contain at least two tips")

    reference_alignment = read_fasta_alignment(arguments.reference_alignment)
    query_alignment = read_fasta_alignment(arguments.query_alignment)
    require_matching_alignment_widths(reference_alignment, query_alignment)
    reference_by_id = {sequence.sample_id: sequence for sequence in reference_alignment.sequences}
    tree_tip_names = set(tree.tip_index)
    reference_names = set(reference_by_id)
    if tree_tip_names != reference_names:
        missing_in_alignment = sorted(tree_tip_names - reference_names)
        missing_in_tree = sorted(reference_names - tree_tip_names)
        details: list[str] = []
        if missing_in_alignment:
            details.append("tips absent from reference alignment: " + ", ".join(missing_in_alignment[:10]))
        if missing_in_tree:
            details.append("alignment IDs absent from tree: " + ", ".join(missing_in_tree[:10]))
        raise ValueError("Reference tree and data IDs differ; " + "; ".join(details))
    queries = query_alignment.sequences
    query_names = {sequence.sample_id for sequence in queries}
    overlap = sorted(query_names & reference_names)
    if overlap:
        raise ValueError(
            "Query sample_id values already occur in the reference tree: " + ", ".join(overlap[:10])
        )

    references_in_tip_order = tuple(
        reference_by_id[tip.name or ""] for tip in tree.tips
    )
    config = placement_config_from_arguments(arguments)
    logger.event(
        1,
        "input",
        "Validated reference inputs",
        input_kind="fixed_column_fasta_alignment",
        reference_tips=len(references_in_tip_order),
        queries=len(queries),
        characters=len(references_in_tip_order[0].characters),
    )

    workers = resolve_threads(arguments.threads, len(queries))
    logger.event(1, "placement", "Running placements", workers=workers, preset=config.preset)
    started_placement = time.perf_counter()
    if workers == 1:
        results = run_serial_placements(
            queries, references_in_tip_order, tree, config, logger, started_placement
        )
    else:
        try:
            results = run_parallel_placements(
                queries,
                tree.newick(),
                references_in_tip_order,
                config,
                arguments.max_cache_megabytes,
                workers,
                logger,
                started_placement,
            )
        except PermissionError as error:
            logger.event(
                1,
                "placement",
                "Parallel workers are unavailable in this environment; falling back to serial placement",
                error=str(error),
            )
            workers = 1
            results = run_serial_placements(
                queries, references_in_tip_order, tree, config, logger, started_placement
            )

    summary = summarize_results(results)
    logger.event(1, "placement", "Completed placement", **summary)
    metadata = RunMetadata(
        version=__version__,
        tree_path=str(Path(arguments.tree).resolve()),
        reference_alignment_path=str(Path(arguments.reference_alignment).resolve()),
        query_alignment_path=str(Path(arguments.query_alignment).resolve()),
        tree_sha256=file_sha256(arguments.tree),
        reference_alignment_sha256=file_sha256(arguments.reference_alignment),
        query_alignment_sha256=file_sha256(arguments.query_alignment),
        config={**asdict(config), "input_kind": "fixed_column_fasta_alignment", "workers": workers},
    )
    logger.event(1, "output", "Writing interoperable outputs", directory=str(output))
    write_run_outputs(output, tree, results, metadata, logger.events)
    (output / "run.log").write_text(
        "".join(
            f"[{event['elapsed_seconds']:8.2f}s] {event['stage']}: {event['message']}"
            + (" " + json.dumps({key: value for key, value in event.items() if key not in {'elapsed_seconds', 'level', 'stage', 'message'}}, sort_keys=True) if len(event) > 4 else "")
            + "\n"
            for event in logger.events
        ),
        encoding="utf-8",
    )
    logger.event(1, "complete", "Run completed", **summary)
    append_final_event(output, logger.events[-1])
    return 0


def run_create_reference_package(arguments: argparse.Namespace) -> int:
    """Create a portable paired-reference manifest without copying reference data."""

    logger = EventLogger(arguments.verbose)
    logger.event(1, "start", "Creating paired reference package", version=__version__)
    output = Path(arguments.output)
    package = create_paired_reference_package(
        output,
        package_name=arguments.name or output.stem,
        classical_tree=arguments.classical_tree,
        classical_alignment=arguments.classical_reference_alignment,
        cryptic_tree=arguments.cryptic_tree,
        cryptic_alignment=arguments.cryptic_reference_alignment,
        minimum_distance_separation=arguments.min_panel_distance_separation,
        classical_tip_clades=arguments.classical_tip_clades,
        cryptic_tip_clades=arguments.cryptic_tip_clades,
        classical_expected_clade_count=arguments.classical_expected_clade_count,
        cryptic_expected_clade_count=arguments.cryptic_expected_clade_count,
    )
    logger.event(
        1,
        "complete",
        "Created paired reference package",
        output=str(package.manifest_path),
        minimum_distance_separation=package.minimum_distance_separation,
    )
    return 0


def run_validate_reference_package(arguments: argparse.Namespace) -> int:
    """Validate a package and print a machine-readable summary to stdout."""

    package = load_paired_reference_package(arguments.reference_package)
    summary = validate_paired_reference_package(package)
    print(
        json.dumps(
            {
                "reference_package": str(package.manifest_path),
                "minimum_distance_separation": package.minimum_distance_separation,
                "panels": summary,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


def run_auto_place(arguments: argparse.Namespace) -> int:
    """Map raw input to both panels, place independently, and select conservatively."""

    logger = EventLogger(arguments.verbose)
    logger.event(1, "start", "Starting automatic classical/cryptic panel placement", version=__version__)
    package, minimum_distance_separation = resolve_auto_reference_inputs(arguments)
    output = prepare_output_directory(arguments.output)
    source_records = read_fasta_records(arguments.query_fasta)
    source_ids = tuple(sample_id for sample_id, _ in source_records)
    logger.event(1, "input", "Validated raw query FASTA", queries=len(source_ids))

    mapping_directory = output / "mapping"
    mapping_directory.mkdir()
    classical_queries = map_queries_to_panel(
        panel_name="classical",
        mafft=arguments.mafft,
        mafft_threads=arguments.mafft_threads,
        query_fasta=arguments.query_fasta,
        query_ids=source_ids,
        reference_alignment=str(package.classical.alignment_path),
        mapping_directory=mapping_directory,
        logger=logger,
    )
    cryptic_queries = map_queries_to_panel(
        panel_name="cryptic",
        mafft=arguments.mafft,
        mafft_threads=arguments.mafft_threads,
        query_fasta=arguments.query_fasta,
        query_ids=source_ids,
        reference_alignment=str(package.cryptic.alignment_path),
        mapping_directory=mapping_directory,
        logger=logger,
    )

    panels_root = output / "panels"
    panels_root.mkdir()
    classical = run_panel_placement(
        panel_name="classical",
        tree_path=str(package.classical.tree_path),
        reference_alignment_path=str(package.classical.alignment_path),
        query_alignment_path=classical_queries,
        output_directory=panels_root / "classical",
        arguments=arguments,
        logger=logger,
    )
    cryptic = run_panel_placement(
        panel_name="cryptic",
        tree_path=str(package.cryptic.tree_path),
        reference_alignment_path=str(package.cryptic.alignment_path),
        query_alignment_path=cryptic_queries,
        output_directory=panels_root / "cryptic",
        arguments=arguments,
        logger=logger,
    )
    decisions = choose_reference_panel(
        classical.results,
        cryptic.results,
        minimum_distance_separation,
    )
    selection_summary = summarize_decisions(decisions)
    logger.event(1, "selection", "Completed conservative panel selection", **selection_summary)
    classical_tip_clades_path = resolve_tip_clades_path(arguments.classical_tip_clades, package.classical)
    cryptic_tip_clades_path = resolve_tip_clades_path(arguments.cryptic_tip_clades, package.cryptic)
    tip_clades = {
        "classical": load_tip_clades(classical_tip_clades_path, classical.tree, "classical"),
        "cryptic": load_tip_clades(cryptic_tip_clades_path, cryptic.tree, "cryptic"),
    }
    clade_assignments = assign_phylogenetic_clades(
        decisions,
        {
            "classical": (classical.tree, classical.results),
            "cryptic": (cryptic.tree, cryptic.results),
        },
        tip_clades,
    )
    clade_summary = summarize_clade_assignments(clade_assignments)
    logger.event(1, "clade_assignment", "Completed fixed-tree clade assignment", **clade_summary)
    selection_metadata = {
        "version": __version__,
        "query_fasta_path": str(Path(arguments.query_fasta).resolve()),
        "query_fasta_sha256": file_sha256(arguments.query_fasta),
        "reference_package_path": str(package.manifest_path) if package.manifest_path else None,
        "classical_tree_path": str(package.classical.tree_path),
        "classical_reference_alignment_path": str(package.classical.alignment_path),
        "cryptic_tree_path": str(package.cryptic.tree_path),
        "cryptic_reference_alignment_path": str(package.cryptic.alignment_path),
        "minimum_distance_separation": minimum_distance_separation,
        "selection_method": "nearest_reference_normalized_hamming_distance",
        "classical_tip_clades_path": (
            str(classical_tip_clades_path)
            if classical_tip_clades_path
            else None
        ),
        "cryptic_tip_clades_path": (
            str(cryptic_tip_clades_path)
            if cryptic_tip_clades_path
            else None
        ),
        "clade_assignment_method": "selected_edge_descendant_tip_labels_with_nearest_reference_fallback",
        "warning": "Distance-separation thresholds require empirical calibration for each pair of panels.",
    }
    write_panel_selection_outputs(output, decisions, selection_metadata)
    write_clade_assignment_outputs(output, clade_assignments, selection_metadata)
    write_selected_panel_trees(output, classical, cryptic, decisions)
    write_auto_report(output, decisions, selection_metadata, clade_assignments, tip_clades)
    (output / "run.log").write_text(format_event_log(logger.events), encoding="utf-8")
    logger.event(1, "complete", "Automatic panel placement completed", **selection_summary)
    append_final_event(output, logger.events[-1])
    return 0


def resolve_tip_clades_path(explicit_path: str | None, panel: PanelReference) -> Path | None:
    """Prefer an explicit metadata path, otherwise use package-bound metadata."""

    return Path(explicit_path).resolve() if explicit_path else panel.tip_clades_path


def resolve_auto_reference_inputs(
    arguments: argparse.Namespace,
) -> tuple[PairedReferencePackage, float]:
    """Resolve either a checked package manifest or a complete explicit file quartet."""

    explicit_values = (
        arguments.classical_tree,
        arguments.classical_reference_alignment,
        arguments.cryptic_tree,
        arguments.cryptic_reference_alignment,
    )
    if arguments.reference_package:
        if any(explicit_values):
            raise ValueError(
                "--reference-package cannot be combined with explicit classical/cryptic tree or alignment paths"
            )
        package = load_paired_reference_package(arguments.reference_package)
        threshold = (
            arguments.min_panel_distance_separation
            if arguments.min_panel_distance_separation is not None
            else package.minimum_distance_separation
        )
        return package, threshold
    if not all(explicit_values):
        raise ValueError(
            "Provide --reference-package or all four explicit classical/cryptic tree and alignment paths"
        )
    package = PairedReferencePackage(
        manifest_path=None,
        classical=PanelReference(
            tree_path=Path(arguments.classical_tree).resolve(),
            alignment_path=Path(arguments.classical_reference_alignment).resolve(),
            tree_sha256=None,
            alignment_sha256=None,
        ),
        cryptic=PanelReference(
            tree_path=Path(arguments.cryptic_tree).resolve(),
            alignment_path=Path(arguments.cryptic_reference_alignment).resolve(),
            tree_sha256=None,
            alignment_sha256=None,
        ),
        minimum_distance_separation=0.01,
    )
    threshold = (
        arguments.min_panel_distance_separation
        if arguments.min_panel_distance_separation is not None
        else package.minimum_distance_separation
    )
    return package, threshold


def build_panel_gate_summary(
    *,
    classical_controls_count: int,
    cryptic_controls_count: int,
    default_threshold: float,
    conflict_count: int,
    sweep: Sequence[dict[str, object]],
    calibration_rows: Sequence[dict[str, object]],
    classical_reference_alignment_sha256: str,
    cryptic_reference_alignment_sha256: str,
) -> dict[str, object]:
    """Build the panel-gate calibration provenance and outcome summary.

    Extracted as a pure function so the outcome-count contract this feeds to
    write_calibration_report can be exercised by a test without requiring
    MAFFT, which run_calibrate_panel_gate otherwise depends on end to end.
    """

    return {
        "version": __version__,
        "method": "leave_one_unique_sequence_out_panel_distance_calibration",
        "warning": "This calibrates cross-panel distance gating only; it is not external validation or topology-placement accuracy.",
        "classical_controls": classical_controls_count,
        "cryptic_controls": cryptic_controls_count,
        "default_threshold": default_threshold,
        "flagged_panel_label_conflicts": conflict_count,
        "default_result": next(
            row for row in sweep if row["threshold"] == default_threshold
        ),
        "threshold_sweep": sweep,
        "summary": summarize_clade_calibration(calibration_rows),
        "classical_reference_alignment_sha256": classical_reference_alignment_sha256,
        "cryptic_reference_alignment_sha256": cryptic_reference_alignment_sha256,
    }


def run_calibrate_panel_gate(arguments: argparse.Namespace) -> int:
    """Calibrate panel assignment on internal controls without topology placement."""

    logger = EventLogger(arguments.verbose)
    logger.event(1, "start", "Starting internal classical/cryptic panel-gate calibration", version=__version__)
    output = prepare_output_directory(arguments.output)
    mapping_directory = output / "mapping"
    mapping_directory.mkdir()
    classical_alignment = read_fasta_alignment(arguments.classical_reference_alignment)
    cryptic_alignment = read_fasta_alignment(arguments.cryptic_reference_alignment)
    classical_controls = select_unique_sequence_controls(
        classical_alignment.sequences, arguments.max_controls_per_panel
    )
    cryptic_controls = select_unique_sequence_controls(
        cryptic_alignment.sequences, arguments.max_controls_per_panel
    )
    logger.event(
        1,
        "input",
        "Selected unique-sequence internal controls",
        classical_controls=len(classical_controls),
        cryptic_controls=len(cryptic_controls),
    )

    classical_queries = calibration_queries("classical", classical_controls)
    cryptic_queries = calibration_queries("cryptic", cryptic_controls)
    classical_raw = mapping_directory / "classical_controls.ungapped.fasta"
    cryptic_raw = mapping_directory / "cryptic_controls.ungapped.fasta"
    write_fasta_records(classical_raw, ungapped_fasta_records(classical_queries))
    write_fasta_records(cryptic_raw, ungapped_fasta_records(cryptic_queries))
    classical_to_cryptic = map_queries_to_panel(
        panel_name="classical_controls_to_cryptic",
        mafft=arguments.mafft,
        mafft_threads=arguments.mafft_threads,
        query_fasta=str(classical_raw),
        query_ids=tuple(query.sample_id for query in classical_queries),
        reference_alignment=arguments.cryptic_reference_alignment,
        mapping_directory=mapping_directory,
        logger=logger,
    )
    cryptic_to_classical = map_queries_to_panel(
        panel_name="cryptic_controls_to_classical",
        mafft=arguments.mafft,
        mafft_threads=arguments.mafft_threads,
        query_fasta=str(cryptic_raw),
        query_ids=tuple(query.sample_id for query in cryptic_queries),
        reference_alignment=arguments.classical_reference_alignment,
        mapping_directory=mapping_directory,
        logger=logger,
    )
    cross_to_cryptic = {
        sequence.sample_id: sequence
        for sequence in read_fasta_alignment(classical_to_cryptic).sequences
    }
    cross_to_classical = {
        sequence.sample_id: sequence
        for sequence in read_fasta_alignment(cryptic_to_classical).sequences
    }

    classical_results: list[PlacementResult] = []
    cryptic_results: list[PlacementResult] = []
    truth_by_id: dict[str, str] = {}
    source_id_by_id: dict[str, str] = {}
    for query in classical_queries:
        classical_results.append(
            qc_only_result(query, references_without_identical_sequence(query, classical_alignment.sequences), arguments)
        )
        cryptic_results.append(qc_only_result(cross_to_cryptic[query.sample_id], cryptic_alignment.sequences, arguments))
        truth_by_id[query.sample_id] = "CLASSICAL"
        source_id_by_id[query.sample_id] = query.sample_id.removeprefix("classical::")
    for query in cryptic_queries:
        classical_results.append(qc_only_result(cross_to_classical[query.sample_id], classical_alignment.sequences, arguments))
        cryptic_results.append(
            qc_only_result(query, references_without_identical_sequence(query, cryptic_alignment.sequences), arguments)
        )
        truth_by_id[query.sample_id] = "CRYPTIC"
        source_id_by_id[query.sample_id] = query.sample_id.removeprefix("cryptic::")

    decisions = choose_reference_panel(
        classical_results, cryptic_results, arguments.min_panel_distance_separation
    )
    conflict_count, calibration_rows, calibration_fields = write_calibration_records(
        output, decisions, truth_by_id, source_id_by_id
    )
    sweep = panel_threshold_sweep(
        classical_results,
        cryptic_results,
        truth_by_id,
        arguments.min_panel_distance_separation,
    )
    summary = build_panel_gate_summary(
        classical_controls_count=len(classical_controls),
        cryptic_controls_count=len(cryptic_controls),
        default_threshold=arguments.min_panel_distance_separation,
        conflict_count=conflict_count,
        sweep=sweep,
        calibration_rows=calibration_rows,
        classical_reference_alignment_sha256=file_sha256(arguments.classical_reference_alignment),
        cryptic_reference_alignment_sha256=file_sha256(arguments.cryptic_reference_alignment),
    )
    (output / "calibration_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_calibration_report(
        output,
        "Panel-gate calibration",
        summary,
        calibration_rows,
        calibration_fields,
        "outcome",
        frozenset({"INCORRECT_CALL"}),
        (
            ("calibration_records.tsv", "Per-control records (TSV)"),
            ("panel_label_conflicts.tsv", "Flagged label conflicts (TSV)"),
            ("calibration_summary.json", "Calibration summary (JSON)"),
        ),
        threshold_sweep=sweep,
    )
    logger.event(1, "calibration", "Completed panel-gate calibration", **summary["default_result"])
    (output / "run.log").write_text(format_event_log(logger.events), encoding="utf-8")
    logger.event(1, "complete", "Panel-gate calibration completed")
    append_final_event(output, logger.events[-1])
    return 0


def run_calibrate_clade_assignment(arguments: argparse.Namespace) -> int:
    """Evaluate discrete calls by holding out labelled unique-sequence tips.

    Each control is removed from both the immutable backbone tree and reference
    alignment before scoring. This is an internal accuracy diagnostic for the
    fixed-tree clade rule, not topology inference or phylogenetic support.
    """

    if arguments.bootstrap:
        raise ValueError("Clade-assignment calibration does not accept bootstrap replicates; use --bootstrap 0")
    logger = EventLogger(arguments.verbose)
    logger.event(1, "start", "Starting fixed-tree clade leave-one-out calibration", version=__version__)
    output = prepare_output_directory(arguments.output)
    tree = Tree.from_path(arguments.tree, arguments.max_cache_megabytes)
    if len(tree.tips) < 4:
        raise ValueError("Clade-assignment calibration requires at least four reference tips")
    alignment = read_fasta_alignment(arguments.reference_alignment)
    reference_by_id = {sequence.sample_id: sequence for sequence in alignment.sequences}
    if set(tree.tip_index) != set(reference_by_id):
        raise ValueError("Reference tree tips do not exactly match reference alignment IDs")
    labels = load_tip_clades(arguments.tip_clades, tree, "reference")
    label_counts = clade_label_counts(labels)
    controls = select_clade_calibration_controls(
        alignment.sequences, labels, arguments.max_controls_per_clade
    )
    if not controls:
        raise ValueError("No labelled unique-sequence controls are available for clade calibration")
    logger.event(
        1,
        "input",
        "Selected labelled unique-sequence leave-one-out controls",
        controls=len(controls),
        labelled_tips=len(labels),
        eligible_clades=sum(count >= 2 for count in label_counts.values()),
        excluded_singleton_clades=sum(count == 1 for count in label_counts.values()),
    )
    config = placement_config_from_arguments(arguments)
    records: list[dict[str, object]] = []
    for number, held_out_id in enumerate(controls, start=1):
        pruned_tree, _ = tree.without_tip(held_out_id)
        held_out = reference_by_id[held_out_id]
        pruned_references = tuple(
            reference_by_id[tip.name or ""] for tip in pruned_tree.tips
        )
        result = place_query(held_out, pruned_references, pruned_tree, config)
        decision = PanelDecision(
            sample_id=held_out_id,
            assignment="REFERENCE",
            decision_status=result.final_status,
            selected_panel="reference",
            selected_placement_status=result.final_status,
            classical_nearest_reference=None,
            classical_nearest_distance=None,
            cryptic_nearest_reference=None,
            cryptic_nearest_distance=None,
            distance_margin=None,
            minimum_distance_separation=0.0,
            findings=("INTERNAL_LEAVE_ONE_OUT_CONTROL",),
        )
        remaining_labels = {tip_id: clade for tip_id, clade in labels.items() if tip_id != held_out_id}
        assignment = assign_phylogenetic_clades(
            [decision], {"reference": (pruned_tree, (result,))}, {"reference": remaining_labels}
        )[0]
        predicted = assignment.phylogenetic_clade
        truth = labels[held_out_id]
        if predicted is None:
            outcome = "UNRESOLVED"
        elif predicted == truth:
            outcome = "CORRECT_PASS" if assignment.clade_assignment_status == "PASS" else "CORRECT_WARN"
        else:
            outcome = "INCORRECT_CALL"
        records.append(
            {
                "source_sample_id": held_out_id,
                "truth_clade": truth,
                "predicted_clade": predicted,
                "clade_assignment_status": assignment.clade_assignment_status,
                "assignment_method": assignment.assignment_method,
                "outcome": outcome,
                "nearest_reference": assignment.nearest_reference,
                "supporting_reference_tips": " | ".join(assignment.supporting_reference_tips),
                "findings": " | ".join(assignment.findings),
                "placement_status": result.final_status,
                "placement_edge_id": result.best_fit.edge_id if result.best_fit else None,
                "coverage": result.qc.coverage,
                "comparable_sites": result.qc.comparable_sites,
                "nearest_reference_distance": result.qc.nearest_distance,
                "relative_fit_weight": (
                    result.best_fit.relative_fit_weight if result.best_fit else None
                ),
                "best_objective": result.best_fit.objective if result.best_fit else None,
                "best_rmse": result.best_fit.rmse if result.best_fit else None,
            }
        )
        if number == len(controls) or number % 10 == 0:
            logger.event(2, "calibration", "Completed leave-one-out controls", completed=number, controls=len(controls))
    write_clade_calibration_outputs(output, records, arguments, tree, alignment, label_counts)
    summary = summarize_clade_calibration(records)
    logger.event(1, "calibration", "Completed fixed-tree clade calibration", **summary)
    (output / "run.log").write_text(format_event_log(logger.events), encoding="utf-8")
    logger.event(1, "complete", "Clade-assignment calibration completed")
    append_final_event(output, logger.events[-1])
    return 0


def select_clade_calibration_controls(
    sequences: Sequence[AlignedSequence], labels: dict[str, str], maximum_per_clade: int
) -> tuple[str, ...]:
    """Select controls only when their clade retains a labelled reference tip."""

    counts: dict[tuple[str | None, ...], int] = {}
    for sequence in sequences:
        counts[sequence.characters] = counts.get(sequence.characters, 0) + 1
    label_counts = clade_label_counts(labels)
    by_clade: dict[str, list[str]] = {}
    for sequence in sequences:
        clade = labels.get(sequence.sample_id)
        if (
            clade is not None
            and label_counts[clade] >= 2
            and counts[sequence.characters] == 1
        ):
            by_clade.setdefault(clade, []).append(sequence.sample_id)
    return tuple(
        sample_id
        for clade in sorted(by_clade)
        for sample_id in sorted(by_clade[clade])[:maximum_per_clade]
    )


def clade_label_counts(labels: dict[str, str]) -> dict[str, int]:
    """Count annotated reference tips per discrete clade."""

    counts: dict[str, int] = {}
    for clade in labels.values():
        counts[clade] = counts.get(clade, 0) + 1
    return counts


def write_clade_calibration_outputs(
    output: Path,
    records: Sequence[dict[str, object]],
    arguments: argparse.Namespace,
    tree: Tree,
    alignment: object,
    label_counts: dict[str, int],
) -> None:
    """Write auditable records, confusion counts, and provenance for clade calibration."""

    fields = [
        "source_sample_id", "truth_clade", "predicted_clade", "clade_assignment_status",
        "assignment_method", "outcome", "nearest_reference", "supporting_reference_tips",
        "findings", "placement_status", "placement_edge_id", "coverage", "comparable_sites",
        "nearest_reference_distance", "relative_fit_weight", "best_objective", "best_rmse",
    ]
    with (output / "clade_calibration_records.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    confusion: dict[tuple[str, str], int] = {}
    for record in records:
        key = (str(record["truth_clade"]), str(record["predicted_clade"] or "UNRESOLVED"))
        confusion[key] = confusion.get(key, 0) + 1
    with (output / "clade_confusion.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["truth_clade", "predicted_clade", "controls"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for (truth, predicted), controls in sorted(confusion.items()):
            writer.writerow({"truth_clade": truth, "predicted_clade": predicted, "controls": controls})
    with (output / "clade_calibration_exclusions.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["clade", "labelled_reference_tips", "reason"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for clade, count in sorted(label_counts.items()):
            if count < 2:
                writer.writerow(
                    {
                        "clade": clade,
                        "labelled_reference_tips": count,
                        "reason": "ONLY_ONE_LABELLED_REFERENCE_TIP_REMAINS_AFTER_HOLDOUT",
                    }
                )
    summary = {
        "version": __version__,
        "method": "fixed_tree_leave_one_unique_sequence_out_discrete_clade_assignment",
        "warning": "Internal validation only: this does not infer a new topology, recalculate phylogenetic support, or estimate external generalization.",
        "tree_path": str(Path(arguments.tree).resolve()),
        "tree_sha256": file_sha256(arguments.tree),
        "reference_alignment_path": str(Path(arguments.reference_alignment).resolve()),
        "reference_alignment_sha256": file_sha256(arguments.reference_alignment),
        "tip_clades_path": str(Path(arguments.tip_clades).resolve()),
        "tip_clades_sha256": file_sha256(arguments.tip_clades),
        "reference_tips": len(tree.tips),
        "reference_columns": getattr(alignment, "site_count"),
        "max_controls_per_clade": arguments.max_controls_per_clade,
        "labelled_reference_tips_per_clade": dict(sorted(label_counts.items())),
        "excluded_singleton_clades": sorted(
            clade for clade, count in label_counts.items() if count < 2
        ),
        "summary": summarize_clade_calibration(records),
        "summary_by_truth_clade": summarize_clade_calibration_by_truth(records),
    }
    (output / "clade_calibration_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_calibration_report(
        output,
        "Clade-assignment calibration",
        summary,
        records,
        fields,
        "outcome",
        frozenset({"INCORRECT_CALL"}),
        (
            ("clade_calibration_records.tsv", "Per-control records (TSV)"),
            ("clade_confusion.tsv", "Clade confusion counts (TSV)"),
            ("clade_calibration_exclusions.tsv", "Excluded singleton clades (TSV)"),
            ("clade_calibration_summary.json", "Calibration summary (JSON)"),
        ),
    )


def summarize_clade_calibration(records: Sequence[dict[str, object]]) -> dict[str, int]:
    """Count calibration outcomes using stable machine-readable keys."""

    counts: dict[str, int] = {"controls": len(records)}
    for record in records:
        outcome = str(record["outcome"]).lower()
        counts[f"outcome_{outcome}"] = counts.get(f"outcome_{outcome}", 0) + 1
    return counts


def summarize_clade_calibration_by_truth(
    records: Sequence[dict[str, object]],
) -> dict[str, dict[str, int]]:
    """Count diagnostic outcomes by held-out reference clade."""

    summary: dict[str, dict[str, int]] = {}
    for record in records:
        clade = str(record["truth_clade"])
        outcome = str(record["outcome"]).lower()
        counts = summary.setdefault(clade, {"controls": 0})
        counts["controls"] += 1
        key = f"outcome_{outcome}"
        counts[key] = counts.get(key, 0) + 1
    return {clade: summary[clade] for clade in sorted(summary)}


def run_calibrate_placement_accuracy(arguments: argparse.Namespace) -> int:
    """Measure whether leave-one-out placement recovers a held-out tip's position.

    Each control is removed from both the immutable backbone tree and the
    reference alignment, placed back with the production placement path, and
    graded in edge hops against its recorded attachment. This is an internal
    topology diagnostic, not phylogenetic support or external validation.
    """

    logger = EventLogger(arguments.verbose)
    logger.event(1, "start", "Starting leave-one-out placement-accuracy calibration", version=__version__)
    output = prepare_output_directory(arguments.output)
    tree = Tree.from_path(arguments.tree, arguments.max_cache_megabytes)
    if len(tree.tips) < 4:
        raise ValueError("Placement-accuracy calibration requires at least four reference tips")
    alignment = read_fasta_alignment(arguments.reference_alignment)
    reference_by_id = {sequence.sample_id: sequence for sequence in alignment.sequences}
    tree_tip_names = set(tree.tip_index)
    alignment_names = set(reference_by_id)
    if tree_tip_names != alignment_names:
        missing_in_alignment = sorted(tree_tip_names - alignment_names)
        missing_in_tree = sorted(alignment_names - tree_tip_names)
        details: list[str] = []
        if missing_in_alignment:
            details.append("tips absent from the alignment: " + ", ".join(missing_in_alignment[:10]))
        if missing_in_tree:
            details.append("alignment IDs absent from the tree: " + ", ".join(missing_in_tree[:10]))
        raise ValueError("Reference tree and alignment IDs differ; " + "; ".join(details))

    controls = tuple(
        sequence.sample_id
        for sequence in select_unique_sequence_controls(alignment.sequences, arguments.max_controls)
    )
    logger.event(
        1,
        "input",
        "Selected unique-sequence leave-one-out controls",
        reference_tips=len(tree.tips),
        controls=len(controls),
        near_hop_threshold=arguments.near_hop_threshold,
    )
    config = placement_config_from_arguments(arguments)
    records: list[dict[str, object]] = []
    for number, held_out_id in enumerate(controls, start=1):
        pruned_tree, attachment = tree.without_tip(held_out_id)
        pruned_references = tuple(
            reference_by_id[tip.name or ""] for tip in pruned_tree.tips
        )
        result = place_query(reference_by_id[held_out_id], pruned_references, pruned_tree, config)
        grade = grade_placement_accuracy(
            result, pruned_tree, attachment, arguments.near_hop_threshold
        )
        records.append(build_accuracy_record(result, grade, attachment, held_out_id))
        if number == len(controls) or number % 10 == 0:
            logger.event(
                2, "calibration", "Completed leave-one-out controls", completed=number, controls=len(controls)
            )

    summary = write_placement_accuracy_outputs(output, records, arguments, tree, alignment)
    logger.event(1, "calibration", "Completed placement-accuracy calibration", **summary["summary"])
    (output / "run.log").write_text(format_event_log(logger.events), encoding="utf-8")
    logger.event(1, "complete", "Placement-accuracy calibration completed")
    append_final_event(output, logger.events[-1])
    return 0


def write_placement_accuracy_outputs(
    output: Path,
    records: Sequence[dict[str, object]],
    arguments: argparse.Namespace,
    tree: Tree,
    alignment: object,
) -> dict[str, object]:
    """Write auditable per-control records and a provenance summary."""

    with (output / "placement_accuracy_records.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(ACCURACY_RECORD_COLUMNS), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for record in records:
            writer.writerow({key: "" if value is None else value for key, value in record.items()})
    summary = {
        "version": __version__,
        "method": "fixed_tree_leave_one_unique_sequence_out_placement_accuracy",
        "warning": (
            "Internal topology diagnostics only: this is not phylogenetic support, "
            "not a branch-length error, and not evidence of external generalization."
        ),
        "tree_path": str(Path(arguments.tree).resolve()),
        "tree_sha256": file_sha256(arguments.tree),
        "reference_alignment_path": str(Path(arguments.reference_alignment).resolve()),
        "reference_alignment_sha256": file_sha256(arguments.reference_alignment),
        "reference_tips": len(tree.tips),
        "reference_columns": getattr(alignment, "site_count"),
        "max_controls": arguments.max_controls,
        "near_hop_threshold": arguments.near_hop_threshold,
        "summary": summarize_placement_accuracy(records),
        "hop_distribution": hop_distribution(records),
    }
    (output / "placement_accuracy_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_calibration_report(
        output,
        "Placement-accuracy calibration",
        summary,
        records,
        ACCURACY_RECORD_COLUMNS,
        "outcome",
        frozenset({"FAR"}),
        (
            ("placement_accuracy_records.tsv", "Per-control records (TSV)"),
            ("placement_accuracy_summary.json", "Calibration summary (JSON)"),
        ),
    )
    return summary


def select_unique_sequence_controls(
    sequences: Sequence[AlignedSequence], maximum: int
) -> tuple[AlignedSequence, ...]:
    """Choose deterministic representatives of distinct aligned sequences."""

    groups: dict[tuple[str | None, ...], list[AlignedSequence]] = {}
    for sequence in sequences:
        groups.setdefault(sequence.characters, []).append(sequence)
    representatives = [
        min(group, key=lambda sequence: sequence.sample_id) for group in groups.values()
    ]
    return tuple(sorted(representatives, key=lambda sequence: sequence.sample_id)[:maximum])


def calibration_queries(panel: str, controls: Sequence[AlignedSequence]) -> tuple[AlignedSequence, ...]:
    """Prefix query identifiers so equal IDs cannot collide between panels."""

    return tuple(
        AlignedSequence(f"{panel}::{sequence.sample_id}", sequence.characters)
        for sequence in controls
    )


def ungapped_fasta_records(
    sequences: Sequence[AlignedSequence],
) -> tuple[tuple[str, str], ...]:
    """Create raw-like records for cross-panel MAFFT mapping from aligned controls."""

    return tuple(
        (sequence.sample_id, "".join(base for base in sequence.characters if base is not None))
        for sequence in sequences
    )


def references_without_identical_sequence(
    query: AlignedSequence, references: Sequence[AlignedSequence]
) -> tuple[AlignedSequence, ...]:
    """Remove every exact sequence duplicate to prevent trivial self matches."""

    retained = tuple(reference for reference in references if reference.characters != query.characters)
    if not retained:
        raise ValueError("Leave-one-unique-sequence-out calibration removed every reference")
    return retained


def qc_only_result(
    query: AlignedSequence, references: Sequence[AlignedSequence], arguments: argparse.Namespace
) -> PlacementResult:
    """Evaluate only the panel gate; no tree placement is implied or performed."""

    qc = evaluate_sequence_qc(
        query,
        references,
        min_coverage=arguments.min_coverage,
        warn_coverage=arguments.warn_coverage,
        min_comparable_sites=arguments.min_comparable_sites,
        outlier_distance=arguments.outlier_distance,
    )
    return PlacementResult(sample_id=query.sample_id, qc=qc, final_status=qc.status)


def panel_threshold_sweep(
    classical_results: Sequence[PlacementResult],
    cryptic_results: Sequence[PlacementResult],
    truth_by_id: dict[str, str],
    default_threshold: float,
) -> list[dict[str, object]]:
    """Report accuracy and retained-call coverage across useful distance margins."""

    thresholds = tuple(sorted({0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.01, 0.02, default_threshold}))
    rows: list[dict[str, object]] = []
    for threshold in thresholds:
        decisions = choose_reference_panel(classical_results, cryptic_results, threshold)
        called = [decision for decision in decisions if decision.assignment in {"CLASSICAL", "CRYPTIC"}]
        correct = sum(decision.assignment == truth_by_id[decision.sample_id] for decision in called)
        rows.append(
            {
                "threshold": threshold,
                "controls": len(decisions),
                "called": len(called),
                "ambiguous_or_unresolved": len(decisions) - len(called),
                "correct_calls": correct,
                "wrong_calls": len(called) - correct,
                "called_accuracy": correct / len(called) if called else None,
                "coverage": len(called) / len(decisions) if decisions else 0.0,
            }
        )
    return rows


def write_calibration_records(
    output_directory: Path,
    decisions: Sequence[object],
    truth_by_id: dict[str, str],
    source_id_by_id: dict[str, str],
) -> tuple[int, list[dict[str, object]], tuple[str, ...]]:
    """Write one auditable calibration record per held-out panel control.

    Returns the flagged-conflict count, the rows written, and the column
    order used, so a report can render them without re-reading the file or
    restating the column list.
    """

    fields = [
        "calibration_id", "source_sample_id", "truth", "assignment", "decision_status",
        "outcome", "classical_nearest_distance", "cryptic_nearest_distance", "distance_margin",
    ]
    conflicts: list[dict[str, object]] = []
    written_rows: list[dict[str, object]] = []
    with (output_directory / "calibration_records.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for decision in decisions:
            sample_id = getattr(decision, "sample_id")
            assignment = getattr(decision, "assignment")
            truth = truth_by_id[sample_id]
            called = assignment in {"CLASSICAL", "CRYPTIC"}
            outcome = "CORRECT_CALL" if assignment == truth else "INCORRECT_CALL" if called else "AMBIGUOUS_OR_UNRESOLVED"
            row: dict[str, object] = {
                "calibration_id": sample_id,
                "source_sample_id": source_id_by_id[sample_id],
                "truth": truth,
                "assignment": assignment,
                "decision_status": getattr(decision, "decision_status"),
                "outcome": outcome,
                "classical_nearest_distance": getattr(decision, "classical_nearest_distance"),
                "cryptic_nearest_distance": getattr(decision, "cryptic_nearest_distance"),
                "distance_margin": getattr(decision, "distance_margin"),
            }
            for field in ("classical_nearest_distance", "cryptic_nearest_distance", "distance_margin"):
                row[field] = "" if row[field] is None else format(row[field], ".12g")
            writer.writerow(row)
            written_rows.append(row)
            if outcome == "INCORRECT_CALL":
                conflicts.append(row)
    with (output_directory / "panel_label_conflicts.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(conflicts)
    return len(conflicts), written_rows, tuple(fields)


def map_queries_to_panel(
    *,
    panel_name: str,
    mafft: str,
    mafft_threads: int,
    query_fasta: str,
    query_ids: Sequence[str],
    reference_alignment: str,
    mapping_directory: Path,
    logger: EventLogger,
) -> Path:
    """Use MAFFT add/keeplength mapping and retain only mapped query records."""

    logger.event(1, "mapping", "Mapping queries to reference panel", panel=panel_name)
    mapped_with_reference = mapping_directory / f"{panel_name}.mapped_with_reference.fasta"
    mafft_stderr = mapping_directory / f"{panel_name}.mafft.stderr.log"
    command = [
        mafft,
        "--add",
        str(Path(query_fasta).resolve()),
        "--keeplength",
        "--thread",
        str(mafft_threads),
        str(Path(reference_alignment).resolve()),
    ]
    try:
        with mapped_with_reference.open("w", encoding="utf-8", newline="\n") as stdout:
            completed = subprocess.run(command, stdout=stdout, stderr=subprocess.PIPE, text=True, check=False)
    except FileNotFoundError as error:
        raise ValueError(f"MAFFT executable was not found: {mafft}") from error
    mafft_stderr.write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise ValueError(
            f"MAFFT mapping failed for the {panel_name} panel (exit {completed.returncode}); "
            f"see {mafft_stderr}"
        )
    mapped_by_id = dict(read_fasta_records(mapped_with_reference))
    missing = [sample_id for sample_id in query_ids if sample_id not in mapped_by_id]
    if missing:
        raise ValueError(
            f"MAFFT output for the {panel_name} panel omitted query IDs: {', '.join(missing[:10])}"
        )
    mapped_queries = mapping_directory / f"{panel_name}.queries_mapped.fasta"
    write_fasta_records(mapped_queries, tuple((sample_id, mapped_by_id[sample_id]) for sample_id in query_ids))
    require_matching_alignment_widths(
        read_fasta_alignment(reference_alignment), read_fasta_alignment(mapped_queries)
    )
    logger.event(
        1,
        "mapping",
        "Completed panel mapping",
        panel=panel_name,
        mapped_queries=len(query_ids),
        columns=read_fasta_alignment(mapped_queries).site_count,
    )
    return mapped_queries


def run_panel_placement(
    *,
    panel_name: str,
    tree_path: str,
    reference_alignment_path: str,
    query_alignment_path: Path,
    output_directory: Path,
    arguments: argparse.Namespace,
    logger: EventLogger,
) -> PanelRun:
    """Run the existing fixed-tree pipeline for one named panel."""

    tree = Tree.from_path(tree_path, arguments.max_cache_megabytes)
    if len(tree.tips) < 2:
        raise ValueError(f"The {panel_name} reference tree must contain at least two tips")
    reference_alignment = read_fasta_alignment(reference_alignment_path)
    query_alignment = read_fasta_alignment(query_alignment_path)
    require_matching_alignment_widths(reference_alignment, query_alignment)
    reference_by_id = {sequence.sample_id: sequence for sequence in reference_alignment.sequences}
    if set(tree.tip_index) != set(reference_by_id):
        raise ValueError(f"The {panel_name} reference tree tips do not exactly match its alignment IDs")
    overlap = set(sequence.sample_id for sequence in query_alignment.sequences) & set(reference_by_id)
    if overlap:
        raise ValueError(
            f"Query sample_id values already occur in the {panel_name} reference tree: "
            + ", ".join(sorted(overlap)[:10])
        )
    references = tuple(reference_by_id[tip.name or ""] for tip in tree.tips)
    config = placement_config_from_arguments(arguments)
    workers = resolve_threads(arguments.threads, len(query_alignment.sequences))
    logger.event(
        1,
        "placement",
        "Running panel placements",
        panel=panel_name,
        reference_tips=len(references),
        queries=len(query_alignment.sequences),
        columns=reference_alignment.site_count,
        workers=workers,
    )
    started = time.perf_counter()
    if workers == 1:
        results = run_serial_placements(query_alignment.sequences, references, tree, config, logger, started)
    else:
        try:
            results = run_parallel_placements(
                query_alignment.sequences,
                tree.newick(),
                references,
                config,
                arguments.max_cache_megabytes,
                workers,
                logger,
                started,
            )
        except PermissionError as error:
            logger.event(1, "placement", "Parallel workers unavailable; using serial placement", panel=panel_name, error=str(error))
            workers = 1
            results = run_serial_placements(query_alignment.sequences, references, tree, config, logger, started)
    summary = summarize_results(results)
    logger.event(1, "placement", "Completed panel placements", panel=panel_name, **summary)
    metadata = RunMetadata(
        version=__version__,
        tree_path=str(Path(tree_path).resolve()),
        reference_alignment_path=str(Path(reference_alignment_path).resolve()),
        query_alignment_path=str(query_alignment_path.resolve()),
        tree_sha256=file_sha256(tree_path),
        reference_alignment_sha256=file_sha256(reference_alignment_path),
        query_alignment_sha256=file_sha256(query_alignment_path),
        config={**asdict(config), "input_kind": "mafft_add_keeplength_mapped_fasta", "panel": panel_name, "workers": workers},
    )
    write_run_outputs(output_directory, tree, results, metadata, logger.events)
    return PanelRun(panel_name, tree, tuple(results))


def write_selected_panel_trees(
    output_directory: Path,
    classical: PanelRun,
    cryptic: PanelRun,
    decisions: Sequence[object],
) -> None:
    """Write one directly drawable placed tree per selected panel assignment."""

    selected_ids = {
        panel: {
            getattr(decision, "sample_id")
            for decision in decisions
            if getattr(decision, "selected_panel") == panel
        }
        for panel in ("classical", "cryptic")
    }
    trees_directory = output_directory / "selected_trees"
    trees_directory.mkdir()
    for panel in (classical, cryptic):
        selected_results = [result for result in panel.results if result.sample_id in selected_ids[panel.name]]
        (trees_directory / f"{panel.name}.selected_placed_tree.nwk").write_text(
            render_placed_tree(panel.tree, selected_results), encoding="utf-8"
        )


def summarize_decisions(decisions: Sequence[object]) -> dict[str, int]:
    """Count the durable cross-panel assignment labels."""

    counts: dict[str, int] = {}
    for decision in decisions:
        assignment = getattr(decision, "assignment")
        counts[assignment] = counts.get(assignment, 0) + 1
    return {f"queries_{label.lower()}": count for label, count in sorted(counts.items())}


def summarize_clade_assignments(assignments: Sequence[object]) -> dict[str, int]:
    """Count durable fixed-tree clade-assignment status labels."""

    counts: dict[str, int] = {}
    for assignment in assignments:
        status = str(getattr(assignment, "clade_assignment_status"))
        counts[status] = counts.get(status, 0) + 1
    return {
        f"clade_assignments_{status.lower()}": count
        for status, count in sorted(counts.items())
    }


def placement_config_from_arguments(arguments: argparse.Namespace) -> PlacementConfig:
    """Build the immutable core configuration shared by placement commands."""

    return PlacementConfig(
        preset=arguments.preset,
        candidate_count=arguments.candidate_count,
        candidate_hops=arguments.candidate_hops,
        bootstrap_replicates=arguments.bootstrap,
        bootstrap_seed=arguments.seed,
        min_coverage=arguments.min_coverage,
        warn_coverage=arguments.warn_coverage,
        min_comparable_sites=arguments.min_comparable_sites,
        outlier_distance=arguments.outlier_distance,
        unresolved_weight=arguments.unresolved_weight,
        force=arguments.force,
    )


def run_parallel_placements(
    queries: Sequence[AlignedSequence],
    tree_newick: str,
    references: tuple[AlignedSequence, ...],
    config: PlacementConfig,
    max_cache_megabytes: int,
    workers: int,
    logger: EventLogger,
    started_at: float,
):
    """Distribute independent queries while retaining deterministic output order."""

    results: list[object | None] = [None] * len(queries)
    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=worker_initialize,
        initargs=(tree_newick, references, config, max_cache_megabytes),
    ) as executor:
        futures = {
            executor.submit(worker_place, query): index for index, query in enumerate(queries)
        }
        completed = 0
        for future in as_completed(futures):
            index = futures[future]
            results[index] = future.result()
            completed += 1
            log_progress(logger, completed, len(queries), started_at)
    return [result for result in results if result is not None]


def run_serial_placements(
    queries: Sequence[AlignedSequence],
    references: tuple[AlignedSequence, ...],
    tree: Tree,
    config: PlacementConfig,
    logger: EventLogger,
    started_at: float,
):
    """Run placement sequentially when one worker is requested or required."""

    results = []
    for position, query in enumerate(queries, start=1):
        results.append(
            place_query(
                query,
                references,
                tree,
                config,
            )
        )
        log_progress(logger, position, len(queries), started_at)
    return results


def log_progress(
    logger: EventLogger,
    completed: int,
    total: int,
    started_at: float,
    *,
    stage: str = "placement",
) -> None:
    """Emit bounded progress events without flooding logs for large batches."""

    interval = max(1, total // 20)
    if completed != total and completed % interval:
        return
    elapsed = time.perf_counter() - started_at
    rate = completed / elapsed if elapsed else 0.0
    remaining = (total - completed) / rate if rate else None
    logger.event(
        1,
        stage,
        "Progress",
        completed=completed,
        total=total,
        percent=round(completed / total * 100.0, 1),
        queries_per_second=round(rate, 3),
        estimated_remaining_seconds=round(remaining, 1) if remaining is not None else None,
    )


def summarize_results(results: Sequence[object]) -> dict[str, int]:
    """Count final statuses without coupling the CLI to result serialization."""

    counts: dict[str, int] = {}
    for result in results:
        status = getattr(result, "final_status")
        counts[status] = counts.get(status, 0) + 1
    return {f"queries_{status.lower()}": count for status, count in sorted(counts.items())}


def prepare_output_directory(path: str) -> Path:
    """Create a new output directory and refuse implicit overwrites."""

    output = Path(path)
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"Output directory already exists and is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    return output


def resolve_threads(value: str, query_count: int) -> int:
    """Resolve the CLI thread policy without creating idle worker processes."""

    if value == "auto":
        return max(1, min(query_count, max(1, (os.cpu_count() or 2) - 1)))
    try:
        threads = int(value)
    except ValueError as error:
        raise ValueError("--threads must be 'auto' or a positive integer") from error
    if threads < 1:
        raise ValueError("--threads must be a positive integer")
    return min(threads, query_count)


def file_sha256(path: str) -> str:
    """Calculate an input checksum for reproducible run provenance."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def append_final_event(output: Path, event: dict[str, object]) -> None:
    """Append the terminal event after all primary output artifacts are durable."""

    with (output / "events.jsonl").open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(event, sort_keys=True) + "\n")
    extra = {key: value for key, value in event.items() if key not in {"elapsed_seconds", "level", "stage", "message"}}
    with (output / "run.log").open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(
            f"[{event['elapsed_seconds']:8.2f}s] {event['stage']}: {event['message']}"
            + (" " + json.dumps(extra, sort_keys=True) if extra else "")
            + "\n"
        )


def format_event_log(events: Sequence[dict[str, object]]) -> str:
    """Render the human-readable counterpart of the structured event stream."""

    return "".join(
        f"[{event['elapsed_seconds']:8.2f}s] {event['stage']}: {event['message']}"
        + (
            " "
            + json.dumps(
                {
                    key: value
                    for key, value in event.items()
                    if key not in {"elapsed_seconds", "level", "stage", "message"}
                },
                sort_keys=True,
            )
            if len(event) > 4
            else ""
        )
        + "\n"
        for event in events
    )


def positive_integer(value: str) -> int:
    """Argparse validator for positive integers."""

    integer = int(value)
    if integer < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return integer


def nonnegative_integer(value: str) -> int:
    """Argparse validator for non-negative integers."""

    integer = int(value)
    if integer < 0:
        raise argparse.ArgumentTypeError("must be a non-negative integer")
    return integer


def probability(value: str) -> float:
    """Argparse validator for values in the closed unit interval."""

    number = float(value)
    if not 0.0 <= number <= 1.0:
        raise argparse.ArgumentTypeError("must be between 0 and 1")
    return number


def nonnegative_float(value: str) -> float:
    """Argparse validator for a non-negative floating-point threshold."""

    number = float(value)
    if number < 0.0:
        raise argparse.ArgumentTypeError("must be non-negative")
    return number


if __name__ == "__main__":
    raise SystemExit(main())
