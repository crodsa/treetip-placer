"""Command-line interface for the TreeTip Placer MVP."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path
from typing import Sequence

from . import __version__
from .alignment import read_fasta_alignment, require_matching_alignment_widths
from .models import PlacementConfig, Profile, RunMetadata
from .outputs import write_run_outputs
from .placement import place_query
from .profiles import (
    ProfileTable,
    project_queries_to_reference,
    read_profile_table,
    require_reference_profiles_complete,
)
from .tree import Tree
from .validation import (
    run_parallel_validation,
    run_serial_validation,
    select_holdout_ids,
    summarize_validation,
    write_validation_outputs,
)


_WORKER_TREE: Tree | None = None
_WORKER_REFERENCES: tuple[Profile, ...] = ()
_WORKER_CONFIG: PlacementConfig | None = None
_WORKER_MISSING_LOCI: tuple[str, ...] = ()


def worker_initialize(
    newick: str,
    references: tuple[Profile, ...],
    config: PlacementConfig,
    missing_loci: tuple[str, ...],
    max_cache_megabytes: int,
) -> None:
    """Create isolated worker state so caches do not need to be pickled per query."""

    global _WORKER_TREE, _WORKER_REFERENCES, _WORKER_CONFIG, _WORKER_MISSING_LOCI
    _WORKER_TREE = Tree.from_newick(newick, max_cache_megabytes)
    _WORKER_REFERENCES = references
    _WORKER_CONFIG = config
    _WORKER_MISSING_LOCI = missing_loci


def worker_place(query: Profile):
    """Place a single query in a process worker."""

    if _WORKER_TREE is None or _WORKER_CONFIG is None:
        raise RuntimeError("Placement worker was not initialized")
    return place_query(
        query,
        _WORKER_REFERENCES,
        _WORKER_TREE,
        _WORKER_CONFIG,
        missing_reference_loci=_WORKER_MISSING_LOCI,
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
        description="Quality-controlled fixed-tree placement for MLST profiles and aligned characters.",
    )
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    place = subparsers.add_parser(
        "place", help="Place query profiles or aligned sequences on a fixed reference tree."
    )
    place.add_argument("--tree", required=True, help="Reference branch-length Newick tree.")
    input_group = place.add_argument_group("input adapters")
    input_group.add_argument(
        "--reference-profiles", help="Reference MLST allele-profile TSV file."
    )
    input_group.add_argument("--query-profiles", help="Query MLST allele-profile TSV file.")
    input_group.add_argument(
        "--reference-alignment", help="Reference fixed-column FASTA alignment."
    )
    input_group.add_argument("--query-alignment", help="Query fixed-column FASTA alignment.")
    place.add_argument("--output", required=True, help="New, empty output directory.")
    place.add_argument(
        "--preset",
        choices=("fast", "balanced", "thorough"),
        default="balanced",
        help="Search scope: fast local, balanced local, or exhaustive thorough search.",
    )
    place.add_argument(
        "--threads",
        default="auto",
        help="Parallel query workers; use 'auto' or a positive integer.",
    )
    place.add_argument(
        "--bootstrap",
        type=nonnegative_integer,
        default=0,
        help="Locus/site bootstrap replicates per successfully placed query.",
    )
    place.add_argument("--seed", type=int, default=1, help="Bootstrap random seed.")
    place.add_argument("--candidate-count", type=positive_integer, default=16)
    place.add_argument("--candidate-hops", type=positive_integer, default=4)
    place.add_argument("--min-coverage", type=probability, default=0.70)
    place.add_argument("--warn-coverage", type=probability, default=0.90)
    place.add_argument("--min-comparable-loci", type=positive_integer, default=3)
    place.add_argument("--outlier-distance", type=probability, default=0.50)
    place.add_argument("--unresolved-weight", type=probability, default=0.60)
    place.add_argument(
        "--force",
        action="store_true",
        help="Attempt placement after a pre-placement QC failure; failures remain annotated.",
    )
    place.add_argument(
        "--max-cache-megabytes",
        type=positive_integer,
        default=256,
        help="Per-process limit for cached node-to-tip distance rows.",
    )
    place.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase stderr progress detail; repeat for more detail.",
    )
    validate = subparsers.add_parser(
        "validate-mlst",
        help="Run pruned leave-one-out validation for an MLST reference panel.",
    )
    validate.add_argument("--tree", required=True, help="Reference branch-length Newick tree.")
    validate.add_argument(
        "--reference-profiles", required=True, help="Reference MLST allele-profile TSV file."
    )
    validate.add_argument("--output", required=True, help="New, empty validation output directory.")
    validate.add_argument(
        "--max-holdouts",
        type=nonnegative_integer,
        default=0,
        help="Maximum profiles to hold out; 0 validates every reference profile.",
    )
    validate.add_argument(
        "--preset",
        choices=("fast", "balanced", "thorough"),
        default="balanced",
        help="Search scope for each held-out placement.",
    )
    validate.add_argument(
        "--threads", default="auto", help="Parallel holdout workers; use 'auto' or a positive integer."
    )
    validate.add_argument(
        "--bootstrap",
        type=nonnegative_integer,
        default=0,
        help="Locus bootstrap replicates per held-out profile.",
    )
    validate.add_argument("--seed", type=int, default=1, help="Sampling and bootstrap random seed.")
    validate.add_argument("--candidate-count", type=positive_integer, default=16)
    validate.add_argument("--candidate-hops", type=positive_integer, default=4)
    validate.add_argument("--min-coverage", type=probability, default=0.70)
    validate.add_argument("--warn-coverage", type=probability, default=0.90)
    validate.add_argument("--min-comparable-loci", type=positive_integer, default=3)
    validate.add_argument("--outlier-distance", type=probability, default=0.50)
    validate.add_argument("--unresolved-weight", type=probability, default=0.60)
    validate.add_argument(
        "--max-cache-megabytes",
        type=positive_integer,
        default=256,
        help="Per-process limit for cached node-to-tip distance rows.",
    )
    validate.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase stderr progress detail; repeat for more detail.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command line and convert expected errors to a clear exit status."""

    parser = build_parser()
    arguments = parser.parse_args(argv)
    try:
        if arguments.command == "place":
            return run_place(arguments)
        if arguments.command == "validate-mlst":
            return run_validate_mlst(arguments)
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

    references, queries, input_kind, absent_loci = load_input_adapter(arguments)
    reference_by_id = references.by_id
    tree_tip_names = set(tree.tip_index)
    reference_names = set(reference_by_id)
    if tree_tip_names != reference_names:
        missing_in_profiles = sorted(tree_tip_names - reference_names)
        missing_in_tree = sorted(reference_names - tree_tip_names)
        details: list[str] = []
        if missing_in_profiles:
            details.append("tips absent from reference profiles: " + ", ".join(missing_in_profiles[:10]))
        if missing_in_tree:
            details.append("profile IDs absent from tree: " + ", ".join(missing_in_tree[:10]))
        raise ValueError("Reference tree and data IDs differ; " + "; ".join(details))
    query_names = {profile.sample_id for profile in queries}
    overlap = sorted(query_names & reference_names)
    if overlap:
        raise ValueError(
            "Query sample_id values already occur in the reference tree: " + ", ".join(overlap[:10])
        )

    references_in_tip_order = tuple(
        reference_by_id[tip.name or ""] for tip in tree.tips
    )
    config = PlacementConfig(
        preset=arguments.preset,
        candidate_count=arguments.candidate_count,
        candidate_hops=arguments.candidate_hops,
        bootstrap_replicates=arguments.bootstrap,
        bootstrap_seed=arguments.seed,
        min_coverage=arguments.min_coverage,
        warn_coverage=arguments.warn_coverage,
        min_comparable_loci=arguments.min_comparable_loci,
        outlier_distance=arguments.outlier_distance,
        unresolved_weight=arguments.unresolved_weight,
        force=arguments.force,
    )
    logger.event(
        1,
        "input",
        "Validated reference inputs",
        input_kind=input_kind,
        reference_tips=len(references_in_tip_order),
        queries=len(queries),
        characters=len(references_in_tip_order[0].alleles),
    )

    workers = resolve_threads(arguments.threads, len(queries))
    logger.event(1, "placement", "Running placements", workers=workers, preset=config.preset)
    started_placement = time.perf_counter()
    if workers == 1:
        results = run_serial_placements(
            queries, references_in_tip_order, tree, config, absent_loci, logger, started_placement
        )
    else:
        try:
            results = run_parallel_placements(
                queries,
                tree.newick(),
                references_in_tip_order,
                config,
                absent_loci,
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
                queries, references_in_tip_order, tree, config, absent_loci, logger, started_placement
            )

    summary = summarize_results(results)
    logger.event(1, "placement", "Completed placement", **summary)
    metadata = RunMetadata(
        version=__version__,
        tree_path=str(Path(arguments.tree).resolve()),
        reference_profiles_path=str(
            Path(arguments.reference_profiles or arguments.reference_alignment).resolve()
        ),
        query_profiles_path=str(Path(arguments.query_profiles or arguments.query_alignment).resolve()),
        tree_sha256=file_sha256(arguments.tree),
        reference_profiles_sha256=file_sha256(
            arguments.reference_profiles or arguments.reference_alignment
        ),
        query_profiles_sha256=file_sha256(arguments.query_profiles or arguments.query_alignment),
        config={**asdict(config), "input_kind": input_kind, "workers": workers},
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


def run_validate_mlst(arguments: argparse.Namespace) -> int:
    """Hold out known MLST references to calibrate placement and QC behavior."""

    logger = EventLogger(arguments.verbose)
    logger.event(1, "start", "Starting MLST leave-one-out validation", version=__version__)
    output = prepare_output_directory(arguments.output)
    tree = Tree.from_path(arguments.tree, arguments.max_cache_megabytes)
    if len(tree.tips) < 4:
        raise ValueError("MLST leave-one-out validation requires at least four reference tips")
    references = read_profile_table(arguments.reference_profiles)
    require_reference_profiles_complete(references)
    tree_tip_names = set(tree.tip_index)
    reference_by_id = references.by_id
    reference_names = set(reference_by_id)
    if tree_tip_names != reference_names:
        missing_in_profiles = sorted(tree_tip_names - reference_names)
        missing_in_tree = sorted(reference_names - tree_tip_names)
        details: list[str] = []
        if missing_in_profiles:
            details.append("tips absent from reference profiles: " + ", ".join(missing_in_profiles[:10]))
        if missing_in_tree:
            details.append("profile IDs absent from tree: " + ", ".join(missing_in_tree[:10]))
        raise ValueError("Reference tree and data IDs differ; " + "; ".join(details))
    references_in_tip_order = tuple(
        reference_by_id[tip.name or ""] for tip in tree.tips
    )
    holdout_ids = select_holdout_ids(
        references_in_tip_order, arguments.max_holdouts, arguments.seed
    )
    config = PlacementConfig(
        preset=arguments.preset,
        candidate_count=arguments.candidate_count,
        candidate_hops=arguments.candidate_hops,
        bootstrap_replicates=arguments.bootstrap,
        bootstrap_seed=arguments.seed,
        min_coverage=arguments.min_coverage,
        warn_coverage=arguments.warn_coverage,
        min_comparable_loci=arguments.min_comparable_loci,
        outlier_distance=arguments.outlier_distance,
        unresolved_weight=arguments.unresolved_weight,
    )
    workers = resolve_threads(arguments.threads, len(holdout_ids))
    logger.event(
        1,
        "validation",
        "Validated MLST reference panel",
        reference_tips=len(references_in_tip_order),
        holdouts=len(holdout_ids),
        loci=len(references.loci),
        workers=workers,
        preset=config.preset,
    )
    started_validation = time.perf_counter()

    def validation_progress(completed: int, total: int) -> None:
        log_progress(logger, completed, total, started_validation, stage="validation")

    if workers == 1:
        records = run_serial_validation(
            holdout_ids,
            tree,
            references_in_tip_order,
            config,
            validation_progress,
        )
    else:
        try:
            records = run_parallel_validation(
                holdout_ids,
                tree.newick(),
                references_in_tip_order,
                config,
                arguments.max_cache_megabytes,
                workers,
                validation_progress,
            )
        except PermissionError as error:
            logger.event(
                1,
                "validation",
                "Parallel workers are unavailable in this environment; falling back to serial validation",
                error=str(error),
            )
            workers = 1
            records = run_serial_validation(
                holdout_ids,
                tree,
                references_in_tip_order,
                config,
                validation_progress,
            )
    summary = summarize_validation(records)
    logger.event(1, "validation", "Completed MLST leave-one-out validation", **summary)
    manifest = {
        "version": __version__,
        "validation_method": "pruned_leave_one_out_mlst",
        "tree_path": str(Path(arguments.tree).resolve()),
        "reference_profiles_path": str(Path(arguments.reference_profiles).resolve()),
        "tree_sha256": file_sha256(arguments.tree),
        "reference_profiles_sha256": file_sha256(arguments.reference_profiles),
        "config": {**asdict(config), "workers": workers, "max_holdouts": arguments.max_holdouts},
    }
    logger.event(1, "output", "Writing validation outputs", directory=str(output))
    logger.event(1, "complete", "Validation run completed", **summary)
    write_validation_outputs(output, records, summary, manifest, logger.events)
    (output / "run.log").write_text(format_event_log(logger.events), encoding="utf-8")
    return 0


def load_input_adapter(
    arguments: argparse.Namespace,
) -> tuple[ProfileTable, tuple[Profile, ...], str, tuple[str, ...]]:
    """Load either categorical MLST tables or fixed-column FASTA alignments."""

    profile_arguments = (arguments.reference_profiles, arguments.query_profiles)
    alignment_arguments = (arguments.reference_alignment, arguments.query_alignment)
    using_profiles = all(profile_arguments)
    using_alignment = all(alignment_arguments)
    if using_profiles == using_alignment:
        raise ValueError(
            "Provide exactly one complete input pair: --reference-profiles/--query-profiles "
            "or --reference-alignment/--query-alignment."
        )
    if using_profiles:
        references = read_profile_table(arguments.reference_profiles)
        require_reference_profiles_complete(references)
        raw_queries = read_profile_table(arguments.query_profiles)
        projected_queries, absent_loci, _ = project_queries_to_reference(raw_queries, references.loci)
        return references, projected_queries.profiles, "mlst_allele_profiles", absent_loci

    reference_alignment = read_fasta_alignment(arguments.reference_alignment)
    query_alignment = read_fasta_alignment(arguments.query_alignment)
    require_matching_alignment_widths(reference_alignment, query_alignment)
    return (
        ProfileTable(loci=(), profiles=reference_alignment.profiles),
        query_alignment.profiles,
        "fixed_column_fasta_alignment",
        (),
    )


def run_parallel_placements(
    queries: Sequence[Profile],
    tree_newick: str,
    references: tuple[Profile, ...],
    config: PlacementConfig,
    absent_loci: tuple[str, ...],
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
        initargs=(tree_newick, references, config, absent_loci, max_cache_megabytes),
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
    queries: Sequence[Profile],
    references: tuple[Profile, ...],
    tree: Tree,
    config: PlacementConfig,
    absent_loci: tuple[str, ...],
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
                missing_reference_loci=absent_loci,
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


if __name__ == "__main__":
    raise SystemExit(main())
