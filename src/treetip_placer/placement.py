"""Fixed-backbone weighted least-squares placement for categorical matrices."""

from __future__ import annotations

import hashlib
import math
import random
from dataclasses import replace
from statistics import median
from typing import Iterable, Sequence

from .models import AlignedSequence, EdgeFit, PlacementConfig, PlacementResult, QcFinding
from .sequence import character_distance, evaluate_sequence_qc
from .tree import Edge, Tree


def place_query(
    query: AlignedSequence,
    references_in_tip_order: Sequence[AlignedSequence],
    tree: Tree,
    config: PlacementConfig,
) -> PlacementResult:
    """Place one query on a fixed tree and report pre- and post-placement quality."""

    qc = evaluate_sequence_qc(
        query,
        references_in_tip_order,
        min_coverage=config.min_coverage,
        warn_coverage=config.warn_coverage,
        min_comparable_sites=config.min_comparable_sites,
        outlier_distance=config.outlier_distance,
    )
    if qc.status == "FAIL" and not config.force:
        return PlacementResult(
            sample_id=query.sample_id,
            qc=qc,
            final_status="FAIL",
            note="Pre-placement QC failed; use --force only after reviewing the findings.",
        )

    observed = distance_vector(query, references_in_tip_order)
    candidates, scope = choose_candidate_edges(observed, references_in_tip_order, tree, config)
    if not candidates:
        qc.findings.append(
            QcFinding(
                code="NO_CANDIDATE_EDGES",
                severity="error",
                message="The reference tree did not yield candidate edges for placement.",
            )
        )
        return PlacementResult(
            sample_id=query.sample_id,
            qc=qc,
            final_status="FAIL",
            search_scope=scope,
        )

    fits = score_edges(tree, candidates, observed)
    if not fits:
        qc.findings.append(
            QcFinding(
                code="NO_SCORABLE_EDGES",
                severity="error",
                message="No candidate edge could be scored from the available alignment sites.",
            )
        )
        return PlacementResult(
            sample_id=query.sample_id,
            qc=qc,
            final_status="FAIL",
            search_scope=scope,
        )

    weighted_fits = relative_fit_weights(fits)
    best = weighted_fits[0]
    alternatives = weighted_fits[: config.top_alternatives]

    bootstrap_support: float | None = None
    if config.bootstrap_replicates:
        bootstrap_support = bootstrap_edge_support(
            query,
            references_in_tip_order,
            tree,
            candidates,
            best.edge_id,
            config,
        )

    final_status = qc.status
    if best.relative_fit_weight < config.unresolved_weight:
        qc.findings.append(
            QcFinding(
                code="DIFFUSE_DISTANCE_WEIGHT",
                severity="warning",
                message=(
                    f"Best relative-fit weight is {best.relative_fit_weight:.3f}, below "
                    f"the configured resolution threshold {config.unresolved_weight:.3f}. "
                    "This distance heuristic is reported as a warning; requested bootstrap "
                    "support can independently make the placement unresolved."
                ),
            )
        )
        if bootstrap_support is None and final_status == "PASS":
            final_status = "WARN"
    if bootstrap_support is not None and bootstrap_support < config.unresolved_weight:
        final_status = "UNRESOLVED"
        qc.findings.append(
            QcFinding(
                code="UNSTABLE_BOOTSTRAP_PLACEMENT",
                severity="warning",
                message=(
                    f"Bootstrap edge support is {bootstrap_support:.3f}, below the configured "
                    f"resolution threshold {config.unresolved_weight:.3f}."
                ),
            )
        )
    if qc.nearest_distance is not None:
        terminal_limit = max(qc.nearest_distance, 1e-12) * config.terminal_outlier_multiplier
        if best.pendant_length > terminal_limit:
            if final_status == "PASS":
                final_status = "WARN"
            qc.findings.append(
                QcFinding(
                    code="LONG_TERMINAL_BRANCH",
                    severity="warning",
                    message=(
                        f"Fitted terminal branch length {best.pendant_length:.6g} exceeds "
                        f"{config.terminal_outlier_multiplier:.1f}x the nearest-reference "
                        "allelic distance."
                    ),
                )
            )

    return PlacementResult(
        sample_id=query.sample_id,
        qc=qc,
        final_status=final_status,
        best_fit=best,
        alternatives=alternatives,
        bootstrap_support=bootstrap_support,
        search_scope=scope,
    )


def distance_vector(
    query: AlignedSequence,
    references: Sequence[AlignedSequence],
    site_indices: Iterable[int] | None = None,
) -> list[float | None]:
    """Compute query-to-reference categorical distances in tree-tip order."""

    return [character_distance(query, reference, site_indices)[0] for reference in references]


def choose_candidate_edges(
    observed: Sequence[float | None],
    references: Sequence[AlignedSequence],
    tree: Tree,
    config: PlacementConfig,
) -> tuple[list[Edge], str]:
    """Choose a local candidate neighbourhood or the complete edge set."""

    if config.preset == "thorough":
        return tree.edges, "EXHAUSTIVE"

    candidate_count = max(1, config.candidate_count)
    hops = max(1, config.candidate_hops)
    if config.preset == "fast":
        candidate_count = min(candidate_count, 8)
        hops = min(hops, 2)

    nearest = sorted(
        (
            (distance, reference.sample_id)
            for distance, reference in zip(observed, references, strict=True)
            if distance is not None
        ),
        key=lambda item: (item[0], item[1]),
    )[:candidate_count]
    edges = tree.nearby_edges((name for _, name in nearest), hops)
    if edges:
        return edges, "LOCAL"
    return tree.edges, "EXHAUSTIVE_FALLBACK"


def score_edges(tree: Tree, edges: Sequence[Edge], observed: Sequence[float | None]) -> list[EdgeFit]:
    """Fit every supplied edge and return fits ordered by the distance objective."""

    valid_indices = [index for index, distance in enumerate(observed) if distance is not None]
    if len(valid_indices) < 2:
        return []
    fits = [fit_edge(tree, edge, observed, valid_indices) for edge in edges]
    return sorted((fit for fit in fits if fit is not None), key=lambda fit: fit.objective)


def fit_edge(
    tree: Tree,
    edge: Edge,
    observed: Sequence[float | None],
    valid_indices: Sequence[int],
) -> EdgeFit | None:
    """Constrained least-squares fit of an insertion on one branch.

    The query joins the branch at distal distance ``x`` from the Newick parent and
    has terminal length ``p``. Both values are constrained to be non-negative and
    x is constrained to the edge length.
    """

    parent_distances = tree.tip_distances(edge.parent_id)
    child_distances = tree.tip_distances(edge.child_id)
    child_tip_start, child_tip_end = tree.tip_span(edge.child_id)

    count = 0
    sum_y = 0.0
    sum_signed_y = 0.0
    sum_y_squared = 0.0
    signed_count = 0.0
    for index in valid_indices:
        distance = observed[index]
        if distance is None:
            continue
        if child_tip_start <= index < child_tip_end:
            base = child_distances[index] + edge.length
            sign = -1.0
        else:
            base = parent_distances[index]
            sign = 1.0
        residual = distance - base
        count += 1
        sum_y += residual
        sum_signed_y += sign * residual
        sum_y_squared += residual * residual
        signed_count += sign

    if count < 2:
        return None
    denominator = count * count - signed_count * signed_count
    if denominator <= 1e-12:
        return None

    unconstrained_pendant = (count * sum_y - signed_count * sum_signed_y) / denominator
    unconstrained_distal = (count * sum_signed_y - signed_count * sum_y) / denominator

    candidates: list[tuple[float, float]] = []
    if unconstrained_pendant >= 0.0 and 0.0 <= unconstrained_distal <= edge.length:
        candidates.append((unconstrained_pendant, unconstrained_distal))

    for distal in (0.0, edge.length):
        pendant = max(0.0, (sum_y - distal * signed_count) / count)
        candidates.append((pendant, distal))
    pendant_zero_distal = min(max(sum_signed_y / count, 0.0), edge.length)
    candidates.append((0.0, pendant_zero_distal))

    pendant, distal = min(
        candidates,
        key=lambda values: squared_error(
            values[0],
            values[1],
            count,
            sum_y,
            sum_signed_y,
            sum_y_squared,
            signed_count,
        ),
    )
    objective = max(
        0.0,
        squared_error(
            pendant,
            distal,
            count,
            sum_y,
            sum_signed_y,
            sum_y_squared,
            signed_count,
        ),
    )
    return EdgeFit(
        edge_id=edge.edge_id,
        parent_id=edge.parent_id,
        child_id=edge.child_id,
        distal_length=distal,
        pendant_length=pendant,
        objective=objective,
        rmse=math.sqrt(objective / count),
    )


def squared_error(
    pendant: float,
    distal: float,
    count: int,
    sum_y: float,
    sum_signed_y: float,
    sum_y_squared: float,
    signed_count: float,
) -> float:
    """Evaluate the least-squares objective from sufficient statistics."""

    return (
        sum_y_squared
        - 2.0 * pendant * sum_y
        - 2.0 * distal * sum_signed_y
        + count * pendant * pendant
        + 2.0 * pendant * distal * signed_count
        + count * distal * distal
    )


def relative_fit_weights(fits: Sequence[EdgeFit]) -> list[EdgeFit]:
    """Transform objective differences into explicitly heuristic relative weights."""

    if not fits:
        return []
    best = fits[0]
    positive_objectives = [fit.objective for fit in fits if fit.objective > 0.0]
    scale = max(best.rmse * best.rmse, median(positive_objectives) if positive_objectives else 0.0, 1e-12)
    raw_weights = [math.exp(-0.5 * (fit.objective - best.objective) / scale) for fit in fits]
    total = sum(raw_weights)
    return [replace(fit, relative_fit_weight=weight / total) for fit, weight in zip(fits, raw_weights, strict=True)]


def bootstrap_edge_support(
    query: AlignedSequence,
    references: Sequence[AlignedSequence],
    tree: Tree,
    candidates: Sequence[Edge],
    expected_best_edge_id: int,
    config: PlacementConfig,
) -> float:
    """Estimate fixed-candidate edge stability by resampling sites with replacement."""

    stable_seed = int.from_bytes(
        hashlib.sha256(f"{config.bootstrap_seed}:{query.sample_id}".encode("utf-8")).digest()[:8],
        "big",
    )
    generator = random.Random(stable_seed)
    site_count = len(query.characters)
    successful = 0
    matching = 0
    for _ in range(config.bootstrap_replicates):
        sampled_sites = [generator.randrange(site_count) for _ in range(site_count)]
        sampled_observed = distance_vector(query, references, sampled_sites)
        fits = score_edges(tree, candidates, sampled_observed)
        if not fits:
            continue
        successful += 1
        if fits[0].edge_id == expected_best_edge_id:
            matching += 1
    return matching / successful if successful else 0.0
