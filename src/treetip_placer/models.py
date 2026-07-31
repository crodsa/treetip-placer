"""Immutable data models shared by the placement pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping


@dataclass(frozen=True)
class Profile:
    """An ordered, categorical MLST allele profile."""

    sample_id: str
    alleles: tuple[str | None, ...]


@dataclass(frozen=True)
class QcFinding:
    """A machine-readable QC finding."""

    code: str
    severity: str
    message: str


@dataclass
class QueryQc:
    """Pre-placement QC result for one query."""

    sample_id: str
    status: str
    coverage: float
    observed_loci: int
    comparable_loci: int
    nearest_reference: str | None
    nearest_distance: float | None
    findings: list[QcFinding] = field(default_factory=list)


@dataclass(frozen=True)
class PlacementConfig:
    """Configuration controlling QC, search scope, and support estimation."""

    preset: str = "balanced"
    candidate_count: int = 16
    candidate_hops: int = 4
    bootstrap_replicates: int = 0
    bootstrap_seed: int = 1
    min_coverage: float = 0.70
    warn_coverage: float = 0.90
    min_comparable_loci: int = 3
    outlier_distance: float = 0.50
    unresolved_weight: float = 0.60
    terminal_outlier_multiplier: float = 3.0
    force: bool = False
    top_alternatives: int = 10


@dataclass(frozen=True)
class EdgeFit:
    """A fitted insertion on one reference edge."""

    edge_id: int
    parent_id: int
    child_id: int
    distal_length: float
    pendant_length: float
    objective: float
    rmse: float
    relative_fit_weight: float = 0.0


@dataclass
class PlacementResult:
    """QC plus all placement information for one query."""

    sample_id: str
    qc: QueryQc
    final_status: str
    best_fit: EdgeFit | None = None
    alternatives: list[EdgeFit] = field(default_factory=list)
    bootstrap_support: float | None = None
    search_scope: str | None = None
    note: str | None = None


@dataclass(frozen=True)
class RunMetadata:
    """Provenance fields serialized to the run manifest."""

    version: str
    tree_path: str
    reference_profiles_path: str
    query_profiles_path: str
    tree_sha256: str
    reference_profiles_sha256: str
    query_profiles_sha256: str
    config: Mapping[str, object]
