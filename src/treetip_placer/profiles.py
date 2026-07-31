"""Reading, projecting, and comparing categorical MLST profiles."""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from .models import Profile, QcFinding, QueryQc


MISSING_VALUES = frozenset({"", ".", "-", "?", "NA", "N", "NULL"})
MIXED_CALL_SEPARATORS = ("/", ";", "|", ",")


@dataclass(frozen=True)
class ProfileTable:
    """A profile table with its declared locus order."""

    loci: tuple[str, ...]
    profiles: tuple[Profile, ...]

    @property
    def by_id(self) -> dict[str, Profile]:
        return {profile.sample_id: profile for profile in self.profiles}


def normalise_allele(value: str | None) -> str | None:
    """Normalize missing tokens while preserving categorical allele identifiers."""

    if value is None:
        return None
    stripped = value.strip()
    return None if stripped.upper() in MISSING_VALUES else stripped


def read_profile_table(path: str | Path, sample_id_column: str = "sample_id") -> ProfileTable:
    """Read a tab-separated MLST profile table and validate its basic schema."""

    table_path = Path(path)
    with table_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Profile table is empty: {table_path}")
        fieldnames = [field.strip() for field in reader.fieldnames]
        if sample_id_column not in fieldnames:
            raise ValueError(
                f"Profile table must contain a '{sample_id_column}' column: {table_path}"
            )
        if len(set(fieldnames)) != len(fieldnames):
            raise ValueError(f"Profile table contains duplicate headers: {table_path}")
        loci = tuple(field for field in fieldnames if field != sample_id_column)
        if not loci:
            raise ValueError(f"Profile table has no locus columns: {table_path}")

        profiles: list[Profile] = []
        seen_ids: set[str] = set()
        for line_number, row in enumerate(reader, start=2):
            sample_id = (row.get(sample_id_column) or "").strip()
            if not sample_id:
                raise ValueError(f"Missing sample_id at {table_path}:{line_number}")
            if sample_id in seen_ids:
                raise ValueError(
                    f"Duplicate sample_id '{sample_id}' at {table_path}:{line_number}"
                )
            seen_ids.add(sample_id)
            alleles = tuple(normalise_allele(row.get(locus)) for locus in loci)
            profiles.append(Profile(sample_id=sample_id, alleles=alleles))

    if not profiles:
        raise ValueError(f"Profile table has no rows: {table_path}")
    return ProfileTable(loci=loci, profiles=tuple(profiles))


def project_queries_to_reference(
    queries: ProfileTable, reference_loci: Sequence[str]
) -> tuple[ProfileTable, tuple[str, ...], tuple[str, ...]]:
    """Project query profiles to reference locus order and reject unknown loci."""

    query_index = {locus: index for index, locus in enumerate(queries.loci)}
    reference_locus_set = set(reference_loci)
    extra_loci = tuple(locus for locus in queries.loci if locus not in reference_locus_set)
    if extra_loci:
        joined = ", ".join(extra_loci)
        raise ValueError(f"Query table contains loci absent from the reference: {joined}")

    absent_loci = tuple(locus for locus in reference_loci if locus not in query_index)
    projected: list[Profile] = []
    for profile in queries.profiles:
        alleles = tuple(
            profile.alleles[query_index[locus]] if locus in query_index else None
            for locus in reference_loci
        )
        projected.append(Profile(sample_id=profile.sample_id, alleles=alleles))
    return ProfileTable(tuple(reference_loci), tuple(projected)), absent_loci, extra_loci


def profile_distance(
    query: Profile,
    reference: Profile,
    locus_indices: Iterable[int] | None = None,
) -> tuple[float | None, int]:
    """Return normalized categorical Hamming distance and comparable locus count."""

    indices = range(len(query.alleles)) if locus_indices is None else locus_indices
    mismatches = 0
    comparable = 0
    for index in indices:
        query_allele = query.alleles[index]
        reference_allele = reference.alleles[index]
        if query_allele is None or reference_allele is None:
            continue
        comparable += 1
        mismatches += query_allele != reference_allele
    if comparable == 0:
        return None, 0
    return mismatches / comparable, comparable


def inspect_profile_tokens(profile: Profile) -> list[QcFinding]:
    """Report allele values that resemble unresolved mixed calls."""

    findings: list[QcFinding] = []
    for locus_index, allele in enumerate(profile.alleles):
        if allele is None:
            continue
        if any(separator in allele for separator in MIXED_CALL_SEPARATORS):
            findings.append(
                QcFinding(
                    code="MIXED_ALLELE_CALL",
                    severity="warning",
                    message=(
                        f"Locus index {locus_index} contains a mixed or composite allele "
                        f"call: '{allele}'."
                    ),
                )
            )
    return findings


def evaluate_query_qc(
    query: Profile,
    reference_profiles: Sequence[Profile],
    *,
    min_coverage: float,
    warn_coverage: float,
    min_comparable_loci: int,
    outlier_distance: float,
    missing_reference_loci: Sequence[str] = (),
) -> QueryQc:
    """Apply deterministic pre-placement QC and identify the nearest reference."""

    total_loci = len(query.alleles)
    observed_loci = sum(allele is not None for allele in query.alleles)
    coverage = observed_loci / total_loci if total_loci else 0.0
    findings = inspect_profile_tokens(query)

    if missing_reference_loci:
        findings.append(
            QcFinding(
                code="REFERENCE_LOCI_ABSENT_FROM_QUERY_TABLE",
                severity="warning",
                message=(
                    "The query table omitted reference loci; those loci are treated as missing: "
                    + ", ".join(missing_reference_loci)
                ),
            )
        )

    nearest_reference: str | None = None
    nearest_distance: float | None = None
    comparable_loci = 0
    for reference in reference_profiles:
        distance, compared = profile_distance(query, reference)
        if distance is None:
            continue
        if nearest_distance is None or distance < nearest_distance:
            nearest_distance = distance
            nearest_reference = reference.sample_id
            comparable_loci = compared

    if coverage < min_coverage:
        findings.append(
            QcFinding(
                code="LOW_COVERAGE",
                severity="error",
                message=(
                    f"Only {observed_loci}/{total_loci} loci are present "
                    f"({coverage:.1%}); the minimum is {min_coverage:.1%}."
                ),
            )
        )
    elif coverage < warn_coverage:
        findings.append(
            QcFinding(
                code="PARTIAL_COVERAGE",
                severity="warning",
                message=(
                    f"Only {observed_loci}/{total_loci} loci are present "
                    f"({coverage:.1%})."
                ),
            )
        )

    if comparable_loci < min_comparable_loci:
        findings.append(
            QcFinding(
                code="INSUFFICIENT_COMPARABLE_LOCI",
                severity="error",
                message=(
                    f"The nearest reference shares {comparable_loci} callable loci; "
                    f"the minimum is {min_comparable_loci}."
                ),
            )
        )
    if nearest_distance is None:
        findings.append(
            QcFinding(
                code="NO_COMPARABLE_REFERENCE",
                severity="error",
                message="No reference profile shares a callable locus with this query.",
            )
        )
    elif nearest_distance > outlier_distance:
        findings.append(
            QcFinding(
                code="DISTANT_FROM_REFERENCE_PANEL",
                severity="warning",
                message=(
                    f"Nearest reference '{nearest_reference}' has allelic distance "
                    f"{nearest_distance:.4f}, above the configured warning threshold "
                    f"{outlier_distance:.4f}."
                ),
            )
        )

    if any(finding.severity == "error" for finding in findings):
        status = "FAIL"
    elif findings:
        status = "WARN"
    else:
        status = "PASS"

    return QueryQc(
        sample_id=query.sample_id,
        status=status,
        coverage=coverage,
        observed_loci=observed_loci,
        comparable_loci=comparable_loci,
        nearest_reference=nearest_reference,
        nearest_distance=nearest_distance,
        findings=findings,
    )


def require_reference_profiles_complete(reference: ProfileTable) -> None:
    """Reject a reference panel that is unusable for a stable MLST metric."""

    incomplete = [
        profile.sample_id
        for profile in reference.profiles
        if any(allele is None for allele in profile.alleles)
    ]
    if incomplete:
        preview = ", ".join(incomplete[:5])
        suffix = "..." if len(incomplete) > 5 else ""
        raise ValueError(
            "Reference profiles must be complete for the MLST MVP. Incomplete samples: "
            f"{preview}{suffix}"
        )


def finite_number(value: float) -> bool:
    """Small helper kept here to make validation intent explicit."""

    return math.isfinite(value)
