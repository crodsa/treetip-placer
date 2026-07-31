"""Distance and quality control for fixed-column nucleotide alignments."""

from __future__ import annotations

from typing import Iterable, Sequence

from .models import AlignedSequence, QcFinding, QueryQc


AMBIGUOUS_BASES = frozenset({"R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"})


def character_distance(
    query: AlignedSequence,
    reference: AlignedSequence,
    site_indices: Iterable[int] | None = None,
) -> tuple[float | None, int]:
    """Return normalized Hamming distance over pairwise callable alignment sites."""

    indices = range(len(query.characters)) if site_indices is None else site_indices
    mismatches = 0
    comparable = 0
    for index in indices:
        query_base = query.characters[index]
        reference_base = reference.characters[index]
        if query_base is None or reference_base is None:
            continue
        comparable += 1
        mismatches += query_base != reference_base
    if comparable == 0:
        return None, 0
    return mismatches / comparable, comparable


def evaluate_sequence_qc(
    query: AlignedSequence,
    references: Sequence[AlignedSequence],
    *,
    min_coverage: float,
    warn_coverage: float,
    min_comparable_sites: int,
    outlier_distance: float,
) -> QueryQc:
    """Apply deterministic QC to one query in a fixed alignment coordinate system."""

    total_sites = len(query.characters)
    observed_sites = sum(base is not None for base in query.characters)
    coverage = observed_sites / total_sites if total_sites else 0.0
    findings: list[QcFinding] = []
    ambiguous_count = sum(base in AMBIGUOUS_BASES for base in query.characters if base is not None)
    if ambiguous_count:
        findings.append(
            QcFinding(
                code="AMBIGUOUS_NUCLEOTIDE",
                severity="warning",
                message=f"Query contains {ambiguous_count} ambiguous IUPAC nucleotide sites.",
            )
        )

    nearest_reference: str | None = None
    nearest_distance: float | None = None
    comparable_sites = 0
    for reference in references:
        distance, compared = character_distance(query, reference)
        if distance is None:
            continue
        if nearest_distance is None or distance < nearest_distance:
            nearest_distance = distance
            nearest_reference = reference.sample_id
            comparable_sites = compared

    if coverage < min_coverage:
        findings.append(QcFinding("LOW_COVERAGE", "error", f"Only {observed_sites}/{total_sites} alignment sites are callable ({coverage:.1%}); the minimum is {min_coverage:.1%}."))
    elif coverage < warn_coverage:
        findings.append(QcFinding("PARTIAL_COVERAGE", "warning", f"Only {observed_sites}/{total_sites} alignment sites are callable ({coverage:.1%})."))
    if comparable_sites < min_comparable_sites:
        findings.append(QcFinding("INSUFFICIENT_COMPARABLE_SITES", "error", f"The nearest reference shares {comparable_sites} callable sites; the minimum is {min_comparable_sites}."))
    if nearest_distance is None:
        findings.append(QcFinding("NO_COMPARABLE_REFERENCE", "error", "No reference sequence shares a callable site with this query."))
    elif nearest_distance > outlier_distance:
        findings.append(QcFinding("DISTANT_FROM_REFERENCE_PANEL", "warning", f"Nearest reference '{nearest_reference}' has normalized distance {nearest_distance:.4f}, above the configured warning threshold {outlier_distance:.4f}."))

    status = "FAIL" if any(finding.severity == "error" for finding in findings) else "WARN" if findings else "PASS"
    return QueryQc(query.sample_id, status, coverage, observed_sites, comparable_sites, nearest_reference, nearest_distance, findings)
