"""Portable paired-reference package loading and checksum validation."""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .alignment import read_fasta_alignment
from .tree import Tree


@dataclass(frozen=True)
class PanelReference:
    """Immutable files that define one placement reference panel."""

    tree_path: Path
    alignment_path: Path
    tree_sha256: str | None
    alignment_sha256: str | None


@dataclass(frozen=True)
class PairedReferencePackage:
    """A classical/cryptic pair and its calibrated selection threshold."""

    manifest_path: Path | None
    classical: PanelReference
    cryptic: PanelReference
    minimum_distance_separation: float


def load_paired_reference_package(path: str | Path) -> PairedReferencePackage:
    """Load a schema-v1 paired package and validate declared file checksums."""

    manifest_path = Path(path).resolve()
    try:
        document = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"Reference package is not valid JSON: {manifest_path}") from error
    if not isinstance(document, dict) or document.get("schema_version") != 1:
        raise ValueError("Reference package must declare schema_version 1")
    try:
        classical = parse_panel_reference(document["classical"], manifest_path.parent, "classical")
        cryptic = parse_panel_reference(document["cryptic"], manifest_path.parent, "cryptic")
        threshold = float(document["minimum_distance_separation"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"Reference package has an invalid required field: {error}") from error
    if threshold < 0.0:
        raise ValueError("Reference package minimum_distance_separation must be non-negative")
    validate_panel_checksums(classical, "classical")
    validate_panel_checksums(cryptic, "cryptic")
    package = PairedReferencePackage(manifest_path, classical, cryptic, threshold)
    validate_paired_reference_package(package)
    return package


def create_paired_reference_package(
    output_path: str | Path,
    *,
    package_name: str,
    classical_tree: str | Path,
    classical_alignment: str | Path,
    cryptic_tree: str | Path,
    cryptic_alignment: str | Path,
    minimum_distance_separation: float,
) -> PairedReferencePackage:
    """Create a new manifest with portable paths when files share a drive."""

    if minimum_distance_separation < 0.0:
        raise ValueError("minimum_distance_separation must be non-negative")
    manifest_path = Path(output_path).resolve()
    if manifest_path.exists():
        raise ValueError(f"Reference package output already exists: {manifest_path}")
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    classical = build_panel_reference(classical_tree, classical_alignment)
    cryptic = build_panel_reference(cryptic_tree, cryptic_alignment)
    validate_panel_structure(classical, "classical")
    validate_panel_structure(cryptic, "cryptic")
    document = {
        "schema_version": 1,
        "package_name": package_name,
        "minimum_distance_separation": minimum_distance_separation,
        "classical": serialize_panel_reference(classical, manifest_path.parent),
        "cryptic": serialize_panel_reference(cryptic, manifest_path.parent),
    }
    manifest_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return load_paired_reference_package(manifest_path)


def build_panel_reference(tree_path: str | Path, alignment_path: str | Path) -> PanelReference:
    """Resolve existing artifacts and calculate their immutable checksums."""

    tree = Path(tree_path).resolve()
    alignment = Path(alignment_path).resolve()
    if not tree.is_file():
        raise ValueError(f"Reference tree does not exist: {tree}")
    if not alignment.is_file():
        raise ValueError(f"Reference alignment does not exist: {alignment}")
    return PanelReference(tree, alignment, sha256(tree), sha256(alignment))


def serialize_panel_reference(panel: PanelReference, manifest_directory: Path) -> dict[str, str]:
    """Serialize one panel with relative paths when the platform permits it."""

    return {
        "tree": portable_path(panel.tree_path, manifest_directory),
        "alignment": portable_path(panel.alignment_path, manifest_directory),
        "tree_sha256": panel.tree_sha256 or sha256(panel.tree_path),
        "alignment_sha256": panel.alignment_sha256 or sha256(panel.alignment_path),
    }


def portable_path(path: Path, manifest_directory: Path) -> str:
    """Return a relative POSIX-style path, falling back to absolute across drives."""

    try:
        return Path(os.path.relpath(path, manifest_directory)).as_posix()
    except ValueError:
        return path.as_posix()


def parse_panel_reference(
    value: Any, base_directory: Path, panel_name: str
) -> PanelReference:
    """Parse one panel object, resolving relative paths beside the manifest."""

    if not isinstance(value, dict):
        raise ValueError(f"{panel_name} must be an object")
    tree_path = resolve_package_path(value["tree"], base_directory, f"{panel_name}.tree")
    alignment_path = resolve_package_path(
        value["alignment"], base_directory, f"{panel_name}.alignment"
    )
    tree_sha256 = optional_sha256(value.get("tree_sha256"), f"{panel_name}.tree_sha256")
    alignment_sha256 = optional_sha256(
        value.get("alignment_sha256"), f"{panel_name}.alignment_sha256"
    )
    return PanelReference(tree_path, alignment_path, tree_sha256, alignment_sha256)


def resolve_package_path(value: Any, base_directory: Path, field_name: str) -> Path:
    """Resolve a package path and ensure the referenced file exists."""

    if not isinstance(value, str) or not value:
        raise ValueError(f"{field_name} must be a non-empty path string")
    candidate = Path(value)
    resolved = (base_directory / candidate).resolve() if not candidate.is_absolute() else candidate.resolve()
    if not resolved.is_file():
        raise ValueError(f"Reference package file does not exist: {resolved}")
    return resolved


def optional_sha256(value: Any, field_name: str) -> str | None:
    """Validate an optional lowercase SHA-256 checksum."""

    if value is None:
        return None
    if not isinstance(value, str) or len(value) != 64 or any(character not in "0123456789abcdef" for character in value):
        raise ValueError(f"{field_name} must be a lowercase SHA-256 hex string")
    return value


def validate_panel_checksums(panel: PanelReference, panel_name: str) -> None:
    """Reject silently modified package files when a checksum is declared."""

    for file_path, expected, artifact in (
        (panel.tree_path, panel.tree_sha256, "tree"),
        (panel.alignment_path, panel.alignment_sha256, "alignment"),
    ):
        if expected is not None and sha256(file_path) != expected:
            raise ValueError(f"{panel_name} {artifact} checksum does not match its reference package")


def validate_paired_reference_package(package: PairedReferencePackage) -> dict[str, dict[str, int]]:
    """Validate each panel's tree/alignment identity contract and return a summary."""

    return {
        "classical": validate_panel_structure(package.classical, "classical"),
        "cryptic": validate_panel_structure(package.cryptic, "cryptic"),
    }


def validate_panel_structure(panel: PanelReference, panel_name: str) -> dict[str, int]:
    """Reject a panel whose immutable tree and alignment do not define the same tips."""

    tree = Tree.from_path(panel.tree_path)
    alignment = read_fasta_alignment(panel.alignment_path)
    if len(tree.tips) < 2:
        raise ValueError(f"{panel_name} reference tree must contain at least two tips")
    tree_ids = set(tree.tip_index)
    alignment_ids = {sequence.sample_id for sequence in alignment.sequences}
    if tree_ids != alignment_ids:
        missing_in_alignment = sorted(tree_ids - alignment_ids)
        missing_in_tree = sorted(alignment_ids - tree_ids)
        details: list[str] = []
        if missing_in_alignment:
            details.append("tips absent from alignment: " + ", ".join(missing_in_alignment[:10]))
        if missing_in_tree:
            details.append("alignment IDs absent from tree: " + ", ".join(missing_in_tree[:10]))
        raise ValueError(f"{panel_name} tree and alignment IDs differ; " + "; ".join(details))
    return {"tips": len(tree.tips), "alignment_records": len(alignment.sequences), "columns": alignment.site_count}


def sha256(path: Path) -> str:
    """Return a streaming SHA-256 digest for a local reference artifact."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
