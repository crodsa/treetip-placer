"""Portable paired-reference package loading and checksum validation."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


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
    return PairedReferencePackage(manifest_path, classical, cryptic, threshold)


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


def sha256(path: Path) -> str:
    """Return a streaming SHA-256 digest for a local reference artifact."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
