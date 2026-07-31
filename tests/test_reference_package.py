"""Regression tests for portable paired reference package loading."""

from __future__ import annotations

import hashlib
import json
import shutil
import unittest
from pathlib import Path

from treetip_placer.reference_package import load_paired_reference_package


def digest(path: Path) -> str:
    """Calculate the test artifact checksum using the package's contract."""

    return hashlib.sha256(path.read_bytes()).hexdigest()


class ReferencePackageTests(unittest.TestCase):
    def test_loads_relative_paths_and_validates_checksums(self) -> None:
        root = Path.cwd() / ".treetip-reference-package-test"
        if root.exists():
            shutil.rmtree(root)
        root.mkdir()
        try:
            classical_tree = root / "classical.nwk"
            classical_alignment = root / "classical.fasta"
            cryptic_tree = root / "cryptic.nwk"
            cryptic_alignment = root / "cryptic.fasta"
            for path, contents in (
                (classical_tree, "(A:1,B:1);\n"),
                (classical_alignment, ">A\nAAAA\n>B\nAAAT\n"),
                (cryptic_tree, "(C:1,D:1);\n"),
                (cryptic_alignment, ">C\nCCCC\n>D\nCCCT\n"),
            ):
                path.write_text(contents, encoding="utf-8")
            manifest = {
                "schema_version": 1,
                "minimum_distance_separation": 0.01,
                "classical": {
                    "tree": classical_tree.name,
                    "alignment": classical_alignment.name,
                    "tree_sha256": digest(classical_tree),
                    "alignment_sha256": digest(classical_alignment),
                },
                "cryptic": {
                    "tree": cryptic_tree.name,
                    "alignment": cryptic_alignment.name,
                    "tree_sha256": digest(cryptic_tree),
                    "alignment_sha256": digest(cryptic_alignment),
                },
            }
            manifest_path = root / "paired_reference.json"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

            package = load_paired_reference_package(manifest_path)

            self.assertEqual(package.minimum_distance_separation, 0.01)
            self.assertEqual(package.classical.tree_path, classical_tree)
            self.assertEqual(package.cryptic.alignment_path, cryptic_alignment)
        finally:
            if root.exists():
                shutil.rmtree(root)

    def test_rejects_checksum_mismatch(self) -> None:
        root = Path.cwd() / ".treetip-reference-package-test"
        if root.exists():
            shutil.rmtree(root)
        root.mkdir()
        try:
            tree = root / "tree.nwk"
            alignment = root / "alignment.fasta"
            tree.write_text("(A:1,B:1);\n", encoding="utf-8")
            alignment.write_text(">A\nAAAA\n>B\nAAAT\n", encoding="utf-8")
            manifest = {
                "schema_version": 1,
                "minimum_distance_separation": 0.01,
                "classical": {"tree": tree.name, "alignment": alignment.name, "tree_sha256": "0" * 64},
                "cryptic": {"tree": tree.name, "alignment": alignment.name},
            }
            manifest_path = root / "paired_reference.json"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "checksum"):
                load_paired_reference_package(manifest_path)
        finally:
            if root.exists():
                shutil.rmtree(root)


if __name__ == "__main__":
    unittest.main()
