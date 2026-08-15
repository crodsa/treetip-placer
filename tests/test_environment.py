"""Tests for the managed MAFFT environment contract."""

from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from treetip_placer.environment import (
    MAFFT_VERSION,
    MafftResolution,
    configuration_path,
    default_config_directory,
    default_data_directory,
    find_environment_manager,
    load_managed_mafft,
    managed_environment_prefix,
    resolve_mafft_for_classification,
    save_managed_mafft,
    setup_mafft_environment,
    verify_mafft,
)


def successful_runner(command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
    """Emulate manager success and MAFFT's conventional version response."""

    if command[-1] == "--version":
        return subprocess.CompletedProcess(command, 0, stdout="v7.526 (2024/Jun/15)\n", stderr="")
    return subprocess.CompletedProcess(command, 0, stdout="", stderr="")


class EnvironmentTests(unittest.TestCase):
    def test_manager_lookup_uses_stable_priority(self) -> None:
        found = find_environment_manager(lambda name: {"mamba": "/tools/mamba", "conda": "/tools/conda"}.get(name))

        self.assertEqual(found, ("mamba", "/tools/mamba"))

    def test_platform_directories_honor_explicit_overrides(self) -> None:
        values = {"TREETIP_CONFIG_DIR": "/tmp/treetip-config", "TREETIP_DATA_DIR": "/tmp/treetip-data"}

        self.assertEqual(default_config_directory(values), Path("/tmp/treetip-config"))
        self.assertEqual(default_data_directory(values), Path("/tmp/treetip-data"))

    def test_verify_mafft_returns_reported_version(self) -> None:
        version = verify_mafft("/tools/mafft", required_version=MAFFT_VERSION, runner=successful_runner)

        self.assertEqual(version, MAFFT_VERSION)

    def test_verify_mafft_rejects_wrong_managed_version(self) -> None:
        def old_version_runner(command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
            return subprocess.CompletedProcess(command, 0, stdout="v7.525\n", stderr="")

        with self.assertRaisesRegex(ValueError, "expected pinned version"):
            verify_mafft("/tools/mafft", required_version=MAFFT_VERSION, runner=old_version_runner)

    @patch("treetip_placer.environment.sys.platform", "linux")
    def test_setup_requires_a_supported_manager(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(ValueError, "micromamba, mamba, or conda"):
                setup_mafft_environment(
                    config_directory=Path(temporary) / "config",
                    data_directory=Path(temporary) / "data",
                    path_lookup=lambda _name: None,
                    runner=successful_runner,
                )

    @patch("treetip_placer.environment.sys.platform", "win32")
    def test_setup_rejects_native_windows(self) -> None:
        with self.assertRaisesRegex(ValueError, "supported on macOS and Linux"):
            setup_mafft_environment(
                path_lookup=lambda _name: "/tools/micromamba",
                runner=successful_runner,
            )

    @patch("treetip_placer.environment.sys.platform", "linux")
    def test_setup_creates_configured_environment(self) -> None:
        commands: list[list[str]] = []

        def recording_runner(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
            commands.append(command)
            return successful_runner(command, **kwargs)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            result = setup_mafft_environment(
                config_directory=root / "config",
                data_directory=root / "data",
                path_lookup=lambda name: "/tools/micromamba" if name == "micromamba" else None,
                runner=recording_runner,
            )

            self.assertTrue(result.created)
            self.assertEqual(result.manager, "micromamba")
            self.assertEqual(result.mafft.version, MAFFT_VERSION)
            self.assertEqual(commands[0][1:4], ["env", "create", "--yes"])
            self.assertTrue(result.config_path.is_file())
            saved = json.loads(result.config_path.read_text(encoding="utf-8"))
            self.assertEqual(saved["mafft_version"], MAFFT_VERSION)
            self.assertEqual(saved["source"], "managed_environment")

    @patch("treetip_placer.environment.sys.platform", "linux")
    def test_existing_environment_is_verified_without_recreation(self) -> None:
        commands: list[list[str]] = []

        def recording_runner(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
            commands.append(command)
            return successful_runner(command, **kwargs)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prefix = managed_environment_prefix(root / "data")
            prefix.mkdir(parents=True)
            result = setup_mafft_environment(
                config_directory=root / "config",
                data_directory=root / "data",
                path_lookup=lambda name: "/tools/conda" if name == "conda" else None,
                runner=recording_runner,
            )

            self.assertFalse(result.created)
            self.assertEqual(len(commands), 1)
            self.assertEqual(commands[0][-1], "--version")

    @patch("treetip_placer.environment.sys.platform", "linux")
    def test_refresh_removes_only_the_managed_prefix_before_recreating(self) -> None:
        commands: list[list[str]] = []

        def recording_runner(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
            commands.append(command)
            return successful_runner(command, **kwargs)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            managed_environment_prefix(root / "data").mkdir(parents=True)
            setup_mafft_environment(
                refresh=True,
                config_directory=root / "config",
                data_directory=root / "data",
                path_lookup=lambda name: "/tools/conda" if name == "conda" else None,
                runner=recording_runner,
            )

            self.assertEqual(commands[0][1:4], ["env", "remove", "--yes"])
            self.assertEqual(commands[1][1:4], ["env", "create", "--yes"])

    def test_classification_prefers_explicit_mafft(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            resolution = resolve_mafft_for_classification(
                "/explicit/mafft",
                config_directory=Path(temporary),
                path_lookup=lambda _name: "/path/mafft",
                runner=successful_runner,
            )

            self.assertEqual(resolution.executable, Path("/explicit/mafft"))
            self.assertEqual(resolution.source, "explicit")

    def test_classification_uses_verified_setup_configuration_before_path(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config_directory = Path(temporary)
            save_managed_mafft(
                manager="conda",
                prefix=Path("/managed"),
                resolution=MafftResolution(Path("/managed/bin/mafft"), MAFFT_VERSION, "managed_environment"),
                config_directory=config_directory,
            )
            resolution = resolve_mafft_for_classification(
                None,
                config_directory=config_directory,
                path_lookup=lambda _name: "/path/mafft",
                runner=successful_runner,
            )

            self.assertEqual(resolution.executable, Path("/managed/bin/mafft"))
            self.assertEqual(resolution.source, "managed_environment")

    def test_stale_setup_configuration_does_not_silently_fall_back_to_path(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config_directory = Path(temporary)
            save_managed_mafft(
                manager="conda",
                prefix=Path("/managed"),
                resolution=MafftResolution(Path("/managed/bin/mafft"), MAFFT_VERSION, "managed_environment"),
                config_directory=config_directory,
            )

            def missing_runner(_command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
                raise FileNotFoundError

            with self.assertRaisesRegex(ValueError, "could not be started"):
                resolve_mafft_for_classification(
                    None,
                    config_directory=config_directory,
                    path_lookup=lambda _name: "/path/mafft",
                    runner=missing_runner,
                )

    def test_classification_uses_path_when_setup_was_not_run(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            resolution = resolve_mafft_for_classification(
                None,
                config_directory=Path(temporary),
                path_lookup=lambda _name: "/path/mafft",
                runner=successful_runner,
            )

            self.assertEqual(resolution.source, "PATH")

    def test_classification_explains_how_to_recover_when_mafft_is_missing(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(ValueError, "treetip setup"):
                resolve_mafft_for_classification(
                    None,
                    config_directory=Path(temporary),
                    path_lookup=lambda _name: None,
                    runner=successful_runner,
                )

    def test_configuration_path_is_inside_given_directory(self) -> None:
        self.assertEqual(configuration_path(Path("/tmp/config")), Path("/tmp/config/mafft.json"))

    def test_load_managed_mafft_rejects_invalid_document(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = configuration_path(Path(temporary))
            path.write_text("not json", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "configuration is invalid"):
                load_managed_mafft(Path(temporary))


if __name__ == "__main__":
    unittest.main()
