"""Managed MAFFT environment discovery, setup, and verification.

This module deliberately owns only the external alignment dependency. It does
not know about reference packages or placement so the installation contract
can be tested without biological inputs.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import sysconfig
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Mapping, Sequence


MAFFT_VERSION = "7.526"
CONFIG_FILENAME = "mafft.json"
ENVIRONMENT_DIRECTORY_NAME = f"mafft-{MAFFT_VERSION}"
_VERSION_PATTERN = re.compile(r"(?:v)?(\d+\.\d+(?:\.\d+)?)")
_REPO_ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class MafftResolution:
    """A verified MAFFT executable and the way it was selected."""

    executable: Path
    version: str
    source: str


@dataclass(frozen=True)
class SetupResult:
    """Details emitted after a successful managed-environment setup."""

    manager: str
    environment_prefix: Path
    config_path: Path
    mafft: MafftResolution
    created: bool


def default_config_directory(
    environ: Mapping[str, str] | None = None,
    *,
    platform: str | None = None,
    home: Path | None = None,
) -> Path:
    """Return the user-scoped directory that stores TreeTip configuration."""

    values = os.environ if environ is None else environ
    if values.get("TREETIP_CONFIG_DIR"):
        return Path(values["TREETIP_CONFIG_DIR"]).expanduser()
    user_home = home or Path.home()
    current_platform = platform or sys.platform
    if current_platform == "darwin":
        return user_home / "Library" / "Application Support" / "treetip-placer"
    if current_platform.startswith("win"):
        return Path(values.get("APPDATA", str(user_home / "AppData" / "Roaming"))) / "treetip-placer"
    return Path(values.get("XDG_CONFIG_HOME", str(user_home / ".config"))) / "treetip-placer"


def default_data_directory(
    environ: Mapping[str, str] | None = None,
    *,
    platform: str | None = None,
    home: Path | None = None,
) -> Path:
    """Return the user-scoped directory that owns managed environments."""

    values = os.environ if environ is None else environ
    if values.get("TREETIP_DATA_DIR"):
        return Path(values["TREETIP_DATA_DIR"]).expanduser()
    user_home = home or Path.home()
    current_platform = platform or sys.platform
    if current_platform == "darwin":
        return user_home / "Library" / "Application Support" / "treetip-placer"
    if current_platform.startswith("win"):
        return Path(values.get("LOCALAPPDATA", str(user_home / "AppData" / "Local"))) / "treetip-placer"
    return Path(values.get("XDG_DATA_HOME", str(user_home / ".local" / "share"))) / "treetip-placer"


def configuration_path(config_directory: Path | None = None) -> Path:
    """Return the MAFFT configuration file location."""

    return (config_directory or default_config_directory()) / CONFIG_FILENAME


def managed_environment_prefix(data_directory: Path | None = None) -> Path:
    """Return the sole prefix managed by TreeTip setup."""

    return (data_directory or default_data_directory()) / "environments" / ENVIRONMENT_DIRECTORY_NAME


def environment_specification_path() -> Path:
    """Return the checked-in or wheel-installed Conda specification used by setup."""

    repo_spec = _REPO_ROOT / "environment.yml"
    if repo_spec.is_file():
        return repo_spec
    installed_spec = Path(sysconfig.get_path("data")) / "share" / "treetip-placer" / "environment.yml"
    if installed_spec.is_file():
        return installed_spec
    return repo_spec


def find_environment_manager(
    path_lookup: Callable[[str], str | None] = shutil.which,
) -> tuple[str, str] | None:
    """Return the first supported manager in a stable priority order."""

    for name in ("micromamba", "mamba", "conda"):
        executable = path_lookup(name)
        if executable:
            return name, executable
    return None


def manager_create_command(manager: str, manager_executable: str, prefix: Path, specification: Path) -> list[str]:
    """Build a non-interactive Conda-compatible create command."""

    if manager not in {"micromamba", "mamba", "conda"}:
        raise ValueError(f"Unsupported environment manager: {manager}")
    return [manager_executable, "env", "create", "--yes", "--prefix", str(prefix), "--file", str(specification)]


def manager_remove_command(manager: str, manager_executable: str, prefix: Path) -> list[str]:
    """Build a non-interactive command limited to the named managed prefix."""

    if manager not in {"micromamba", "mamba", "conda"}:
        raise ValueError(f"Unsupported environment manager: {manager}")
    return [manager_executable, "env", "remove", "--yes", "--prefix", str(prefix)]


def verify_mafft(
    executable: str | Path,
    *,
    required_version: str | None = None,
    runner: Callable[..., subprocess.CompletedProcess[str]] = subprocess.run,
) -> str:
    """Run MAFFT's version command and return its reported semantic version."""

    path = Path(executable).expanduser()
    try:
        completed = runner(
            [str(path), "--version"],
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError as error:
        raise ValueError(f"MAFFT executable could not be started: {path}") from error
    output = "\n".join(part for part in (completed.stdout, completed.stderr) if part).strip()
    if completed.returncode != 0:
        detail = f": {output}" if output else ""
        raise ValueError(f"MAFFT verification failed for {path} (exit {completed.returncode}){detail}")
    version_match = _VERSION_PATTERN.search(output)
    if version_match is None:
        raise ValueError(f"MAFFT verification did not report a version for {path}")
    version = version_match.group(1)
    if required_version is not None and version != required_version:
        raise ValueError(
            f"Managed MAFFT environment reported version {version}; expected pinned version {required_version}. "
            "Run 'treetip setup --refresh'."
        )
    return version


def load_managed_mafft(config_directory: Path | None = None) -> MafftResolution | None:
    """Read the last verified managed MAFFT selection, if one was saved."""

    path = configuration_path(config_directory)
    if not path.is_file():
        return None
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
        executable = Path(document["mafft_executable"])
        version = str(document["mafft_version"])
        source = str(document["source"])
    except (OSError, KeyError, TypeError, json.JSONDecodeError) as error:
        raise ValueError(f"TreeTip MAFFT configuration is invalid: {path}. Run 'treetip setup --refresh'.") from error
    if source != "managed_environment":
        raise ValueError(
            f"TreeTip MAFFT configuration has an unsupported source: {path}. Run 'treetip setup --refresh'."
        )
    return MafftResolution(executable=executable, version=version, source="managed_environment")


def save_managed_mafft(
    *,
    manager: str,
    prefix: Path,
    resolution: MafftResolution,
    config_directory: Path | None = None,
) -> Path:
    """Atomically save a verified TreeTip-managed MAFFT selection."""

    path = configuration_path(config_directory)
    path.parent.mkdir(parents=True, exist_ok=True)
    document = {
        "schema_version": 1,
        "manager": manager,
        "environment_prefix": str(prefix),
        "mafft_executable": str(resolution.executable),
        "mafft_version": resolution.version,
        "source": resolution.source,
    }
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="\n",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as temporary:
        json.dump(document, temporary, indent=2, sort_keys=True)
        temporary.write("\n")
        temporary_path = Path(temporary.name)
    temporary_path.replace(path)
    return path


def _mafft_in_prefix(prefix: Path) -> Path:
    """Resolve the conventional MAFFT executable location inside a Conda prefix."""

    executable_name = "mafft.exe" if sys.platform.startswith("win") else "mafft"
    return prefix / ("Library" / "bin" if sys.platform.startswith("win") else "bin") / executable_name


def _run_manager(command: Sequence[str], runner: Callable[..., subprocess.CompletedProcess[str]]) -> None:
    """Run a manager command and preserve its diagnostic output on failure."""

    try:
        completed = runner(list(command), check=False, capture_output=True, text=True)
    except OSError as error:
        raise ValueError(f"Environment manager could not be started: {command[0]}") from error
    if completed.returncode != 0:
        output = "\n".join(part for part in (completed.stdout, completed.stderr) if part).strip()
        detail = f"\n{output}" if output else ""
        raise ValueError(f"Environment manager failed (exit {completed.returncode}){detail}")


def setup_mafft_environment(
    *,
    refresh: bool = False,
    config_directory: Path | None = None,
    data_directory: Path | None = None,
    path_lookup: Callable[[str], str | None] = shutil.which,
    runner: Callable[..., subprocess.CompletedProcess[str]] = subprocess.run,
) -> SetupResult:
    """Create or verify the TreeTip-owned, version-pinned MAFFT environment."""

    if sys.platform.startswith("win"):
        raise ValueError(
            "Managed MAFFT setup is currently supported on macOS and Linux. "
            "Use a Linux environment (for example WSL) or pass a verified --mafft path explicitly."
        )
    prefix = managed_environment_prefix(data_directory)
    specification = environment_specification_path()
    if not specification.is_file():
        raise ValueError(f"TreeTip environment specification is missing: {specification}")
    manager = find_environment_manager(path_lookup)
    if manager is None:
        raise ValueError(
            "No supported environment manager was found. Install micromamba, mamba, or conda, then run 'treetip setup'."
        )
    manager_name, manager_executable = manager
    mafft = _mafft_in_prefix(prefix)
    created = False
    if prefix.exists() and not refresh:
        version = verify_mafft(mafft, required_version=MAFFT_VERSION, runner=runner)
    else:
        if refresh and prefix.exists():
            _run_manager(manager_remove_command(manager_name, manager_executable, prefix), runner)
        _run_manager(manager_create_command(manager_name, manager_executable, prefix, specification), runner)
        version = verify_mafft(mafft, required_version=MAFFT_VERSION, runner=runner)
        created = True
    resolution = MafftResolution(executable=mafft, version=version, source="managed_environment")
    config_path = save_managed_mafft(
        manager=manager_name,
        prefix=prefix,
        resolution=resolution,
        config_directory=config_directory,
    )
    return SetupResult(
        manager=manager_name,
        environment_prefix=prefix,
        config_path=config_path,
        mafft=resolution,
        created=created,
    )


def resolve_mafft_for_classification(
    explicit: str | None,
    *,
    config_directory: Path | None = None,
    path_lookup: Callable[[str], str | None] = shutil.which,
    runner: Callable[..., subprocess.CompletedProcess[str]] = subprocess.run,
) -> MafftResolution:
    """Resolve and verify MAFFT using explicit, managed, then PATH precedence."""

    if explicit:
        executable = Path(explicit).expanduser()
        return MafftResolution(
            executable=executable, version=verify_mafft(executable, runner=runner), source="explicit"
        )
    configured = load_managed_mafft(config_directory)
    if configured is not None:
        version = verify_mafft(configured.executable, required_version=MAFFT_VERSION, runner=runner)
        return MafftResolution(executable=configured.executable, version=version, source="managed_environment")
    discovered = path_lookup("mafft")
    if discovered:
        executable = Path(discovered)
        return MafftResolution(executable=executable, version=verify_mafft(executable, runner=runner), source="PATH")
    raise ValueError(
        "MAFFT is unavailable. Run 'treetip setup' to install the managed dependency, "
        "install MAFFT independently, or pass --mafft /path/to/mafft."
    )
