# User-first README and managed-environment design

## Goal

Make the routine TreeTip Placer experience understandable and executable by a
scientific end user without requiring them to understand reference-package
construction, calibration workflows, or MAFFT installation. The routine
workflow starts with a raw FASTA and produces a named result directory and a
human-readable HTML report.

## Scope

This change has two connected parts:

1. Replace the repository README's maintainer-first narrative with a concise
   end-user guide ordered as: problem, method, installation, one-command run,
   result interpretation, other applications, and links to technical and
   maintenance documentation.
2. Add a `treetip setup` command that creates a dedicated, reproducible Conda
   environment containing MAFFT and records the exact executable selected for
   future `treetip classify` runs.

The supported routine command remains explicit about the selected bundled
reference package:

```bash
treetip classify genome.fasta --organism cdifficile-mlst
```

It accepts raw FASTA, resolves its reference package and output directory, and
uses the MAFFT configured by `treetip setup` when no `--mafft` override is
provided.

## User-facing behaviour

### Installation and setup

Users install the Python package in the normal supported way, then run:

```bash
treetip setup
```

`setup` finds one supported Conda-compatible environment manager in this order:
`micromamba`, `mamba`, then `conda`. It creates a TreeTip-owned prefix under
the user's TreeTip data directory (not in an existing shared environment),
installs `mafft=7.526` from `conda-forge` as declared in the repository-owned
`environment.yml`, resolves the environment's `mafft` executable, and writes a
small user configuration file that records the manager, environment prefix,
executable path, and observed MAFFT version. Version 7.526 is selected because
conda-forge supplies it for macOS Intel, macOS Apple Silicon, and Linux x86_64
and aarch64 at design time.

The command prints the selected manager, the installed MAFFT version, and the
configuration-file location. It must never edit an existing user environment,
modify `PATH`, or silently choose an arbitrary existing MAFFT installation.

If an existing TreeTip-managed environment already satisfies the pinned
specification, setup reports the verified version and exits successfully. A
dedicated `--refresh` option may recreate or update only the TreeTip-owned
prefix; the command prints that exact prefix immediately before performing the
operation. The initial implementation does not perform network package-manager
installation unless the user explicitly runs `treetip setup` themselves.

### Classification

MAFFT resolution for `classify` follows this precedence:

1. An explicit `--mafft PATH`.
2. A previously verified executable recorded by `treetip setup`.
3. An executable named `mafft` found on `PATH`.

The selected path is checked by running `mafft --version` before processing
queries. A missing, non-executable, or non-responsive MAFFT produces one
actionable error: either run `treetip setup`, install MAFFT through a supported
manager, or supply `--mafft`.

The existing fixed-tree, dual-panel mapping, QC, conservative panel choice,
clade call, output naming, and `report.html` behavior remain unchanged.

## Files and boundaries

- `src/treetip_placer/environment.py`: isolated discovery, environment setup,
  configuration read/write, and MAFFT validation. It must have no placement
  logic.
- `src/treetip_placer/cli.py`: exposes `treetip setup` and delegates MAFFT
  discovery to the environment module. It remains responsible for CLI parsing
  and messages.
- `environment.yml`: declarative environment input with only
  `conda-forge` and `mafft=7.526`.
- `README.md`: end-user landing page, deliberately concise.
- `docs/MANUAL.md`: expanded end-user instructions and troubleshooting.
- `docs/architecture.md` and `docs/HANDOFF.md`: technical implementation,
  provenance, and maintenance details.

The bundled reference package keeps its separate redistribution restriction.
This design does not alter repository visibility, licensing, package contents,
or scientific interpretation.

## README structure

The README will contain these sections only:

1. What question the tool answers, and its fixed-reference limitation.
2. How it works at a practical level (reference mapping, QC, conservative
   panel/clade assignment, two complementary report views).
3. Installation and one-time `treetip setup`.
4. One normal `treetip classify` invocation.
5. How to read `report.html`, including `CLASSICAL`, `CRYPTIC`, `AMBIGUOUS`,
   `UNRESOLVED`, `PASS`, and `WARN`.
6. Other appropriate uses and the boundaries of those uses.
7. Short links to advanced commands, reference-package maintenance,
   contribution/security material, and the licence boundary.

Reference-package creation, calibration, schema definitions, QA details, and
developer test commands move behind these links rather than appearing in the
normal reading path.

## Error handling and reproducibility

- Setup has a single, named environment and a recorded versioned configuration.
- Configuration writes are atomic and use a user-scoped platform-appropriate
  directory (`$XDG_CONFIG_HOME`/`$XDG_DATA_HOME` when set; platform defaults
  otherwise); exact paths are reported. `TREETIP_CONFIG_DIR` and
  `TREETIP_DATA_DIR` are documented test and administrator overrides.
- Setup validates both the environment-manager exit status and `mafft --version`
  before saving configuration.
- Classification validates the configured executable at use time, because the
  environment may have been removed or changed after setup.
- All setup and classification decisions are written into a machine-readable
  run/provenance record where one exists; no implicit fallback is hidden.
- `--mafft` remains an auditable override for constrained HPC or offline
  environments.

## Tests

Tests use mocked subprocess calls and temporary configuration directories to
cover manager priority, manager absence, successful setup, failed environment
creation, failed MAFFT verification, configuration precedence, stale
configuration, explicit override, and `PATH` fallback. Existing classification
tests continue to exercise the full synthetic placement path. Documentation
examples are checked against the parser/help text where practical.

## Acceptance criteria

- A newcomer can follow the README from install to `report.html` without
  needing to construct a reference manifest or locate MAFFT manually.
- `treetip setup` creates and verifies only a dedicated TreeTip environment
  using a pinned MAFFT specification.
- `treetip classify genome.fasta --organism cdifficile-mlst` uses the setup
  configuration automatically, while preserving explicit `--mafft` and PATH
  alternatives.
- Missing prerequisites produce clear recovery commands and no ambiguous run
  directory.
- The README contains no maintainer workflow details beyond compact links.
- Full unit tests, static checks, and a synthetic end-to-end CLI smoke test
  pass before merge.
