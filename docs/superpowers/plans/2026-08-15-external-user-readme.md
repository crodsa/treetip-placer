# External-user README implementation plan

_Documentation implementation plan for TreeTip Placer 0.7.0, approved 2026-08-15_

---

## 🎯 Goal and constraints

Implement the approved
[external-user README design](../specs/2026-08-15-external-user-readme-design.md)
without changing placement behavior, CLI interfaces, reference artifacts,
repository visibility, licences, version, or release status.

The repository remains English-only. Discussion with the owner may remain in
Spanish. All commands and claims must be checked against the current `main`
checkout rather than copied from an earlier TreeTip or ML-classifier task.

The system `/usr/bin/python3` is Python 3.8.2 and is outside the declared
Python 3.11+ support range. Verification must use
`/Users/crodsa/miniforge3/bin/python` (Python 3.12.10) with `PYTHONPATH=src`.

## 📋 Intended file scope

| File | Planned change |
| --- | --- |
| `README.md` | Replace the compressed landing page with the approved layered guide |
| `docs/MANUAL.md` | Correct only statements that conflict with the verified README contract |
| `docs/HANDOFF.md` | Remove only the stale current-state contradiction identified during review |
| `CHANGELOG.md` | Record the external-user documentation redesign under `Unreleased` |
| `docs/superpowers/plans/2026-08-15-external-user-readme.md` | Preserve this implementation plan |

Do not modify source, tests, package data, workflow configuration, version
metadata, licence files, or generated report assets unless a newly discovered
documentation defect makes that necessary and is reported before expansion.

## 🔍 Task 1: Build the evidence register

- [ ] Reconfirm the worktree, branch, remote, and divergence from `origin/main`
- [ ] Record `treetip --version` and help for `classify`, `setup`, `place-auto`,
  and `place` with Python 3.12
- [ ] Validate the bundled paired-reference manifest and retain its counts,
  coordinate width, clade counts, and distance threshold
- [ ] Reconfirm the output directory from `run_auto_place`,
  `write_run_outputs`, panel/clade writers, and report renderer
- [ ] Reconcile every validation number against the current reference-package
  README and `docs/HANDOFF.md`
- [ ] Verify any scientific background or upstream-service links against
  primary or authoritative sources before adding them
- [ ] Confirm that commit `0128a09` is absent from this repository and exclude
  the prior ML-classifier result from the implementation

## ✍️ Task 2: Rewrite the beginner path

- [ ] Add one clear outcome statement and a short audience description
- [ ] Define isolate, FASTA, MLST, locus, reference panel, clade, fixed tree,
  and placement before relying on those terms
- [ ] Explain suitable and unsuitable inputs, including the absence of a
  genome-to-MLST extraction command
- [ ] Add a synthetic two-record FASTA example with the required seven-locus
  order and no biological identifiers
- [ ] Explain what a terminal is and present installation commands without
  shell prompt characters
- [ ] State repository-access, Python, environment-manager, MAFFT, operating
  system, and `treetip setup` requirements truthfully
- [ ] Present the parser-verified `classify` command as the routine path
- [ ] Explain automatic reference, MAFFT, and output-directory resolution

## 🔄 Task 3: Document the workflow and decisions

- [ ] Add one accessible Mermaid flowchart with `accTitle` and `accDescr`
- [ ] Show separate mapping to classical and cryptic panels
- [ ] Show pre-placement QC and the route to unresolved output
- [ ] Show the distance-margin gate and the route to ambiguous output
- [ ] Show selected-panel fixed-tree placement and clade assignment
- [ ] Keep panel assignment, placement/QC status, and clade status as three
  separate interpretation layers
- [ ] Define `CLASSICAL`, `CRYPTIC`, `AMBIGUOUS`, `UNRESOLVED`, `PASS`, `WARN`,
  and `FAIL` only in contexts where the implementation can emit them
- [ ] Add synthetic interpretation examples for pass, warning, and no-call
  outcomes without assigning probabilities

## 📊 Task 4: Document outputs and methods

- [ ] Identify `report.html` as the first human-readable file and explain its
  two display-only tree views
- [ ] Inventory root panel/clade assignments, mappings, selected trees, both
  complete panel runs, placement alternatives, QC, manifests, and logs
- [ ] State that TSV/JSON carry authoritative decisions and QC
- [ ] State that Newick/Nexus are drawing-oriented and the distance-JPlace file
  is experimental rather than likelihood JPlace
- [ ] Explain normalized Hamming distance, the calibrated 0.01 panel margin,
  constrained WLS placement, optional bootstrap, and clade-label arbitration
- [ ] Keep nearest distance, distance margin, relative-fit weight, RMSE, and
  bootstrap support conceptually separate
- [ ] Explain checksums, parameters, software versions, paths, and event logs as
  the reproducibility trail

## 📚 Task 5: State package evidence and limitations

- [ ] Document the validated 201-reference classical panel, 68-reference
  cryptic panel, 3501 columns, five classical clades, and seven cryptic clades
- [ ] Describe the package's build tools only at the level supported by its
  versioned provenance record
- [ ] Summarize the panel and clade leave-one-out results with exact
  denominators, exclusions, and an explicit internal-validation label
- [ ] Avoid claiming external, clinical, diagnostic, or taxonomic validation
- [ ] Explain MLST resolution, reference dependence, unrepresented diversity,
  weakly resolved branches, and the need for core-genome or whole-genome work
  when finer resolution is required
- [ ] Preserve the Apache-2.0 software boundary and the separate restricted
  reference-package boundary
- [ ] State that the GitHub repository remains private and do not change that
  visibility

## 🔗 Task 6: Preserve advanced access and reconcile linked docs

- [ ] Add a compact table of contents and documentation index
- [ ] Retain links to `docs/MANUAL.md`, `docs/data-contract.md`,
  `docs/architecture.md`, the bundled package README, `docs/HANDOFF.md`,
  `docs/public-release.md`, `docs/release-checklist.md`, `CONTRIBUTING.md`,
  `SECURITY.md`, `CITATION.cff`, `CHANGELOG.md`, and `LICENSE`
- [ ] Introduce `place-auto`, `place`, package creation/validation, and the three
  calibration commands without placing maintainer flags in the quick start
- [ ] Correct only contradictory claims in the manual and handoff document
- [ ] Add an `Unreleased` changelog entry without changing the project version

## ✅ Task 7: Verify the complete documentation change

- [ ] Re-run all documented CLI help commands with Python 3.12
- [ ] Re-run bundled reference-package validation
- [ ] Programmatically check every local Markdown link and image target
- [ ] Check Markdown fence balance, heading hierarchy, one H1, and Mermaid
  accessibility fields
- [ ] Check the Mermaid source for supported flowchart syntax, semantic node
  identifiers, one flow direction, and no inline theme/style directives
- [ ] Scan for placeholder text, stale Windows paths, unsupported promises,
  accidental real identifiers, and references to the wrong repository
- [ ] Run the full supported-Python suite:
  `PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src /Users/crodsa/miniforge3/bin/python -m unittest discover -s tests -v`
- [ ] Run `git diff --check`
- [ ] Inspect the complete diff and confirm that only intended documentation
  files changed

## 🚀 Task 8: Commit and synchronize private main

- [ ] Stage only the approved documentation files
- [ ] Commit with an English message that describes the complete README change
- [ ] Re-run the final verification gate on the committed tree
- [ ] Reconfirm that local `main` is based on current `origin/main`
- [ ] Push the verified commits directly to the confirmed private
  `crodsa/treetip-placer` `main`, as explicitly requested by the owner
- [ ] Verify the remote `main` commit and GitHub Actions result when available
- [ ] Report exact commits, tests, documentation checks, remote state, and any
  release blockers that remain

Direct synchronization to private `main` does not authorize a pull request,
public release, tag, licence change, reference-data redistribution, or
repository-visibility change.
