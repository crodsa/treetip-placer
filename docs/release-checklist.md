# Repository publication checklist

## Private development checkpoint

These gates allow synchronization only to the confirmed private GitHub backup.
They do not authorize a public release.

- [x] The complete interface and documentation are English-only.
- [x] The complete local regression suite passes.
- [x] At least one redistributable synthetic fixture exercises the complete workflow.
- [x] Real aligned-sequence validation has been performed locally with permitted data; data and outputs remain untracked.
- [x] Critical and high-severity defects identified for this checkpoint are resolved or explicitly recorded in `docs/HANDOFF.md`.
- [x] Dependencies and licensing boundaries are reviewed for the private checkpoint.
- [x] No private genomes, alignments, identifiers, credentials, or restricted database material are tracked.
- [x] The destination `crodsa/treetip-placer` and private visibility are confirmed by the repository owner.

## Public release

Do not publish publicly until every applicable item is complete.

- [ ] Unit, integration, and reproducibility tests pass in continuous integration on the candidate public-release revision.
- [ ] Placement accuracy, support calibration, runtime, and memory benchmarks are documented.
- [ ] Newick, Nexus, TSV, JSON, and JPlace outputs have been inspected with compatible downstream software.
- [ ] Source code is licensed under Apache-2.0 and every bundled reference-data artifact has separate, documented redistribution permission or has been removed from the public release.
- [x] `README.md`, `CITATION.cff`, `CONTRIBUTING.md`, `CHANGELOG.md`, security/community files, and issue templates are present; complete a final release-specific review.
- [x] GitHub Actions is configured for linting, formatting, tests, package installation, and a synthetic smoke placement; confirm a passing run on the candidate release revision.
- [ ] A version tag and release notes have been prepared.
