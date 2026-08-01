# Repository publication checklist

## Private development checkpoint

These gates allow synchronization only to the confirmed private GitHub backup.
They do not authorize a public release.

- [x] The complete interface and documentation are English-only.
- [x] The complete local regression suite passes.
- [x] At least one redistributable synthetic fixture exercises the complete workflow.
- [x] Real aligned-sequence validation has been performed locally with permitted data; data and outputs remain untracked.
- [x] Critical and high-severity defects identified for this checkpoint are resolved or explicitly recorded in `docs/HANDOFF.md`.
- [x] Dependencies and the private proprietary license are reviewed.
- [x] No private genomes, alignments, identifiers, credentials, or restricted database material are tracked.
- [x] The destination `crodsa/treetip-placer` and private visibility are confirmed by the repository owner.

## Public release

Do not publish publicly until every applicable item is complete.

- [ ] Unit, integration, and reproducibility tests pass in continuous integration.
- [ ] Placement accuracy, support calibration, runtime, and memory benchmarks are documented.
- [ ] Newick, Nexus, TSV, JSON, and JPlace outputs have been inspected with compatible downstream software.
- [ ] A public license has been selected and added as `LICENSE`.
- [ ] `README.md`, `CITATION.cff`, `CONTRIBUTING.md`, `CHANGELOG.md`, and issue templates are complete.
- [ ] GitHub Actions runs linting, tests, and a smoke placement.
- [ ] A version tag and release notes have been prepared.
