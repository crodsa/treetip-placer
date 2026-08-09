# Public release readiness

## Current decision

The repository must remain private. The source code is Apache-2.0, but the
bundled `mlst-cdifficile-20260806` reference package is separately restricted
until its source-by-source redistribution rights are documented. A code licence
does not license reference sequences, alignments, tree labels, or clade metadata.

## Required gates before changing visibility

- [ ] Record the original source, version/date, licence or terms, and explicit
  redistribution basis for every artifact in the bundled reference package.
- [ ] Remove the package from the public release or add a maintainer-approved
  public package licence once every source is cleared.
- [ ] Audit tracked files and Git history for sensitive identifiers, restricted
  data, local operational paths, and credentials; rewrite history if needed.
- [ ] Run the pinned CI workflow on the candidate release revision.
- [ ] Record reproducible runtime and memory benchmarks for supported panel
  sizes and thread counts.
- [ ] Inspect Newick, Nexus, TSV, JSON, and experimental JPlace outputs with
  named compatible downstream software and retain the inspection record.
- [ ] Create a version tag and release notes that state the biological scope,
  calibration boundaries, and non-claims.

## Completed software scaffolding

- [x] Apache-2.0 licence and NOTICE for source code
- [x] CITATION.cff, CONTRIBUTING.md, SECURITY.md, CODE_OF_CONDUCT.md, and
  CHANGELOG.md
- [x] Issue and pull-request templates
- [x] Dependabot configuration
- [x] CI checks for formatting, linting, package installation, unit tests, and
  a deterministic synthetic smoke placement
