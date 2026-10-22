# Changelog

All notable changes are documented here. The project follows semantic versioning
once public releases begin.

## Unreleased

### Changed

- Expanded the README into a layered external-user guide covering first use,
  input preparation, decision logic, complete outputs, methods, package
  evidence, interpretation, limitations, and advanced documentation; corrected
  linked manual and handoff statements that contradicted the verified contract.
- Reworked automatic-placement HTML reports into two complementary per-query
  views: a reference-name-free global clade context and a topology-derived local
  placement branch with named neighbours, isolating each query from other placed
  queries in multi-query batch runs.
- Enforced `final_status != "FAIL"` during panel selection and prioritized
  references meeting `min_comparable_sites` during nearest-reference QC.
- Relicensed software and documentation under Apache-2.0 while explicitly
  retaining the current private reference package outside that licence.

### Added

- Public-release governance, citation, contribution, security, and community
  files.
- CI checks for Ruff, package installation, unit tests, and a synthetic CLI
  smoke placement.
