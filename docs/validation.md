# MLST validation protocol

## Purpose

This protocol calibrates TreeTip Placer for a specific MLST reference panel without adding restricted data to the repository. It evaluates the categorical-distance backend only; it does not validate a future maximum-likelihood backend.

## Required inputs

- A branch-length Newick tree generated from the intended complete reference panel.
- A complete MLST allele-profile TSV with exactly the same sample IDs as the tree tips.
- Permission to use the profiles locally for validation.

Keep those inputs outside version control unless their redistribution rights are explicit.

## Run design

Start with all profiles when the panel is modest:

```powershell
treetip validate-mlst `
  --tree reference_tree.nwk `
  --reference-profiles reference_profiles.tsv `
  --output validation_all `
  --preset thorough `
  --bootstrap 100 `
  --threads auto `
  --verbose
```

For a large panel, select a reproducible subset and retain the command line, seed, checksums, and software commit with the results:

```powershell
treetip validate-mlst `
  --tree reference_tree.nwk `
  --reference-profiles reference_profiles.tsv `
  --output validation_subset `
  --max-holdouts 250 `
  --seed 20260731 `
  --preset balanced `
  --threads auto `
  --verbose
```

Use `thorough` at least once for the final calibration dataset. A balanced or fast run may be used to explore runtime, but its candidate restriction must be reported.

## What is measured

Each validation task removes the held-out sample and its matching tip from a copied backbone before placement. `leave_one_out.tsv` reports the known attachment representation, selected edge, QC status, relative-fit weight, optional bootstrap support, and `topological_distance_edges`.

A recovery has topological distance zero. This means the selected edge touches the known retained attachment node or exactly matches the merged edge containing the original attachment after pruning. It is intentionally a local attachment-neighbourhood criterion, not a claim of full-tree reconstruction accuracy.

Review at least:

1. Recovery rate among records with a selected edge.
2. Distribution of topological distances for non-recovered placements.
3. Recovery stratified by nearest-reference distance, coverage, and number of comparable loci.
4. The rate of `WARN`, `FAIL`, and `UNRESOLVED` status by the same strata.
5. Calibration of `relative_fit_weight` and bootstrap support against the observed recovery rate.
6. Runtime and peak memory at the chosen thread count.

## Calibration record

Create a local results note containing:

- Tool version and Git commit.
- Input SHA-256 values from `validation_manifest.json`.
- Tree construction method, alignment/profile provenance, and mask or allele-calling version.
- Command line, hardware, operating system, worker count, and elapsed time.
- Summary metrics and stratified tables or plots.
- Chosen QC thresholds with a scientific justification.
- Known failure modes and the decision on the intended operational range.

Do not commit patient identifiers, unreleased isolates, restricted allele profiles, or other non-redistributable material.
