# Project handoff

Last updated: 2026-07-31

## Purpose

TreeTip Placer adds one or more new genomes or profiles to an immutable reference tree without rebuilding the reference alignment or backbone topology. It reports placement uncertainty, quality-control decisions, and portable output files.

## Current implementation

Version `0.1.0` is a tested MVP with a categorical-distance backend.

- MLST allele-profile input: tab-separated files with `sample_id` plus categorical locus columns.
- Fixed-column FASTA input: pre-aligned or reference-mapped concatenated MLST alleles, core SNPs, or concatenated core genes.
- Fixed-tree weighted least-squares placement with non-negative pendant length and in-edge attachment optimization.
- `fast`, `balanced`, and exhaustive `thorough` candidate-search presets.
- Pre-placement QC and post-placement resolution status.
- Optional bootstrap resampling of loci or alignment sites.
- Parallel query execution, with explicit serial fallback when process workers are unavailable.
- Newick, Nexus, TSV, JSON, JSONL, and experimental distance-JPlace outputs.

## Important limitations

- The current FASTA adapter uses categorical-distance placement; it is not a maximum-likelihood implementation.
- The current implementation does not yet use EPA-ng, pplacer, UShER, or any other external placement engine.
- The current implementation does not yet perform conditional local NNI/SPR refinement.
- Relative-fit weights are heuristic distance-objective weights, not likelihood weight ratios or posterior probabilities.
- QC thresholds are defaults. They must be calibrated with real, held-out MLST data before scientific use.

## Verified commands

```powershell
python -m pip install --editable . --no-deps
python -m unittest discover -s tests -v
```

Example profile placement:

```powershell
treetip place `
  --tree reference_tree.nwk `
  --reference-profiles reference_profiles.tsv `
  --query-profiles query_profiles.tsv `
  --output placement_run `
  --preset thorough `
  --bootstrap 100 `
  --verbose
```

Example aligned-sequence placement:

```powershell
treetip place `
  --tree reference_tree.nwk `
  --reference-alignment reference_alignment.fasta `
  --query-alignment query_alignment.fasta `
  --output placement_run `
  --preset thorough
```

## Required real-data validation

When real MLST data are available:

1. Verify that tree tips, reference profile IDs, and locus names are consistent.
2. Withhold known reference profiles, place them, and measure edge accuracy and topological distance.
3. Stratify validation by nearest-reference distance, missingness, and novelty.
4. Calibrate QC thresholds and assess bootstrap/relative-fit calibration.
5. Compare against a reconstruction with held-out samples restored and, where appropriate, established placement software.
6. Record results in `docs/validation.md` without committing restricted data.

## Next development milestones

1. Real-data MLST validation and QC calibration.
2. Robust reference-package manifest with immutable masks, partitions, and provenance.
3. Maximum-likelihood aligned-sequence backend with explicit SNP ascertainment metadata.
4. Conditional local topology refinement and local support recomputation.
5. Benchmarking, release hardening, selected license, and GitHub publication.

## Portability protocol

1. Read this document and `AGENTS.md` at the start of a new Codex chat.
2. Use the same repository revision and run the verification commands above.
3. Keep real validation data outside the repository unless sharing rights are explicit.
4. Update this document after each completed milestone, major decision, or validation result.
5. Before using Codex cloud or another local computer, place the repository in a private remote Git repository or another approved synchronized storage location. A Git remote transports files and history; it does not automatically transport private local data or credentials.

## Repository state

The project is prepared for a private GitHub backup. The initial local Git history is intentionally separate from any future public release; public publication remains subject to the release checklist and explicit user confirmation. Confirm the active remote with `git remote -v` before continuing work on another computer.
