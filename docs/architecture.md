# Architecture

## Design principle

TreeTip Placer treats a branch-length reference tree and its fixed-column FASTA alignment as an immutable package. Query sequences must already be mapped to the same columns; the tool does not perform an alignment during placement.

## Current engine

The FASTA adapter represents callable nucleotide columns as an ordered sequence. For every query, the engine computes pairwise normalized Hamming distances to reference tips and scores candidate insertion edges by constrained least squares. `fast` and `balanced` search a neighbourhood around nearest references; `thorough` scores every edge.

The QC gate reports coverage, comparable sites, ambiguous bases, panel distance, diffuse placements, and long terminal branches. Query placements run independently in process workers, with a recorded serial fallback when workers are unavailable. Bootstrap resampling occurs inside each worker to avoid nested parallelism.

`place-auto` maps an uploaded FASTA independently to a classical and a cryptic immutable package using MAFFT `--add --keeplength`. Each mapped query receives a full independent placement. The selector compares only nearest-reference normalized Hamming distances because WLS edge objectives and heuristic weights are not comparable across differently sized trees. It emits a panel label only when the winning distance margin exceeds a recorded calibration threshold; otherwise it emits `AMBIGUOUS` or `UNRESOLVED`.

Paired-reference manifests separate package identity from a specific workstation. They resolve panel tree/alignment paths relative to the manifest and may declare checksums. The automatic command validates declared checksums before using either reference, so an edited tree or alignment cannot silently inherit a prior calibration threshold.

`calibrate-panel-gate` validates that distance gate with deterministic internal leave-one-unique-sequence-out controls. Exact duplicate own-panel sequences are removed before a held-out control is scored. Controls are mapped only to the opposing panel; the procedure deliberately measures panel discrimination rather than WLS insertion accuracy or tree topology recovery.

## Reference packages

A reference package must preserve the original alignment, tree, sample-ID mapping, filtering/masking rules, sequence partitions when available, software versions, and checksums. A cryptic-only package and a mixed cryptic–classical package require separate calibration because their distance distributions differ.

## Planned work

1. Maximum-likelihood placement for aligned nucleotide data with explicit substitution model, partitions, and SNP ascertainment metadata.
2. Alignment leave-one-out validation and calibration reports.
3. Conditional local topology refinement with support recomputation only in the affected region.

## Outputs

Newick and Nexus are drawing-oriented outputs. TSV and JSON sidecars contain QC and uncertainty. The experimental JPlace-shaped output is distance-based and must not be interpreted as maximum-likelihood JPlace.
