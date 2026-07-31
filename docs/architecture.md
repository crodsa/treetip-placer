# Architecture

## Design principle

TreeTip Placer treats a reference tree as an immutable backbone. A query is projected onto the reference feature space, screened by quality-control gates, scored on candidate edges, and reported with explicit uncertainty. The program does not imply that a fixed-backbone result is identical to a full de novo tree reconstruction.

## Input adapters

`mlst_allele_profiles` reads tab-separated categorical allele calls. It uses normalized allelic Hamming distances and reports missing or composite calls.

`fixed_column_fasta_alignment` reads aligned FASTA records. Each fixed column is a categorical character. It is appropriate for concatenated MLST alleles, core SNP alignments, and concatenated core genes only when queries are already in the exact reference coordinate system.

Both adapters produce an ordered character profile consumed by the same placement engine. This isolates format-specific validation from placement mathematics.

## Current placement engine

For every query, the engine obtains distances to reference tips and scores an insertion on a branch by constrained weighted least squares. It estimates the attachment position along the edge and a non-negative pendant length. The `fast` and `balanced` presets use a neighbourhood around the nearest reference profiles; `thorough` evaluates every edge.

The reported `relative_fit_weight` is a heuristic normalization of distance-objective differences. It is not a likelihood weight ratio and is labelled accordingly in every output.

Bootstrap resampling is available over categorical loci or alignment columns. It measures stability among the original candidate edges; it is not a replacement for a whole-tree bootstrap.

## Planned backends

1. Maximum-likelihood placement adapter for aligned nucleotide or amino-acid data, preserving model parameters, partitions, and SNP ascertainment metadata from the reference analysis.
2. Mutation-annotated parsimony adapter for exceptionally large, low-divergence genome panels.
3. Conditional local refinement: extract the uncertain local region, keep the exterior fixed through boundary likelihood information, optimize local NNI/SPR alternatives, and report newly computed support only for affected branches.

## Parallel execution

Queries are independent against a reference snapshot and are distributed across process workers. A worker creates its own reference cache, so the requested worker count must respect available memory. The command falls back to serial execution if the host forbids process workers and records that decision in `events.jsonl` and `run.log`.

Bootstrap is intentionally performed inside each query worker. This avoids nested parallelism and oversubscription when many queries are supplied.

## Output contract

`placed_tree.nwk` is a simple drawing-oriented tree. `placed_tree.nexus` offers the same topology in Nexus form. `placements.tsv` is the authoritative tabular result, `alternative_edges.tsv` holds ambiguity information, and `query_qc.tsv` describes all gate findings. `run_manifest.json` and `events.jsonl` make a run reproducible.

The optional JPlace-shaped JSON document is explicitly marked experimental for the distance backend because its objective is not a maximum-likelihood score.
