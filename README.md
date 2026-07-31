# TreeTip Placer

TreeTip Placer adds query sequences to an immutable, branch-length reference tree without rebuilding its topology. It is designed exclusively for fixed-column FASTA alignments: concatenated MLST allele sequences, concatenated core genes, or core-genome SNP matrices.

## Scope

The current backend performs quality-controlled, distance-based placement on a fixed reference tree. It does not replace a de novo phylogenetic reconstruction or a maximum-likelihood placement engine. The reference alignment, coordinate system, masking policy, and tree must be retained as one versioned reference package.

Tabular MLST allele-profile inputs are intentionally unsupported.

## Input contract

- A branch-length Newick reference tree.
- A fixed-column reference FASTA alignment whose IDs exactly match tree tips.
- A query FASTA alignment mapped to the identical reference columns.

The program will reject unequal alignment widths rather than silently realigning, trimming, padding, or dropping query sites. For SNP-only matrices, retain the original SNP-selection, masking, and ascertainment metadata beside the reference package.

## Quick start

```bash
treetip place \
  --tree reference_tree.nwk \
  --reference-alignment reference_aligned.fasta \
  --query-alignment batch_aligned.fasta \
  --output placement_run \
  --preset balanced \
  --threads auto \
  --bootstrap 100 \
  --verbose
```

`placed_tree.nwk` is directly usable in FigTree, iTOL, Dendroscope, ETE, and most Newick-compatible programs. `placements.tsv`, `alternative_edges.tsv`, `query_qc.tsv`, JSON, and event logs provide machine-readable placement and QC detail.

## Quality control

Each query receives `PASS`, `WARN`, `FAIL`, or `UNRESOLVED`. The gate checks alignment coverage, callable-site count, ambiguous nucleotide symbols, nearest-reference distance, and post-placement ambiguity. Thresholds are panel-specific: a cryptic-only reference panel must not reuse distance thresholds calibrated on a mixed cryptic–classical panel.

## Current limitations

- The scorer is categorical fixed-column distance, not maximum likelihood.
- Relative-fit weights are heuristic distance-objective weights, not likelihood weight ratios or posterior probabilities.
- Local NNI/SPR refinement and recomputed local support are not yet implemented.

## Development

The test suite uses only synthetic fixtures:

```powershell
python -m unittest discover -s tests -v
```

Do not commit private alignments, genomes, isolates, identifiers, or credentials. Project continuity notes are in [docs/HANDOFF.md](docs/HANDOFF.md).
