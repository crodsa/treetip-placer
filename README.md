# TreeTip Placer

TreeTip Placer is a quality-controlled, distance-based placement tool for adding new MLST allele profiles to a fixed reference tree without recalculating the backbone topology.

The first implementation evaluates query quality, finds a locally restricted or exhaustive set of candidate edges, fits a weighted least-squares insertion for each candidate, and writes interoperable tree and tabular outputs. It accepts both MLST allele-profile tables and fixed-column FASTA alignments. The same aligned-character adapter is the foundation for concatenated MLST alleles, core SNPs, and concatenated core genes; maximum-likelihood scoring is the next planned backend.

## Current capabilities

- Validates a reference tree and either an MLST table or a fixed-column alignment.
- Applies pre-placement quality control to every query profile.
- Places valid profiles on a fixed tree using normalized allelic Hamming distances.
- Supports `fast`, `balanced`, and `thorough` candidate-search presets.
- Runs query placements in parallel process workers.
- Optionally estimates edge stability with locus bootstrap resampling.
- Writes Newick, Nexus, TSV, JSON, JSONL event logs, and an experimental distance-placement JPlace document.
- Keeps the reference topology unchanged and records every QC or placement warning.

## Installation

The MVP has no runtime dependencies outside the Python standard library.

```bash
python -m pip install -e .
```

## Required input

Reference and query profile files must be tab-separated tables. The first column is `sample_id`; every other column is an MLST locus. Allele identifiers are categorical strings.

```text
sample_id	abcZ	adk	aroE	fumC	gdh	pdhC	pgm
isolate_01	4	7	2	1	5	3	8
```

The Newick tree must contain every reference `sample_id` exactly once and every branch must have a non-negative numeric length.

### Concatenated alleles, core SNPs, and core genes

Use the FASTA adapter when the reference and query sequences have already been mapped or aligned to the same immutable columns. The program never modifies the reference alignment or attempts a de novo alignment.

```bash
treetip place \
  --tree reference_tree.nwk \
  --reference-alignment reference_concatenated.fasta \
  --query-alignment new_genomes_in_reference_coordinates.fasta \
  --output placement_run
```

For SNP-only matrices, preserve the original SNP selection, masking, and ascertainment-correction metadata alongside the reference package. The current distance backend does not replace a likelihood model; it supplies a reproducible placement baseline while the ML backend is integrated.

## Quick start

```bash
treetip place \
  --tree reference_tree.nwk \
  --reference-profiles reference_profiles.tsv \
  --query-profiles new_profiles.tsv \
  --output placement_run \
  --preset balanced \
  --threads auto \
  --bootstrap 100 \
  --verbose
```

`placed_tree.nwk` is directly usable in FigTree, iTOL, Dendroscope, ETE, and most Newick-compatible programs. Detailed placement uncertainty belongs in `placements.tsv` and `placements.json` rather than in overloaded Newick node labels.

## MLST leave-one-out validation

Before applying the tool to new samples, use known reference profiles as held-out queries. `validate-mlst` removes both the selected profile and its original tip from a copied backbone, suppresses the resulting unary node when needed, then compares the selected placement edge with the known attachment neighbourhood. This avoids the information leak of leaving the held-out tip in the tree.

```bash
treetip validate-mlst \
  --tree reference_tree.nwk \
  --reference-profiles reference_profiles.tsv \
  --output mlst_leave_one_out \
  --preset thorough \
  --threads auto \
  --bootstrap 100 \
  --verbose
```

`leave_one_out.tsv` is directly importable into R, Python, or spreadsheet software. `validation_summary.json` reports recovery rate and graph-edge distance to the known attachment; these are topology diagnostics, not likelihood supports. For large panels, use `--max-holdouts` with `--seed` for a reproducible subset. The validation protocol is in [docs/validation.md](docs/validation.md).

## Scientific scope

This command answers: *where should a query attach if the reference topology is held fixed?* It does not claim to reproduce a de novo phylogenetic reconstruction with the query added. Existing branch supports are preserved only as backbone metadata; they are not recalculated supports for the expanded taxon set.

Allele identifier numbers are never treated as numeric distances. A mismatch at a locus is one categorical difference; missing loci are excluded pairwise and reported by QC.

## Status policy

Every query receives one of the following statuses:

- `PASS`: adequate input and a resolved placement.
- `WARN`: placement is produced, with an input quality warning.
- `FAIL`: structurally invalid or insufficient input; no placement is produced unless `--force` is used.
- `UNRESOLVED`: input passed QC but placement alternatives remain too diffuse.

## Development and release

All source code, interface text, logs, schemas, tests, and documentation are English-only. Local real-world validation data are never included in the repository without explicit redistribution permission. The GitHub publication checklist is in [docs/release-checklist.md](docs/release-checklist.md).

The component design and planned likelihood/local-refinement extensions are described in [docs/architecture.md](docs/architecture.md).

For continuity across Codex chats or computers, read and update [docs/HANDOFF.md](docs/HANDOFF.md).
