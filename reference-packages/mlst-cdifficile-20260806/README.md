# mlst-cdifficile-20260806

A paired classical/cryptic *Clostridioides difficile* reference package for `treetip classify --organism cdifficile-mlst` and `treetip place-auto`, built from concatenated 7-locus MLST alleles (`adk`, `atpA`, `dxr`, `glyA`, `recA`, `sodA`, `tpi`; 3,501 bp).

## Usage

```bash
treetip classify your_isolates.fasta --organism cdifficile-mlst
```

Or with explicit paths via `place-auto`:

```bash
treetip place-auto \
  --reference-package reference-packages/mlst-cdifficile-20260806/paired_reference.json \
  --query-fasta your_isolates.fasta \
  --mafft /path/to/mafft \
  --output your_results
```

See `docs/MANUAL.md` for a step-by-step guide to running the workflow and interpreting the outputs.

## Contents

- `classical/` — 201 sequence types (STs) across all 5 classical MLST clades (`Clade-1`: 49, `Clade-2`: 49, `Clade-3`: 19, `Clade-4`: 49, `Clade-5`: 35), curated from PubMLST.
- `cryptic/` — 68 isolates across all 7 cryptic clades (`C-I`: 19, `C-II`: 10, `C-III`: 30 [`C-IIIa`: 20, `C-IIIb`: 10], `C-IV`: 2, `C-V`: 5, `C-VI`: 1, `C-VII`: 1), curated from EnteroBase and PubMLST genome assemblies.
- Both panels: 3,501-bp `mafft --auto` alignment and `IQ-TREE 3 -m MFP -B 1000 --alrt 1000 --bnni` maximum-likelihood reference tree.
- `paired_reference.json` — checksum-verified manifest with `minimum_distance_separation: 0.01`.

## Leave-one-out calibration

- **Panel gate** (`calibrate-panel-gate`): 201/201 classical (100.0%) and 68/68 cryptic (100.0%) leave-one-out controls separated cleanly across thresholds `0.00` to `0.01` (observed LOO margins `0.0323–0.0400` classical and `0.0200–0.0591` cryptic). The `0.01` margin in this package is the calibrated default for this paired reference.
- **Clade assignment** (`calibrate-clade-assignment`): 200/201 classical (99.50%) and 66/66 cryptic (100.00%, excluding the two singleton clades `C-VI` and `C-VII`, which cannot be leave-one-out evaluated with one reference tip each), corroborated by independent FastBAPS clustering (`99.50%` classical at `k = 5`, `100.00%` cryptic cluster purity at `k = 7`).
- **Placement accuracy** (`calibrate-placement-accuracy`): 201/201 classical and 68/68 cryptic leave-one-out placements succeeded with small branch-length-weighted patristic distances from the true pre-pruned attachment (classical mean `0.0030`, median `0.0024`, max `0.0100` on a `0.141`-length tree; cryptic mean `0.0251`, median `0.0273`, max `0.0365` on a `0.324`-length tree) and small nearest-reference p-distances (classical mean `0.0006`, max `0.0031`, `<= 0.002` in `199/201` controls; cryptic mean `0.0024`, median `0.0011`, max `0.0237`, `<= 0.002` in `52/68` controls).

## Scope and limitations

Seven MLST loci reliably separate classical from cryptic *Clostridioides* isolates and assign them to classical (`Clade-1`–`Clade-5`) and cryptic (`C-I`–`C-VII`) clades. They do not resolve fine-scale transmission or whole-genome differences within clonal complexes; see `docs/MANUAL.md` and `README.md`.
