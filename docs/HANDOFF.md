# Project handoff

Last updated: 2026-07-31

## Purpose

TreeTip Placer places aligned nucleotide sequences on an immutable reference tree without recalculating the backbone topology. Supported inputs are fixed-column FASTA alignments for concatenated MLST allele sequences, concatenated core genes, and core-genome SNP matrices.

## Current implementation

- Fixed-column FASTA input only; tabular allele-profile input is intentionally unsupported.
- Distance-based constrained least-squares placement with local or exhaustive candidate search.
- Pre-placement QC, optional site bootstrap, parallel query execution, and explicit progress logs.
- Newick, Nexus, TSV, JSON, JSONL, and experimental distance-JPlace outputs.
- `place-auto` maps an uploaded FASTA independently to paired classical and cryptic references with MAFFT `--add --keeplength`, runs both placements, and writes an auditable conservative panel assignment plus selected-panel trees.

## Important limitations

- The current scorer is distance-based, not maximum likelihood.
- Queries must already use the immutable reference alignment coordinates.
- Relative-fit weights are heuristic, not likelihood weight ratios or posterior probabilities.
- Reference-panel QC thresholds require local calibration and cannot be transferred blindly between cryptic-only and mixed cryptic–classical panels.

## Verified command

```powershell
python -m unittest discover -s tests -v
```

## Local cryptic reference build

A local, non-versioned 135-sequence cryptic reference alignment and IQ-TREE 3 inference were created at `C:\Users\crodriguez\Desktop\treetip-placer-local-runs\cryptic_reference_iqtree3_20260731`. The source sequences remain intact. MAFFT produced 3501 columns and retained the twelve 3500-nt records with one gap each. IQ-TREE 3 used `-m MFP -B 1000 --alrt 1000 --bnni`.

## Automatic-panel stress test

The original cryptic reference panel contained a user-curated classical control. A derived local `cryptic_strict` package was created by pruning that tip without modifying the original data. The remaining six user-curated classical controls were not present in this source alignment and require their exact source panel or identifiers before they can be excluded.

The local stress test at `C:\Users\crodriguez\Desktop\treetip-placer-local-runs\auto_panel_strict_stress_test_20260731` confirmed two known cryptic queries as `CRYPTIC` and one known clade-4 classical control as `CLASSICAL` against the strict package.

The expanded leave-one-unique-sequence-out calibration at `C:\Users\crodriguez\Desktop\treetip-placer-local-runs\panel_gate_calibration_strict_20260731` supports a distance margin of `0.01` for this exact classical/strict-cryptic pair: 194/195 controls were called correctly and one was retained as ambiguous. This validates conservative gating on these internal panels only, not external generalization or placement-topology accuracy.

## Next milestones

1. Build a labelled classical/cryptic calibration set to estimate the automatic panel-selection threshold and boundary error rate.
2. Add alignment leave-one-out validation and panel-pair calibration reports.
3. Implement maximum-likelihood placement and partition-aware reference manifests.

## Repository state

Private backup: `https://github.com/crodsa/treetip-placer` on `main`. Real data and local analysis outputs must remain outside Git.
