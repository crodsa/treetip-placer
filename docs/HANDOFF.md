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

The local run at `C:\Users\crodriguez\Desktop\treetip-placer-local-runs\auto_panel_stress_test_20260731` confirmed two known cryptic queries as `CRYPTIC` and one known clade-4 classical control as `CLASSICAL`. The provisional default distance margin is `0.005`; it requires broader labelled calibration before a scientific release.

A local 80-control leave-one-unique-sequence-out gate calibration is available at `C:\Users\crodriguez\Desktop\treetip-placer-local-runs\panel_gate_calibration_command_fixed_20260731`. At margin `0.005`, 58/80 controls received a panel call, all 58 were correct, and 22 were retained as ambiguous; at margin `0`, one control was misclassified. The reproducible `calibrate-panel-gate` command now produces this report. This validates conservative gating on these internal panels only, not external generalization or placement-topology accuracy.

The expanded 196-control local run at `C:\Users\crodriguez\Desktop\treetip-placer-local-runs\panel_gate_calibration_98_per_panel_20260731` found one high-confidence label conflict at margin `0.005`. The command now writes `panel_label_conflicts.tsv` to surface such cases for provenance review without changing labels or panel membership.

## Next milestones

1. Build a labelled classical/cryptic calibration set to estimate the automatic panel-selection threshold and boundary error rate.
2. Add alignment leave-one-out validation and panel-pair calibration reports.
3. Implement maximum-likelihood placement and partition-aware reference manifests.

## Repository state

Private backup: `https://github.com/crodsa/treetip-placer` on `main`. Real data and local analysis outputs must remain outside Git.
