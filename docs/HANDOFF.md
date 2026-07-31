# Project handoff

Last updated: 2026-07-31

## Purpose

TreeTip Placer places aligned nucleotide sequences on an immutable reference tree without recalculating the backbone topology. Supported inputs are fixed-column FASTA alignments for concatenated MLST allele sequences, concatenated core genes, and core-genome SNP matrices.

## Current implementation

- Fixed-column FASTA input only; tabular allele-profile input is intentionally unsupported.
- Distance-based constrained least-squares placement with local or exhaustive candidate search.
- Pre-placement QC, optional site bootstrap, parallel query execution, and explicit progress logs.
- Newick, Nexus, TSV, JSON, JSONL, and experimental distance-JPlace outputs.

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

## Next milestones

1. Map new batch sequences to the cryptic reference alignment and run placement.
2. Add alignment leave-one-out validation and calibrate QC thresholds for the cryptic panel.
3. Implement maximum-likelihood placement and partition-aware reference manifests.

## Repository state

Private backup: `https://github.com/crodsa/treetip-placer` on `main`. Real data and local analysis outputs must remain outside Git.
