# Project handoff

Last updated: 2026-08-01

## Purpose

TreeTip Placer places aligned nucleotide sequences on an immutable reference tree without recalculating the backbone topology. Supported inputs are fixed-column FASTA alignments for concatenated MLST allele sequences, concatenated core genes, and core-genome SNP matrices.

## Current implementation

- Fixed-column FASTA input only; tabular allele-profile input is intentionally unsupported.
- Distance-based constrained least-squares placement with local or exhaustive candidate search.
- Pre-placement QC, optional site bootstrap, parallel query execution, and explicit progress logs.
- Newick, Nexus, TSV, JSON, JSONL, and experimental distance-JPlace outputs.
- `place-auto` maps an uploaded FASTA independently to paired classical and cryptic references with MAFFT `--add --keeplength`, runs both placements, and writes an auditable conservative panel assignment plus selected-panel trees.
- Tip-clade metadata can be checksummed inside a paired-reference package; package-bound metadata loads automatically unless a run explicitly supplies another TSV.
- `calibrate-clade-assignment` performs deterministic labelled leave-one-unique-sequence-out controls against a fixed tree, writing record-level outcomes and a clade confusion table without recalculating topology or phylogenetic support.

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

The production local reference pair is a separate classical package and the derived `cryptic_strict` package. No further record-level curation is planned for this workflow because the selected classical and cryptic panels are already non-overlapping.

The local stress test at `C:\Users\crodriguez\Desktop\treetip-placer-local-runs\auto_panel_strict_stress_test_20260731` confirmed two known cryptic queries as `CRYPTIC` and one known clade-4 classical control as `CLASSICAL` against the strict package.

The expanded leave-one-unique-sequence-out calibration at `C:\Users\crodriguez\Desktop\treetip-placer-local-runs\panel_gate_calibration_strict_20260731` supports a distance margin of `0.01` for this exact classical/strict-cryptic pair: 194/195 controls were called correctly and one was retained as ambiguous. This validates conservative gating on these internal panels only, not external generalization or placement-topology accuracy.

The portable paired manifest for this local production pair is `C:\Users\crodriguez\Desktop\treetip-placer-local-runs\paired_reference_created_20260731.json`. It was produced by `create-reference-package`, uses relative paths and checksums, and was smoke-tested through `place-auto --reference-package`.

`validate-reference-package` also passed for this manifest: the classical panel has 200 matching tree/alignment IDs across 3537 columns, and the cryptic panel has 134 matching IDs across 3501 columns.

Automatic runs now include a self-contained `report.html` index beside the authoritative tables and trees. For a selected panel, it shows a complete-tree context dendrogram and a branch zoom for each newly placed query; no unselected panel is shown. Report rendering is covered by the complete unit-test suite and was verified in a fresh permitted-data end-to-end smoke run at `C:\Users\crodriguez\Desktop\treetip-placer-local-runs\single_cryptic_report_context_smoke_20260731`: the query was assigned `CRYPTIC`, and the report contained the selected cryptic full-tree view and one placement zoom with no classical full-tree view.

The fixed-tree clade workflow was validated locally with the permitted supplementary `Table S2` metadata. The resulting non-versioned cryptic tip-clade TSV contains 105 included labels, 28 unlabelled reference tips, and one excluded non-cryptic table label. In the one-query end-to-end smoke run at `C:\Users\crodriguez\Desktop\treetip-placer-local-runs\single_cryptic_clade_warn_smoke_20260801`, the selected edge contained mixed clade labels, so the tool emitted an explicit `WARN` nearest-labelled-reference fallback instead of an unsupported `PASS` clade call. The output contained `clade_assignment.tsv`/JSON and an HTML report with only the selected cryptic tree.

The known classical `ST|clade` tree-tip convention was used to create the non-versioned `classical_tip_clades_reference_20260801.tsv` with 200 labels. The new package `C:\Users\crodriguez\Desktop\treetip-placer-local-runs\paired_reference_with_tip_clades_20260801.json` checksums both the classical and cryptic metadata TSVs and validated successfully: 200 classical tips across 3537 columns and 134 cryptic tips across 3501 columns.

`calibrate-clade-assignment` was run locally with five deterministic unique-sequence controls per represented clade. `classical_clade_calibration_20260801` yielded 25/25 clade-correct calls, all explicitly `CORRECT_WARN`; this requires interpretation before any confidence promotion. `cryptic_clade_calibration_20260801` yielded 22 `CORRECT_WARN`, four `UNRESOLVED`, and one `INCORRECT_CALL` across 27 controls. The unresolved and incorrect controls are retained in the local records for focused inspection; no labels were changed automatically.

## Next milestones

1. Build a labelled classical/cryptic calibration set to estimate the automatic panel-selection threshold and boundary error rate.
2. Review the cryptic `UNRESOLVED` and `INCORRECT_CALL` leave-one-out records, especially C-V/C-VI/C-VII boundary neighborhoods; keep unannotated tips unlabelled unless independently curated.
3. Decide whether systematic `WARN` calls in the classical calibration reflect appropriate conservative QC thresholds or require a separately validated confidence policy.
4. Add alignment leave-one-out validation and panel-pair calibration reports.

## Repository state

Private development backup: `https://github.com/crodsa/treetip-placer` on `main`. Real data, local analysis outputs, and supplementary input tables must remain outside Git. The fixed-tree clade-assignment interface and its internal validation command are implemented, but local clade-specific calibration must be reviewed before it can be considered production-ready.
