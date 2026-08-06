# Project handoff

Last updated: 2026-08-06

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
- `calibrate-placement-accuracy` grades leave-one-out placement geometry against each control's recorded attachment in edge hops and in branch-length-weighted (patristic) distance, writing per-control records, a hop distribution, a patristic-distance summary, and an HTML report. `Tree.patristic_distance_to_attachment` was added after a real cryptic-panel run (68 tips, MLST-scheme resolution) graded nearly every control `FAR` or `NEAR` by hop count alone while `nearest_reference_distance` stayed under 0.002 throughout: the hop count and the genetic distance recovered can disagree sharply on a tree with many near-zero-length branches, and the patristic distance is what makes that visible instead of reading as a placement failure.
- All three calibration commands write a `report.html` index beside their authoritative TSV and JSON sidecars.

## Important limitations

- The current scorer is distance-based, not maximum likelihood.
- Queries must already use the immutable reference alignment coordinates.
- Relative-fit weights are heuristic, not likelihood weight ratios or posterior probabilities.
- Reference-panel QC thresholds require local calibration and cannot be transferred blindly between cryptic-only and mixed cryptic–classical panels.
- The placement-accuracy grade is an unweighted topology diagnostic on a fixed backbone, not phylogenetic support or branch-length error. The `--near-hop-threshold` default of 2 is provisional and has not been calibrated against real panels. On a poorly resolved reference tree (many near-zero-length branches, typical of a short MLST scheme rather than core-genome or SNP-based alignment), the hop-based grade alone can look poor while `patristic_distance_from_truth` and `nearest_reference_distance` both stay small; read both before concluding placement itself is inaccurate.

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

`calibrate-clade-assignment` was run locally with five deterministic unique-sequence controls per represented clade. `classical_clade_calibration_20260801` yielded 25/25 clade-correct calls, all explicitly `CORRECT_WARN`; this requires interpretation before any confidence promotion. The first cryptic run showed one apparent incorrect C-VI call, but investigation established that C-VI has only one labelled reference tip. Version 0.6.3 therefore excludes singleton-labelled clades from leave-one-out evaluation and writes `clade_calibration_exclusions.tsv`. The corrected diagnostic at `cryptic_clade_calibration_edge_diagnostic_v2_20260801` evaluated 26 controls: 22 `CORRECT_WARN`, four `UNRESOLVED`, and zero incorrect calls; C-VI is explicitly excluded. The unresolved C-V/C-VII controls remain evidence of insufficient labelled edge resolution, and no labels were changed automatically.

Version 0.6.3 also prevents any underlying `UNRESOLVED` placement from emitting a discrete clade label. This fixes an overly permissive fallback that could otherwise present a nearest labelled tip despite unresolved placement geometry.

Version 0.6.4 changes the no-bootstrap interpretation of a diffuse WLS relative-fit weight from `UNRESOLVED` to explicit `WARN`. The weight is a heuristic, not phylogenetic support; when no bootstrap was requested, it must not silently make the tool unusable. Requested bootstrap support below the configured threshold remains `UNRESOLVED`. This policy was smoke-tested with the package-bound metadata at `single_cryptic_package_metadata_warn_smoke_20260801`: the query was selected as `CRYPTIC` (nearest-distance margin 0.04313) and received `C-VII` only as a `WARN` nearest-labelled-reference fallback on a mixed edge. The report and sidecar retain all reasons.

Version 0.6.5 extends `clade_calibration_records.tsv` with per-control coverage, comparable sites, nearest-reference distance, relative fit weight, objective, and RMSE; its JSON summary now also reports outcomes by truth clade. Exhaustive local leave-one-out calibration was retained outside Git under `classical_clade_calibration_full_metrics_20260801` and `cryptic_clade_calibration_full_metrics_20260801`. The classical panel evaluated all 200 controls: 194 `CORRECT_WARN` and six `INCORRECT_CALL` (97.0% retained-call agreement), all involving clades 1/2, with two clade-1 controls placing in clade 5. The cryptic panel evaluated 53 eligible unique-sequence controls: 48 `CORRECT_WARN` and five `UNRESOLVED`, with no incorrect calls. C-VI remains excluded because it has only one labelled tip; the 28 unlabelled cryptic tips and exact-sequence duplicates reduce the eligible control count. Distance, relative-fit weight, and RMSE did not supply a defensible simple threshold for suppressing only the six classical conflicts, so no unvalidated rule was added.

Version 0.6.6 improves `report.html` without changing any placement or clade result. The complete selected-tree SVG now displays every reference tip name. Each query additionally has a topology-preserving clade-context overview: the query clade stays expanded, while a non-target subtree is collapsed only when its annotated descendant tips all share one clade; it is labelled `Clade (n tips)`. Unlabelled pure subtrees are displayed as `Unlabelled (n tips)`. The local run `single_cryptic_clade_context_smoke_20260801` verified the full-tip labels, clade context, collapsed labels, and target C-VII expansion.

Version 0.6.7 inverts the report views following user review: the full tree is now the labelled clade-context overview, and the local placement zoom displays every visible reference tip name. The smoke run `single_cryptic_collapsed_fulltree_smoke_20260801` verified the target C-VII expansion, collapsed C-I labels, and visible reference names in the zoom.

Version 0.6.8 adds `expected_clade_count` to each optional panel definition. A complete ANI-defined seven-clade partition can now be required by passing `--classical-expected-clade-count 7 --cryptic-expected-clade-count 7` during reference-package creation. Validation rejects missing tip labels and a label cardinality other than seven; the exact-partition regression tests passed locally.

The target reference definition is seven ANI-defined clades in each panel: seven classical and seven cryptic. FastBAPS can provide a reproducible segmentation of each existing reference alignment/tree, but its automatic optimum must not be misrepresented as an ANI truth. The workflow therefore fixes the intended cardinality at seven, reconciles each FastBAPS cluster to independently defined ANI labels through reference anchors, and rejects a non-bijective or incomplete result for review. The package schema now supports `expected_clade_count`; with `7`, it requires a label for every reference tip and exactly seven distinct labels. Query placement then inherits this fixed reference partition, while `WARN` and `UNRESOLVED` continue to expose a QC or edge-boundary limitation rather than being erased. The current workstation has neither R nor the R `fastbaps` package, and WSL has no installed distribution; the available Python port documents Linux/macOS rather than Windows. Run FastBAPS later in a reproducible supported Linux/macOS environment on the private alignment, retain its output outside Git, and import only the validated tip-clade TSV into the private reference package.

## Next milestones

1. Build a labelled classical/cryptic calibration set to estimate the automatic panel-selection threshold and boundary error rate.
2. Review the cryptic `UNRESOLVED` leave-one-out records, especially C-V/C-VII boundary neighborhoods; keep unannotated tips unlabelled unless independently curated. C-VI requires at least one additional independently curated labelled reference before internal leave-one-out evaluation is possible.
3. Decide whether systematic `WARN` calls in the classical calibration reflect appropriate conservative QC thresholds or require a separately validated confidence policy.
4. Obtain an independent classical test set, especially clade-1/2 boundary isolates, before introducing any clade-specific confidence promotion or suppression rule.
5. Run broader local clade calibration after independently adding labels for currently singleton/unlabelled cryptic neighborhoods, particularly C-VI.
6. Run `calibrate-placement-accuracy` on the private classical and cryptic panels, then set `--near-hop-threshold` from the observed hop distribution instead of the provisional default of 2.
7. Run FastBAPS in a supported reproducible environment for each private reference alignment, constrain/reconcile the result to the seven ANI-defined clades per panel, and create a complete, checksummed reference-tip TSV with `expected_clade_count: 7`.
8. Enrich the production `place-auto` and `place` `report.html` with the per-query quantitative metrics that currently live only in the sidecar TSV — coverage, comparable sites, RMSE, objective, and bootstrap support — so an end user reads one document instead of opening tables alongside it.

## Repository state

Private development backup: `https://github.com/crodsa/treetip-placer` on `main`. Real data, local analysis outputs, and supplementary input tables must remain outside Git. The fixed-tree clade-assignment interface and its internal validation command are implemented, but local clade-specific calibration must be reviewed before it can be considered production-ready.
