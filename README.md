# TreeTip Placer

TreeTip Placer adds query sequences to an immutable, branch-length reference tree without rebuilding its topology. It is designed exclusively for fixed-column FASTA alignments: concatenated MLST allele sequences, concatenated core genes, or core-genome SNP matrices.

**Running this on your own isolates and not building or maintaining the tool?** See [`docs/MANUAL.md`](docs/MANUAL.md) for a plain-language guide to installing it, running `place-auto`, and reading the results. Everything below this point is written for whoever builds and maintains the reference packages and the tool itself.

## Bundled reference package

A ready-to-use paired classical/cryptic *Clostridioides difficile* reference package, built from concatenated 7-locus MLST alleles, ships in this repository at [`reference-packages/mlst-cdifficile-20260806/`](reference-packages/mlst-cdifficile-20260806/). It is committed with the repository owner's explicit permission — see `AGENTS.md`'s data-in-version-control rule and `docs/HANDOFF.md`'s "Current reference panel" section for full provenance, methodology, and calibration numbers. The package's own `README.md` has the quick-start command.

`treetip classify <query.fasta> --organism cdifficile-mlst` resolves this bundled package, an on-`PATH` MAFFT, and an auto-named output directory without any other flags — see `docs/MANUAL.md`. It is a thin wrapper around `place-auto` (below): `ORGANISM_REFERENCE_PACKAGES` in `cli.py` maps `--organism` keys to `reference-packages/` subdirectories, and adding a new bundled package only requires a new registry entry plus its directory, not a CLI change.

## Scope

The current backend performs quality-controlled, distance-based placement on a fixed reference tree. It does not replace a de novo phylogenetic reconstruction or a maximum-likelihood placement engine. The reference alignment, coordinate system, masking policy, and tree must be retained as one versioned reference package.

Tabular MLST allele-profile inputs are intentionally unsupported.

## Input contract

- A branch-length Newick reference tree.
- A fixed-column reference FASTA alignment whose IDs exactly match tree tips.
- A query FASTA alignment mapped to the identical reference columns.

The program will reject unequal alignment widths rather than silently realigning, trimming, padding, or dropping query sites. For SNP-only matrices, retain the original SNP-selection, masking, and ascertainment metadata beside the reference package.

## Automatic classical/cryptic selection

When the origin of an uploaded isolate is unknown, use `place-auto`. It maps the raw query FASTA independently to the immutable classical and cryptic alignments with MAFFT `--add --keeplength`, runs the complete QC/placement workflow for both panels, and then compares the nearest-reference normalized distances.

For routine use, provide one portable paired-reference manifest instead of four individual file paths. Relative paths are resolved beside the manifest and declared SHA-256 checksums are validated before mapping:

```json
{
  "schema_version": 1,
  "minimum_distance_separation": 0.01,
  "classical": {"tree": "classical/tree.nwk", "alignment": "classical/aligned.fasta"},
  "cryptic": {"tree": "cryptic/tree.nwk", "alignment": "cryptic/aligned.fasta"}
}
```

Create the manifest without hand-editing JSON:

```bash
treetip create-reference-package \
  --classical-tree classical/tree.nwk \
  --classical-reference-alignment classical/aligned.fasta \
  --cryptic-tree cryptic/tree.nwk \
  --cryptic-reference-alignment cryptic/aligned.fasta \
  --classical-tip-clades classical_tip_clades.tsv \
  --cryptic-tip-clades cryptic_tip_clades.tsv \
  --min-panel-distance-separation 0.01 \
  --output paired_reference.json
```

Validate a package before sharing or analyzing it:

```bash
treetip validate-reference-package --reference-package paired_reference.json
```

```bash
treetip place-auto \
  --reference-package paired_reference.json \
  --query-fasta uploaded_queries.fasta \
  --mafft /path/to/mafft \
  --output automatic_run \
  --bootstrap 100 \
  --verbose
```

The result is deliberately conservative. `panel_assignment.tsv` contains one `CLASSICAL`, `CRYPTIC`, `AMBIGUOUS`, or `UNRESOLVED` call per query. A call is made only when its nearest-reference distance is better than the other panel by `--min-panel-distance-separation` (default `0.01`). This provisional default was calibrated for the current paired panels and must be recalibrated with known classical, cryptic, and boundary samples whenever either reference package changes. `report.html` is a human-readable index with direct artifact links. For every query, it shows a full clade-context tree and a local placement zoom. The full view expands the query clade and collapses pure non-target clade subtrees as labelled `Clade (n tips)` branches; the zoom displays the query and its visible reference-tip names. Both are display-only and do not modify the fixed tree. `selected_trees/` contains one directly drawable Newick tree per selected panel; `panels/` retains both complete, auditable placement runs.

Synthetic report layout preview (no biological data or real identifiers):

![Synthetic placement report preview](docs/assets/placement-report-preview-synthetic.png)

Tip-clade metadata is optional and uses a two-column UTF-8 TSV with the exact headers `reference_tip_id` and `clade`. It can be checksummed into the paired reference package and then loads automatically; command-line TSV paths can explicitly override it for a run. If supplied, `place-auto` writes `clade_assignment.tsv` and JSON. A `PASS` clade call has one label among the annotated descendant reference tips on the selected placement edge, and that label agrees with the query's nearest labelled reference by direct genetic distance. On a weakly resolved reference tree (many near-zero-length branches, typical of a short MLST scheme) a selected edge's full descendant subtree can span tips that are topologically nearby but genetically distant; when the edge consensus and the nearest labelled reference disagree, the nearest-reference label is used instead as an explicit `WARN` (`SELECTED_EDGE_CONSENSUS_DISAGREES_WITH_NEAREST_REFERENCE`). Mixed or unlabelled edge labels also use a nearest-labelled-reference fallback as explicit `WARN`; unavailable labels remain `UNRESOLVED`. An `UNRESOLVED` placement never emits a clade label. These are discrete calls on a fixed reference tree, not a new phylogenetic reconstruction or bootstrap calculation.

For the intended ANI-defined seven-clade classical and cryptic panels, derive a complete reference-tip partition with FastBAPS in a supported environment, map its seven clusters to the independently defined ANI labels, and package the resulting TSV with `--classical-expected-clade-count 7 --cryptic-expected-clade-count 7`. The package validator then rejects any missing reference tip or any partition other than exactly seven labels. A query inherits the label from its fixed-tree placement; `WARN` or `UNRESOLVED` still records a placement/QC limitation rather than silently hiding it.

Validate the discrete clade rule without modifying the backbone with deterministic labelled leave-one-unique-sequence-out controls:

```bash
treetip calibrate-clade-assignment \
  --tree cryptic/tree.nwk \
  --reference-alignment cryptic/aligned.fasta \
  --tip-clades cryptic_tip_clades.tsv \
  --max-controls-per-clade 20 \
  --output cryptic_clade_calibration \
  --preset thorough \
  --verbose
```

For each control, the command removes its tip from the fixed tree and its sequence from the reference alignment, then applies the production placement and discrete clade rule. A clade with only one labelled reference tip is excluded because its label cannot remain after holdout; `clade_calibration_exclusions.tsv` records that fact. `clade_calibration_records.tsv` includes coverage, comparable sites, nearest-reference distance, relative fit weight, objective, and RMSE for every control; `clade_confusion.tsv` and `clade_calibration_summary.json` summarize `CORRECT_PASS`, `CORRECT_WARN`, `INCORRECT_CALL`, and `UNRESOLVED` overall and by truth clade. This is internal validation only; it does not recalculate a topology or phylogenetic support.

Calibrate the paired-panel distance gate with internal leave-one-unique-sequence-out controls:

```bash
treetip calibrate-panel-gate \
  --classical-reference-alignment classical_aligned.fasta \
  --cryptic-reference-alignment cryptic_aligned.fasta \
  --mafft /path/to/mafft \
  --max-controls-per-panel 50 \
  --output gate_calibration \
  --verbose
```

This writes `calibration_records.tsv`, `panel_label_conflicts.tsv`, and `calibration_summary.json`, including a fixed threshold sweep. The conflict table flags high-confidence calls that disagree with the control panel label for provenance review; it never changes source labels or reference membership. The command evaluates the panel-assignment gate only; it does not estimate external generalization, tree-placement accuracy, or phylogenetic support.

Measure whether leave-one-out placement recovers a held-out reference tip's own position:

```bash
treetip calibrate-placement-accuracy \
  --tree cryptic/tree.nwk \
  --reference-alignment cryptic/aligned.fasta \
  --max-controls 50 \
  --near-hop-threshold 2 \
  --output cryptic_placement_accuracy \
  --preset thorough \
  --verbose
```

For each control, the command removes its tip from the fixed tree and its sequence from the reference alignment, places it back through the production placement path, and compares the selected edge with the pruned tip's recorded attachment. `placement_accuracy_records.tsv` reports `CORRECT` for an exact recovery, `NEAR` within `--near-hop-threshold` edge hops, `FAR` beyond it, and `UNRESOLVED` when placement QC produced no defensible call. `placement_accuracy_summary.json` adds the full `hop_distribution`, so the threshold can be re-evaluated without repeating the run. The `--near-hop-threshold` default of `2` is provisional and uncalibrated; treat it as a starting point and set it from an observed hop distribution. This is an unweighted topology diagnostic on a fixed backbone: it is not phylogenetic support, not a branch-length error, and not evidence of external generalization.

Each record also carries `patristic_distance_from_truth`, the branch-length-weighted distance between the placed point and the true attachment, and the summary adds `patristic_distance_summary` (mean, median, and max over graded controls). Edge-hop count and patristic distance can disagree on a poorly resolved reference tree: a run with many near-zero-length branches can grade most controls `FAR` by hop count while `nearest_reference_distance` and the patristic distance both stay small, because a single true nearest-neighbor swap can cross several hops along a chain of branches that are genetically almost indistinguishable. Read the hop-based outcome and the patristic distance together before concluding a reference tree's topology is poorly recovered; a low-resolution locus set (a short MLST scheme, for example) is a more likely explanation than a placement defect when nearest-reference distances stay small.

Every calibration command also writes a `report.html` index carrying the run's provenance, scope warning, outcome counts, and per-control table; the TSV and JSON sidecars remain authoritative.

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
- Without requested bootstrap replicates, a diffuse relative-fit weight is reported as `WARN`, not a phylogenetic support value; low requested bootstrap support remains `UNRESOLVED`.
- Local NNI/SPR refinement and recomputed local support are not yet implemented.
- Cross-panel assignment is a distance gate, not a taxonomic classifier or a phylogenetic likelihood-ratio test.

## Development

The test suite uses only synthetic fixtures:

```powershell
python -m unittest discover -s tests -v
```

Do not commit private alignments, genomes, isolates, identifiers, or credentials. Project continuity notes are in [docs/HANDOFF.md](docs/HANDOFF.md).

## License

This private development repository is proprietary and research/development-only. Redistribution or public release requires explicit written permission; see [LICENSE](LICENSE).
