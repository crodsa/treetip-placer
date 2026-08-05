# Placement-accuracy calibration and calibration reports

Date: 2026-08-05
Status: approved design, not yet implemented
Target version: 0.7.0

## Problem

`docs/HANDOFF.md` milestone 6 reads "Add alignment leave-one-out validation and
panel-pair calibration reports." It bundles two separable gaps.

**Gap A — placement accuracy is never measured.** Two calibration commands
exist. `calibrate-panel-gate` measures whether a control is routed to the
correct panel. `calibrate-clade-assignment` measures whether a control receives
the correct discrete clade label. Neither asks the prior question: when a
reference tip is removed and placed back, does the placement land where the tip
actually was? A clade call can be correct while the placement geometry is wrong,
because a clade spans many edges. Nothing in the repository currently
distinguishes those cases.

**Gap B — calibration output is not readable.** Each calibration command writes
TSV and JSON only. The panel-gate threshold sweep is the evidence behind the
`0.01` distance margin baked into every reference package, and reading it
requires opening three files and cross-referencing them by hand. The production
`place-auto` command has had an HTML index since 0.6.6; the calibration
commands have no equivalent.

## Existing infrastructure

The pruning and scoring machinery for Gap A is already built and tested; only
the driving command is missing.

- `Tree.without_tip(tip_name)` (`src/treetip_placer/tree.py:262`) returns a
  pruned tree plus a `PrunedTipAttachment` recording where the removed tip was
  attached, as either a surviving node or a point inside an edge. Its docstring
  states the record is intended for leave-one-out validation.
- `Tree.edge_distance_to_attachment(edge, attachment)`
  (`src/treetip_placer/tree.py:317`) returns the unweighted graph-edge distance
  between a selected edge and a true attachment, via the BFS in
  `Tree._node_distance` (`src/treetip_placer/tree.py:381`). It is covered by
  `tests/test_placement.py:97`.
- `select_unique_sequence_controls` (`src/treetip_placer/cli.py:1004`) chooses
  deterministic representatives of distinct aligned sequences.
- `place_query` (`src/treetip_placer/placement.py:17`) is the production
  placement entry point used unchanged by both existing calibration commands.

`calibrate-clade-assignment` already calls `without_tip` but discards the
attachment with `_` (`src/treetip_placer/cli.py:807`); it only checks the clade
label, never the geometry.

## Scope

In scope:

1. A new `calibrate-placement-accuracy` command (Gap A).
2. A shared HTML calibration report used by all three calibration commands
   (Gap B).
3. Deleting the dead `run_validate_mlst` function.

Out of scope, recorded here so it is not lost:

- Enriching the production `place-auto` / `place` `report.html` with the
  quantitative per-query metrics that today live only in the sidecar TSV
  (coverage, comparable sites, RMSE, objective, bootstrap support). This serves
  the end user of the laboratory workflow rather than internal calibration, and
  is the intended next spec after this one.
- Patristic (branch-length) placement error. The metric here is a topology
  diagnostic in edge hops. A branch-length error measure remains future work,
  in the same way `AGENTS.md` holds likelihood placement separate until
  implemented and validated.
- Any change to placement scoring, QC thresholds, panel gating, or clade rules.
  This design measures and displays; it does not alter a single call.

## Dead code removal

`run_validate_mlst` (`src/treetip_placer/cli.py:1338`) is a pre-pivot remnant.
It is unreachable — no subparser dispatches to it — and it cannot run: it calls
`read_profile_table`, `require_reference_profiles_complete`,
`select_holdout_ids`, `run_serial_validation`, `run_parallel_validation`,
`summarize_validation`, and `write_validation_outputs`, none of which are
imported or defined in the module, and it passes a `min_comparable_loci`
argument that `PlacementConfig` does not accept. It targets the tabular
allele-profile input that `README.md` documents as intentionally unsupported.

Delete the function. Delete `load_input_adapter`
(`src/treetip_placer/cli.py:1447`) and any other helper only reachable from it,
verified by grep before removal. Leaving a broken function whose name promises
exactly the validation this spec adds would mislead the next reader.

## Feature A: `calibrate-placement-accuracy`

### Command

```
treetip calibrate-placement-accuracy \
  --tree reference/tree.nwk \
  --reference-alignment reference/aligned.fasta \
  --max-controls 50 \
  --near-hop-threshold 2 \
  --output placement_accuracy_calibration \
  --preset thorough \
  --verbose
```

Arguments:

- `--tree` (required): reference branch-length Newick tree.
- `--reference-alignment` (required): reference fixed-column FASTA alignment.
- `--max-controls` (positive integer, default 50): maximum deterministic
  unique-sequence controls.
- `--near-hop-threshold` (positive integer, default 2): the largest edge-hop
  distance still reported as `NEAR`.
- `--output` (required): new, empty output directory.
- Everything `add_placement_controls()` supplies, matching
  `calibrate-clade-assignment`.

Unlike `calibrate-clade-assignment`, `--bootstrap` is accepted rather than
rejected. Correlating bootstrap support against measured placement accuracy is
the evidence HANDOFF milestone 3 needs before any confidence-promotion policy
can be justified.

Tree and alignment identifiers must match exactly, reusing the existing
mismatch-reporting style. Fewer than four reference tips is an error, matching
the existing `without_tip` guard.

### Per-control procedure

For each control, in deterministic order:

1. `pruned_tree, attachment = tree.without_tip(control_id)`.
2. Build `pruned_references` in pruned-tip order, exactly as
   `calibrate-clade-assignment` does.
3. `result = place_query(held_out, pruned_references, pruned_tree, config)`.
4. Check `result.final_status` **before** touching `result.best_fit`. A status
   of `FAIL` or `UNRESOLVED` grades the control `UNRESOLVED` with an empty
   `edge_hops_from_truth`, and so does a `None` `best_fit`. The status check
   must come first because `place_query` populates `best_fit` even when it
   returns `UNRESOLVED` (`src/treetip_placer/placement.py:142`); grading on
   `best_fit is None` alone would silently score placements the project's own
   0.6.3 policy says must not receive a discrete call.
5. A `WARN` status is graded normally. `WARN` marks a heuristic concern, not an
   unusable placement, and `placement_status` is carried into every record so
   `PASS` and `WARN` controls can be separated during analysis.
6. Otherwise compute
   `hops = pruned_tree.edge_distance_to_attachment(selected_edge, attachment)`
   and grade: `hops == 0` is `CORRECT`; `hops <= near_hop_threshold` is `NEAR`;
   otherwise `FAR`.

The selected edge is recovered from `pruned_tree.edges` by matching
`result.best_fit.edge_id`.

### Outputs

`placement_accuracy_records.tsv`, one row per control:

```
source_sample_id, outcome, edge_hops_from_truth, true_attachment_kind,
placement_status, placement_edge_id, coverage, comparable_sites,
nearest_reference_distance, relative_fit_weight, best_objective, best_rmse,
bootstrap_support
```

`placement_accuracy_summary.json`:

- `version`, `method`
  (`fixed_tree_leave_one_unique_sequence_out_placement_accuracy`), and an
  explicit `warning` that this is internal topology diagnostics only — not
  external generalization, not phylogenetic support, not branch-length error.
- Resolved paths and SHA-256 for tree and alignment; `reference_tips`;
  `reference_columns`; `max_controls`; `near_hop_threshold`.
- `summary`: control count and counts per outcome, keyed as
  `outcome_correct`, `outcome_near`, `outcome_far`, `outcome_unresolved`,
  following `summarize_clade_calibration`.
- `hop_distribution`: control counts keyed by observed hop count, so a
  `--near-hop-threshold` change can be re-evaluated without re-running.

Plus `run.log` and `report.html`, matching the existing commands.

## Feature B: shared calibration reports

### Module

A new `src/treetip_placer/calibration_reports.py`. Not `reports.py`: that module
is 472 lines of tree-SVG rendering for the production report, shares no
rendering logic with tabular calibration output, and serves a different reader.

### Interface

```python
def write_calibration_report(
    output_directory: str | Path,
    title: str,
    summary: Mapping[str, object],
    records: Sequence[Mapping[str, object]],
    record_columns: Sequence[str],
    outcome_column: str,
    problem_outcomes: frozenset[str],
    artifacts: Sequence[tuple[str, str]],
    threshold_sweep: Sequence[Mapping[str, object]] = (),
) -> None:
```

All three calibration commands call it. `threshold_sweep` is populated only by
`calibrate-panel-gate`; the sweep section is omitted when empty.

### Content

- Provenance header: version, method, control counts, and every path with its
  SHA-256, read from the `summary` mapping the command already builds.
- The `warning` string from the summary, rendered prominently rather than
  buried in JSON.
- Outcome counts as horizontal bars in plain HTML and CSS. No JavaScript, no
  external assets, consistent with the existing dependency-free report.
- For panel-gate only: the threshold sweep as a table with the default-threshold
  row highlighted, showing the accuracy-versus-coverage trade-off that justifies
  the packaged margin.
- The full per-control record table, with rows whose outcome is in
  `problem_outcomes` marked by CSS class — the pattern
  `render_decision_row` already uses for `AMBIGUOUS` and `UNRESOLVED`.
- Links to the TSV and JSON artifacts, stating they remain authoritative.

Per-command configuration:

| Command | Problem outcomes | Sweep |
| --- | --- | --- |
| `calibrate-panel-gate` | `INCORRECT_CALL` | yes |
| `calibrate-clade-assignment` | `INCORRECT_CALL` | no |
| `calibrate-placement-accuracy` | `FAR` | no |

All sample identifiers and free-text findings pass through `html.escape`, as
`test_report_escapes_sample_identifiers_and_links_artifacts` requires of the
production report.

## Testing

Test-driven, following `AGENTS.md`. New `tests/test_placement_accuracy.py` and
`tests/test_calibration_reports.py`; extend `tests/test_cli.py`.

Feature A:

- A control whose placement returns to its original attachment grades
  `CORRECT` with zero hops. Built on a small synthetic tree where the expected
  edge is known, in the style of
  `test_pruning_records_the_original_tip_attachment`.
- A placement one edge away grades `NEAR` under the default threshold and `FAR`
  under `--near-hop-threshold 0`, proving the threshold is honoured rather than
  hard-coded.
- A control failing pre-placement QC grades `UNRESOLVED` with an empty
  `edge_hops_from_truth`, and never `FAR`.
- A control whose placement returns `UNRESOLVED` **while carrying a non-`None`
  `best_fit`** grades `UNRESOLVED`. This is the regression that guards step 4;
  a status-blind implementation passes every other test and fails only this one.
- Both `PrunedTipAttachment` kinds — `node` and `edge` — are exercised, since
  `edge_distance_to_attachment` branches on kind.
- CLI end-to-end: the command writes both output files with the expected
  headers, and `hop_distribution` sums to the control count.

Feature B:

- The report escapes hostile sample identifiers.
- Rows with a problem outcome carry the marking CSS class; ordinary rows do not.
- The sweep section appears for a supplied sweep and is absent for an empty one.
- The `warning` string reaches the rendered document.
- Each of the three commands writes a `report.html`.

Full suite before handoff:

```bash
python -m unittest discover -s tests -v
```

## Documentation

- `README.md`: document `calibrate-placement-accuracy`, and note that
  calibration runs now emit `report.html`.
- `docs/HANDOFF.md`: record what milestone 6 delivered, and restate the
  remaining half — production-report metric enrichment — as its own milestone.
- `pyproject.toml`: version to 0.7.0.
- No real data, alignments, or local run outputs enter version control.

## Risks

**The accuracy grade could be misread as phylogenetic support.** It is an
unweighted topology diagnostic on a fixed backbone. Mitigation: the `warning`
field, the report header, and the README all state this, matching how the
project already frames relative-fit weight.

**`--near-hop-threshold` has no empirical basis yet.** The default of 2 is a
starting point, not a calibrated value. Mitigation: `hop_distribution` in the
summary lets the threshold be re-evaluated from a completed run without
re-computation, and the README labels the default provisional — the same
treatment `--min-panel-distance-separation` receives.

**Exhaustive controls on a 200-tip panel may be slow under `--preset
thorough`.** Mitigation: `--max-controls` defaults to 50, and the existing
progress logging reports every ten controls.
