# Placement Accuracy and Calibration Reports Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `calibrate-placement-accuracy` command that measures whether leave-one-out placement recovers a removed reference tip's original position, and give all three calibration commands a readable HTML report.

**Architecture:** Two new focused modules hold pure logic — `placement_accuracy.py` grades a placement against a recorded attachment in edge hops, and `calibration_reports.py` renders a dependency-free HTML report shared by every calibration command. `cli.py` gains only orchestration, and sheds two unreachable pre-pivot functions. No placement, QC, gating, or clade logic changes.

**Tech Stack:** Python 3.11+, standard library only (`unittest`, `csv`, `json`, `html.escape`, `dataclasses`). No new dependencies.

## Global Constraints

- Design source of truth: `docs/superpowers/specs/2026-08-05-placement-accuracy-and-calibration-reports-design.md`.
- All source, identifiers, comments, docstrings, CLI output, logs, schemas, tests, documentation, and Git history in English (`AGENTS.md`). User-facing chat may be Spanish.
- Never describe a placement result as equivalent to de novo phylogenetic inference. The accuracy grade is an unweighted topology diagnostic — not phylogenetic support, not branch-length error, not external generalization.
- No new runtime dependencies. Reports are self-contained HTML and CSS with no JavaScript and no external assets.
- No real genomic data, alignments, identifiers, or local run outputs enter version control.
- Full suite must pass before handoff: `python -m unittest discover -s tests -v`.
- Development environment: `source .venv/bin/activate` from the repository root.
- Target version: `0.7.0`.

## File Structure

| File | Responsibility |
| --- | --- |
| `src/treetip_placer/placement_accuracy.py` (create) | Pure grading and summarizing: `AccuracyGrade`, `grade_placement_accuracy`, `build_accuracy_record`, `summarize_placement_accuracy`, `hop_distribution`. No I/O, no argparse. |
| `src/treetip_placer/calibration_reports.py` (create) | `write_calibration_report` — dependency-free HTML shared by all three calibration commands. |
| `src/treetip_placer/cli.py` (modify) | Delete dead code; add the `calibrate-placement-accuracy` subparser and its orchestration; call the report writer from three commands. |
| `tests/test_placement_accuracy.py` (create) | Grading unit tests and the CLI end-to-end test for the new command. |
| `tests/test_calibration_reports.py` (create) | Report rendering tests. |
| `tests/test_cli.py` (modify) | Assert the two pre-existing calibration commands now emit `report.html`. |
| `README.md`, `docs/HANDOFF.md`, `pyproject.toml` (modify) | Documentation and version. |

`placement_accuracy.py` exists as its own module rather than growing `cli.py`, which is already 1679 lines — by far the largest file in a codebase whose next largest is 472. Keeping the grading pure also makes the `UNRESOLVED` regression in Task 2 testable without constructing a CLI run.

---

### Task 1: Remove unreachable pre-pivot code

`run_validate_mlst` and `load_input_adapter` target the tabular allele-profile input that `README.md` documents as intentionally unsupported. Neither is reachable: no subparser dispatches to `run_validate_mlst`, and nothing calls `load_input_adapter`. Neither could run — both call `read_profile_table`, `require_reference_profiles_complete`, and friends, none of which are defined or imported anywhere, and they reference `ProfileTable` and `Profile` types absent from `models.py`.

**Files:**
- Modify: `src/treetip_placer/cli.py` (delete `run_validate_mlst` at 1338-1444 and `load_input_adapter` at 1447-1477)

**Interfaces:**
- Consumes: nothing.
- Produces: nothing. This task only removes code.

- [ ] **Step 1: Confirm both functions are unreachable**

```bash
source .venv/bin/activate
grep -rn "run_validate_mlst\|load_input_adapter" src/ tests/ --include="*.py"
```

Expected: exactly two lines, both the `def` lines in `src/treetip_placer/cli.py`. No callers, no test references. If anything else appears, stop and report it rather than deleting.

- [ ] **Step 2: Confirm the helpers they call are genuinely undefined**

```bash
grep -n "def read_profile_table\|def select_holdout_ids\|def summarize_validation\|def write_validation_outputs\|def run_serial_validation\|def run_parallel_validation\|def require_reference_profiles_complete" src/treetip_placer/*.py
```

Expected: no output. These names appear only as call sites inside the two dead functions.

- [ ] **Step 3: Record the current test count as a baseline**

```bash
python -m unittest discover -s tests 2>&1 | tail -3
```

Expected: `Ran 30 tests` and `OK`. This is the number to preserve in Step 5.

- [ ] **Step 4: Delete both functions**

Delete the whole of `def run_validate_mlst(arguments: argparse.Namespace) -> int:` through its closing `return 0`, and the whole of `def load_input_adapter(` through the end of its body. Delete nothing else — the functions that follow `load_input_adapter` (`run_parallel_placements`, `run_serial_placements`, and the rest) are live and must stay.

- [ ] **Step 5: Verify nothing broke**

```bash
python -m unittest discover -s tests 2>&1 | tail -3
python -c "import treetip_placer.cli"
treetip --help
```

Expected: `Ran 30 tests` and `OK`; a clean import; and `--help` listing the six existing subcommands.

- [ ] **Step 6: Commit**

```bash
git add src/treetip_placer/cli.py
git commit -m "Remove unreachable pre-pivot profile-input functions

run_validate_mlst and load_input_adapter target the tabular allele-profile
input the README documents as intentionally unsupported. Neither is
dispatched or called, and neither could execute: both call helpers that are
defined nowhere and reference ProfileTable and Profile types absent from
models.py.

Leaving a broken function whose name promises leave-one-out validation would
mislead the next reader, since a working one is being added alongside it."
```

---

### Task 2: Grade a placement against its true attachment

The metric itself already exists: `Tree.edge_distance_to_attachment` (`src/treetip_placer/tree.py:317`) returns unweighted graph-edge distance and is covered by `tests/test_placement.py:97`. This task adds the grading policy around it.

The subtle requirement: `place_query` populates `best_fit` even when it returns `UNRESOLVED` (`src/treetip_placer/placement.py:142`). Status must be checked *before* `best_fit`, or unresolved placements get silently graded — violating the policy 0.6.3 established for clade calls.

**Files:**
- Create: `src/treetip_placer/placement_accuracy.py`
- Test: `tests/test_placement_accuracy.py`

**Interfaces:**
- Consumes: `PlacementResult` and `PlacementConfig` from `.models`; `Tree`, `Edge`, `PrunedTipAttachment` from `.tree`.
- Produces, relied on by Task 3:
  - `AccuracyGrade` — frozen dataclass with `outcome: str` and `edge_hops_from_truth: int | None`
  - `grade_placement_accuracy(result: PlacementResult, pruned_tree: Tree, attachment: PrunedTipAttachment, near_hop_threshold: int) -> AccuracyGrade`
  - `build_accuracy_record(result: PlacementResult, grade: AccuracyGrade, attachment: PrunedTipAttachment, source_sample_id: str) -> dict[str, object]`
  - `summarize_placement_accuracy(records: Sequence[Mapping[str, object]]) -> dict[str, int]`
  - `hop_distribution(records: Sequence[Mapping[str, object]]) -> dict[str, int]`
  - `ACCURACY_RECORD_COLUMNS: tuple[str, ...]`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_placement_accuracy.py`. Every fixture below was verified against the real placement code — the concordant tree grades all four controls `CORRECT`, and the discordant tree places held-out `A` exactly two hops from truth.

```python
"""Tests for leave-one-out placement-accuracy grading."""

from __future__ import annotations

import unittest

from treetip_placer.models import AlignedSequence, PlacementConfig, PlacementResult
from treetip_placer.placement import place_query
from treetip_placer.placement_accuracy import (
    grade_placement_accuracy,
    hop_distribution,
    summarize_placement_accuracy,
)
from treetip_placer.tree import Tree


CONCORDANT_SEQUENCES = {
    "A": "AAAAAA",
    "B": "AAAAAT",
    "C": "CCCCCC",
    "D": "CCCCCG",
}

# 'A' carries a C-clade sequence while the topology keeps it beside 'B',
# so its leave-one-out placement lands two edges from its true attachment.
DISCORDANT_SEQUENCES = {
    "A": "CCCCCC",
    "B": "AAAAAA",
    "C": "CCCCCC",
    "D": "CCCCCG",
}


def grade_holdout(newick: str, sequences: dict[str, object], held_out: str, near_hop_threshold: int):
    """Run one leave-one-out control and grade it.

    A sequence may be given as a string of characters or, to model absent data,
    as a tuple already containing None. Reading FASTA maps '-' and 'N' to None
    (`alignment.SEQUENCE_MISSING`), but an AlignedSequence built directly in a
    test keeps whatever characters it is handed — so a literal '-' string here
    would be scored as real data and grade WARN rather than failing QC.
    """

    tree = Tree.from_newick(newick)
    by_id = {
        name: AlignedSequence(name, tuple(text) if isinstance(text, str) else text)
        for name, text in sequences.items()
    }
    pruned, attachment = tree.without_tip(held_out)
    references = tuple(by_id[tip.name or ""] for tip in pruned.tips)
    config = PlacementConfig(
        preset="thorough", min_coverage=0.5, min_comparable_sites=2, unresolved_weight=0.0
    )
    result = place_query(by_id[held_out], references, pruned, config)
    return grade_placement_accuracy(result, pruned, attachment, near_hop_threshold), attachment


class PlacementAccuracyTests(unittest.TestCase):
    def test_recovered_attachment_grades_correct_with_zero_hops(self) -> None:
        grade, attachment = grade_holdout(
            "((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);", CONCORDANT_SEQUENCES, "A", 2
        )

        self.assertEqual(grade.outcome, "CORRECT")
        self.assertEqual(grade.edge_hops_from_truth, 0)
        self.assertEqual(attachment.kind, "edge")

    def test_node_kind_attachment_is_graded(self) -> None:
        # A polytomy leaves a node with two surviving children after pruning,
        # exercising the 'node' branch of edge_distance_to_attachment.
        grade, attachment = grade_holdout(
            "((A:0.1,B:0.1,E:0.1):0.3,(C:0.1,D:0.1):0.3);",
            {**CONCORDANT_SEQUENCES, "E": "AAAATT"},
            "A",
            2,
        )

        self.assertEqual(attachment.kind, "node")
        self.assertEqual(grade.outcome, "CORRECT")
        self.assertEqual(grade.edge_hops_from_truth, 0)

    def test_displaced_placement_grades_near_within_the_threshold(self) -> None:
        grade, _ = grade_holdout(
            "((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);", DISCORDANT_SEQUENCES, "A", 2
        )

        self.assertEqual(grade.edge_hops_from_truth, 2)
        self.assertEqual(grade.outcome, "NEAR")

    def test_the_same_placement_grades_far_under_a_stricter_threshold(self) -> None:
        grade, _ = grade_holdout(
            "((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);", DISCORDANT_SEQUENCES, "A", 1
        )

        self.assertEqual(grade.edge_hops_from_truth, 2)
        self.assertEqual(grade.outcome, "FAR")

    def test_failed_qc_grades_unresolved_without_hops(self) -> None:
        grade, _ = grade_holdout(
            "((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);",
            {**CONCORDANT_SEQUENCES, "A": (None,) * 6},
            "A",
            2,
        )

        self.assertEqual(grade.outcome, "UNRESOLVED")
        self.assertIsNone(grade.edge_hops_from_truth)

    def test_unresolved_status_is_not_graded_even_with_a_best_fit(self) -> None:
        """The regression guarding the status-before-best_fit ordering.

        place_query populates best_fit even when it returns UNRESOLVED, so an
        implementation testing only `best_fit is None` passes every other test
        here and fails exactly this one.
        """

        tree = Tree.from_newick("((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);")
        pruned, attachment = tree.without_tip("A")
        placed = place_query(
            AlignedSequence("A", tuple("AAAAAA")),
            tuple(AlignedSequence(tip.name or "", tuple(CONCORDANT_SEQUENCES[tip.name or ""])) for tip in pruned.tips),
            pruned,
            PlacementConfig(preset="thorough", min_coverage=0.5, min_comparable_sites=2, unresolved_weight=0.0),
        )
        self.assertIsNotNone(placed.best_fit)
        unresolved = PlacementResult(
            sample_id=placed.sample_id,
            qc=placed.qc,
            final_status="UNRESOLVED",
            best_fit=placed.best_fit,
            bootstrap_support=0.4,
        )

        grade = grade_placement_accuracy(unresolved, pruned, attachment, 2)

        self.assertEqual(grade.outcome, "UNRESOLVED")
        self.assertIsNone(grade.edge_hops_from_truth)

    def test_summaries_count_outcomes_and_hops(self) -> None:
        records = [
            {"outcome": "CORRECT", "edge_hops_from_truth": 0},
            {"outcome": "CORRECT", "edge_hops_from_truth": 0},
            {"outcome": "NEAR", "edge_hops_from_truth": 2},
            {"outcome": "UNRESOLVED", "edge_hops_from_truth": None},
        ]

        self.assertEqual(
            summarize_placement_accuracy(records),
            {"controls": 4, "outcome_correct": 2, "outcome_near": 1, "outcome_unresolved": 1},
        )
        self.assertEqual(hop_distribution(records), {"0": 2, "2": 1, "unresolved": 1})


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
source .venv/bin/activate
python -m unittest tests.test_placement_accuracy -v
```

Expected: collection fails with `ModuleNotFoundError: No module named 'treetip_placer.placement_accuracy'`.

- [ ] **Step 3: Write the implementation**

Create `src/treetip_placer/placement_accuracy.py`:

```python
"""Grade leave-one-out placements against their recorded true attachment.

The grade is an unweighted topology diagnostic on a fixed backbone. It is not
phylogenetic support, not a branch-length error, and not evidence of external
generalization.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

from .models import PlacementResult
from .tree import PrunedTipAttachment, Tree


ACCURACY_RECORD_COLUMNS: tuple[str, ...] = (
    "source_sample_id",
    "outcome",
    "edge_hops_from_truth",
    "true_attachment_kind",
    "placement_status",
    "placement_edge_id",
    "coverage",
    "comparable_sites",
    "nearest_reference_distance",
    "relative_fit_weight",
    "best_objective",
    "best_rmse",
    "bootstrap_support",
)

# A placement with either status must never receive an accuracy grade, matching
# the 0.6.3 policy that an unresolved placement cannot emit a discrete call.
UNGRADABLE_STATUSES = frozenset({"FAIL", "UNRESOLVED"})


@dataclass(frozen=True)
class AccuracyGrade:
    """One control's discrete accuracy outcome and its edge-hop distance."""

    outcome: str
    edge_hops_from_truth: int | None


def grade_placement_accuracy(
    result: PlacementResult,
    pruned_tree: Tree,
    attachment: PrunedTipAttachment,
    near_hop_threshold: int,
) -> AccuracyGrade:
    """Compare a placed edge against the held-out tip's true attachment.

    Status is checked before best_fit: place_query populates best_fit even when
    it returns UNRESOLVED, so grading on a missing best_fit alone would score
    placements the project treats as having no defensible discrete call.
    """

    if result.final_status in UNGRADABLE_STATUSES or result.best_fit is None:
        return AccuracyGrade(outcome="UNRESOLVED", edge_hops_from_truth=None)
    selected = pruned_tree.edge_by_child.get(result.best_fit.child_id)
    if selected is None:
        return AccuracyGrade(outcome="UNRESOLVED", edge_hops_from_truth=None)
    hops = pruned_tree.edge_distance_to_attachment(selected, attachment)
    if hops == 0:
        outcome = "CORRECT"
    elif hops <= near_hop_threshold:
        outcome = "NEAR"
    else:
        outcome = "FAR"
    return AccuracyGrade(outcome=outcome, edge_hops_from_truth=hops)


def build_accuracy_record(
    result: PlacementResult,
    grade: AccuracyGrade,
    attachment: PrunedTipAttachment,
    source_sample_id: str,
) -> dict[str, object]:
    """Build one auditable record row for a graded control."""

    best = result.best_fit
    return {
        "source_sample_id": source_sample_id,
        "outcome": grade.outcome,
        "edge_hops_from_truth": grade.edge_hops_from_truth,
        "true_attachment_kind": attachment.kind,
        "placement_status": result.final_status,
        "placement_edge_id": best.edge_id if best else None,
        "coverage": result.qc.coverage,
        "comparable_sites": result.qc.comparable_sites,
        "nearest_reference_distance": result.qc.nearest_distance,
        "relative_fit_weight": best.relative_fit_weight if best else None,
        "best_objective": best.objective if best else None,
        "best_rmse": best.rmse if best else None,
        "bootstrap_support": result.bootstrap_support,
    }


def summarize_placement_accuracy(
    records: Sequence[Mapping[str, object]],
) -> dict[str, int]:
    """Count controls per outcome using stable machine-readable keys."""

    counts: dict[str, int] = {"controls": len(records)}
    for record in records:
        key = f"outcome_{str(record['outcome']).lower()}"
        counts[key] = counts.get(key, 0) + 1
    return counts


def hop_distribution(records: Sequence[Mapping[str, object]]) -> dict[str, int]:
    """Count controls per observed hop distance.

    Retaining the raw distribution lets a near-hop threshold be re-evaluated
    from a completed run without repeating the placements.
    """

    counts: dict[str, int] = {}
    for record in records:
        hops = record["edge_hops_from_truth"]
        key = "unresolved" if hops is None else str(hops)
        counts[key] = counts.get(key, 0) + 1
    # Sort numerically, not lexically, so 10 hops follows 2 rather than preceding it.
    return dict(
        sorted(
            counts.items(),
            key=lambda item: (item[0] == "unresolved", int(item[0]) if item[0] != "unresolved" else 0),
        )
    )
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
python -m unittest tests.test_placement_accuracy -v
```

Expected: `Ran 7 tests` and `OK`.

- [ ] **Step 5: Run the full suite**

```bash
python -m unittest discover -s tests 2>&1 | tail -3
```

Expected: `Ran 37 tests` and `OK`.

- [ ] **Step 6: Commit**

```bash
git add src/treetip_placer/placement_accuracy.py tests/test_placement_accuracy.py
git commit -m "Grade leave-one-out placements against their true attachment

Tree.edge_distance_to_attachment already measured the edge-hop distance
between a selected edge and a pruned tip's recorded attachment, but nothing
consumed it. This adds the grading policy around it: zero hops is CORRECT,
within the near-hop threshold is NEAR, beyond it is FAR.

Status is checked before best_fit because place_query populates best_fit even
when it returns UNRESOLVED; grading on a missing best_fit alone would score
placements that the 0.6.3 policy says must not receive a discrete call."
```

---

### Task 3: Add the `calibrate-placement-accuracy` command

**Files:**
- Modify: `src/treetip_placer/cli.py` (import block; new subparser after the `calibrate-clade-assignment` block ending at line 188; dispatch in `main`; new `run_calibrate_placement_accuracy` and `write_placement_accuracy_outputs`)
- Test: `tests/test_placement_accuracy.py` (append a CLI class)

**Interfaces:**
- Consumes from Task 2: `ACCURACY_RECORD_COLUMNS`, `build_accuracy_record`, `grade_placement_accuracy`, `hop_distribution`, `summarize_placement_accuracy`.
- Consumes existing `cli.py` helpers: `prepare_output_directory`, `file_sha256`, `placement_config_from_arguments`, `select_unique_sequence_controls`, `add_placement_controls`, `EventLogger`, `format_event_log`, `append_final_event`, `positive_integer`.
- Produces, relied on by Task 5: a `summary` mapping containing `version`, `method`, `warning`, and paths; written to `placement_accuracy_summary.json`.

- [ ] **Step 1: Write the failing CLI test**

Append to `tests/test_placement_accuracy.py`:

```python
class PlacementAccuracyCliTests(unittest.TestCase):
    def test_command_writes_records_and_summary(self) -> None:
        import json
        import shutil
        from pathlib import Path

        from treetip_placer.cli import main

        root = Path.cwd() / ".treetip-placement-accuracy-test"
        if root.exists():
            shutil.rmtree(root)
        root.mkdir()
        try:
            tree = root / "tree.nwk"
            alignment = root / "reference.fasta"
            output = root / "output"
            tree.write_text("((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);\n", encoding="utf-8")
            alignment.write_text(
                ">A\nAAAAAA\n>B\nAAAAAT\n>C\nCCCCCC\n>D\nCCCCCG\n", encoding="utf-8"
            )

            exit_code = main(
                [
                    "calibrate-placement-accuracy",
                    "--tree", str(tree),
                    "--reference-alignment", str(alignment),
                    "--output", str(output),
                    "--preset", "thorough",
                    "--min-comparable-sites", "2",
                    "--unresolved-weight", "0",
                ]
            )

            self.assertEqual(exit_code, 0)
            records = (output / "placement_accuracy_records.tsv").read_text(encoding="utf-8")
            self.assertIn("edge_hops_from_truth", records)
            self.assertIn("true_attachment_kind", records)
            summary = json.loads((output / "placement_accuracy_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["summary"]["controls"], 4)
            self.assertEqual(summary["summary"]["outcome_correct"], 4)
            self.assertEqual(summary["near_hop_threshold"], 2)
            self.assertEqual(sum(summary["hop_distribution"].values()), 4)
            self.assertIn("not phylogenetic support", summary["warning"])
        finally:
            if root.exists():
                shutil.rmtree(root)
```

- [ ] **Step 2: Run it to verify it fails**

```bash
python -m unittest tests.test_placement_accuracy.PlacementAccuracyCliTests -v
```

Expected: FAIL — argparse rejects `invalid choice: 'calibrate-placement-accuracy'` and exits non-zero.

- [ ] **Step 3: Add the import**

In `src/treetip_placer/cli.py`, after the `from .placement import place_query` line:

```python
from .placement_accuracy import (
    ACCURACY_RECORD_COLUMNS,
    build_accuracy_record,
    grade_placement_accuracy,
    hop_distribution,
    summarize_placement_accuracy,
)
```

- [ ] **Step 4: Add the subparser**

In `build_parser`, immediately after `add_placement_controls(clade_calibration)` (line 188) and before the `calibrate.add_argument("--mafft", ...)` line:

```python
    accuracy = subparsers.add_parser(
        "calibrate-placement-accuracy",
        help="Measure whether leave-one-out placement recovers a held-out reference tip's position.",
    )
    accuracy.add_argument("--tree", required=True, help="Reference branch-length Newick tree.")
    accuracy.add_argument(
        "--reference-alignment", required=True, help="Reference fixed-column FASTA alignment."
    )
    accuracy.add_argument(
        "--max-controls",
        type=positive_integer,
        default=50,
        help="Maximum deterministic unique-sequence controls evaluated.",
    )
    accuracy.add_argument(
        "--near-hop-threshold",
        type=positive_integer,
        default=2,
        help="Largest edge-hop distance still reported as NEAR; provisional and uncalibrated.",
    )
    accuracy.add_argument("--output", required=True, help="New, empty output directory.")
    add_placement_controls(accuracy)
```

- [ ] **Step 5: Add the dispatch**

In `main`, after the `calibrate-clade-assignment` dispatch (line 318-319):

```python
        if arguments.command == "calibrate-placement-accuracy":
            return run_calibrate_placement_accuracy(arguments)
```

- [ ] **Step 6: Add the command implementation**

Insert after `run_calibrate_clade_assignment` and its helpers — placing it immediately before `def select_unique_sequence_controls(` keeps the calibration commands together:

```python
def run_calibrate_placement_accuracy(arguments: argparse.Namespace) -> int:
    """Measure whether leave-one-out placement recovers a held-out tip's position.

    Each control is removed from both the immutable backbone tree and the
    reference alignment, placed back with the production placement path, and
    graded in edge hops against its recorded attachment. This is an internal
    topology diagnostic, not phylogenetic support or external validation.
    """

    logger = EventLogger(arguments.verbose)
    logger.event(1, "start", "Starting leave-one-out placement-accuracy calibration", version=__version__)
    output = prepare_output_directory(arguments.output)
    tree = Tree.from_path(arguments.tree, arguments.max_cache_megabytes)
    if len(tree.tips) < 4:
        raise ValueError("Placement-accuracy calibration requires at least four reference tips")
    alignment = read_fasta_alignment(arguments.reference_alignment)
    reference_by_id = {sequence.sample_id: sequence for sequence in alignment.sequences}
    tree_tip_names = set(tree.tip_index)
    alignment_names = set(reference_by_id)
    if tree_tip_names != alignment_names:
        missing_in_alignment = sorted(tree_tip_names - alignment_names)
        missing_in_tree = sorted(alignment_names - tree_tip_names)
        details: list[str] = []
        if missing_in_alignment:
            details.append("tips absent from the alignment: " + ", ".join(missing_in_alignment[:10]))
        if missing_in_tree:
            details.append("alignment IDs absent from the tree: " + ", ".join(missing_in_tree[:10]))
        raise ValueError("Reference tree and alignment IDs differ; " + "; ".join(details))

    controls = tuple(
        sequence.sample_id
        for sequence in select_unique_sequence_controls(alignment.sequences, arguments.max_controls)
    )
    logger.event(
        1,
        "input",
        "Selected unique-sequence leave-one-out controls",
        reference_tips=len(tree.tips),
        controls=len(controls),
        near_hop_threshold=arguments.near_hop_threshold,
    )
    config = placement_config_from_arguments(arguments)
    records: list[dict[str, object]] = []
    for number, held_out_id in enumerate(controls, start=1):
        pruned_tree, attachment = tree.without_tip(held_out_id)
        pruned_references = tuple(
            reference_by_id[tip.name or ""] for tip in pruned_tree.tips
        )
        result = place_query(reference_by_id[held_out_id], pruned_references, pruned_tree, config)
        grade = grade_placement_accuracy(
            result, pruned_tree, attachment, arguments.near_hop_threshold
        )
        records.append(build_accuracy_record(result, grade, attachment, held_out_id))
        if number == len(controls) or number % 10 == 0:
            logger.event(
                2, "calibration", "Completed leave-one-out controls", completed=number, controls=len(controls)
            )

    summary = write_placement_accuracy_outputs(output, records, arguments, tree, alignment)
    logger.event(1, "calibration", "Completed placement-accuracy calibration", **summary["summary"])
    (output / "run.log").write_text(format_event_log(logger.events), encoding="utf-8")
    logger.event(1, "complete", "Placement-accuracy calibration completed")
    append_final_event(output, logger.events[-1])
    return 0


def write_placement_accuracy_outputs(
    output: Path,
    records: Sequence[dict[str, object]],
    arguments: argparse.Namespace,
    tree: Tree,
    alignment: object,
) -> dict[str, object]:
    """Write auditable per-control records and a provenance summary."""

    with (output / "placement_accuracy_records.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(ACCURACY_RECORD_COLUMNS), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for record in records:
            writer.writerow({key: "" if value is None else value for key, value in record.items()})
    summary = {
        "version": __version__,
        "method": "fixed_tree_leave_one_unique_sequence_out_placement_accuracy",
        "warning": (
            "Internal topology diagnostics only: this is not phylogenetic support, "
            "not a branch-length error, and not evidence of external generalization."
        ),
        "tree_path": str(Path(arguments.tree).resolve()),
        "tree_sha256": file_sha256(arguments.tree),
        "reference_alignment_path": str(Path(arguments.reference_alignment).resolve()),
        "reference_alignment_sha256": file_sha256(arguments.reference_alignment),
        "reference_tips": len(tree.tips),
        "reference_columns": getattr(alignment, "site_count"),
        "max_controls": arguments.max_controls,
        "near_hop_threshold": arguments.near_hop_threshold,
        "summary": summarize_placement_accuracy(records),
        "hop_distribution": hop_distribution(records),
    }
    (output / "placement_accuracy_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary
```

- [ ] **Step 7: Run the CLI test to verify it passes**

```bash
python -m unittest tests.test_placement_accuracy -v
```

Expected: `Ran 8 tests` and `OK`.

- [ ] **Step 8: Run the full suite**

```bash
python -m unittest discover -s tests 2>&1 | tail -3
```

Expected: `Ran 38 tests` and `OK`.

- [ ] **Step 9: Commit**

```bash
git add src/treetip_placer/cli.py tests/test_placement_accuracy.py
git commit -m "Add calibrate-placement-accuracy command

Prunes each unique-sequence control from the backbone tree and alignment,
places it back through the production placement path, and grades the
recovered edge against the recorded attachment.

This is the question neither existing calibration command asks. A clade call
can be correct while the placement geometry is wrong, because a clade spans
many edges. The summary retains the full hop distribution so the provisional
near-hop threshold can be re-evaluated without repeating the run."
```

---

### Task 4: Build the shared calibration report

**Files:**
- Create: `src/treetip_placer/calibration_reports.py`
- Test: `tests/test_calibration_reports.py`

**Interfaces:**
- Consumes: nothing from earlier tasks — this module is standalone by design.
- Produces, relied on by Task 5:
  `write_calibration_report(output_directory: str | Path, title: str, summary: Mapping[str, object], records: Sequence[Mapping[str, object]], record_columns: Sequence[str], outcome_column: str, problem_outcomes: frozenset[str], artifacts: Sequence[tuple[str, str]], threshold_sweep: Sequence[Mapping[str, object]] = ()) -> None`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_calibration_reports.py`:

```python
"""Tests for the shared calibration HTML report."""

from __future__ import annotations

import shutil
import unittest
from pathlib import Path

from treetip_placer.calibration_reports import write_calibration_report


SUMMARY = {
    "version": "0.7.0",
    "method": "fixed_tree_leave_one_unique_sequence_out_placement_accuracy",
    "warning": "Internal topology diagnostics only: this is not phylogenetic support.",
    "tree_path": "/private/tree.nwk",
    "tree_sha256": "abc123",
    "summary": {"controls": 2, "outcome_correct": 1, "outcome_far": 1},
}

COLUMNS = ("source_sample_id", "outcome", "edge_hops_from_truth")

RECORDS = (
    {"source_sample_id": "safe_sample", "outcome": "CORRECT", "edge_hops_from_truth": 0},
    {"source_sample_id": "sample<unsafe>", "outcome": "FAR", "edge_hops_from_truth": 7},
)

ARTIFACTS = (("placement_accuracy_records.tsv", "Per-control records (TSV)"),)


class CalibrationReportTests(unittest.TestCase):
    def setUp(self) -> None:
        self.root = Path.cwd() / ".treetip-calibration-report-test"
        if self.root.exists():
            shutil.rmtree(self.root)
        self.root.mkdir()

    def tearDown(self) -> None:
        if self.root.exists():
            shutil.rmtree(self.root)

    def render(self, **overrides) -> str:
        parameters = {
            "output_directory": self.root,
            "title": "Placement-accuracy calibration",
            "summary": SUMMARY,
            "records": RECORDS,
            "record_columns": COLUMNS,
            "outcome_column": "outcome",
            "problem_outcomes": frozenset({"FAR"}),
            "artifacts": ARTIFACTS,
        }
        parameters.update(overrides)
        write_calibration_report(**parameters)
        return (self.root / "report.html").read_text(encoding="utf-8")

    def test_escapes_sample_identifiers(self) -> None:
        document = self.render()

        self.assertIn("sample&lt;unsafe&gt;", document)
        self.assertNotIn("sample<unsafe>", document)

    def test_marks_only_problem_rows(self) -> None:
        document = self.render()

        self.assertEqual(document.count('class="problem"'), 1)

    def test_renders_the_scope_warning_and_provenance(self) -> None:
        document = self.render()

        self.assertIn("not phylogenetic support", document)
        self.assertIn("abc123", document)
        self.assertIn("0.7.0", document)

    def test_threshold_sweep_is_omitted_when_empty(self) -> None:
        document = self.render()

        self.assertNotIn("Threshold sweep", document)

    def test_threshold_sweep_highlights_the_default_row(self) -> None:
        document = self.render(
            summary={**SUMMARY, "default_threshold": 0.01},
            threshold_sweep=(
                {"threshold": 0.0, "controls": 2, "called": 2, "called_accuracy": 0.5},
                {"threshold": 0.01, "controls": 2, "called": 1, "called_accuracy": 1.0},
            ),
        )

        self.assertIn("Threshold sweep", document)
        self.assertEqual(document.count('class="default-threshold"'), 1)

    def test_links_every_artifact(self) -> None:
        document = self.render()

        self.assertIn('href="placement_accuracy_records.tsv"', document)


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
python -m unittest tests.test_calibration_reports -v
```

Expected: `ModuleNotFoundError: No module named 'treetip_placer.calibration_reports'`.

- [ ] **Step 3: Write the implementation**

Create `src/treetip_placer/calibration_reports.py`:

```python
"""Dependency-free human-readable reports for internal calibration runs.

These reports serve the maintainer calibrating the tool, not the end user of a
placement run; `reports.py` renders the production placement report and shares
no logic with this module.
"""

from __future__ import annotations

from html import escape
from pathlib import Path
from typing import Mapping, Sequence


PROVENANCE_KEYS = (
    "version",
    "method",
    "reference_tips",
    "reference_columns",
    "max_controls",
    "max_controls_per_clade",
    "max_controls_per_panel",
    "near_hop_threshold",
    "default_threshold",
    "tree_path",
    "tree_sha256",
    "reference_alignment_path",
    "reference_alignment_sha256",
    "classical_reference_alignment_sha256",
    "cryptic_reference_alignment_sha256",
    "tip_clades_path",
    "tip_clades_sha256",
)


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
    """Write a self-contained HTML index beside the authoritative sidecars."""

    root = Path(output_directory)
    document = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>TreeTip Placer calibration report</title>
  <style>
    body {{ font-family: system-ui, sans-serif; margin: 2rem; color: #1f2937; }}
    h1, h2 {{ color: #111827; }}
    table {{ border-collapse: collapse; width: 100%; font-size: 0.9rem; }}
    th, td {{ border: 1px solid #d1d5db; padding: 0.45rem; text-align: left; vertical-align: top; }}
    th {{ background: #f3f4f6; }}
    tr.problem td {{ background: #fef3c7; font-weight: 600; }}
    tr.default-threshold td {{ background: #e0f2fe; font-weight: 600; }}
    .warning {{ border-left: 4px solid #b45309; background: #fffbeb; padding: 0.75rem 1rem; margin: 1rem 0; }}
    .bar-track {{ background: #f3f4f6; border: 1px solid #d1d5db; height: 1.1rem; width: 100%; }}
    .bar-fill {{ background: #2563eb; height: 100%; }}
    code {{ word-break: break-all; }}
    dt {{ font-weight: 600; margin-top: 0.4rem; }}
    dd {{ margin: 0 0 0 1rem; }}
  </style>
</head>
<body>
  <h1>{escape(title)}</h1>
  <p>This report is a human-readable index. The TSV and JSON sidecars remain authoritative.</p>
  {render_warning(summary)}
  <h2>Provenance</h2>
  <dl>{render_provenance(summary)}</dl>
  <h2>Outcomes</h2>
  {render_outcome_bars(summary)}
  {render_threshold_sweep(threshold_sweep, summary)}
  <h2>Controls</h2>
  <table>
    <thead><tr>{"".join(f"<th>{escape(column)}</th>" for column in record_columns)}</tr></thead>
    <tbody>{render_record_rows(records, record_columns, outcome_column, problem_outcomes)}</tbody>
  </table>
  <h2>Artifacts</h2>
  <ul>{"".join(f'<li><a href="{escape(name)}">{escape(label)}</a></li>' for name, label in artifacts)}</ul>
</body>
</html>
"""
    (root / "report.html").write_text(document, encoding="utf-8")


def render_warning(summary: Mapping[str, object]) -> str:
    """Render the run's scope warning prominently rather than only in JSON."""

    warning = summary.get("warning")
    if not warning:
        return ""
    return f'<p class="warning">{escape(str(warning))}</p>'


def render_provenance(summary: Mapping[str, object]) -> str:
    """Render the provenance fields the command already recorded."""

    items = []
    for key in PROVENANCE_KEYS:
        if key not in summary:
            continue
        items.append(
            f"<dt>{escape(key)}</dt><dd><code>{escape(format_scalar(summary[key]))}</code></dd>"
        )
    return "".join(items)


def render_outcome_bars(summary: Mapping[str, object]) -> str:
    """Render outcome counts as plain HTML and CSS bars."""

    counts = summary.get("summary")
    if not isinstance(counts, Mapping):
        return "<p>No outcome counts were recorded.</p>"
    total = int(counts.get("controls", 0) or 0)
    rows = []
    for key in sorted(counts):
        if key == "controls":
            continue
        value = int(counts[key] or 0)
        percentage = (value / total * 100) if total else 0.0
        rows.append(
            f"<tr><td>{escape(key)}</td><td>{value}</td>"
            f'<td><div class="bar-track"><div class="bar-fill" style="width: {percentage:.1f}%"></div></div></td></tr>'
        )
    return (
        f"<p>Controls: <strong>{total}</strong></p>"
        f"<table><thead><tr><th>Outcome</th><th>Controls</th><th>Share</th></tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table>"
    )


def render_threshold_sweep(
    sweep: Sequence[Mapping[str, object]], summary: Mapping[str, object]
) -> str:
    """Render the accuracy-versus-coverage trade-off behind a packaged margin."""

    if not sweep:
        return ""
    default_threshold = summary.get("default_threshold")
    columns = list(sweep[0])
    rows = []
    for row in sweep:
        is_default = default_threshold is not None and row.get("threshold") == default_threshold
        css = ' class="default-threshold"' if is_default else ""
        cells = "".join(f"<td>{escape(format_scalar(row.get(column)))}</td>" for column in columns)
        rows.append(f"<tr{css}>{cells}</tr>")
    header = "".join(f"<th>{escape(column)}</th>" for column in columns)
    return (
        "<h2>Threshold sweep</h2>"
        "<p>The highlighted row is the threshold this run treated as the default.</p>"
        f"<table><thead><tr>{header}</tr></thead><tbody>{''.join(rows)}</tbody></table>"
    )


def render_record_rows(
    records: Sequence[Mapping[str, object]],
    record_columns: Sequence[str],
    outcome_column: str,
    problem_outcomes: frozenset[str],
) -> str:
    """Render one escaped row per control, marking problem outcomes."""

    rows = []
    for record in records:
        outcome = str(record.get(outcome_column, ""))
        css = ' class="problem"' if outcome in problem_outcomes else ""
        cells = "".join(
            f"<td>{escape(format_scalar(record.get(column)))}</td>" for column in record_columns
        )
        rows.append(f"<tr{css}>{cells}</tr>")
    return "".join(rows)


def format_scalar(value: object) -> str:
    """Format one cell with stable, locale-independent numbers."""

    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, float):
        return format(value, ".12g")
    return str(value)
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
python -m unittest tests.test_calibration_reports -v
```

Expected: `Ran 6 tests` and `OK`.

- [ ] **Step 5: Run the full suite**

```bash
python -m unittest discover -s tests 2>&1 | tail -3
```

Expected: `Ran 44 tests` and `OK`.

- [ ] **Step 6: Commit**

```bash
git add src/treetip_placer/calibration_reports.py tests/test_calibration_reports.py
git commit -m "Add a shared dependency-free calibration report

Calibration runs wrote TSV and JSON only. Reading the panel-gate threshold
sweep that justifies the packaged distance margin meant opening three files
and cross-referencing them by hand.

One renderer serves all three calibration commands: provenance, the scope
warning shown rather than buried, outcome counts as CSS bars, an optional
threshold sweep with the default row highlighted, and per-control rows with
problem outcomes marked."
```

---

### Task 5: Emit the report from all three calibration commands

**Files:**
- Modify: `src/treetip_placer/cli.py` (import; end of `run_calibrate_panel_gate`; end of `write_clade_calibration_outputs`; end of `write_placement_accuracy_outputs`)
- Modify: `tests/test_cli.py` (assert `report.html` for the two pre-existing commands)
- Modify: `tests/test_placement_accuracy.py` (assert `report.html` for the new command)

**Interfaces:**
- Consumes from Task 4: `write_calibration_report`.
- Consumes from Task 2: `ACCURACY_RECORD_COLUMNS`.
- Produces: nothing new; each command gains a `report.html`.

- [ ] **Step 1: Write the failing assertions**

In `tests/test_cli.py`, inside `test_clade_calibration_writes_auditable_outputs`, after the existing `clade_calibration_summary.json` assertion:

```python
            self.assertTrue((output / "report.html").is_file())
            report = (output / "report.html").read_text(encoding="utf-8")
            self.assertIn("Clade-assignment calibration", report)
```

In `tests/test_placement_accuracy.py`, inside `test_command_writes_records_and_summary`, after the `warning` assertion:

```python
            report = (output / "report.html").read_text(encoding="utf-8")
            self.assertIn("Placement-accuracy calibration", report)
            self.assertIn("not phylogenetic support", report)
```

- [ ] **Step 2: Run them to verify they fail**

```bash
python -m unittest tests.test_cli tests.test_placement_accuracy 2>&1 | tail -5
```

Expected: two failures — `FileNotFoundError` for the missing `report.html`.

- [ ] **Step 3: Add the import**

In `src/treetip_placer/cli.py`, after the `from .alignment import (...)` block:

```python
from .calibration_reports import write_calibration_report
```

- [ ] **Step 4: Return the rows that `write_calibration_records` already builds**

`write_calibration_records` currently builds each row, writes it, and returns
only a conflict count — so the report would otherwise have to re-read the TSV
it just wrote. Have it return the rows alongside the count instead.

Change its signature and its `return`:

```python
def write_calibration_records(
    output_directory: Path,
    decisions: Sequence[object],
    truth_by_id: dict[str, str],
    source_id_by_id: dict[str, str],
) -> tuple[int, list[dict[str, object]]]:
    """Write one auditable calibration record per held-out panel control.

    Returns the flagged-conflict count and the rows written, so a report can
    render them without re-reading the file.
    """
```

Inside the function, collect each row as it is written by adding
`written_rows: list[dict[str, object]] = []` beside the existing
`conflicts: list[dict[str, object]] = []`, appending `row` to `written_rows`
immediately after `writer.writerow(row)`, and changing the final
`return len(conflicts)` to `return len(conflicts), written_rows`.

Then update the one call site in `run_calibrate_panel_gate`:

```python
    conflict_count, calibration_rows = write_calibration_records(
        output, decisions, truth_by_id, source_id_by_id
    )
```

- [ ] **Step 5: Wire the panel-gate report**

In `run_calibrate_panel_gate`, immediately after the `calibration_summary.json` write and before the `logger.event(1, "calibration", ...)` call:

```python
    write_calibration_report(
        output,
        "Panel-gate calibration",
        summary,
        calibration_rows,
        (
            "calibration_id", "source_sample_id", "truth", "assignment",
            "decision_status", "outcome", "classical_nearest_distance",
            "cryptic_nearest_distance", "distance_margin",
        ),
        "outcome",
        frozenset({"INCORRECT_CALL"}),
        (
            ("calibration_records.tsv", "Per-control records (TSV)"),
            ("panel_label_conflicts.tsv", "Flagged label conflicts (TSV)"),
            ("calibration_summary.json", "Calibration summary (JSON)"),
        ),
        threshold_sweep=sweep,
    )
```

- [ ] **Step 6: Wire the clade-assignment report**

At the end of `write_clade_calibration_outputs`, after the `clade_calibration_summary.json` write:

```python
    write_calibration_report(
        output,
        "Clade-assignment calibration",
        summary,
        records,
        fields,
        "outcome",
        frozenset({"INCORRECT_CALL"}),
        (
            ("clade_calibration_records.tsv", "Per-control records (TSV)"),
            ("clade_confusion.tsv", "Clade confusion counts (TSV)"),
            ("clade_calibration_exclusions.tsv", "Excluded singleton clades (TSV)"),
            ("clade_calibration_summary.json", "Calibration summary (JSON)"),
        ),
    )
```

- [ ] **Step 7: Wire the placement-accuracy report**

At the end of `write_placement_accuracy_outputs`, between the `placement_accuracy_summary.json` write and `return summary`:

```python
    write_calibration_report(
        output,
        "Placement-accuracy calibration",
        summary,
        records,
        ACCURACY_RECORD_COLUMNS,
        "outcome",
        frozenset({"FAR"}),
        (
            ("placement_accuracy_records.tsv", "Per-control records (TSV)"),
            ("placement_accuracy_summary.json", "Calibration summary (JSON)"),
        ),
    )
```

- [ ] **Step 8: Run the full suite**

```bash
python -m unittest discover -s tests 2>&1 | tail -3
```

Expected: `Ran 44 tests` and `OK`.

- [ ] **Step 9: Inspect one rendered report by eye**

```bash
rm -rf /tmp/ttp-report-smoke && mkdir -p /tmp/ttp-report-smoke
printf '((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);\n' > /tmp/ttp-report-smoke/tree.nwk
printf '>A\nAAAAAA\n>B\nAAAAAT\n>C\nCCCCCC\n>D\nCCCCCG\n' > /tmp/ttp-report-smoke/ref.fasta
treetip calibrate-placement-accuracy --tree /tmp/ttp-report-smoke/tree.nwk --reference-alignment /tmp/ttp-report-smoke/ref.fasta --output /tmp/ttp-report-smoke/out --preset thorough --min-comparable-sites 2 --unresolved-weight 0
```

Expected: exit 0, and `/tmp/ttp-report-smoke/out/report.html` opens in a browser showing the warning box, provenance with SHA-256 values, four `CORRECT` controls with a full-width bar, and no threshold-sweep section. Confirm no unescaped markup and no broken artifact links.

- [ ] **Step 10: Commit**

```bash
git add src/treetip_placer/cli.py tests/test_cli.py tests/test_placement_accuracy.py
git commit -m "Emit a calibration report from all three calibration commands

Panel-gate, clade-assignment, and placement-accuracy runs now write
report.html beside their sidecars. Only panel-gate supplies a threshold
sweep, so the sweep section is omitted elsewhere.

No calculation, threshold, or call changes: this makes existing output
readable, as 0.6.6 did for the production report."
```

---

### Task 6: Document the release

**Files:**
- Modify: `README.md`, `docs/HANDOFF.md`, `pyproject.toml`

**Interfaces:**
- Consumes: the finished behavior of Tasks 1-5.
- Produces: nothing consumed by later tasks.

- [ ] **Step 1: Bump the version**

In `pyproject.toml`, change `version = "0.6.8"` to `version = "0.7.0"`.

- [ ] **Step 2: Document the command in `README.md`**

After the `calibrate-panel-gate` section:

````markdown
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

Every calibration command also writes a `report.html` index carrying the run's provenance, scope warning, outcome counts, and per-control table; the TSV and JSON sidecars remain authoritative.
````

- [ ] **Step 3: Update `docs/HANDOFF.md`**

Set `Last updated: 2026-08-05`. Add to **Current implementation**:

```markdown
- `calibrate-placement-accuracy` grades leave-one-out placement geometry against each control's recorded attachment in edge hops, writing per-control records, a hop distribution, and an HTML report.
- All three calibration commands write a `report.html` index beside their authoritative TSV and JSON sidecars.
```

Add to **Important limitations**:

```markdown
- The placement-accuracy grade is an unweighted topology diagnostic on a fixed backbone, not phylogenetic support or branch-length error. The `--near-hop-threshold` default of 2 is provisional and has not been calibrated against real panels.
```

`docs/HANDOFF.md` currently lists milestones 1 through 7. Replace the text of milestone 6 with the first item below, keep milestones 1-5 and 7 unchanged, and append the second item as a new milestone 8:

```markdown
6. Run `calibrate-placement-accuracy` on the private classical and cryptic panels, then set `--near-hop-threshold` from the observed hop distribution instead of the provisional default of 2.
8. Enrich the production `place-auto` and `place` `report.html` with the per-query quantitative metrics that currently live only in the sidecar TSV — coverage, comparable sites, RMSE, objective, and bootstrap support — so an end user reads one document instead of opening tables alongside it.
```

- [ ] **Step 4: Verify the documented command runs exactly as written**

```bash
source .venv/bin/activate
python -m unittest discover -s tests 2>&1 | tail -3
treetip calibrate-placement-accuracy --help
python -c "import tomllib, pathlib; print(tomllib.loads(pathlib.Path('pyproject.toml').read_text())['project']['version'])"
```

Expected: `Ran 44 tests` and `OK`; help text listing `--near-hop-threshold` and `--max-controls`; and `0.7.0`.

- [ ] **Step 5: Confirm no data was staged**

```bash
git status --short
```

Expected: only `README.md`, `docs/HANDOFF.md`, and `pyproject.toml`. No FASTA, Newick, TSV, or run output. If any appears, do not commit it — it is covered by the constraint against tracking real data.

- [ ] **Step 6: Commit**

```bash
git add README.md docs/HANDOFF.md pyproject.toml
git commit -m "Document placement-accuracy calibration for 0.7.0

Records what HANDOFF milestone 6 delivered and restates its remaining half —
enriching the production report with per-query quantitative metrics — as its
own milestone rather than leaving it implied.

States plainly that the near-hop threshold default is provisional and that
the accuracy grade is a topology diagnostic, not phylogenetic support."
```

---

## Verification

After Task 6, the repository should show:

- `python -m unittest discover -s tests -v` reporting `Ran 44 tests` and `OK`.
- Seven CLI subcommands, `calibrate-placement-accuracy` among them.
- No reference anywhere to `run_validate_mlst` or `load_input_adapter`.
- `report.html` written by all three calibration commands.
- Version `0.7.0`, with no real data tracked.

Pushing to `crodsa/treetip-placer` requires the private-checkpoint gates in `docs/release-checklist.md` and the repository owner's confirmation; ask before pushing.
