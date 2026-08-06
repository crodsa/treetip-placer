"""Conservative selection between independently calibrated reference panels."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping, Sequence

from .models import PlacementResult


@dataclass(frozen=True)
class PanelDecision:
    """A cross-panel classification decision for one query sequence."""

    sample_id: str
    assignment: str
    decision_status: str
    selected_panel: str | None
    selected_placement_status: str | None
    classical_nearest_reference: str | None
    classical_nearest_distance: float | None
    cryptic_nearest_reference: str | None
    cryptic_nearest_distance: float | None
    distance_margin: float | None
    minimum_distance_separation: float
    findings: tuple[str, ...]


def choose_reference_panel(
    classical_results: Sequence[PlacementResult],
    cryptic_results: Sequence[PlacementResult],
    minimum_distance_separation: float,
) -> list[PanelDecision]:
    """Choose a panel only when its nearest-reference distance is clearly better.

    Panel selection deliberately uses the directly comparable normalized distance to
    the nearest reference, not the least-squares objective or relative fit weight.
    Those quantities depend on a panel's tree size and topology. A query failing a
    panel's pre-placement QC cannot win selection even when ``--force`` was used to
    obtain a diagnostic placement.
    """

    classical_by_id = {result.sample_id: result for result in classical_results}
    cryptic_by_id = {result.sample_id: result for result in cryptic_results}
    if set(classical_by_id) != set(cryptic_by_id):
        raise ValueError("Panel placement results do not contain the same query identifiers")

    decisions: list[PanelDecision] = []
    for sample_id in sorted(classical_by_id):
        classical = classical_by_id[sample_id]
        cryptic = cryptic_by_id[sample_id]
        candidates = [
            ("classical", classical),
            ("cryptic", cryptic),
        ]
        eligible = [
            (panel, result)
            for panel, result in candidates
            if result.qc.status != "FAIL" and result.qc.nearest_distance is not None
        ]
        common = dict(
            sample_id=sample_id,
            classical_nearest_reference=classical.qc.nearest_reference,
            classical_nearest_distance=classical.qc.nearest_distance,
            cryptic_nearest_reference=cryptic.qc.nearest_reference,
            cryptic_nearest_distance=cryptic.qc.nearest_distance,
            minimum_distance_separation=minimum_distance_separation,
        )
        if not eligible:
            decisions.append(
                PanelDecision(
                    assignment="UNRESOLVED",
                    decision_status="FAIL",
                    selected_panel=None,
                    selected_placement_status=None,
                    distance_margin=None,
                    findings=("NO_PANEL_PASSED_QC",),
                    **common,
                )
            )
            continue
        if len(eligible) == 1:
            panel, result = eligible[0]
            decisions.append(
                PanelDecision(
                    assignment=panel.upper(),
                    decision_status="WARN",
                    selected_panel=panel,
                    selected_placement_status=result.final_status,
                    distance_margin=None,
                    findings=("ONLY_ONE_PANEL_PASSED_QC",),
                    **common,
                )
            )
            continue

        ordered = sorted(
            eligible,
            key=lambda item: (
                item[1].qc.nearest_distance
                if item[1].qc.nearest_distance is not None
                else float("inf")
            ),
        )
        chosen_panel, chosen = ordered[0]
        _, runner_up = ordered[1]
        assert chosen.qc.nearest_distance is not None
        assert runner_up.qc.nearest_distance is not None
        margin = runner_up.qc.nearest_distance - chosen.qc.nearest_distance
        if margin < minimum_distance_separation:
            decisions.append(
                PanelDecision(
                    assignment="AMBIGUOUS",
                    decision_status="UNRESOLVED",
                    selected_panel=None,
                    selected_placement_status=None,
                    distance_margin=margin,
                    findings=("PANEL_DISTANCE_MARGIN_TOO_SMALL",),
                    **common,
                )
            )
            continue
        # A WARN here means something about the chosen placement itself was
        # flagged (diffuse fit weight, a long terminal branch, ...); carry its
        # own QC finding codes onto the decision so the report can say why
        # instead of showing an unexplained WARN.
        decisions.append(
            PanelDecision(
                assignment=chosen_panel.upper(),
                decision_status="PASS" if chosen.final_status == "PASS" else "WARN",
                selected_panel=chosen_panel,
                selected_placement_status=chosen.final_status,
                distance_margin=margin,
                findings=tuple(finding.code for finding in chosen.qc.findings),
                **common,
            )
        )
    return decisions


def write_panel_selection_outputs(
    output_directory: str | Path,
    decisions: Sequence[PanelDecision],
    metadata: Mapping[str, object],
) -> None:
    """Write stable TSV and JSON cross-panel decision sidecars."""

    root = Path(output_directory)
    fields = [
        "sample_id",
        "assignment",
        "decision_status",
        "selected_panel",
        "selected_placement_status",
        "classical_nearest_reference",
        "classical_nearest_distance",
        "cryptic_nearest_reference",
        "cryptic_nearest_distance",
        "distance_margin",
        "minimum_distance_separation",
        "findings",
    ]
    with (root / "panel_assignment.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for decision in decisions:
            row = asdict(decision)
            row["findings"] = " | ".join(decision.findings)
            for field in (
                "classical_nearest_distance",
                "cryptic_nearest_distance",
                "distance_margin",
                "minimum_distance_separation",
            ):
                value = row[field]
                row[field] = "" if value is None else format(value, ".12g")
            writer.writerow(row)
    (root / "panel_assignment.json").write_text(
        json.dumps({"metadata": dict(metadata), "decisions": [asdict(item) for item in decisions]}, indent=2, sort_keys=True)
        + "\n",
        encoding="utf-8",
    )
