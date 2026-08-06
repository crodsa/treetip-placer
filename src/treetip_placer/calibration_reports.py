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
    "classical_controls",
    "cryptic_controls",
    "flagged_panel_label_conflicts",
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
  <title>{escape(title)}</title>
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
