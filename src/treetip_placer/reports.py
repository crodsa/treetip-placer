"""Dependency-free human-readable reports for automatic panel placement."""

from __future__ import annotations

from html import escape
from pathlib import Path
from typing import Mapping, Sequence

from .panel_selection import PanelDecision
from .tree import Node, Tree


def write_auto_report(
    output_directory: str | Path,
    decisions: Sequence[PanelDecision],
    metadata: Mapping[str, object],
) -> None:
    """Write a self-contained HTML summary without replacing tabular outputs."""

    root = Path(output_directory)
    assignments: dict[str, int] = {}
    for decision in decisions:
        assignments[decision.assignment] = assignments.get(decision.assignment, 0) + 1
    summary_items = "".join(
        f"<li><strong>{escape(label)}</strong>: {count}</li>"
        for label, count in sorted(assignments.items())
    )
    rows = "\n".join(render_decision_row(decision) for decision in decisions)
    package = metadata.get("reference_package_path") or "Explicit command-line reference paths"
    threshold = format(float(metadata["minimum_distance_separation"]), ".12g")
    document = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>TreeTip Placer automatic panel report</title>
  <style>
    body {{ font-family: system-ui, sans-serif; margin: 2rem; color: #1f2937; }}
    h1, h2 {{ color: #111827; }}
    table {{ border-collapse: collapse; width: 100%; font-size: 0.9rem; }}
    th, td {{ border: 1px solid #d1d5db; padding: 0.45rem; text-align: left; vertical-align: top; }}
    th {{ background: #f3f4f6; }}
    .CLASSICAL, .CRYPTIC {{ color: #065f46; font-weight: 700; }}
    .AMBIGUOUS, .UNRESOLVED {{ color: #92400e; font-weight: 700; }}
    code {{ word-break: break-all; }}
    figure {{ margin: 1.25rem 0; }}
    svg.dendrogram {{ display: block; width: 100%; height: auto; margin-top: 0.4rem; border: 1px solid #d1d5db; background: #fff; }}
  </style>
</head>
<body>
  <h1>Automatic classical/cryptic panel placement report</h1>
  <p>This report is a human-readable index. TSV, JSON, Newick, Nexus, and JPlace sidecars remain authoritative.</p>
  <h2>Run summary</h2>
  <ul>{summary_items}</ul>
  <p>Minimum distance separation: <strong>{escape(threshold)}</strong></p>
  <p>Reference package: <code>{escape(str(package))}</code></p>
  <h2>Selected-tree dendrograms</h2>
  <p>Only panels assigned to at least one query are shown. Each assigned query has a full-tree context view and a zoom around its placement branch. These diagrams visualize immutable selected trees. Orange labels mark newly placed queries; reference-tip labels are omitted to remain readable.</p>
  {render_selected_tree_figures(root, decisions)}
  <h2>Queries</h2>
  <table>
    <thead><tr><th>Sample</th><th>Assignment</th><th>Decision</th><th>Classical nearest distance</th><th>Cryptic nearest distance</th><th>Margin</th><th>Selected placement</th><th>Findings</th></tr></thead>
    <tbody>{rows}</tbody>
  </table>
  <h2>Artifacts</h2>
  <ul>
    <li><a href="panel_assignment.tsv">Panel assignments (TSV)</a></li>
    <li><a href="panel_assignment.json">Panel assignments (JSON)</a></li>
    <li><a href="selected_trees/classical.selected_placed_tree.nwk">Selected classical placed tree (Newick)</a></li>
    <li><a href="selected_trees/cryptic.selected_placed_tree.nwk">Selected cryptic placed tree (Newick)</a></li>
    <li><a href="panels/classical/placements/placements.tsv">Classical placement details (TSV)</a></li>
    <li><a href="panels/cryptic/placements/placements.tsv">Cryptic placement details (TSV)</a></li>
  </ul>
</body>
</html>
"""
    (root / "report.html").write_text(document, encoding="utf-8")


def render_decision_row(decision: PanelDecision) -> str:
    """Render one escaped result row with stable numeric formatting."""

    values = (
        decision.sample_id,
        decision.assignment,
        decision.decision_status,
        format_value(decision.classical_nearest_distance),
        format_value(decision.cryptic_nearest_distance),
        format_value(decision.distance_margin),
        decision.selected_placement_status or "",
        " | ".join(decision.findings),
    )
    cells = "".join(f"<td>{escape(value)}</td>" for value in values)
    return f"<tr class=\"{escape(decision.assignment)}\">{cells}</tr>"


def render_selected_tree_figures(root: Path, decisions: Sequence[PanelDecision]) -> str:
    """Render each selected fixed-backbone tree and a placement zoom per query."""

    panels = ("classical", "cryptic")
    highlighted = {
        panel: {
            decision.sample_id
            for decision in decisions
            if decision.selected_panel == panel
        }
        for panel in panels
    }
    figures = []
    for panel in panels:
        if not highlighted[panel]:
            continue
        tree_path = root / "selected_trees" / f"{panel}.selected_placed_tree.nwk"
        if not tree_path.is_file():
            figures.append(f"<p><strong>{escape(panel.title())}</strong>: selected tree is unavailable.</p>")
            continue
        tree = Tree.from_path(tree_path)
        figures.append(render_tree_figure(panel, tree, highlighted[panel]))
        for sample_id in sorted(highlighted[panel]):
            figures.append(render_tree_zoom(panel, tree, sample_id))
    return "\n".join(figures)


def render_tree_figure(panel: str, tree: Tree, highlighted_tips: set[str]) -> str:
    """Render a compact rooted phylogram with only placed-query tip labels."""

    width = 1100
    left = 28
    right = 230
    top = 28
    row_height = max(3.0, min(9.0, 780.0 / max(len(tree.tips), 1)))
    height = int(top * 2 + max(len(tree.tips) - 1, 1) * row_height + 24)
    distances: dict[int, float] = {}
    y_positions: dict[int, float] = {}

    def assign(node: Node, distance: float) -> None:
        distances[node.node_id] = distance
        for child in node.children:
            assign(child, distance + (child.length or 0.0))

    assign(tree.root, 0.0)
    maximum_distance = max(distances[node.node_id] for node in tree.tips) or 1.0
    for index, tip in enumerate(tree.tips):
        y_positions[tip.node_id] = top + index * row_height

    def assign_internal_y(node: Node) -> float:
        if node.is_leaf:
            return y_positions[node.node_id]
        children_y = [assign_internal_y(child) for child in node.children]
        y_positions[node.node_id] = sum(children_y) / len(children_y)
        return y_positions[node.node_id]

    assign_internal_y(tree.root)

    def x_position(node: Node) -> float:
        return left + distances[node.node_id] / maximum_distance * (width - left - right)

    lines: list[str] = []
    labels: list[str] = []
    for node in tree.nodes.values():
        if not node.children:
            continue
        child_y = [y_positions[child.node_id] for child in node.children]
        lines.append(
            f'<line x1="{x_position(node):.3f}" y1="{min(child_y):.3f}" '
            f'x2="{x_position(node):.3f}" y2="{max(child_y):.3f}" class="branch"/>'
        )
        for child in node.children:
            lines.append(
                f'<line x1="{x_position(node):.3f}" y1="{y_positions[child.node_id]:.3f}" '
                f'x2="{x_position(child):.3f}" y2="{y_positions[child.node_id]:.3f}" class="branch"/>'
            )
    for tip in tree.tips:
        if tip.name not in highlighted_tips:
            continue
        x = x_position(tip)
        y = y_positions[tip.node_id]
        labels.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="3" class="query-tip"/>')
        labels.append(
            f'<text x="{x + 6:.3f}" y="{y + 3:.3f}" class="query-label">{escape(tip.name or "")}</text>'
        )
    return f"""<figure>
  <figcaption><strong>{escape(panel.title())} selected tree</strong> ({len(tree.tips)} tips)</figcaption>
  <svg class="dendrogram" viewBox="0 0 {width} {height}" role="img" aria-label="{escape(panel)} selected-tree dendrogram">
    <style>.branch {{ stroke: #64748b; stroke-width: 0.8; fill: none; }} .query-tip {{ fill: #ea580c; }} .query-label {{ fill: #9a3412; font: 12px system-ui, sans-serif; }}</style>
    {''.join(lines)}
    {''.join(labels)}
  </svg>
</figure>"""


def render_tree_zoom(panel: str, tree: Tree, sample_id: str) -> str:
    """Render a cropped neighborhood around one newly placed tip.

    The SVG retains the original tree coordinates, so it is a genuine zoom of the
    selected placement tree rather than an independently rearranged schematic.
    """

    focus_tip = next((tip for tip in tree.tips if tip.name == sample_id), None)
    if focus_tip is None:
        return (
            f"<p><strong>{escape(sample_id)}</strong>: placement zoom is unavailable "
            "because the selected tree does not contain the query tip.</p>"
        )

    width = 1100
    left = 28
    right = 230
    top = 28
    row_height = max(3.0, min(9.0, 780.0 / max(len(tree.tips), 1)))
    height = int(top * 2 + max(len(tree.tips) - 1, 1) * row_height + 24)
    distances: dict[int, float] = {}
    y_positions: dict[int, float] = {}

    def assign(node: Node, distance: float) -> None:
        distances[node.node_id] = distance
        for child in node.children:
            assign(child, distance + (child.length or 0.0))

    assign(tree.root, 0.0)
    maximum_distance = max(distances[tip.node_id] for tip in tree.tips) or 1.0
    for index, tip in enumerate(tree.tips):
        y_positions[tip.node_id] = top + index * row_height

    def assign_internal_y(node: Node) -> float:
        if node.is_leaf:
            return y_positions[node.node_id]
        children_y = [assign_internal_y(child) for child in node.children]
        y_positions[node.node_id] = sum(children_y) / len(children_y)
        return y_positions[node.node_id]

    assign_internal_y(tree.root)

    def x_position(node: Node) -> float:
        return left + distances[node.node_id] / maximum_distance * (width - left - right)

    lines: list[str] = []
    for node in tree.nodes.values():
        if not node.children:
            continue
        child_y = [y_positions[child.node_id] for child in node.children]
        lines.append(
            f'<line x1="{x_position(node):.3f}" y1="{min(child_y):.3f}" '
            f'x2="{x_position(node):.3f}" y2="{max(child_y):.3f}" class="branch"/>'
        )
        for child in node.children:
            lines.append(
                f'<line x1="{x_position(node):.3f}" y1="{y_positions[child.node_id]:.3f}" '
                f'x2="{x_position(child):.3f}" y2="{y_positions[child.node_id]:.3f}" class="branch"/>'
            )

    focus_x = x_position(focus_tip)
    focus_y = y_positions[focus_tip.node_id]
    crop_width = min(width, 650)
    crop_left = max(0.0, min(focus_x - crop_width + 180, width - crop_width))
    crop_height = min(float(height), 180.0)
    crop_top = max(0.0, min(focus_y - crop_height / 2, height - crop_height))
    label = escape(sample_id)
    return f"""<figure>
  <figcaption><strong>{label} placement zoom</strong> ({escape(panel.title())} panel)</figcaption>
  <svg class="dendrogram" viewBox="{crop_left:.3f} {crop_top:.3f} {crop_width:.3f} {crop_height:.3f}" role="img" aria-label="{label} placement zoom">
    <style>.branch {{ stroke: #64748b; stroke-width: 0.8; fill: none; }} .query-tip {{ fill: #ea580c; }} .query-label {{ fill: #9a3412; font: 12px system-ui, sans-serif; font-weight: 700; }}</style>
    {''.join(lines)}
    <circle cx="{focus_x:.3f}" cy="{focus_y:.3f}" r="4" class="query-tip"/>
    <text x="{focus_x + 7:.3f}" y="{focus_y + 4:.3f}" class="query-label">{label}</text>
  </svg>
</figure>"""


def format_value(value: float | None) -> str:
    """Format optional distance values without locale-dependent output."""

    return "" if value is None else format(value, ".12g")
