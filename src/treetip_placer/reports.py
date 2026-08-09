"""Dependency-free human-readable reports for automatic panel placement."""

from __future__ import annotations

from dataclasses import dataclass, field
from html import escape
from pathlib import Path
from typing import Mapping, Sequence

from .clade_assignment import CladeAssignment
from .panel_selection import PanelDecision
from .tree import Node, Tree


@dataclass
class ContextNode:
    """A display-only tree node used by the clade-context SVG."""

    source: Node
    children: list["ContextNode"] = field(default_factory=list)
    label: str | None = None
    is_query: bool = False
    collapsed: bool = False


def write_auto_report(
    output_directory: str | Path,
    decisions: Sequence[PanelDecision],
    metadata: Mapping[str, object],
    clade_assignments: Sequence[CladeAssignment] = (),
    tip_clades: Mapping[str, Mapping[str, str]] = {},
) -> None:
    """Write a self-contained HTML summary without replacing tabular outputs."""

    root = Path(output_directory)
    assignments: dict[str, int] = {}
    for decision in decisions:
        assignments[decision.assignment] = assignments.get(decision.assignment, 0) + 1
    summary_items = "".join(
        f"<li><strong>{escape(label)}</strong>: {count}</li>" for label, count in sorted(assignments.items())
    )
    clade_by_sample = {assignment.sample_id: assignment for assignment in clade_assignments}
    rows = "\n".join(render_decision_row(decision, clade_by_sample.get(decision.sample_id)) for decision in decisions)
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
  <p>Only panels assigned to at least one query are shown. Each assigned query has exactly two display-only diagrams of its immutable selected tree: a global clade context without reference-tip names and a topology-derived local placement branch with named neighbouring references. Orange marks newly placed queries; the blue branch in the local view is the query attachment branch. Neither diagram changes the fixed tree or clade call.</p>
  {render_selected_tree_figures(root, decisions, clade_by_sample, tip_clades)}
  <h2>Queries</h2>
  <table>
    <thead><tr><th>Sample</th><th>Assignment</th><th>Phylogenetic clade</th><th>Clade status</th><th>Decision</th><th>Classical nearest distance</th><th>Cryptic nearest distance</th><th>Margin</th><th>Selected placement</th><th>Findings</th></tr></thead>
    <tbody>{rows}</tbody>
  </table>
  <h2>Artifacts</h2>
  <ul>
    <li><a href="panel_assignment.tsv">Panel assignments (TSV)</a></li>
    <li><a href="panel_assignment.json">Panel assignments (JSON)</a></li>
    <li><a href="clade_assignment.tsv">Fixed-tree clade assignments (TSV)</a></li>
    <li><a href="clade_assignment.json">Fixed-tree clade assignments (JSON)</a></li>
    <li><a href="selected_trees/classical.selected_placed_tree.nwk">Selected classical placed tree (Newick)</a></li>
    <li><a href="selected_trees/cryptic.selected_placed_tree.nwk">Selected cryptic placed tree (Newick)</a></li>
    <li><a href="panels/classical/placements/placements.tsv">Classical placement details (TSV)</a></li>
    <li><a href="panels/cryptic/placements/placements.tsv">Cryptic placement details (TSV)</a></li>
  </ul>
</body>
</html>
"""
    (root / "report.html").write_text(document, encoding="utf-8")


def render_decision_row(
    decision: PanelDecision,
    clade_assignment: CladeAssignment | None = None,
) -> str:
    """Render one escaped result row with stable numeric formatting.

    Findings combines both sources that can independently raise a WARN --
    the panel decision (the selected placement's own QC findings) and the
    clade assignment -- so a reader sees every reason, not just whichever
    layer happened to have one.
    """

    all_findings = decision.findings + (clade_assignment.findings if clade_assignment else ())
    values = (
        decision.sample_id,
        decision.assignment,
        clade_assignment.phylogenetic_clade if clade_assignment and clade_assignment.phylogenetic_clade else "",
        clade_assignment.clade_assignment_status if clade_assignment else "",
        decision.decision_status,
        format_value(decision.classical_nearest_distance),
        format_value(decision.cryptic_nearest_distance),
        format_value(decision.distance_margin),
        decision.selected_placement_status or "",
        " | ".join(all_findings),
    )
    cells = "".join(f"<td>{escape(value)}</td>" for value in values)
    return f'<tr class="{escape(decision.assignment)}">{cells}</tr>'


def render_selected_tree_figures(
    root: Path,
    decisions: Sequence[PanelDecision],
    clade_by_sample: Mapping[str, CladeAssignment],
    tip_clades: Mapping[str, Mapping[str, str]],
) -> str:
    """Render global clade context and local placement branch per selected query."""

    panels = ("classical", "cryptic")
    highlighted = {
        panel: {decision.sample_id for decision in decisions if decision.selected_panel == panel} for panel in panels
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
        for sample_id in sorted(highlighted[panel]):
            assignment = clade_by_sample.get(sample_id)
            figures.append(
                render_clade_context(
                    panel,
                    tree,
                    sample_id,
                    assignment.phylogenetic_clade if assignment else None,
                    tip_clades.get(panel, {}),
                )
            )
            figures.append(render_tree_zoom(panel, tree, sample_id, assignment))
    return "\n".join(figures)


def render_tree_figure(
    panel: str,
    tree: Tree,
    highlighted_tips: set[str],
    *,
    show_tip_names: bool = True,
) -> str:
    """Render a compact rooted phylogram with readable reference-tip labels."""

    width = 1100
    left = 28
    right = 330
    top = 28
    row_height = max(7.0, min(12.0, 1500.0 / max(len(tree.tips), 1)))
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
        x = x_position(tip)
        y = y_positions[tip.node_id]
        if tip.name in highlighted_tips:
            labels.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="3" class="query-tip"/>')
            labels.append(f'<text x="{x + 6:.3f}" y="{y + 3:.3f}" class="query-label">{escape(tip.name or "")}</text>')
        elif show_tip_names:
            labels.append(
                f'<text x="{x + 4:.3f}" y="{y + 2.8:.3f}" class="reference-label">{escape(tip.name or "")}</text>'
            )
    return f"""<figure>
  <figcaption><strong>{escape(panel.title())} selected tree</strong> ({len(tree.tips)} tips; all tip names shown)</figcaption>
  <svg class="dendrogram" viewBox="0 0 {width} {height}" role="img" aria-label="{escape(panel)} selected-tree dendrogram">
    <style>.branch {{ stroke: #64748b; stroke-width: 0.8; fill: none; }} .query-tip {{ fill: #ea580c; }} .query-label {{ fill: #9a3412; font: 12px system-ui, sans-serif; font-weight: 700; }} .reference-label {{ fill: #334155; font: 8px system-ui, sans-serif; }}</style>
    {"".join(lines)}
    {"".join(labels)}
  </svg>
</figure>"""


def render_tree_zoom(
    panel: str,
    tree: Tree,
    sample_id: str,
    clade_assignment: CladeAssignment | None = None,
    *,
    max_reference_tips: int = 12,
) -> str:
    """Render a topology-derived local subtree around one newly placed tip.

    This deliberately avoids an SVG-coordinate crop: display neighbours are
    selected by graph-hop proximity to the placement branch and all excluded
    reference-only subtrees remain explicit collapsed leaves.  It is therefore
    a faithful local view of the immutable selected tree, not a re-layout or
    local reconstruction.
    """

    focus_tip = next((tip for tip in tree.tips if tip.name == sample_id), None)
    if focus_tip is None:
        return (
            f"<p><strong>{escape(sample_id)}</strong>: placement branch is unavailable "
            "because the selected tree does not contain the query tip.</p>"
        )

    if focus_tip.parent is None:
        return f"<p><strong>{escape(sample_id)}</strong>: placement branch is unavailable because the query is a tree root.</p>"
    if max_reference_tips < 2:
        raise ValueError("max_reference_tips must be at least 2")

    distances: dict[int, float] = {}

    def assign(node: Node, distance: float) -> None:
        distances[node.node_id] = distance
        for child in node.children:
            assign(child, distance + (child.length or 0.0))

    assign(tree.root, 0.0)

    def descendant_reference_names(node: Node) -> tuple[str, ...]:
        if node.is_leaf:
            return () if node.name == sample_id else (node.name or "",)
        return tuple(name for child in node.children for name in descendant_reference_names(child))

    def contains_query(node: Node) -> bool:
        return node.name == sample_id if node.is_leaf else any(contains_query(child) for child in node.children)

    ancestors: list[Node] = []
    current = focus_tip.parent
    while current is not None:
        ancestors.append(current)
        current = current.parent
    display_root = next(
        (node for node in ancestors if 2 <= len(descendant_reference_names(node)) <= max_reference_tips),
        next(
            (node for node in ancestors if len(descendant_reference_names(node)) >= 2),
            ancestors[-1],
        ),
    )

    allowed_node_ids: set[int] = set()

    def collect_node_ids(node: Node) -> None:
        allowed_node_ids.add(node.node_id)
        for child in node.children:
            collect_node_ids(child)

    collect_node_ids(display_root)

    # Graph-hop ordering is topology-based and deterministic because the tree's
    # node identifiers preserve Newick traversal order for tied tips. Restrict
    # the walk to the chosen display subtree so an external cousin cannot use up
    # the local-label budget.
    visible_reference_names: set[str] = set()
    queue: list[tuple[int, int]] = [(focus_tip.parent.node_id, 0)]
    seen = {focus_tip.parent.node_id}
    cursor = 0
    while cursor < len(queue) and len(visible_reference_names) < max_reference_tips:
        node_id, hops = queue[cursor]
        cursor += 1
        node = tree.nodes[node_id]
        if node.is_leaf and node.name != sample_id:
            visible_reference_names.add(node.name or "")
            continue
        for neighbour_id, _ in sorted(tree.adjacency[node_id]):
            if neighbour_id in allowed_node_ids and neighbour_id not in seen:
                seen.add(neighbour_id)
                queue.append((neighbour_id, hops + 1))

    def build(node: Node) -> ContextNode:
        references = descendant_reference_names(node)
        has_query = contains_query(node)
        if not has_query and not any(name in visible_reference_names for name in references):
            return ContextNode(source=node, label=f"… ({len(references)} tips)", collapsed=True)
        if node.is_leaf:
            if node.name == sample_id:
                return ContextNode(source=node, label=node.name, is_query=True)
            return ContextNode(source=node, label=node.name)
        return ContextNode(source=node, children=[build(child) for child in node.children])

    local_root = build(display_root)
    leaves: list[ContextNode] = []

    def collect_leaves(node: ContextNode) -> None:
        if not node.children:
            leaves.append(node)
            return
        for child in node.children:
            collect_leaves(child)

    collect_leaves(local_root)
    width = 1100
    left = 28
    right = 330
    top = 28
    row_height = max(18.0, min(28.0, 620.0 / max(len(leaves), 1)))
    height = int(top * 2 + max(len(leaves) - 1, 1) * row_height + 28)
    y_positions = {id(leaf): top + index * row_height for index, leaf in enumerate(leaves)}

    def assign_y(node: ContextNode) -> float:
        if not node.children:
            return y_positions[id(node)]
        values = [assign_y(child) for child in node.children]
        y_positions[id(node)] = sum(values) / len(values)
        return y_positions[id(node)]

    assign_y(local_root)
    root_distance = distances[display_root.node_id]
    visible_distances = [distances[leaf.source.node_id] for leaf in leaves]
    distance_span = max(visible_distances, default=root_distance) - root_distance or 1.0

    def x_position(node: ContextNode) -> float:
        return left + (distances[node.source.node_id] - root_distance) / distance_span * (width - left - right)

    lines: list[str] = []
    labels: list[str] = []

    def render(node: ContextNode) -> None:
        x = x_position(node)
        if node.children:
            child_y = [y_positions[id(child)] for child in node.children]
            lines.append(
                f'<line x1="{x:.3f}" y1="{min(child_y):.3f}" x2="{x:.3f}" y2="{max(child_y):.3f}" class="branch"/>'
            )
            for child in node.children:
                child_x = x_position(child)
                branch_class = "placement-branch" if child.is_query else "branch"
                lines.append(
                    f'<line x1="{x:.3f}" y1="{y_positions[id(child)]:.3f}" '
                    f'x2="{child_x:.3f}" y2="{y_positions[id(child)]:.3f}" class="{branch_class}"/>'
                )
                render(child)
            return
        y = y_positions[id(node)]
        label = escape(node.label or "")
        if node.is_query:
            labels.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="4" class="query-tip"/>')
            labels.append(f'<text x="{x + 7:.3f}" y="{y + 4:.3f}" class="query-label">{label}</text>')
        elif node.collapsed:
            labels.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="3" class="collapsed-tip"/>')
            labels.append(f'<text x="{x + 6:.3f}" y="{y + 3.5:.3f}" class="collapsed-label">{label}</text>')
        else:
            labels.append(f'<text x="{x + 4:.3f}" y="{y + 3:.3f}" class="reference-label">{label}</text>')

    render(local_root)
    label = escape(sample_id)
    clade = escape(clade_assignment.phylogenetic_clade or "unresolved") if clade_assignment else "unavailable"
    status = escape(clade_assignment.clade_assignment_status) if clade_assignment else "unavailable"
    return f"""<figure>
  <figcaption><strong>{label} placement branch</strong> ({escape(panel.title())} panel; clade {clade}; status {status}; up to {max_reference_tips} named reference tips)</figcaption>
  <svg class="dendrogram" viewBox="0 0 {width} {height}" role="img" aria-label="{label} topology-derived placement branch">
    <style>.branch {{ stroke: #64748b; stroke-width: 0.9; fill: none; }} .placement-branch {{ stroke: #2563eb; stroke-width: 2.8; fill: none; }} .query-tip {{ fill: #ea580c; }} .query-label {{ fill: #9a3412; font: 13px system-ui, sans-serif; font-weight: 700; }} .reference-label {{ fill: #334155; font: 11px system-ui, sans-serif; }} .collapsed-tip {{ fill: #64748b; }} .collapsed-label {{ fill: #475569; font: 11px system-ui, sans-serif; }}</style>
    {"".join(lines)}
    {"".join(labels)}
  </svg>
</figure>"""


def render_clade_context(
    panel: str,
    tree: Tree,
    sample_id: str,
    target_clade: str | None,
    tip_clades: Mapping[str, str],
) -> str:
    """Render the full known reference universe without reference-tip names.

    Every reference-only subtree becomes a single clade, mixed, or unlabelled
    display leaf.  This gives the query global context while keeping all
    reference identifiers out of the zoom-out figure.
    """

    if sample_id not in tree.tip_index:
        return (
            f"<p><strong>{escape(sample_id)}</strong>: global clade context is unavailable "
            "because the selected tree does not contain the query tip.</p>"
        )

    def descendant_tip_names(node: Node) -> tuple[str, ...]:
        if node.is_leaf:
            return (node.name or "",)
        return tuple(name for child in node.children for name in descendant_tip_names(child))

    def build(node: Node) -> ContextNode:
        names = descendant_tip_names(node)
        contains_query = sample_id in names
        reference_names = tuple(name for name in names if name != sample_id)
        clades = {tip_clades[name] for name in reference_names if name in tip_clades}
        if not contains_query:
            if not clades:
                label = f"Unlabelled ({len(reference_names)} tips)"
            elif len(clades) == 1 and len(clades) == len({tip_clades.get(name) for name in reference_names}):
                label = f"{next(iter(clades))} ({len(reference_names)} tips)"
            elif len(clades) == 1:
                label = f"{next(iter(clades))} + unlabelled ({len(reference_names)} tips)"
            else:
                label = f"Mixed ({len(reference_names)} tips)"
            return ContextNode(source=node, label=label, collapsed=True)
        if node.is_leaf:
            name = node.name or ""
            if name == sample_id:
                return ContextNode(source=node, label=name, is_query=True)
            return ContextNode(source=node, label="Unlabelled (1 tips)", collapsed=True)
        return ContextNode(source=node, children=[build(child) for child in node.children])

    display_root = build(tree.root)
    distances: dict[int, float] = {}

    def assign_distances(node: Node, distance: float) -> None:
        distances[node.node_id] = distance
        for child in node.children:
            assign_distances(child, distance + (child.length or 0.0))

    assign_distances(tree.root, 0.0)
    maximum_distance = max(distances[tip.node_id] for tip in tree.tips) or 1.0
    leaves: list[ContextNode] = []

    def collect_leaves(node: ContextNode) -> None:
        if not node.children:
            leaves.append(node)
            return
        for child in node.children:
            collect_leaves(child)

    collect_leaves(display_root)
    width = 1100
    left = 28
    right = 300
    top = 28
    row_height = max(18.0, min(28.0, 900.0 / max(len(leaves), 1)))
    height = int(top * 2 + max(len(leaves) - 1, 1) * row_height + 28)
    y_positions = {id(leaf): top + index * row_height for index, leaf in enumerate(leaves)}

    def assign_y(node: ContextNode) -> float:
        if not node.children:
            return y_positions[id(node)]
        values = [assign_y(child) for child in node.children]
        y_positions[id(node)] = sum(values) / len(values)
        return y_positions[id(node)]

    assign_y(display_root)

    def x_position(node: ContextNode) -> float:
        return left + distances[node.source.node_id] / maximum_distance * (width - left - right)

    lines: list[str] = []
    labels: list[str] = []

    def render(node: ContextNode) -> None:
        x = x_position(node)
        y = y_positions[id(node)]
        if node.children:
            child_y = [y_positions[id(child)] for child in node.children]
            lines.append(
                f'<line x1="{x:.3f}" y1="{min(child_y):.3f}" x2="{x:.3f}" y2="{max(child_y):.3f}" class="branch"/>'
            )
            for child in node.children:
                child_x = x_position(child)
                lines.append(
                    f'<line x1="{x:.3f}" y1="{y_positions[id(child)]:.3f}" '
                    f'x2="{child_x:.3f}" y2="{y_positions[id(child)]:.3f}" class="branch"/>'
                )
                render(child)
            return
        label = escape(node.label or "")
        if node.is_query:
            labels.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="4" class="query-tip"/>')
            labels.append(f'<text x="{x + 7:.3f}" y="{y + 4:.3f}" class="query-label">{label}</text>')
        elif node.collapsed:
            labels.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="3" class="collapsed-tip"/>')
            labels.append(f'<text x="{x + 6:.3f}" y="{y + 3.5:.3f}" class="collapsed-label">{label}</text>')
        else:
            labels.append(f'<text x="{x + 4:.3f}" y="{y + 3:.3f}" class="context-tip-label">{label}</text>')

    render(display_root)
    label = escape(sample_id)
    clade = escape(target_clade or "unresolved")
    return f"""<figure>
  <figcaption><strong>{label} global clade context</strong> (clade {clade}; all reference-only subtrees collapsed; reference-tip names omitted)</figcaption>
  <svg class="dendrogram" viewBox="0 0 {width} {height}" role="img" aria-label="{label} global clade context">
    <style>.branch {{ stroke: #64748b; stroke-width: 0.9; fill: none; }} .query-tip {{ fill: #ea580c; }} .query-label {{ fill: #9a3412; font: 13px system-ui, sans-serif; font-weight: 700; }} .collapsed-tip {{ fill: #2563eb; }} .collapsed-label {{ fill: #1e3a8a; font: 12px system-ui, sans-serif; font-weight: 650; }} .context-tip-label {{ fill: #334155; font: 11px system-ui, sans-serif; }}</style>
    {"".join(lines)}
    {"".join(labels)}
  </svg>
</figure>"""


def format_value(value: float | None) -> str:
    """Format optional distance values without locale-dependent output."""

    return "" if value is None else format(value, ".12g")
