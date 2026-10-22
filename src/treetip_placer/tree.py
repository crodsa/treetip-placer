"""Minimal Newick handling and fixed-tree distance utilities."""

from __future__ import annotations

import math
from array import array
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


@dataclass
class Node:
    """A rooted Newick node; the underlying metric is treated as unrooted."""

    node_id: int
    name: str | None = None
    length: float | None = None
    parent: "Node | None" = None
    children: list["Node"] = field(default_factory=list)

    @property
    def is_leaf(self) -> bool:
        return not self.children


@dataclass(frozen=True)
class Edge:
    """A parent-to-child Newick edge."""

    edge_id: int
    parent_id: int
    child_id: int
    length: float


@dataclass(frozen=True)
class PrunedTipAttachment:
    """Where a removed reference tip attached to its pruned backbone.

    An attachment can remain an internal node or become a point on an edge when
    pruning suppresses a unary internal node.  The record is intended for
    leave-one-out validation and is not a phylogenetic support value.
    """

    kind: str
    node_id: int | None = None
    edge_parent_id: int | None = None
    edge_child_id: int | None = None
    distal_length: float | None = None


class NewickParser:
    """A deliberately small parser for branch-length Newick files."""

    def __init__(self, text: str) -> None:
        self.text = text.strip()
        self.position = 0
        self.next_node_id = 0

    def parse(self) -> Node:
        if not self.text:
            raise ValueError("Newick input is empty")
        root = self._parse_subtree()
        self._skip_whitespace()
        if self._peek() == ";":
            self.position += 1
        self._skip_whitespace()
        if self.position != len(self.text):
            raise ValueError(f"Unexpected Newick content at position {self.position}")
        root.length = None
        return root

    def _parse_subtree(self) -> Node:
        self._skip_whitespace()
        node = self._new_node()
        if self._peek() == "(":
            self.position += 1
            while True:
                child = self._parse_subtree()
                child.parent = node
                node.children.append(child)
                self._skip_whitespace()
                token = self._peek()
                if token == ",":
                    self.position += 1
                    continue
                if token == ")":
                    self.position += 1
                    break
                raise ValueError(f"Expected ',' or ')' at position {self.position}")
            node.name = self._parse_label(optional=True)
        else:
            node.name = self._parse_label(optional=False)
        self._skip_whitespace()
        if self._peek() == ":":
            self.position += 1
            number = self._read_until(",();[] \t\r\n")
            try:
                node.length = float(number)
            except ValueError as error:
                raise ValueError(f"Invalid branch length '{number}' at position {self.position}") from error
        return node

    def _parse_label(self, optional: bool) -> str | None:
        self._skip_whitespace()
        if self._peek() == "'":
            self.position += 1
            chunks: list[str] = []
            while self.position < len(self.text):
                char = self.text[self.position]
                self.position += 1
                if char == "'":
                    if self._peek() == "'":
                        chunks.append("'")
                        self.position += 1
                        continue
                    return "".join(chunks)
                chunks.append(char)
            raise ValueError("Unterminated quoted Newick label")
        label = self._read_until(":,();[] \t\r\n")
        if not label and not optional:
            raise ValueError(f"Expected a tip label at position {self.position}")
        return label or None

    def _read_until(self, stop_characters: str) -> str:
        start = self.position
        while self.position < len(self.text) and self.text[self.position] not in stop_characters:
            self.position += 1
        return self.text[start : self.position]

    def _new_node(self) -> Node:
        node = Node(node_id=self.next_node_id)
        self.next_node_id += 1
        return node

    def _skip_whitespace(self) -> None:
        while self.position < len(self.text) and self.text[self.position].isspace():
            self.position += 1

    def _peek(self) -> str | None:
        return self.text[self.position] if self.position < len(self.text) else None


class Tree:
    """A branch-length tree with cached node-to-tip distance rows."""

    def __init__(self, root: Node, max_cache_megabytes: int = 256) -> None:
        self.root = root
        self.nodes = {node.node_id: node for node in self._walk(root)}
        self.tips = [node for node in self.nodes.values() if node.is_leaf]
        self.tip_index: dict[str, int] = {}
        for index, tip in enumerate(self.tips):
            if not tip.name:
                raise ValueError("Every Newick tip must have a label")
            if tip.name in self.tip_index:
                raise ValueError(f"Duplicate Newick tip label: {tip.name}")
            self.tip_index[tip.name] = index
        self._tip_spans: dict[int, tuple[int, int]] = {}
        self._record_tip_spans(self.root)
        self._validate_lengths()
        self.edges = [
            Edge(edge_id=index, parent_id=node.parent.node_id, child_id=node.node_id, length=node.length or 0.0)
            for index, node in enumerate(node for node in self.nodes.values() if node.parent is not None)
        ]
        self.edge_by_child = {edge.child_id: edge for edge in self.edges}
        self.edge_by_nodes = {frozenset((edge.parent_id, edge.child_id)): edge for edge in self.edges}
        self.adjacency: dict[int, list[tuple[int, float]]] = {node_id: [] for node_id in self.nodes}
        for edge in self.edges:
            self.adjacency[edge.parent_id].append((edge.child_id, edge.length))
            self.adjacency[edge.child_id].append((edge.parent_id, edge.length))
        self._distance_cache: dict[int, array] = {}
        self._cache_limit_bytes = max_cache_megabytes * 1024 * 1024
        self._cached_bytes = 0

    @classmethod
    def from_path(cls, path: str | Path, max_cache_megabytes: int = 256) -> "Tree":
        return cls(NewickParser(Path(path).read_text(encoding="utf-8")).parse(), max_cache_megabytes)

    @classmethod
    def from_newick(cls, text: str, max_cache_megabytes: int = 256) -> "Tree":
        return cls(NewickParser(text).parse(), max_cache_megabytes)

    def newick(self) -> str:
        """Render the reference tree with normalized Newick formatting."""

        return self._render_node(self.root) + ";\n"

    def tip_distances(self, node_id: int) -> array:
        """Return distances from one node to tips in stable tip-index order."""

        cached = self._distance_cache.get(node_id)
        if cached is not None:
            return cached
        values = array("d", [0.0]) * len(self.tips)
        stack: list[tuple[int, int | None, float]] = [(node_id, None, 0.0)]
        while stack:
            current_id, previous_id, distance = stack.pop()
            current = self.nodes[current_id]
            if current.is_leaf:
                values[self.tip_index[current.name or ""]] = distance
            for neighbor_id, branch_length in self.adjacency[current_id]:
                if neighbor_id != previous_id:
                    stack.append((neighbor_id, current_id, distance + branch_length))
        required_bytes = len(values) * values.itemsize
        if required_bytes <= self._cache_limit_bytes - self._cached_bytes:
            self._distance_cache[node_id] = values
            self._cached_bytes += required_bytes
        return values

    def nearby_edges(self, tip_names: Iterable[str], hops: int) -> list[Edge]:
        """Return all edges within a graph-radius around selected reference tips."""

        visited: dict[int, int] = {}
        queue: list[int] = []
        for tip_name in tip_names:
            tip_id = self.tips[self.tip_index[tip_name]].node_id
            if tip_id not in visited:
                visited[tip_id] = 0
                queue.append(tip_id)
        cursor = 0
        while cursor < len(queue):
            node_id = queue[cursor]
            cursor += 1
            depth = visited[node_id]
            if depth >= hops:
                continue
            for neighbor_id, _ in self.adjacency[node_id]:
                if neighbor_id not in visited:
                    visited[neighbor_id] = depth + 1
                    queue.append(neighbor_id)
        return [edge for edge in self.edges if edge.parent_id in visited and edge.child_id in visited]

    def edge_subtree_tip_indices(self, edge: Edge) -> set[int]:
        """Return tip indexes on the child side of an oriented Newick edge."""

        indices: set[int] = set()
        stack = [edge.child_id]
        while stack:
            node_id = stack.pop()
            node = self.nodes[node_id]
            if node.is_leaf:
                indices.add(self.tip_index[node.name or ""])
            else:
                stack.extend(child.node_id for child in node.children)
        return indices

    def without_tip(self, tip_name: str) -> tuple["Tree", PrunedTipAttachment]:
        """Return a copied tree with one tip removed and unary nodes suppressed.

        This operation produces the backbone required for an honest
        leave-one-out placement: the held-out profile and its matching tree tip
        are both absent before the placement is scored.
        """

        if tip_name not in self.tip_index:
            raise ValueError(f"Cannot prune unknown tip '{tip_name}'")
        if len(self.tips) <= 3:
            raise ValueError("Leave-one-out validation requires at least four reference tips")

        def clone_subtree(node: Node) -> Node:
            clone = Node(node_id=node.node_id, name=node.name, length=node.length)
            for child in node.children:
                cloned_child = clone_subtree(child)
                cloned_child.parent = clone
                clone.children.append(cloned_child)
            return clone

        root = clone_subtree(self.root)
        cloned_nodes = {node.node_id: node for node in self._walk(root)}
        removed_tip = cloned_nodes[self.tips[self.tip_index[tip_name]].node_id]
        parent = removed_tip.parent
        if parent is None:
            raise ValueError("Cannot prune a tree root tip")
        parent.children.remove(removed_tip)
        if not parent.children:
            raise ValueError("Pruning produced an empty internal node")

        if len(parent.children) >= 2:
            attachment = PrunedTipAttachment(kind="node", node_id=parent.node_id)
        else:
            survivor = parent.children[0]
            if parent.parent is None:
                survivor.parent = None
                survivor.length = None
                root = survivor
                attachment = PrunedTipAttachment(kind="node", node_id=survivor.node_id)
            else:
                grandparent = parent.parent
                parent_length = parent.length
                survivor_length = survivor.length
                if parent_length is None or survivor_length is None:
                    raise ValueError("Cannot suppress a node without branch lengths")
                survivor.parent = grandparent
                survivor.length = parent_length + survivor_length
                child_index = grandparent.children.index(parent)
                grandparent.children[child_index] = survivor
                attachment = PrunedTipAttachment(
                    kind="edge",
                    edge_parent_id=grandparent.node_id,
                    edge_child_id=survivor.node_id,
                    distal_length=parent_length,
                )

        return Tree(root, self._cache_limit_bytes // (1024 * 1024)), attachment

    def edge_distance_to_attachment(self, edge: Edge, attachment: PrunedTipAttachment) -> int:
        """Return a graph-edge distance from a selected edge to a true attachment.

        Zero means the selected edge touches a node attachment or exactly matches
        an attachment that lies within an edge.  It is a topology diagnostic,
        not a branch-length error or a support statistic.
        """

        selected_nodes = (edge.parent_id, edge.child_id)
        if attachment.kind == "node":
            if attachment.node_id is None:
                raise ValueError("Node attachment is missing a node ID")
            return min(self._node_distance(attachment.node_id, selected_node) for selected_node in selected_nodes)
        if attachment.kind == "edge":
            if attachment.edge_parent_id is None or attachment.edge_child_id is None:
                raise ValueError("Edge attachment is missing edge node IDs")
            expected_nodes = (attachment.edge_parent_id, attachment.edge_child_id)
            if frozenset(selected_nodes) == frozenset(expected_nodes):
                return 0
            return 1 + min(
                self._node_distance(expected_node, selected_node)
                for expected_node in expected_nodes
                for selected_node in selected_nodes
            )
        raise ValueError(f"Unknown attachment kind '{attachment.kind}'")

    def patristic_distance_to_attachment(
        self, edge: Edge, distal_length: float, attachment: PrunedTipAttachment
    ) -> float:
        """Return the branch-length-weighted distance from a placed point to a true attachment.

        A point on ``edge`` at ``distal_length`` from its parent, and the
        attachment point recorded by ``without_tip``, can each only be reached
        from the rest of the tree through the two ends of their own edge (or,
        for a node attachment, through the node itself). Distance is therefore
        the minimum over those endpoint combinations; the tree's single-path
        property makes every other combination strictly longer, so the
        minimum always selects the one true path. This complements
        ``edge_distance_to_attachment``'s topology-only hop count with a
        continuous measure: an edge many hops away along a chain of
        near-zero-length branches is not scored as far as one of equal hop
        count on long branches. It remains a topology diagnostic, not a
        phylogenetic support value.

        When both points sit on the same edge, routing through an endpoint
        would retrace that shared segment and double-count it, exactly the
        degenerate case ``edge_distance_to_attachment`` also special-cases;
        the direct distal-length difference is used instead.
        """

        if (
            attachment.kind == "edge"
            and attachment.edge_parent_id == edge.parent_id
            and attachment.edge_child_id == edge.child_id
        ):
            return abs(distal_length - (attachment.distal_length or 0.0))
        placed_from_parent = distal_length + self._distance_to_attachment(edge.parent_id, attachment)
        placed_from_child = (edge.length - distal_length) + self._distance_to_attachment(edge.child_id, attachment)
        return min(placed_from_parent, placed_from_child)

    def _distance_to_attachment(self, node_id: int, attachment: PrunedTipAttachment) -> float:
        """Return the weighted distance from one node to a true attachment."""

        if attachment.kind == "node":
            if attachment.node_id is None:
                raise ValueError("Node attachment is missing a node ID")
            return self._weighted_node_distance(node_id, attachment.node_id)
        if attachment.kind == "edge":
            if attachment.edge_parent_id is None or attachment.edge_child_id is None:
                raise ValueError("Edge attachment is missing edge node IDs")
            if attachment.distal_length is None:
                raise ValueError("Edge attachment is missing a distal length")
            attachment_edge = self.edge_by_nodes[frozenset((attachment.edge_parent_id, attachment.edge_child_id))]
            via_parent = attachment.distal_length + self._weighted_node_distance(node_id, attachment.edge_parent_id)
            via_child = (attachment_edge.length - attachment.distal_length) + self._weighted_node_distance(
                node_id, attachment.edge_child_id
            )
            return min(via_parent, via_child)
        raise ValueError(f"Unknown attachment kind '{attachment.kind}'")

    def tip_span(self, node_id: int) -> tuple[int, int]:
        """Return the contiguous DFS tip-index span below one node."""

        return self._tip_spans[node_id]

    def _validate_lengths(self) -> None:
        for node in self.nodes.values():
            if node.parent is None:
                continue
            if node.length is None:
                raise ValueError(f"Tree edge ending at node '{node.name or node.node_id}' has no branch length")
            if not math.isfinite(node.length) or node.length < 0:
                raise ValueError(
                    f"Tree edge ending at node '{node.name or node.node_id}' has an invalid branch length"
                )

    def _walk(self, root: Node) -> Iterable[Node]:
        stack = [root]
        while stack:
            node = stack.pop()
            yield node
            stack.extend(reversed(node.children))

    def _record_tip_spans(self, node: Node) -> tuple[int, int]:
        if node.is_leaf:
            index = self.tip_index[node.name or ""]
            span = (index, index + 1)
        else:
            child_spans = [self._record_tip_spans(child) for child in node.children]
            span = (child_spans[0][0], child_spans[-1][1])
        self._tip_spans[node.node_id] = span
        return span

    def _node_distance(self, start_id: int, end_id: int) -> int:
        """Return unweighted graph distance between two nodes in this tree."""

        if start_id == end_id:
            return 0
        if start_id not in self.adjacency or end_id not in self.adjacency:
            raise ValueError("Attachment node is not present in the pruned tree")
        queue = [start_id]
        distances = {start_id: 0}
        for node_id in queue:
            next_distance = distances[node_id] + 1
            for neighbor_id, _ in self.adjacency[node_id]:
                if neighbor_id == end_id:
                    return next_distance
                if neighbor_id not in distances:
                    distances[neighbor_id] = next_distance
                    queue.append(neighbor_id)
        raise ValueError("Tree nodes are disconnected")

    def _weighted_node_distance(self, start_id: int, end_id: int) -> float:
        """Return branch-length-weighted distance between two nodes in this tree.

        A tree has exactly one path between any two nodes, so a plain BFS
        accumulating branch length (rather than Dijkstra's priority queue)
        already finds the unique shortest path, mirroring ``_node_distance``.
        """

        if start_id == end_id:
            return 0.0
        if start_id not in self.adjacency or end_id not in self.adjacency:
            raise ValueError("Attachment node is not present in the pruned tree")
        queue = [start_id]
        distances = {start_id: 0.0}
        for node_id in queue:
            for neighbor_id, branch_length in self.adjacency[node_id]:
                next_distance = distances[node_id] + branch_length
                if neighbor_id == end_id:
                    return next_distance
                if neighbor_id not in distances:
                    distances[neighbor_id] = next_distance
                    queue.append(neighbor_id)
        raise ValueError("Tree nodes are disconnected")

    def _render_node(self, node: Node) -> str:
        label = quote_label(node.name) if node.name else ""
        if node.is_leaf:
            return label
        return (
            "("
            + ",".join(self._render_node(child) + ":" + format_length(child.length or 0.0) for child in node.children)
            + ")"
            + label
        )


def quote_label(label: str) -> str:
    """Quote a Newick label only when the syntax requires it."""

    if label and all(character not in " \t\r\n():,;[]'" for character in label):
        return label
    return "'" + label.replace("'", "''") + "'"


def format_length(value: float) -> str:
    """Produce compact, locale-independent Newick branch lengths."""

    return format(value, ".12g")
