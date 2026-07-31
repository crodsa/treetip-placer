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
                raise ValueError(
                    f"Invalid branch length '{number}' at position {self.position}"
                ) from error
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
            for index, node in enumerate(
                node for node in self.nodes.values() if node.parent is not None
            )
        ]
        self.edge_by_child = {edge.child_id: edge for edge in self.edges}
        self.edge_by_nodes = {
            frozenset((edge.parent_id, edge.child_id)): edge for edge in self.edges
        }
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
        return [
            edge
            for edge in self.edges
            if edge.parent_id in visited and edge.child_id in visited
        ]

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

    def tip_span(self, node_id: int) -> tuple[int, int]:
        """Return the contiguous DFS tip-index span below one node."""

        return self._tip_spans[node_id]

    def _validate_lengths(self) -> None:
        for node in self.nodes.values():
            if node.parent is None:
                continue
            if node.length is None:
                raise ValueError(
                    f"Tree edge ending at node '{node.name or node.node_id}' has no branch length"
                )
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

    def _render_node(self, node: Node) -> str:
        label = quote_label(node.name) if node.name else ""
        if node.is_leaf:
            return label
        return "(" + ",".join(self._render_node(child) + ":" + format_length(child.length or 0.0) for child in node.children) + ")" + label


def quote_label(label: str) -> str:
    """Quote a Newick label only when the syntax requires it."""

    if label and all(character not in " \t\r\n():,;[]'" for character in label):
        return label
    return "'" + label.replace("'", "''") + "'"


def format_length(value: float) -> str:
    """Produce compact, locale-independent Newick branch lengths."""

    return format(value, ".12g")
