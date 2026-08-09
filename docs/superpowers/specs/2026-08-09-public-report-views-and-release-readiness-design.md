# Public report views and release readiness design

## Status

Approved by the repository owner on 2026-08-09. This specification covers the
two per-query report trees and public-release scaffolding. It does not change
the fixed-backbone placement algorithm, QC semantics, reference topology, or
authoritative TSV/JSON/Newick artifacts.

## Report views

`place-auto` renders exactly two SVG figures for every query with a selected
panel.

### Global clade context

The global view answers where the query sits within the complete known
reference universe. It retains branch paths and relative branch lengths from
the selected immutable tree. No reference-tip identifier may be emitted in
this SVG. Every maximal reference-only subtree is represented as one display
leaf labelled as one of:

- `<clade> (<n> tips)` when all annotated descendants have one clade label;
- `Mixed (<n> tips)` when annotated descendants have multiple labels; or
- `Unlabelled (<n> tips)` when it has no annotated descendants.

The query remains an orange named leaf. This rule applies to the target clade
as well as non-target clades, so the overview never expands reference labels.
When clade metadata is unavailable or the query clade assignment is unresolved,
the overview still renders using `Unlabelled` or `Mixed` blocks and declares
the assignment status in its caption. The display must not manufacture a clade
label from the placement tree.

### Local placement branch

The local view answers where the query attaches relative to named neighbouring
reference tips. It is a topology-derived display subtree, not an SVG coordinate
crop. Starting at the query parent, the renderer ascends through ancestors until
the smallest subtree with at least two and at most `max_reference_tips` visible
reference tips is found. If no such ancestor exists, it uses the closest
ancestor with at least two references and contracts excess reference-only
subtrees deterministically as `… (<n> tips)` blocks.

Visible reference tips retain their names. The query is orange and the branch
from its parent to the query is blue and heavier than ordinary branches. The
caption names the selected panel and the clade-assignment status. Default
`max_reference_tips` is 12 and is supplied internally by the report renderer;
it is not a public CLI flag in this release because it affects only HTML
readability and needs no user-facing persistence.

## Tests and documentation

Tests must prove that the global SVG contains no reference identifiers, that it
persists without metadata, and that a non-monophyletic clade is represented in
more than one block when topology requires it. Local-view tests must demonstrate
that a vertical-neighbour distractor outside the selected display subtree is not
shown, while a genuine sister reference is shown. The synthetic HTML screenshot,
README, and data contract must describe only these two views.

## Public release boundaries

The Python source is relicensed under Apache-2.0. The current C. difficile
reference package is explicitly excluded from that code licence and remains a
release blocker until source-by-source redistribution rights are documented.
This change does not make the repository public and does not claim that the
reference package is redistributable. Release documentation must state the
blocker and require either documented permission or package removal/replacement
with redistributable synthetic data before a visibility change.

The release scaffold includes citation, contribution, security, changelog,
issue/PR templates, Dependabot, and CI checks for formatting, tests, package
build, and a deterministic CLI smoke placement. Version tags and empirical
benchmark records remain manual release gates.

## Acceptance criteria

1. The complete supported Python test suite passes after the renderer change.
2. Every selected query in `report.html` has one global and one local figure.
3. The global SVG contains the query label but no reference-tip label.
4. The local SVG includes the query, its highlighted attachment branch, and
   topology-selected named references.
5. Public-release documentation rejects publication while reference-data
   redistribution rights are unresolved.
6. No real reference sequence is relicensed or represented as public-ready by
   this change.
