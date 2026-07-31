# Data contract

TreeTip Placer accepts only fixed-column FASTA alignments. Reference and query records must have identical widths and use stable, unique identifiers. Reference IDs must exactly match branch-length Newick tree tips.

Input sequences must already be mapped to the immutable reference coordinate system. Gaps, `N`, `.`, `?`, and `X` are treated as non-callable. The tool rejects unequal widths instead of changing an alignment.

`treetip place-auto` is the exception at the command boundary, not in the placement engine: it accepts a raw query FASTA, maps it independently to the supplied classical and cryptic reference alignments through MAFFT `--add --keeplength`, and preserves both mapped query alignments under `mapping/`. Raw query records may be incomplete or have unequal raw lengths; the mapped alignments are then subjected to the same fixed-column QC gate. The reference alignments themselves are never changed.

The QC gate reports coverage, comparable sites, ambiguous IUPAC bases, nearest-reference distance, and placement ambiguity. Thresholds are reference-panel specific; do not apply cryptic-only calibrations to a mixed cryptic–classical dataset.

Outputs include drawing-oriented Newick and Nexus trees, authoritative TSV placement/QC tables, JSON, JSONL event logs, and an experimental distance-based JPlace-shaped document. The JPlace objective is not a likelihood. Automatic runs also contain `panel_assignment.tsv`/JSON, two complete panel-specific output directories, MAFFT mapping logs, and selected-panel Newick trees.
