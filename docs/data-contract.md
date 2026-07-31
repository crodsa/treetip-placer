# Data contract

TreeTip Placer accepts only fixed-column FASTA alignments. Reference and query records must have identical widths and use stable, unique identifiers. Reference IDs must exactly match branch-length Newick tree tips.

Input sequences must already be mapped to the immutable reference coordinate system. Gaps, `N`, `.`, `?`, and `X` are treated as non-callable. The tool rejects unequal widths instead of changing an alignment.

`treetip place-auto` is the exception at the command boundary, not in the placement engine: it accepts a raw query FASTA, maps it independently to the supplied classical and cryptic reference alignments through MAFFT `--add --keeplength`, and preserves both mapped query alignments under `mapping/`. Raw query records may be incomplete or have unequal raw lengths; the mapped alignments are then subjected to the same fixed-column QC gate. The reference alignments themselves are never changed.

A paired reference package is a JSON document with `schema_version: 1`, `classical` and `cryptic` panel objects, and `minimum_distance_separation`. Each panel declares `tree` and `alignment` paths, optionally with lowercase SHA-256 checksums. Relative paths are portable relative to the manifest location. When checksums are present, `place-auto` rejects modified package artifacts before analysis.

The QC gate reports coverage, comparable sites, ambiguous IUPAC bases, nearest-reference distance, and placement ambiguity. Thresholds are reference-panel specific; do not apply cryptic-only calibrations to a mixed cryptic–classical dataset.

Outputs include drawing-oriented Newick and Nexus trees, authoritative TSV placement/QC tables, JSON, JSONL event logs, and an experimental distance-based JPlace-shaped document. The JPlace objective is not a likelihood. Automatic runs also contain `panel_assignment.tsv`/JSON, two complete panel-specific output directories, MAFFT mapping logs, and selected-panel Newick trees.

`calibrate-panel-gate` selects deterministic unique-sequence representatives from each labelled panel, removes every exact duplicate of an own-panel control before calculating its nearest distance, and maps the controls to the opposing panel. It writes calibration records, a threshold sweep, and `panel_label_conflicts.tsv` for high-confidence calls that disagree with the supplied control label. It does not alter source labels or panel membership, does not perform topology placement during calibration, and must not be described as an external validation.
