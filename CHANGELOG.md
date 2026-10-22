# Changelog

All notable changes are documented here.

## 1.0 — 2026-10-22

Initial public release of TreeTip Placer (`v1.0`).

- Dual-panel coordinate projection (`MAFFT --add --keeplength`), calibrated `0.01` nearest-reference p-distance separation gate, and constrained non-negative weighted least-squares (`WLS`) fixed-tree placement for concatenated seven-locus *Clostridioides difficile* MLST sequences (`3,501 bp`).
- Routine one-command classification workflow (`treetip classify --organism cdifficile-mlst`) and explicit paired-reference workflow (`treetip place-auto`) with managed MAFFT environment setup (`treetip setup`).
- Descendant-tip consensus and nearest-labelled-reference clade arbitration resolving classical (`Clade-1`–`Clade-5`) and cryptic (`C-I`–`C-VII`) lineages with explicit `PASS`, `WARN`, `AMBIGUOUS`, and `UNRESOLVED` quality flags.
- Self-contained interactive HTML reports with per-query global clade-context and local placement-branch SVG dendrograms, authoritative TSV/JSON sidecars, Newick/Nexus trees, and JPlace-shaped distance-fit outputs.
- Built-in leave-one-out calibration workflows (`calibrate-panel-gate`, `calibrate-clade-assignment`, and `calibrate-placement-accuracy`) and checksummed reference-package validation (`mlst-cdifficile-20260806`).
