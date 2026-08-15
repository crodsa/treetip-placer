# External-user README redesign

## Status

Approved in conversation on 2026-08-15. This specification governs a
documentation-only change for the existing TreeTip Placer 0.7.0 implementation.

## Goal

Replace the current compressed landing page with a layered README that a reader
without prior *Clostridioides difficile* or command-line experience can follow
from first principles, while retaining enough method, output, validation, and
provenance detail for advanced users to judge what the software actually does.

The README must be exhaustive by progressive disclosure rather than by mixing
routine use, reference construction, calibration, and maintenance into one
undifferentiated narrative.

## Evidence basis

Every command, default, method description, output path, and interpretation rule
in the README must be checked against the current repository rather than copied
from historical documentation. The implementation and package validation
currently establish the following facts:

- TreeTip Placer is version 0.7.0 and requires Python 3.11 or newer.
- Routine classification uses the bundled `cdifficile-mlst` paired reference
  package and a MAFFT 7.526 environment managed by `treetip setup`.
- The package contains 201 classical references and 68 cryptic references,
  aligned over 3501 columns, with five classical and seven cryptic clade labels.
- The package declares a calibrated minimum cross-panel distance separation of
  0.01.
- Raw query records are mapped independently to both fixed reference alignments
  with MAFFT `--add --keeplength`.
- Panel selection compares normalized Hamming distance to the nearest reference;
  it does not compare WLS objectives across differently sized trees.
- Eligible queries are inserted onto an immutable reference topology by
  categorical-distance weighted least squares.
- Discrete clade labels come from the selected fixed-tree placement and may use
  an explicitly warned nearest-labelled-reference fallback when edge labels are
  mixed, unavailable, or conflict with the nearest labelled reference.
- `report.html` is a human-readable index. TSV and JSON preserve authoritative
  decisions and QC, while manifests and event logs preserve provenance. Newick
  and Nexus are drawing-oriented, and the distance-JPlace document is
  explicitly experimental.

Historical calibration results may be summarized only as internal validation
of this exact package. They must not be presented as external validation,
clinical validation, a probability of correctness, or evidence of performance
on unrepresented diversity.

## Audience and reading model

The primary reader is a laboratory or genomic-surveillance user who has a
FASTA file but may not know what MLST, a clade, a terminal, or phylogenetic
placement means. The secondary reader is a bioinformatician evaluating the
algorithm, outputs, calibration, and reference package.

The README will therefore use three layers:

1. **Orientation and first successful run** — plain-language concepts,
   prerequisites, input example, installation, one routine command, and where
   to find the result.
2. **Interpretation and audit trail** — decision layers, result labels, output
   directory, findings, reproducibility records, and worked interpretation.
3. **Method and advanced operation** — exact algorithm, defaults, calibration
   boundaries, generic reference-package capabilities, advanced commands, and
   links to specialist documents.

A reader must be able to stop after layer 1 and run the routine workflow, or
continue through layers 2 and 3 without encountering a different or
contradictory description of the same process.

## Scope

### Included

- Rewrite `README.md` in English.
- Preserve and validate links to the user manual, data contract, architecture,
  reference-package documentation, maintainer handoff, public-release status,
  contributing guide, security policy, citation file, changelog, and licence.
- Correct a linked document only when its current wording directly contradicts
  the verified README contract or current repository state.
- Add a concise changelog entry for the documentation redesign.
- Use a small GitHub-rendered Mermaid flowchart when it makes the routine data
  flow easier to understand than prose.

### Excluded

- No placement, QC, report-generation, setup, or reference-package code changes.
- No change to command names, flags, defaults, supported inputs, or outputs.
- No change to repository visibility, source licence, reference-data licence,
  version number, release tag, or release status.
- No new biological data, identifiers, calibration outputs, or restricted
  reference artifacts.
- No claim that the software is clinically validated, externally validated,
  maximum-likelihood placement, de novo phylogenetics, MLST calling, or
  whole-genome analysis.

## README structure

The README will use the following order.

1. **Title and one-sentence outcome.** State that the tool places concatenated
   seven-locus *C. difficile* MLST sequences against a fixed classical/cryptic
   reference package and reports conservative panel and clade assignments.
2. **Who this is for.** Give concrete examples of suitable users and distinguish
   screening from fine-scale genomic investigation.
3. **Concepts in plain language.** Define isolate, FASTA, seven-locus MLST
   sequence, reference panel, clade, fixed tree, and placement. Scientific
   background claims will use traceable authoritative links.
4. **What the tool does and does not do.** Present supported routine behavior
   beside explicit non-capabilities, including the fact that the program does
   not extract MLST loci from a genome assembly or accept allele-number tables.
5. **What input is required.** Show a synthetic two-record FASTA and an input
   checklist: nucleotide sequences, one concatenation per isolate, required
   locus order, unique identifiers, and no whole-genome assemblies. State what
   the user must obtain upstream when only assemblies or allele profiles exist.
6. **Installation.** Define a terminal briefly, show commands individually,
   distinguish repository access from software prerequisites, state supported
   operating systems, and explain what `treetip setup` creates and records.
7. **Quick start.** Use the parser-verified routine command:
   `treetip classify my_isolates.fasta --organism cdifficile-mlst`.
   Explain automatic reference, MAFFT, and timestamped-output resolution.
8. **Workflow diagram and numbered logic.** Show query mapping to both panels,
   pre-placement QC, conservative panel selection, fixed-tree placement, clade
   assignment, and output generation. The diagram must display unresolved and
   ambiguous branches rather than implying that every query receives a call.
9. **First result to open.** Explain `report.html`, its two display-only tree
   views, and its non-authoritative status.
10. **Three decision layers.** Keep panel assignment, placement/QC status, and
    clade-assignment status separate. Define `CLASSICAL`, `CRYPTIC`,
    `AMBIGUOUS`, `UNRESOLVED`, `PASS`, `WARN`, and `FAIL` in their actual
    contexts.
11. **Worked interpretation.** Interpret one synthetic `PASS`, one `WARN`, and
    one unresolved or ambiguous outcome without inventing biological data.
12. **Complete output map.** Group the root assignments, mapping evidence,
    selected trees, both panel-specific runs, placement alternatives, QC,
    manifests, checksums, and logs. For each group, state the intended reader,
    format, and authority.
13. **Methods.** Describe MAFFT mapping, pairwise callable-site normalized
    Hamming distance, the 0.01 panel gate, constrained WLS edge search, optional
    bootstrap, discrete clade assignment, and provenance recording. Link to the
    architecture and data contract for schemas and implementation details.
14. **Current bundled reference package.** Give the verified package counts,
    coordinate width, clade-label counts, build provenance at an appropriate
    level, checksum validation, and the separate reference-data licence.
15. **Validation and evidence limits.** Summarize internal leave-one-out results
    with their denominators and exclusions, clearly label them internal, and
    link to the detailed package/handoff record. Explain that internal controls
    do not establish external generalization.
16. **Limitations and escalation.** Explain MLST resolution, reference
    dependence, unrepresented diversity, heuristic fit weights, distance-based
    JPlace, and when core-genome or whole-genome analysis is needed.
17. **Advanced workflows.** Introduce `place-auto`, `place`, package creation and
    validation, and the three calibration commands without embedding every
    maintainer flag in the beginner path.
18. **Documentation index, citation, licence, and availability.** Retain all
    specialist links and state that the repository remains private while the
    bundled reference package is not cleared for public redistribution.

## Input truth contract

The README must not imply that an arbitrary genome FASTA is accepted. Routine
input is a raw, unaligned concatenated nucleotide sequence of the seven MLST
loci in this order: `adk`, `atpA`, `dxr`, `glyA`, `recA`, `sodA`, `tpi`.
`classify` maps those query sequences to each fixed panel. The lower-level
`place` command instead requires a query already represented in the immutable
reference coordinates.

The README will not invent a locus-extraction command. If an upstream extractor
is named, its current interface and citation must be independently verified and
clearly separated from TreeTip Placer itself.

## Interpretation contract

The routine result has three related but non-interchangeable layers:

- **Panel assignment:** `CLASSICAL`, `CRYPTIC`, `AMBIGUOUS`, or `UNRESOLVED`.
- **Panel/placement decision status:** `PASS`, `WARN`, `FAIL`, or
  `UNRESOLVED`, depending on QC and placement evidence.
- **Clade assignment:** a label plus `PASS`, `WARN`, or `UNRESOLVED`.

`WARN` means that a call was retained with a named review reason; it is neither
a probability nor a failure. `AMBIGUOUS` means that the cross-panel distance
margin did not reach the package threshold. `UNRESOLVED` means that the
available input or evidence did not support the relevant decision. `FAIL`
records a blocking per-panel QC failure and must not be disguised as an
ordinary ambiguous comparison.

Nearest-reference distance, distance margin, relative-fit weight, RMSE, and
bootstrap support must be described as different quantities. In particular,
relative-fit weight is not a likelihood weight ratio or posterior probability,
and the experimental JPlace-shaped output is not maximum-likelihood JPlace.

## Output contract

The complete routine output will be documented in functional groups:

- `report.html`: human-readable, self-contained index and display-only tree
  diagrams.
- `panel_assignment.tsv` and `.json`: authoritative cross-panel decisions and
  run metadata.
- `clade_assignment.tsv` and `.json`: authoritative fixed-tree clade calls,
  methods, supporting tips, and findings.
- `mapping/`: mapped-with-reference FASTA files, query-only mapped FASTA files,
  and MAFFT stderr logs for both panels.
- `selected_trees/`: directly drawable Newick trees containing queries selected
  for each panel.
- `panels/classical/` and `panels/cryptic/`: complete independent placement
  runs, even when the other panel wins the final decision.
- Each panel run's `trees/`, `placements/`, and `qc/` directories: placed and
  normalized trees, Nexus, placement summaries and alternatives, JSON,
  experimental distance-JPlace, and query-level QC.
- Each panel run's `run_manifest.json` and `events.jsonl`, plus the root
  `run.log`: paths, hashes, parameters, worker decisions, progress, and
  provenance.

The README will identify TSV/JSON sidecars and manifests as the basis for audit
and downstream analysis. Tree graphics will never be described as changing the
fixed topology or as carrying support values that actually live in sidecars.

## Writing and presentation rules

- Repository text remains English; the prose uses short sentences and defines
  each technical term before using it as shorthand.
- Commands contain no shell prompt character and identify where they should be
  run.
- Examples use only synthetic identifiers and sequences.
- The first routine path contains no calibration or package-maintenance flags.
- Repeated cautions are consolidated into the relevant interpretation or
  limitation section instead of appearing before every action.
- Advanced detail remains visible through a table of contents and stable links,
  rather than being deleted or placed ahead of the quick start.
- Claims about biology, software requirements, algorithms, and validation carry
  either a repository evidence link, an authoritative external source, or an
  explicit qualification.

## Verification

Before the documentation change is committed and synchronized:

1. Compare every documented command and option with `treetip --help` and the
   relevant subcommand help generated by the current parser.
2. Run `validate-reference-package` on the bundled manifest and compare all
   documented package counts and thresholds with its output.
3. Check every local Markdown link and image target and reject missing targets.
4. Check fenced code-block balance and scan for stale paths, commands,
   unsupported promises, and placeholder text.
5. Run the complete repository suite with a supported interpreter:
   `PYTHONPATH=src /Users/crodsa/miniforge3/bin/python -m unittest discover -s tests -v`.
6. Run `git diff --check` and inspect the final README diff as rendered Markdown.
7. Confirm that the worktree contains no generated outputs, real query data,
   credentials, or unrelated edits.

## Publication path

The owner has confirmed the destination `crodsa/treetip-placer`, the present
private visibility, and a direct update to `main`. The private-development gates
in `docs/release-checklist.md` are currently recorded as complete. After final
verification, the documentation commit may be pushed to the private `main`
branch. This does not authorize changing repository visibility, redistributing
the bundled reference package, tagging a release, or describing the project as
publicly released.

## Acceptance criteria

- A reader unfamiliar with *C. difficile* and CLI tools can identify the exact
  accepted input, install the software, run `classify`, find the output, and
  distinguish a call from a warning or unresolved result.
- A technical reader can trace panel selection, placement, clade assignment,
  defaults, package identity, validation evidence, and every output family to
  current implementation or versioned documentation.
- The README is materially more complete than the current 107-line version but
  keeps routine instructions ahead of calibration and maintenance detail.
- No command, capability, statistic, citation, licence statement, or scientific
  claim is invented or overstated.
- All existing specialist documentation links remain available and valid.
- The complete supported-Python test suite and documentation checks pass before
  synchronization with private `main`.
