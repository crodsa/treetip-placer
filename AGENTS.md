# TreeTip Placer project guidance

## Language

Write all source code, identifiers, comments, docstrings, CLI output, logs, schemas, tests, documentation, Git history, and GitHub materials in English. User-facing discussion may follow the user's language.

## Scientific scope

The repository develops a quality-controlled fixed-backbone placement tool. Do not describe a placement result as equivalent to a full de novo phylogenetic inference. Preserve the reference coordinate system, reference topology, masking, model metadata, and provenance.

The current backend uses categorical-distance weighted least squares for MLST profiles and fixed-column alignments. Keep likelihood placement, SNP ascertainment handling, and conditional local topology refinement clearly separated as future work until implemented and validated.

## Input and output contract

- Treat MLST allele IDs as categorical values, never as numeric distances.
- Do not silently drop invalid or incomplete queries; report `PASS`, `WARN`, `FAIL`, or `UNRESOLVED` with explicit reasons.
- Keep drawing-oriented Newick free of overloaded confidence labels. Store uncertainty and QC in sidecar tables and JSON.
- Do not add real genomic or MLST data to version control unless redistribution permission is explicit.

## Verification

Run the complete test suite before handoff or release:

```powershell
python -m unittest discover -s tests -v
```

Update `docs/HANDOFF.md` whenever a milestone, validation result, dependency, input contract, or known limitation changes.

## Release policy

Do not publish to GitHub until the release gates in `docs/release-checklist.md` are met and the user has confirmed repository destination, visibility, and license.
