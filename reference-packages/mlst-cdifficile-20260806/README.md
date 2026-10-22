# mlst-cdifficile-20260806

A paired classical/cryptic *Clostridioides difficile* reference package for `treetip place-auto`, built from concatenated 7-locus MLST alleles (adk, atpA, dxr, glyA, recA, sodA, tpi). Committed to this repository with the repository owner's explicit permission (see `AGENTS.md` and `docs/HANDOFF.md` "Current reference panel").

## Use it

```bash
treetip place-auto \
  --reference-package reference-packages/mlst-cdifficile-20260806/paired_reference.json \
  --query-fasta your_isolates.fasta \
  --mafft /path/to/mafft \
  --output your_results
```

See `docs/MANUAL.md` for a plain-language guide to running this and reading the results.

## What's in it

- `classical/` — 201 STs across all 5 classical MLST clades, from a PubMLST/BIGSdb export of 1319 STs.
- `cryptic/` — 68 isolates across all 7 cryptic clades (C-I–C-VII), MLST-typed with `fastMLST` from local genome assemblies.
- Both panels: `mafft --auto` alignment (3501 columns), `IQ-TREE 3 -m MFP -B 1000 --alrt 1000 --bnni` tree (support labels stripped from the tree TreeTip Placer reads, per `AGENTS.md`'s rule against overloaded confidence labels).
- `paired_reference.json` — the manifest, `minimum_distance_separation: 0.01`.

## Calibration

- **Panel gate** (`calibrate-panel-gate`): 100/100 leave-one-out controls correct at every threshold from 0.0 to 0.01. The `0.01` margin in this package is a calibrated value for this exact pair, not an inherited default.
- **Clade assignment** (`calibrate-clade-assignment`): 200/201 classical (99.5%), 66/66 cryptic (100%, excluding the two singleton clades C-VI/C-VII, which cannot be leave-one-out evaluated with one reference each). The one remaining classical case is a boundary isolate independently confirmed by FastBAPS clustering, correctly surfaced as `WARN` rather than a silent pass.
- **Placement accuracy** (`calibrate-placement-accuracy`): edge-hop grading looks poor on both panels (a known artifact of MLST-resolution trees with many near-zero-length branches — see `docs/HANDOFF.md`), but the branch-length-weighted (patristic) distance stays small throughout (classical mean 0.005 / max 0.010 on a 0.141-length tree; cryptic mean 0.031 / max 0.037 on a 0.324-length tree), and `nearest_reference_distance` remains small (classical mean 0.0006, max 0.0031, <=0.002 in 199/201 controls; cryptic mean 0.0024, median 0.0009, max 0.0237, <=0.002 in 52/68 controls). Read the patristic summary, not the hop-based outcome column, for this package.

Full provenance, exact numbers, and the FastBAPS cross-validation that motivated the clade-assignment fix are in `docs/HANDOFF.md`. A copy of every calibration run's full output (`report.html` included) is kept outside this repository by the owner.

## Limitations

Seven MLST loci separate classical from cryptic isolates and place them within the right clade reliably. They do not resolve fine genetic differences between isolates that are close to each other — this is not a substitute for core-genome or whole-genome analysis. See `docs/MANUAL.md` for the plain-language version of this caveat.
