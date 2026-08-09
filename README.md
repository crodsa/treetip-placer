# TreeTip Placer

TreeTip Placer answers a practical screening question: given one or more
*Clostridioides difficile* 7-locus MLST sequences, do they fit the known
classical or cryptic diversity, and which reference clade is their closest
defensible placement?

It is for placing sequences against a curated, fixed reference collection. It
does not build a new phylogenetic tree and it is not a replacement for
core-genome or whole-genome analysis when fine-scale resolution is required.

## How it works

Give the tool a raw FASTA containing one concatenated 7-locus MLST sequence per
isolate (`adk`, `atpA`, `dxr`, `glyA`, `recA`, `sodA`, `tpi`). TreeTip Placer:

1. maps each sequence to the bundled classical and cryptic reference panels;
2. checks sequence coverage and comparability before calling anything;
3. selects a panel only when the distance difference is clear enough; and
4. places eligible sequences on that panel's fixed tree and assigns a clade
   when its evidence is sufficient.

The result is a self-contained `report.html`, plus TSV and JSON files for
spreadsheets or downstream workflows.

## Install once

TreeTip Placer currently runs from this repository checkout and needs Python
3.11+ plus one Conda-compatible environment manager: `micromamba`, `mamba`, or
`conda`.

```bash
git clone https://github.com/crodsa/treetip-placer.git
cd treetip-placer
python3 -m venv .venv
source .venv/bin/activate
pip install .
treetip setup
```

`treetip setup` automatically finds a supported manager, creates a dedicated
TreeTip environment, installs the pinned MAFFT version, verifies it, and saves
the configuration. It never modifies your current Conda environment or `PATH`.
Run `treetip setup --refresh` only when you intend to recreate that dedicated
environment.

## Run your samples

```bash
treetip classify my_isolates.fasta --organism cdifficile-mlst
```

That is the routine command. It finds the bundled reference package, the MAFFT
installed by setup, and a timestamped output directory automatically. Use
`--output results_name` only when you want to choose the result-directory name.

Your FASTA must contain the 7-locus MLST concatenation, not an unprocessed
whole-genome assembly. Use short, unique isolate IDs in the FASTA headers.

## Read the report

Open `report.html` in the output directory first.

- `CLASSICAL` or `CRYPTIC`: the isolate is distinctly closer to one panel.
- `AMBIGUOUS`: both panels are too similar; the tool deliberately does not
  guess.
- `UNRESOLVED`: input quality or overlap was insufficient for a defensible
  call.
- `PASS`: no placement concern was detected.
- `WARN`: a call was made, but the report records a reason to inspect it.

For each placed isolate, the report shows two complementary trees:

- a **clade overview**, which shows the placement in the full known reference
  universe using collapsed clades rather than reference-tip names; and
- a **placement branch**, which zooms into the local branch, marks the query
  attachment, and names only nearby reference tips.

The report's Findings column explains every `WARN`, `AMBIGUOUS`, or
`UNRESOLVED` result. Treat such calls as review targets, not as failures to be
hidden.

## Other uses

The same fixed-tree workflow can support another organism, marker scheme, or
core-genome/SNP reference collection when a curated, versioned reference
package is available. Use `place-auto` for an explicitly supplied paired
reference package, or `place` when inputs are already in the immutable
reference alignment coordinates. These are advanced workflows; package
construction and calibration are documented separately.

## Technical and maintenance details

- [User manual](docs/MANUAL.md): requirements, troubleshooting, and advanced
  options.
- [Reference package](reference-packages/mlst-cdifficile-20260806/README.md):
  current panel scope and its separate redistribution restriction.
- [Data contract](docs/data-contract.md) and [architecture](docs/architecture.md):
  method and file-format details.
- [Maintainer handoff](docs/HANDOFF.md), [contributing](CONTRIBUTING.md), and
  [security policy](SECURITY.md): maintenance and reporting procedures.

TreeTip Placer source code and documentation are licensed under
[Apache-2.0](LICENSE). The bundled reference package has a separate licence
and is not cleared for public redistribution; see its
[package licence](reference-packages/mlst-cdifficile-20260806/LICENSE) and the
[public-release gates](docs/public-release.md).
