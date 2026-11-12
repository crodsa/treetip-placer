# TreeTip Placer — User Manual

This manual is for lab staff running TreeTip Placer on their own isolates. See `README.md` and `docs/architecture.md` for mathematical and algorithmic details.

## What this tool does

You give it one FASTA file with the sequences of your isolates. It tells you, for each one:

- Whether it looks like a **classical** or a **cryptic** *Clostridioides difficile* isolate.
- Which known **clade** it is closest to within that group.
- How confident that call is, in plain terms (`PASS`, `WARN`, `AMBIGUOUS`, or `UNRESOLVED` — explained below).

It does this by comparing your sequence to a fixed reference collection of already-classified isolates. It does not build a new evolutionary tree from your data, and it does not replace expert review — it is a triage and screening aid.

## Before you start: what you need

1. **The tool installed** (see below).
2. **A Conda-compatible environment manager**: `micromamba`, `mamba`, or `conda`. The one-time `treetip setup` command installs and verifies the MAFFT version TreeTip Placer uses internally. You do not need to install MAFFT separately.
3. **Your own sequences as a FASTA file.** Each sequence should be the concatenation of the same 7 MLST loci (adk, atpA, dxr, glyA, recA, sodA, tpi) the bundled reference was built from, one entry per isolate, with a short unique name after the `>` on each header line. TreeTip validates FASTA structure and mapped-sequence QC, but it cannot infer whether an upstream workflow extracted the correct loci or joined them in the correct order. Confirm that preparation before running.

You do **not** need to separately locate a reference package — one ships with this repository and `classify` finds it on its own; see below.

## Installing

You need Python 3.11 or newer. From a terminal, inside the project folder:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
treetip setup
```

`treetip setup` creates a dedicated TreeTip environment, installs the pinned
MAFFT version through the first available manager (`micromamba`, then `mamba`,
then `conda`), verifies it, and records its exact path. It does not alter your
current environment or `PATH`. If this dedicated environment needs to be
recreated, use `treetip setup --refresh` deliberately.
Managed setup supports macOS and Linux; on Windows, run TreeTip Placer in a
Linux environment such as WSL or provide a verified MAFFT path explicitly.

Confirm the command-line installation:

```bash
treetip --version
```

This should print `1.0.0`. If it doesn't, something went wrong with the steps above — ask for help before continuing.

## Running it on your samples

The command you will use almost every time is `classify`. Give it your FASTA file and say which organism/reference package you want — it figures out the rest (finds MAFFT, names the output folder, and decides classical vs. cryptic per isolate) without any other setup:

```bash
treetip classify my_isolates.fasta --organism cdifficile-mlst
```

That's the whole command for routine use. It prints the three things it resolved automatically before it starts:

```
Using organism 'cdifficile-mlst' -> reference-packages/mlst-cdifficile-20260806/paired_reference.json
Using MAFFT 7.526 (managed_environment): /path/to/treetip-placer/environments/mafft-7.526/bin/mafft
Writing output to: my_isolates_classify_20260806T142530Z
```

`--organism` is the one thing you always state explicitly — the tool will never guess which reference collection you mean, even when only one is bundled today. Run `treetip classify --help` to see every currently available `--organism` value.

Override anything it would have picked automatically if you need to: `--output my_results` to name the folder yourself, `--mafft /path/to/mafft` if you have more than one copy installed, or `--min-panel-distance-separation` to override the package's own calibrated threshold. Everything not given a value falls back to the automatic choice above.

Runtime depends on the number of isolates, hardware, search preset, worker
count, and whether bootstrap replicates were requested; use the completion event
and output manifest to confirm that a run finished.

### If you need more control: `place-auto`

`place-auto` is the same underlying command with nothing automatic — every input is an explicit flag, including the reference package path, in case you're pointing at a reference collection other than the one bundled with this repository:

```bash
treetip place-auto \
  --reference-package reference-packages/mlst-cdifficile-20260806/paired_reference.json \
  --query-fasta my_isolates.fasta \
  --mafft /path/to/mafft \
  --output my_results
```

Use this if `classify`'s auto-detection picks something you don't want and no override flag covers it, or if you maintain your own reference package outside this repository.

## Reading the results

Open **`report.html`** inside your output folder in any web browser. This is the file to look at first — it summarizes the run in one place:

- A summary count of how many isolates were called `CLASSICAL`, `CRYPTIC`, `AMBIGUOUS`, or `UNRESOLVED`.
- A table with one row per isolate: its call, its clade, and the distances that call was based on.
- A tree picture for each isolate that was placed, showing where it landed relative to the known reference collection.

The underlying data behind that report also lives in the same folder as plain tables, in case you need to load them into a spreadsheet or another program:

- `panel_assignment.tsv` — the classical/cryptic call, per isolate.
- `clade_assignment.tsv` — the clade call, per isolate.

The report is a human-readable index. The TSV and JSON sidecars retain the
authoritative decisions and QC. See the
[complete output map](../README.md#-complete-output-map) for mappings,
panel-specific placements, alternatives, trees, manifests, and logs.

### What the labels mean

**Panel call** (classical vs. cryptic):

| Label | Meaning |
| --- | --- |
| `CLASSICAL` / `CRYPTIC` | Clearly closer to one reference collection than the other. |
| `AMBIGUOUS` | Close to both — the tool declines to guess. Treat as needing manual review. |
| `UNRESOLVED` | The sequence itself had a quality problem (too much missing data, too little overlap with the reference, etc.) before a call could even be attempted. |

**Decision and clade status**, shown per isolate in the tables:

- `PASS` means no rule changed that layer to a warning or no-call. Still inspect
  the Findings column, especially when bootstrap was requested.
- `WARN` means the call was still made, but something about it deserves a second look — for example, nearby branches carry different clade labels or candidate edge fits were diffuse. `WARN` is not a failure or a probability; read its named findings.
- `UNRESOLVED` means the available evidence did not support that decision
  layer. A per-panel `FAIL` records blocking QC; the other panel can still be
  considered if it remains eligible.

If you see a lot of `WARN` or `AMBIGUOUS` results on isolates you have independent reason to trust, that is worth flagging to whoever maintains the reference package — it may mean the reference collection needs to grow to cover that part of the diversity better.

### Why a call is `WARN`: reading the Findings column

Every `report.html` row and every line of `panel_assignment.tsv` / `clade_assignment.tsv` has a **Findings** column. Whenever a status is `WARN`, `AMBIGUOUS`, or `UNRESOLVED`, this column names the specific reason — it is never left blank. One isolate can show more than one reason at once, separated by `|`.

| Code | In plain terms |
| --- | --- |
| `SELECTED_EDGE_CONSENSUS_DISAGREES_WITH_NEAREST_REFERENCE` | The selected edge's descendant labels and the nearest labelled reference disagree. The configured rule retains the nearest-reference label as `WARN`; review the conflict rather than treating it as a clean consensus. |
| `MIXED_CLADE_LABELS_ON_SELECTED_EDGE` | The isolate landed in a part of the reference tree where neighboring reference isolates don't all agree on their clade label, so the closest single match was used instead of a broader consensus. |
| `DIFFUSE_DISTANCE_WEIGHT` | The heuristic weight was spread across alternative tree edges — the selected edge was still the best available fit, but it was not sharply separated. This weight is not a probability. |
| `LONG_TERMINAL_BRANCH` | The fitted terminal branch is more than three times the nearest-reference distance. Review the input, nearest match, and placement alternatives. |
| `UNSTABLE_BOOTSTRAP_PLACEMENT` | (Only appears if bootstrap replicates were requested.) Repeated resampling of the alignment did not consistently support the same placement. |
| `PANEL_DISTANCE_MARGIN_TOO_SMALL` | The isolate is nearly equally close to both the classical and cryptic collections — this is what produces an `AMBIGUOUS` call. |
| `PARTIAL_COVERAGE` / `LOW_COVERAGE` | Part of the isolate's sequence didn't align to usable reference positions. `LOW_COVERAGE` is severe enough to block a call (`UNRESOLVED`); `PARTIAL_COVERAGE` is not. |
| `DISTANT_FROM_REFERENCE_PANEL` | Even the closest reference isolate is unusually far away — treat the call with extra caution. |
| `NO_COMPARABLE_REFERENCE` / `INSUFFICIENT_COMPARABLE_SITES` | Too little of the isolate's sequence overlaps with the reference collection to compare at all. |

### A note on precision

The current reference collection is built from 7 MLST genes, not full genomes.
Internal leave-one-out controls support panel discrimination and clade-label
recovery for this exact package; they are not an external performance estimate.
Seven loci are **not** enough to resolve fine differences between two isolates
that are genetically very close to each other — do not over-interpret small
distance differences as if they had whole-genome resolution. Use an
appropriately validated core-genome or whole-genome workflow when finer
resolution is required.

## If something looks wrong

- **The command reports duplicate IDs, reference-ID overlap, or invalid FASTA structure**: make each first header token unique, keep sequence content after a header, and do not reuse a bundled reference-tip ID.
- **Every one of your isolates comes back `UNRESOLVED`**: inspect the Findings column for low coverage or insufficient comparable sites, then confirm that the input contains the correct seven nucleotide loci in the required order and is not a genome assembly.
- **A result surprises you**: open `report.html` and look at the tree picture and the distance numbers for that isolate before assuming the tool is wrong — the numbers are there specifically so you can check the reasoning, not just the label.

For anything else, contact whoever maintains the reference package for your lab.

## Administrator details

`treetip setup` stores its configuration and managed environment in standard
user-scoped platform directories. Administrators and test systems can redirect
them with `TREETIP_CONFIG_DIR` and `TREETIP_DATA_DIR`; this makes managed setup
independent of shared Conda environments.
