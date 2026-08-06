# TreeTip Placer — User Manual

This manual is for lab staff running TreeTip Placer on their own isolates. It does not explain how the tool is built internally; see `README.md` and `docs/HANDOFF.md` for that.

## What this tool does

You give it one FASTA file with the sequences of your isolates. It tells you, for each one:

- Whether it looks like a **classical** or a **cryptic** *Clostridioides difficile* isolate.
- Which known **clade** it is closest to within that group.
- How confident that call is, in plain terms (`PASS`, `WARN`, `AMBIGUOUS`, or `UNRESOLVED` — explained below).

It does this by comparing your sequence to a fixed reference collection of already-classified isolates. It does not build a new evolutionary tree from your data, and it does not replace expert review — it is a triage and screening aid.

## Before you start: what you need

1. **The tool installed** (see below).
2. **MAFFT** installed and on your system (a free, widely used sequence-alignment program: <https://mafft.cbrc.jp/alignment/software/>). TreeTip Placer uses it internally to line up your sequences against the reference. `classify` (below) finds it automatically if it's on your system `PATH`.
3. **Your own sequences as a FASTA file.** Each sequence should be the concatenation of the same 7 MLST loci (adk, atpA, dxr, glyA, recA, sodA, tpi) the bundled reference was built from, one entry per isolate, with a short unique name after the `>` on each header line. If you are not sure your sequences are in the right format, placement will refuse to run rather than silently guess.

You do **not** need to separately locate a reference package — one ships with this repository and `classify` finds it on its own; see below.

## Installing

You need Python 3.11 or newer. From a terminal, inside the project folder:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

Confirm it installed correctly:

```bash
treetip --version
```

This should print a version number (for example `0.7.0`). If it doesn't, something went wrong with the steps above — ask for help before continuing.

## Running it on your samples

The command you will use almost every time is `classify`. Give it your FASTA file and say which organism/reference package you want — it figures out the rest (finds MAFFT, names the output folder, and decides classical vs. cryptic per isolate) without any other setup:

```bash
treetip classify my_isolates.fasta --organism cdifficile-mlst
```

That's the whole command for routine use. It prints the three things it resolved automatically before it starts:

```
Using organism 'cdifficile-mlst' -> reference-packages/mlst-cdifficile-20260806/paired_reference.json
Using MAFFT: /usr/local/bin/mafft
Writing output to: my_isolates_classify_20260806T142530Z
```

`--organism` is the one thing you always state explicitly — the tool will never guess which reference collection you mean, even when only one is bundled today. Run `treetip classify --help` to see every currently available `--organism` value.

Override anything it would have picked automatically if you need to: `--output my_results` to name the folder yourself, `--mafft /path/to/mafft` if you have more than one copy installed, or `--min-panel-distance-separation` to override the package's own calibrated threshold. Everything not given a value falls back to the automatic choice above.

A run over a handful of isolates typically finishes in well under a minute.

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

Open **`report.html`** inside your output folder in any web browser. This is the file to look at first — it has everything in one place:

- A summary count of how many isolates were called `CLASSICAL`, `CRYPTIC`, `AMBIGUOUS`, or `UNRESOLVED`.
- A table with one row per isolate: its call, its clade, and the distances that call was based on.
- A tree picture for each isolate that was placed, showing where it landed relative to the known reference collection.

The underlying data behind that report also lives in the same folder as plain tables, in case you need to load them into a spreadsheet or another program:

- `panel_assignment.tsv` — the classical/cryptic call, per isolate.
- `clade_assignment.tsv` — the clade call, per isolate.

### What the labels mean

**Panel call** (classical vs. cryptic):

| Label | Meaning |
| --- | --- |
| `CLASSICAL` / `CRYPTIC` | Clearly closer to one reference collection than the other. |
| `AMBIGUOUS` | Close to both — the tool declines to guess. Treat as needing manual review. |
| `UNRESOLVED` | The sequence itself had a quality problem (too much missing data, too little overlap with the reference, etc.) before a call could even be attempted. |

**Confidence status** (`PASS` / `WARN`), shown per isolate in the tables:

- `PASS` means nothing about the placement looked unusual.
- `WARN` means the call was still made, but something about it deserves a second look — for example, the isolate landed in a part of the reference tree where nearby branches carry different clade labels, or the fit to its nearest neighbor was a little looser than typical. `WARN` is not a failure; it is the tool being honest about a borderline case instead of hiding it.

If you see a lot of `WARN` or `AMBIGUOUS` results on isolates you have independent reason to trust, that is worth flagging to whoever maintains the reference package — it may mean the reference collection needs to grow to cover that part of the diversity better.

### Why a call is `WARN`: reading the Findings column

Every `report.html` row and every line of `panel_assignment.tsv` / `clade_assignment.tsv` has a **Findings** column. Whenever a status is `WARN`, `AMBIGUOUS`, or `UNRESOLVED`, this column names the specific reason — it is never left blank. One isolate can show more than one reason at once, separated by `|`.

| Code | In plain terms |
| --- | --- |
| `SELECTED_EDGE_CONSENSUS_DISAGREES_WITH_NEAREST_REFERENCE` | The clade call was decided by the isolate's single closest match in the reference collection, because a broader look at nearby branches would have suggested a different (and less reliable) answer. This is the tool correctly resolving a conflict, not a defect. |
| `MIXED_CLADE_LABELS_ON_SELECTED_EDGE` | The isolate landed in a part of the reference tree where neighboring reference isolates don't all agree on their clade label, so the closest single match was used instead of a broader consensus. |
| `DIFFUSE_DISTANCE_WEIGHT` | The statistical fit to the nearest reference was looser than typical — still the best available match, just not a sharp, unambiguous one. |
| `LONG_TERMINAL_BRANCH` | The isolate is noticeably more distant from its nearest reference than is typical for this collection. Worth a second look if that distance matters for your purpose. |
| `UNSTABLE_BOOTSTRAP_PLACEMENT` | (Only appears if bootstrap replicates were requested.) Repeated resampling of the alignment did not consistently support the same placement. |
| `PANEL_DISTANCE_MARGIN_TOO_SMALL` | The isolate is nearly equally close to both the classical and cryptic collections — this is what produces an `AMBIGUOUS` call. |
| `PARTIAL_COVERAGE` / `LOW_COVERAGE` | Part of the isolate's sequence didn't align to usable reference positions. `LOW_COVERAGE` is severe enough to block a call (`UNRESOLVED`); `PARTIAL_COVERAGE` is not. |
| `DISTANT_FROM_REFERENCE_PANEL` | Even the closest reference isolate is unusually far away — treat the call with extra caution. |
| `NO_COMPARABLE_REFERENCE` / `INSUFFICIENT_COMPARABLE_SITES` | Too little of the isolate's sequence overlaps with the reference collection to compare at all. |

### A note on precision

The current reference collection is built from 7 MLST genes (not full genomes). That is enough to separate classical from cryptic isolates very reliably, and to place an isolate within the right clade. It is **not** enough to resolve fine differences between two isolates that are genetically very close to each other — don't over-interpret small distance differences as if they were whole-genome resolution. If you need finer resolution, ask whoever maintains the reference package whether a core-genome or whole-genome version is planned.

## If something looks wrong

- **The command fails immediately with a message about mismatched columns or IDs**: your FASTA file's sequences are probably not the same 7-locus MLST format the reference was built from. Confirm the format with whoever built the reference package.
- **Every one of your isolates comes back `UNRESOLVED`**: check your FASTA file isn't empty or truncated, and that sequence names don't contain unusual characters.
- **A result surprises you**: open `report.html` and look at the tree picture and the distance numbers for that isolate before assuming the tool is wrong — the numbers are there specifically so you can check the reasoning, not just the label.

For anything else, contact whoever maintains the reference package for your lab.
