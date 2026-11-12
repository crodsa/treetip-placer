# TreeTip Placer (`v1.0.0`)

**Dual-panel fixed-tree phylogenetic placement and fine-clade assignment (`Clade-1`–`Clade-5` & `C-I`–`C-VII`) for *Clostridioides difficile* seven-locus MLST sequences.**

`TreeTip Placer` takes concatenated 7-locus *C. difficile* MLST nucleotide sequences (`FASTA`, ~3,501 bp), aligns each isolate independently against bundled **classical** (`201` STs) and **cryptic** (`68` isolates) reference panels using `MAFFT --add --keeplength`, enforces a calibrated `0.01` p-distance separation gate, places queries onto fixed reference trees via constrained non-negative least squares, and reports the supported clade (`Clade-1`–`Clade-5` or `C-I`–`C-VII`) in a self-contained `report.html` and machine-readable `TSV`/`JSON`/Newick outputs.

> **Part of a Two-Tier *C. difficile* MLST Suite:**
> - **Tier 1 — [`CdiffCladeML`](https://github.com/crodsa/ML-classifier):** Screens tabular FastMLST/PubMLST allele-number tables in seconds (`Clade1`–`Clade5` vs pooled `Cryptic`).
> - **Tier 2 — `TreeTip Placer` (this repository):** Resolves nucleotide sequences into individual cryptic clades (`C-I` to `C-VII`), rescues `Clade 1`/`Clade 2` boundary STs, types novel-allele (`~n`) environmental strains, and explicitly abstains (`AMBIGUOUS`) on equidistant transitional lineages (such as clade `C7` prototype `ST881`).

---

## Quick Start

### Requirements
- **Python** `3.11+` (zero third-party Python library dependencies—uses pure Python standard library).
- **MAFFT** (`7.x`):
  - **Already have `mafft` on your `PATH`** (e.g., via `brew install mafft`, `sudo apt install mafft`, or Bioconda)? `treetip classify` finds and verifies it automatically—no extra setup needed.
  - **Don't have `mafft` yet?** Run `treetip setup` once after installing (uses `micromamba`, `mamba`, or `conda` on Linux/macOS to create an isolated `MAFFT 7.526` environment without modifying your `PATH`).

### 1. Install `TreeTip Placer`

**From GitHub (includes the runnable 4-isolate demo FASTA):**

```bash
git clone https://github.com/crodsa/treetip-placer.git
cd treetip-placer
python3 -m venv .venv
source .venv/bin/activate
python -m pip install .

# Only needed if `mafft` is not already installed on your PATH:
treetip setup

treetip --version    # Prints: 1.0.0
```

**Or install the standalone release wheel offline** (already bundles the `mlst-cdifficile-20260806` reference package):

```bash
python -m pip install treetip_placer-1.0.0-py3-none-any.whl
```

### 2. Run the bundled 4-isolate demo

```bash
treetip classify examples/demo_isolates.fasta --organism cdifficile-mlst
```

In ~1–2 seconds, this creates `demo_isolates_classify_<UTC_TIMESTAMP>/` demonstrating all core outcomes on real isolates:
- `soil_G15S1_Clade4` → **`CLASSICAL` · `Clade-4`** (`0 SNPs` to `ST109`, margin `0.0363`)
- `soil_G3S6_PB_Cryptic_CIII` → **`CRYPTIC` · `C-III`** (`3 SNPs` to `CLO_CA6460AA`, margin `0.0608`)
- `soil_G9S1_Cryptic_CV` → **`CRYPTIC` · `C-V`** (toxigenic *tcdB*7+ tropical soil isolate, `7 SNPs` to `CLO_AA9252AA`, margin `0.0431`)
- `pubmlst_ST881_transitional` → **`AMBIGUOUS` (`UNRESOLVED`)** (margin `0.0080 < 0.01`; genuine transitional clade `C7` prototype)

Open **`report.html`** inside the output directory in any web browser to explore the summary table and interactive local placement trees.

---

## Preparing Your Input FASTA

`treetip classify` expects an **unaligned nucleotide FASTA** with **one record per isolate** (do **not** pass a whole-genome assembly, raw reads, or an allele-number CSV). Each sequence must concatenate the **7 *C. difficile* housekeeping MLST loci in this exact order** (~3,501 bp total):

```text
adk (501 bp) -> atpA (543 bp) -> dxr (411 bp) -> glyA (516 bp) -> recA (564 bp) -> sodA (450 bp) -> tpi (516 bp)
```

```fasta
>isolate_01
ATG...3501_bp_concatenated_adk_atpA_dxr_glyA_recA_sodA_tpi...TAA
>isolate_02
ATG...3501_bp_concatenated_adk_atpA_dxr_glyA_recA_sodA_tpi...TAA
```

- **Header IDs:** Must be unique and must not duplicate a reference tip ID (`ST1`, `CLO_AA0175AA`, etc.).
- **Partial sequences:** Allowed (`MAFFT --add --keeplength` maps available bases onto the fixed 3,501 reference columns; coverage `< 0.70` or `< 100` comparable sites is caught automatically by QC).

<details>
<summary><strong>Helper snippet: Concatenate per-locus FASTA sequences in the required order</strong></summary>

If your MLST caller exports the 7 locus sequences for each isolate (e.g., `adk`, `atpA`, `dxr`, `glyA`, `recA`, `sodA`, `tpi`), you can join them in Python:

```python
from pathlib import Path

LOCI = ["adk", "atpA", "dxr", "glyA", "recA", "sodA", "tpi"]

def concat_isolate(isolate_id: str, locus_seqs: dict[str, str]) -> str:
    joined = "".join(locus_seqs[locus].strip().upper().replace("-", "") for locus in LOCI)
    return f">{isolate_id}\n{joined}\n"
```

</details>

### Classify your FASTA
```bash
# Auto-creates a timestamped output directory:
treetip classify my_isolates.fasta --organism cdifficile-mlst

# Or specify an explicit output directory and parallel worker count:
treetip classify my_isolates.fasta --organism cdifficile-mlst --output results/run_01 --threads 4
```

---

## Output Files & How to Interpret Results

### Key files in your output directory
| File / Folder | Purpose |
| :--- | :--- |
| **`report.html`** | **Start here.** Self-contained browser report with summary counts, distances, findings, and local placement subtrees per isolate. |
| **`panel_assignment.tsv` / `.json`** | Authoritative `CLASSICAL` / `CRYPTIC` / `AMBIGUOUS` / `UNRESOLVED` call, nearest reference in both panels, and cross-panel `distance_margin`. |
| **`clade_assignment.tsv` / `.json`** | Authoritative clade label (`Clade-1`–`Clade-5` or `C-I`–`C-VII`), `clade_assignment_status` (`PASS`, `WARN`, `UNRESOLVED`), arbitration method, and `findings`. |
| **`selected_trees/*.nwk`** | Ready-to-view Newick trees (`FigTree`, `iTOL`, `ggtree`) with placed queries grafted onto the selected reference tree. |
| **`panels/{classical,cryptic}/`** | Full audit trail for both panels (`placements.tsv`, `alternative_edges.tsv`, `query_qc.tsv`, `run_manifest.json` with SHA-256 hashes). |

### Understanding statuses and common `findings` codes
- **`CLASSICAL` or `CRYPTIC`:** Winning panel was closer by at least `0.01` p-distance (`>= 35` SNPs across 3,501 bp).
- **`AMBIGUOUS` (`PANEL_DISTANCE_MARGIN_TOO_SMALL`):** Isolate passed QC on both panels, but `|d_cryptic - d_classical| < 0.01`. Treat as a genuine borderline/mosaic or novel transitional lineage (e.g., `C7`) and confirm by whole-genome `ANI`/phylogenomics.
- **`PASS` vs `WARN`:**
  - `PASS` means clean QC, concentrated edge fit, and unanimous clade agreement below the placement branch.
  - `WARN` retains a valid clade call while surfacing useful context in `findings`:
    - `DIFFUSE_DISTANCE_WEIGHT`: Very common inside tight MLST clonal complexes where multiple near-zero-length branches differ by 0–2 SNPs and fit equally well (`relative_fit_weight < 0.60`). Check `supporting_reference_tips`—if all neighbouring tips share the same clade, the clade call itself is unambiguous.
    - `MIXED_CLADE_LABELS_ON_SELECTED_EDGE` or `SELECTED_EDGE_CONSENSUS_DISAGREES_WITH_NEAREST_REFERENCE`: Query attached to a basal/polytomous internal branch; `TreeTip Placer` arbitrated to the sequence-nearest labelled reference tip (`nearest_reference`) and flagged it for transparency.
- **`UNRESOLVED` (`LOW_COVERAGE`, `INSUFFICIENT_COMPARABLE_SITES`):** Input sequence had `< 70%` callable coverage across the 3,501-bp marker or `< 100` comparable sites; verify that the input FASTA contains the 7 MLST loci rather than another gene or whole-genome contigs.

---

## Bundled Reference Package & Validation (`v1.0.0`)

`TreeTip Placer v1.0.0` ships with the checksum-verified paired package **[`mlst-cdifficile-20260806`](reference-packages/mlst-cdifficile-20260806/README.md)**:
- **Classical panel (`classical/`):** `201` curated PubMLST STs (`3,501` bp, `IQ-TREE 3` ML tree) across `Clade-1` (`49`), `Clade-2` (`49`), `Clade-3` (`19`), `Clade-4` (`49`), and `Clade-5` (`35`).
- **Cryptic panel (`cryptic/`):** `68` curated EnteroBase/PubMLST isolates (`3,501` bp, `IQ-TREE 3` ML tree) across `C-I` (`19`), `C-II` (`10`), `C-III` (`30`: `C-IIIa` `20`, `C-IIIb` `10`), `C-IV` (`2`), `C-V` (`5`), `C-VI` (`1`), and `C-VII` (`1`).

### Validated benchmark highlights
- **Leave-one-out (`LOO`) internal calibration (`269` reference tips):**
  - **Panel gate (`margin = 0.01`):** `100.0%` (`269/269`: `201/201` classical, `68/68` cryptic; `0` misassignments).
  - **Clade recovery (with nearest-reference arbitration):** `99.63%` (`266/267`: `200/201 = 99.50%` classical; `66/66 = 100.00%` non-singleton cryptic controls), matching independent `FastBAPS` clustering (`99.50%` and `100.00%`).
- **External held-out PubMLST benchmark (`319` STs / `4,650` isolates):** `100.0%` (`319/319`) panel separation, `98.43%` (`314/319`) ST accuracy (`99.81%` [`4,641/4,650`] isolate accuracy; `97.70%` [`212/217`] on STs absent from the classical reference panel), rescuing `14/16` Clade 2 boundary STs missed by allele-table screening and resolving all `25/25` Cryptic STs into `C-I`–`C-VII`.
- **Prospective tropical soil surveillance (`48` Costa Rican *Clostridioides* genomes, ENA BioProject [`PRJEB103842`](https://www.ebi.ac.uk/ena/browser/view/PRJEB103842)):** `100.0%` (`48/48`) panel accuracy and `97.92%` (`47/48`) fine-clade concordance with whole-genome `FastANI` species clusters (plus `1` informative `C-I`/`C-VII` boundary `WARN` on deep-branching `ST997` isolate `G10S3`).

---

## Advanced CLI Subcommands & Documentation

Most users only ever need `treetip classify`. Advanced subcommands support custom paired reference packages and internal calibration suites:

| Subcommand | Purpose |
| :--- | :--- |
| `treetip classify` | One-command routine classification (`--organism cdifficile-mlst`, `--preset {fast,balanced,thorough}`, `--bootstrap N`, `--threads N`). |
| `treetip setup` | Create or refresh (`--refresh`) a dedicated Conda/Mamba `MAFFT 7.526` environment. |
| `treetip place-auto` | Run dual-panel mapping and placement against any custom `paired_reference.json` manifest. |
| `treetip place` | Place pre-aligned fixed-column query FASTAs onto a single reference tree + alignment. |
| `treetip create-reference-package` / `validate-reference-package` | Build and verify SHA-256 manifests for new marker schemes or organisms. |
| `treetip calibrate-panel-gate` / `calibrate-clade-assignment` / `calibrate-placement-accuracy` | Run deterministic leave-one-out calibration suites on a reference package. |

- **Additional documentation:** [Step-by-Step Lab Manual (`docs/MANUAL.md`)](docs/MANUAL.md) · [Architecture (`docs/architecture.md`)](docs/architecture.md) · [Data Contract (`docs/data-contract.md`)](docs/data-contract.md) · [Reference Package Provenance (`reference-packages/mlst-cdifficile-20260806/README.md`)](reference-packages/mlst-cdifficile-20260806/README.md) · [`SECURITY.md`](SECURITY.md)

---

## Citation & License

If you use `TreeTip Placer` in your research or surveillance workflow, please cite [`CITATION.cff`](CITATION.cff) and the accompanying manuscript:

> Rodríguez C. *Two-tier multilocus sequence typing clade classification and fixed-tree phylogenetic placement resolve classical and cryptic Clostridioides difficile lineages.* (2026). Companion Tier-1 repository: https://github.com/crodsa/ML-classifier · Environmental validation genomes: ENA BioProject [`PRJEB103842`](https://www.ebi.ac.uk/ena/browser/view/PRJEB103842).

Source code and documentation are licensed under **[Apache-2.0](LICENSE)**; bundled reference package files are governed by [`reference-packages/mlst-cdifficile-20260806/LICENSE`](reference-packages/mlst-cdifficile-20260806/LICENSE).
