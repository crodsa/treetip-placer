# Data contract

## Reference package

The MLST MVP consumes two coupled artifacts:

1. A branch-length Newick tree.
2. A tab-separated profile table with one `sample_id` column and one column per locus.

The tree tips and reference `sample_id` values must be identical sets. Tree branch lengths must be finite and non-negative. The program accepts rooted or arbitrarily rooted Newick representations; distance placement is performed on the underlying unrooted tree metric.

## Query profiles

Query tables use the same `sample_id` column and locus names as the reference. Columns may appear in a different order. Unknown extra columns are rejected. The following values are treated as missing by default: an empty string, `.`, `-`, `?`, `NA`, `N`, and `NULL`.

An allele identifier is a categorical state. For example, allele `18` has no implied numerical proximity to allele `17`.

## QC policy

The pre-placement gate checks schema compatibility, missingness, comparable locus count, ambiguous or mixed allele calls, duplicate profiles, and nearest-reference distance. It emits `PASS`, `WARN`, or `FAIL`; it does not silently discard a query.

The post-placement gate can mark a placed query as `UNRESOLVED` when its relative-fit distribution is diffuse or its fitted terminal branch is an outlier.

Thresholds are configuration values, not biological constants. They should be calibrated on leave-one-out placement experiments using the intended reference panel.

## Outputs

`placed_tree.nwk` contains all successfully placed queries, independently scored against the original reference tree. `placed_tree.nexus` contains the same tree in Nexus form. `placements.tsv` is the authoritative human-readable table. `run_manifest.json` records input checksums, configuration, timing, and software version.

The experimental `placements.jplace` file uses a distance-objective field rather than likelihood values and must not be interpreted as a maximum-likelihood JPlace result.
