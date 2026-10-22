# Contributing to TreeTip Placer

## Scope and scientific guardrails

TreeTip Placer is a fixed-backbone, distance-based placement tool. Do not
describe it as de novo phylogenetic reconstruction, maximum-likelihood
placement, or a source of posterior probabilities. Preserve the input coordinate
system and report uncertainty through the existing status and sidecar fields.

Do not submit identifiable clinical metadata, restricted database exports, or
credentials in issues or pull requests.

## Development workflow

1. Open an issue describing the biological and software impact.
2. Create a focused branch from `main`; do not mix algorithm, reference-data,
   and documentation changes.
3. Add deterministic synthetic tests before changing behaviour.
4. Run the checks below and include the command output in the pull request.
5. Update `README.md`, `docs/data-contract.md`, and `CHANGELOG.md` when the
   interface, limitation, or validation claim changes.

```bash
python -m pip install -e '.[dev]'
ruff check src tests
ruff format --check src tests
python -m unittest discover -s tests -v
```

## Reference packages

Reference packages are versioned scientific inputs. A change must include the
source provenance, licence or redistribution permission, checksums, alignment
coordinate-system description, masking policy, tree method, validation plan,
and panel-specific QC calibration. Never overwrite an existing package or
replace a calibration result silently.

## Contributions

By submitting a contribution, you agree that it may be distributed under the
Apache-2.0 licence. This does not grant rights over third-party data you include;
only submit material you are authorised to license.
