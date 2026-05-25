# Contributing to UnCorrupt

Thanks for considering a contribution. This project aims to be the most rigorous Excel-corruption detector in genomics — every change needs to defend that bar.

## Where to start

- **Issues with a real reproduction case** are the most welcome. If you can ship a tiny xlsx fixture that the current detector fails on, we can almost always turn it into a test + fix.
- Code style: ruff defaults (line-length 100), mypy strict, pytest. Run `ruff check`, `mypy src`, `pytest -q` before opening a PR.
- The detector adheres to Ziemann's "Five Pillars of Computational Reproducibility" (Briefings in Bioinformatics 2023). Proposals that compromise reproducibility, validation rigor, or audit trail honesty will not land — see [`docs/methods.md`](docs/methods.md) for the project's working stance.

## Development workflow

```bash
git clone <fork-url> uncorrupt && cd uncorrupt

# Install pinned deps locally (same exact versions as the production container)
python -m venv .venv && source .venv/bin/activate
pip install --upgrade pip==24.3.1
pip install --require-hashes -r requirements.txt
pip install -e .

# Run the test suite (must pass before you push)
pytest -q

# Or run inside the production container exactly as CI does
docker build -t uncorrupt:dev .
docker run --rm uncorrupt:dev
```

## Pull-request checklist

- [ ] The full test suite (`pytest -q`) passes locally, in the Docker container, and in CI
- [ ] New behavior is covered by a unit test (or by a fixture in `tests/fixtures/`)
- [ ] If you fixed a detector miss, add a regression test pinned to the exact `(pmc_id, file, sheet, column, row, value)` case — see existing tests in `tests/test_detector.py` for shape
- [ ] If you changed the detector decision rules, update `docs/methods.md` with the rationale and the recall/precision delta against Ziemann S2
- [ ] If you bumped a dependency, regenerate `requirements.txt` via `pip-compile --generate-hashes --output-file=requirements.txt requirements.lock`
- [ ] `CHANGELOG.md` has an entry under `[Unreleased]`
- [ ] No commits add `# noqa` or `# type: ignore` without an inline comment explaining why
- [ ] No commit relaxes a test assertion or removes a failing test without escalating to the maintainers first

## What we will not accept

- "Improvements" that silently rewrite data. The detector emits **flags** with confidence scores; the human reviewer or downstream code does the rewrite. We never auto-correct.
- New detection rules without ground-truth validation. Every rule needs to point to a Ziemann S2 entry (or equivalently documented corruption) it newly catches AND a negative-control sample showing it doesn't false-positive.
- Removing the column classifier or placeholder filtering to "catch more". Those layers cost real precision; relaxing them needs a thorough false-positive count before/after.
- Skipping the Docker container build. The published numbers are reproducible only because the environment is bit-for-bit pinned — patches that say "works on my machine" need to also say "works in the container."

## Reporting security issues

Do not file security issues in the public tracker. Open a GitHub security advisory on the repository (private to maintainers). We'll triage within 48 h.

## Citation

If your change ends up in a release, you'll appear in the contributors list of the corresponding `CITATION.cff` entry. Real names + ORCID welcomed.
