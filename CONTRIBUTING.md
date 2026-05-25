# Contributing to UnCorrupt

Thanks for considering a contribution. The fastest way to land a change is to attach the file that motivated it.

## The single most useful thing you can do

If you found a corruption pattern UnCorrupt missed, or a false positive on a real file, **attach the file** (or the smallest redacted version that still reproduces it) on a bug report. Almost every fix in this codebase started as a one-cell fixture committed alongside a regression test. No file, no test, no fix.

## Set up a dev environment

```bash
git clone https://github.com/<your-fork>/UnCorrupt.git
cd UnCorrupt

python -m venv .venv
source .venv/bin/activate    # Windows: .venv\Scripts\activate

pip install -e ".[dev]"

# Sanity check
pytest -q
ruff check src tests
mypy src
```

That installs the package in editable mode plus pytest, ruff, mypy, hypothesis, build, and pytest-cov.

To run the UI locally while you work:

```bash
uncorrupt-app
```

To run the container build the way CI does:

```bash
docker build -t uncorrupt:dev .
docker run --rm uncorrupt:dev python -c "from uncorrupt.detector import detect_file; print('ok')"
```

## Pull-request checklist

- [ ] `pytest -q` passes locally
- [ ] `ruff check src tests` is clean
- [ ] `mypy src` is clean
- [ ] New behaviour has a unit test, or a fixture in `tests/fixtures/`
- [ ] A missed detection has a regression test pinned to the exact corruption pattern
- [ ] If you changed detector decision rules, `docs/methods.md` shows the recall and precision delta on the published validation corpora
- [ ] `CHANGELOG.md` has a bullet under `[Unreleased]`
- [ ] No new `# noqa` or `# type: ignore` without an inline reason

## What we will not accept

- **Silent data rewrites.** The detector emits flags with confidence scores. The reviewer or the downstream tool commits the change. We never auto-overwrite a cell.
- **New detection rules without evidence.** Every new rule needs a real-world positive example it catches AND a real-world negative example showing it does not false-positive. Synthetic-only test fixtures are not enough.
- **Loosening the column classifier or the placeholder filter to "catch more."** Those layers cost real precision on the published validation walks. Relaxing them needs a full false-positive count, before and after, on the Ziemann S2 and the EPMC-expanded corpus.
- **Skipping the container build.** The published accuracy numbers are reproducible because the container is pinned. Patches that say "works on my machine" need to also say "works in the container."

## Reporting security issues

Do not file security problems in the public issue tracker. Open a private GitHub Security Advisory at https://github.com/shitcoinsherpa/UnCorrupt/security/advisories/new. We respond within 72 hours. Full policy is in [`SECURITY.md`](SECURITY.md).

## Code of conduct

By participating you agree to abide by [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md) (Contributor Covenant 2.1).

## Citation

If your change ships in a release, you appear in the contributors list of the corresponding `CITATION.cff` entry. Real names plus ORCID welcomed.
