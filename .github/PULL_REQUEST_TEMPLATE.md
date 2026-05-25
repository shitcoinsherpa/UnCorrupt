<!-- Thanks for the PR. The faster checklist clears, the faster this merges. -->

## What this changes

<!-- One short paragraph. What did you change and why? -->

## How to verify it

<!-- Concrete steps a reviewer can run. Include the test command and any fixture path. -->

```
# example
pytest -q tests/test_detector.py::test_new_thing
```

## Checklist

- [ ] Tests pass locally (`pytest -q`)
- [ ] Lint passes (`ruff check src tests`)
- [ ] Type check passes (`mypy src`)
- [ ] New behaviour has a test (or a fixture in `tests/fixtures/`)
- [ ] If this fixes a missed detection, the regression test pins the exact corruption pattern
- [ ] If this changes detector rules, `docs/methods.md` is updated with the recall / precision delta
- [ ] `CHANGELOG.md` has an entry under `[Unreleased]`
- [ ] No silent data rewrites: detector emits flags, the human commits the change

## Related issues

<!-- Closes #123, related to #456, etc. -->
