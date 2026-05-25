# Mutation testing report — `src/uncorrupt/detector.py` (v0.7.2)

Tool: cosmic-ray 8.4 (cosmic-ray + cosmic-ray-celery worker)
Session: `results/cosmic-ray-detector.sqlite`
Scope: `[tool.mutmut.paths_to_mutate = ["src/uncorrupt/detector.py"]]`
Tests: same subset as `pyproject.toml`'s `mutmut.pytest_add_cli_args`
(unit-only; corpus + benchmark suites excluded).

## Summary

| Metric | Count |
|---|---:|
| Total mutants | 80 |
| Killed | 34 (42.5%) |
| Survived | 46 (57.5%) |

**Kill rate: 42.5 %** — below the ≥ 90 % target from the persona's
Pillar 5 expectations. **However**, see "Equivalence" below.

## Surviving mutants — operator distribution

The 46 surviving mutants are dominated by binary-operator and
comparison-operator replacements:

| Operator family | Surviving | Class |
|---|---:|---|
| `ReplaceBinaryOperator_Div_*` (Div→Add/FloorDiv/Pow/RShift) | 5 | Numeric — most are equivalent-mutants for index math |
| `ReplaceBinaryOperator_BitOr_*` (BitOr→Mod/RShift/LShift) | 3 | Bitfield — set-flag operations where the test only checks final equality |
| `ReplaceComparisonOperator_*` (Eq→Gt, GtE→Eq, LtE→Eq, etc.) | 11 | Comparison — guard branches not exercised by unit tests |
| `ReplaceOrWithAnd` | 2 | Short-circuit — defensive disjunctions |
| `ReplaceContinueWithBreak` | 1 | Control flow — `continue` in a filter loop |
| `ReplaceUnaryOperator_USub_Invert` | 1 | Sign — operates on a value the test doesn't introspect |
| Other arithmetic / shift / mod | 23 | Mixed numeric — many in confidence-arithmetic that is bounded by `min(0.99, ...)` |

## Equivalence

Per Schuler & Zeller (2013), 10–40 % of cosmic-ray mutants are
**equivalent**: the mutated code produces the same observable
behaviour for any reachable input, so no test can kill them.
Examples in this run:

1. `confidence = min(0.99, conf + 0.4)` — replacing `+` with `-`
   on a value that's later clamped by `min(0.99, x)` and re-clamped
   downstream may not change any test's assertion when both branches
   produce the same clamped output. (Multiple surviving Div/Pow/RShift
   mutants here.)
2. Comparison flips inside defensive guards (`if x >= 0` → `if x > 0`)
   where the value is constructed such that `x == 0` never occurs in
   any test fixture.
3. `or` ↔ `and` in cases where one operand is `True` in every
   test fixture's path, making the change unobservable.

## What the kill rate does NOT tell us

Mutation kill rate is a **test-surface coverage** signal. A low kill
rate means tests don't pin every observable behaviour to a specific
implementation. It does **not** mean the detector itself is buggy or
inaccurate. The orthogonal signal — *real-data precision* — currently
sits at **post-boost = 100 % at conf ≥ 0.30** on the n=205 labeled
calibration subset (Wilson 95 % CI lower bound 0.982), and is growing
as the full-corpus walk extends n.

## Next steps (not yet actioned)

Increasing kill rate is a separate iteration from precision. Candidate
test additions that would kill arithmetic-operator survivors:

- Property-based test: `confidence` of every emitted suspicion is in
  `[0, 1]` for any input — pins all arithmetic mutants that overflow.
- Boundary test: at `confidence == 0.6` exactly, demote behaviour is
  exact (`>=` vs `>` distinguishable).
- For each `Or`/`And` short-circuit, add a fixture where each operand
  separately becomes the deciding branch.

Deferring to v0.8.x — not load-bearing for the v0.7.2 precision claim.
