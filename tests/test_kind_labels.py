"""Coverage guard for plain-English `kind` labels.

The detector emits internal enum values like `gene-date-serial` and
`long-int-precision-loss`. Three places turn those into something a
bench scientist reads:

  - src/uncorrupt/cli.py        KIND_LABELS  (detect subcommand)
  - src/uncorrupt/cli.py        KIND_LABELS  (audit subcommand, separate dict)
  - src/uncorrupt/app.py        _KIND_LABELS (Gradio UI)

A regression where a new detector kind ships without an entry in one of
these dicts surfaces the raw enum to the user. These tests fail in CI
the moment that happens.
"""
from __future__ import annotations

import ast
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src" / "uncorrupt"


def _detector_kinds() -> set[str]:
    """Every `kind="..."` literal passed to a Suspicion in detector.py."""
    text = (SRC / "detector.py").read_text()
    return set(re.findall(r'kind="([a-z][a-z0-9-]+)"', text))


def _label_dicts_in(path: Path) -> list[dict[str, str]]:
    """Return every top-level / module-level dict literal in `path` whose
    name contains 'KIND' and whose values are str-typed (i.e. label maps).
    Multiple such dicts can live in one file (cli.py has two)."""
    tree = ast.parse(path.read_text())
    out: list[dict[str, str]] = []

    class Visitor(ast.NodeVisitor):
        def visit_Assign(self, node: ast.Assign) -> None:
            for t in node.targets:
                if isinstance(t, ast.Name) and "KIND" in t.id.upper():
                    if isinstance(node.value, ast.Dict):
                        d: dict[str, str] = {}
                        for k, v in zip(node.value.keys, node.value.values,
                                         strict=False):
                            if (isinstance(k, ast.Constant)
                                    and isinstance(v, ast.Constant)
                                    and isinstance(k.value, str)
                                    and isinstance(v.value, str)):
                                d[k.value] = v.value
                        if d:
                            out.append(d)
            self.generic_visit(node)

    Visitor().visit(tree)
    return out


def test_detector_kinds_are_nonempty() -> None:
    kinds = _detector_kinds()
    assert len(kinds) >= 10, f"detector.py kind set looks suspiciously small: {kinds}"


def test_cli_detect_dict_covers_every_detector_kind() -> None:
    kinds = _detector_kinds()
    dicts = _label_dicts_in(SRC / "cli.py")
    assert dicts, "cli.py has no KIND_LABELS dict"
    detect_dict = max(dicts, key=len)
    missing = kinds - set(detect_dict)
    assert not missing, (
        f"cli.py detect KIND_LABELS missing entries for: {sorted(missing)}. "
        f"Add a plain-English label for each to avoid leaking the raw enum "
        f"to users."
    )


def test_cli_audit_dict_covers_every_detector_kind() -> None:
    kinds = _detector_kinds()
    dicts = _label_dicts_in(SRC / "cli.py")
    assert len(dicts) >= 2, (
        "cli.py is expected to have two KIND_LABELS dicts (detect + audit)"
    )
    audit_dict = min(dicts, key=len) if len(dicts[0]) != len(dicts[1]) else dicts[1]
    # Both dicts cover the same key set; require it from the smaller one too.
    missing = kinds - set(audit_dict)
    assert not missing, (
        f"cli.py audit KIND_LABELS missing entries for: {sorted(missing)}."
    )


def test_app_dict_covers_every_detector_kind() -> None:
    kinds = _detector_kinds()
    dicts = _label_dicts_in(SRC / "app.py")
    assert dicts, "app.py has no _KIND_LABELS dict"
    ui_dict = dicts[0]
    missing = kinds - set(ui_dict)
    assert not missing, (
        f"app.py _KIND_LABELS missing entries for: {sorted(missing)}. "
        f"Add a plain-English label for each to avoid leaking the raw enum "
        f"to users of the Gradio UI."
    )


def test_no_dict_has_stale_keys() -> None:
    """The reverse direction: catch entries for kinds the detector no longer
    emits. Either the kind was renamed or the dict is stale."""
    kinds = _detector_kinds()
    for path in [SRC / "cli.py", SRC / "app.py"]:
        for d in _label_dicts_in(path):
            stale = set(d) - kinds
            assert not stale, (
                f"{path.name} KIND label dict has entries for kinds not "
                f"emitted by detector.py: {sorted(stale)}. Either remove "
                f"the stale entries or restore the detector code that "
                f"emitted them."
            )
