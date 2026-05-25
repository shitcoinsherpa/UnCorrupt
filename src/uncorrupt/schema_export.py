"""Emit Frictionless Table Schema sidecars from detector reports.

The detection side tells the scientist *what's wrong now*. The schema sidecar
tells *every downstream consumer of this data* what each column is supposed
to be : so corruption introduced later (by Excel, pandas-with-default-args,
or any other re-handling) gets flagged at the next import boundary instead
of silently propagating.

Frictionless Table Schema spec: https://specs.frictionlessdata.io/table-schema/
The sidecar is a JSON file conventionally named `<datafile>.schema.json` or
embedded in a Data Package `datapackage.json`.

Per Pillar 4 (FAIR data sharing, Ziemann et al. 2023): identifier columns
get a schema entry pinning the type to `string` and a `constraints.pattern`
regex matching the registry. Downstream readers honoring the schema will
refuse to coerce : preventing the same corruption at import time.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from .detector import Report
from .registries import (
    GENE_LIKE_PATTERN,
    HGNC_NEW_SYMBOLS,
    RIKEN_PATTERN,
)


@dataclass
class SchemaField:
    """One field (column) in a Frictionless Table Schema."""
    name: str
    type: str  # "string", "number", "integer", "date", "datetime", "boolean"
    format: str | None = None
    constraints: dict[str, Any] | None = None
    description: str | None = None
    rdfType: str | None = None  # ontology URL for semantic typing

    def to_dict(self) -> dict[str, Any]:
        d: dict[str, Any] = {"name": self.name, "type": self.type}
        if self.format:
            d["format"] = self.format
        if self.constraints:
            d["constraints"] = self.constraints
        if self.description:
            d["description"] = self.description
        if self.rdfType:
            d["rdfType"] = self.rdfType
        return d


def _infer_field_type(series: pd.Series) -> str:
    """Coarse type inference for a column with no detector signal."""
    non_null = series.dropna()
    if len(non_null) == 0:
        return "string"
    types = {type(v).__name__ for v in non_null.head(50)}
    if types <= {"int", "float"} or types == {"int"}:
        return "integer" if types == {"int"} else "number"
    if types == {"bool"}:
        return "boolean"
    if "date" in types or "datetime" in types or "Timestamp" in types:
        return "date"
    return "string"


def _identifier_constraints(column_name: str, series: pd.Series) -> tuple[str, dict[str, Any], str]:
    """Return (description, constraints, registry_label) for an identifier column.

    Decides which registry pattern best matches the column's content and
    pins both type and a regex constraint accordingly.
    """
    non_null = series.dropna()
    sample_strings = [v for v in non_null.head(200) if isinstance(v, str)]

    n_known_hgnc = sum(1 for v in sample_strings if v in HGNC_NEW_SYMBOLS)
    n_riken = sum(1 for v in sample_strings if RIKEN_PATTERN.fullmatch(v))
    n_gene_like = sum(1 for v in sample_strings if GENE_LIKE_PATTERN.fullmatch(v))

    if n_known_hgnc >= 1 or (n_gene_like > 0 and n_gene_like >= n_riken):
        return (
            "HGNC gene symbol : pinned to text to prevent date / numeric coercion",
            {"pattern": "^[A-Z][A-Z0-9-]{1,14}$"},
            "HGNC",
        )
    if n_riken >= 1:
        return (
            "RIKEN-style identifier : pinned to text to prevent scientific-notation coercion",
            {"pattern": "^\\d{7}[A-Z]\\d{2}$"},
            "RIKEN",
        )
    return (
        "Identifier column : pinned to text to prevent any auto-coercion on re-import",
        {},
        "generic",
    )


def emit_table_schema(
    df: pd.DataFrame,
    report: Report,
    name: str = "table",
) -> dict[str, Any]:
    """Build a Frictionless Table Schema dict from a DataFrame + detector report.

    Identifier columns flagged in `report.identifier_columns` get type="string"
    plus a regex pattern constraint matching the inferred registry. Other
    columns get coarse type inference.
    """
    # report.identifier_columns may be in "Sheet!Col" form for multi-sheet
    id_cols: set[str] = set()
    for c in report.identifier_columns:
        # Take the part after "!" if present (per detect_file convention)
        id_cols.add(c.rsplit("!", 1)[-1])

    fields: list[dict[str, Any]] = []
    for col in df.columns:
        col_str = str(col)
        if col_str in id_cols:
            description, constraints, _registry = _identifier_constraints(col_str, df[col])
            f = SchemaField(
                name=col_str, type="string",
                constraints=constraints or None,
                description=description,
            )
        else:
            f = SchemaField(name=col_str, type=_infer_field_type(df[col]))
        fields.append(f.to_dict())

    return {
        "$schema": "https://specs.frictionlessdata.io/schemas/table-schema.json",
        "name": name,
        "fields": fields,
        "missingValues": ["", "NA", "N/A", "null", "None"],
    }


def write_schema_sidecar(
    df: pd.DataFrame,
    report: Report,
    sidecar_path: Path,
    name: str | None = None,
) -> Path:
    """Write a Frictionless Table Schema JSON next to a data file.

    Convention: if `data.csv` exists, write `data.csv.schema.json`.
    """
    schema = emit_table_schema(df, report, name=name or sidecar_path.stem)
    sidecar_path.parent.mkdir(parents=True, exist_ok=True)
    with sidecar_path.open("w") as fh:
        json.dump(schema, fh, indent=2)
    return sidecar_path
