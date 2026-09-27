"""Input validation and molecular identifier handling for the query core."""

from __future__ import annotations

import csv
import gzip
import math
import re
from pathlib import Path

DISEASE_COLUMNS = (
    "study_id", "species", "tissue", "layer", "feature_id", "id_namespace",
    "log2_fc", "statistic", "p_value", "q_value", "contrast",
)
REFERENCE_COLUMNS = DISEASE_COLUMNS + ("timepoint", "contrast_category")
OPTIONAL_COLUMNS = ("gene_symbol", "uniprot", "refmet_id", "refmet_name",
                    "pathway_collection", "effect_scale", "statistic_type")
LAYERS = frozenset(("rna", "protein", "phosphosite", "metabolite", "pathway"))
SPECIES = frozenset(("human", "rat"))
MISSING_TOKENS = frozenset(("na", "n/a", "nan", "null", "none", "-"))


def clean_optional(value: str | None) -> str:
    value = (value or "").strip()
    return "" if value.casefold() in MISSING_TOKENS else value


def open_text(path: Path, mode: str):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, mode, encoding="utf-8-sig", newline="")
    return path.open(mode, encoding="utf-8-sig", newline="")


def load_csv(path: Path, *, reference: bool) -> list[dict[str, str]]:
    """Load a signature and reject absent columns or invalid biological axes."""
    path = Path(path)
    required = REFERENCE_COLUMNS if reference else DISEASE_COLUMNS
    with open_text(path, "rt") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"{path}: empty CSV")
        missing = set(required) - set(reader.fieldnames)
        if missing:
            raise ValueError(f"{path}: missing required columns: {', '.join(sorted(missing))}")
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"{path}: duplicate column names")
        rows: list[dict[str, str]] = []
        for line, source in enumerate(reader, start=2):
            if None in source:
                raise ValueError(f"{path}:{line}: more fields than column names")
            row = {key: (value or "").strip() for key, value in source.items()}
            for key in OPTIONAL_COLUMNS:
                row[key] = clean_optional(row.get(key))
            for key in ("study_id", "species", "tissue", "layer", "feature_id",
                        "id_namespace", "contrast") + (("timepoint", "contrast_category")
                                                   if reference else ()):
                if not row[key]:
                    raise ValueError(f"{path}:{line}: {key} is empty")
            row["species"] = row["species"].lower()
            row["tissue"] = row["tissue"].lower()
            row["layer"] = row["layer"].lower()
            if row["species"] not in SPECIES:
                raise ValueError(f"{path}:{line}: unsupported species {row['species']!r}")
            if row["layer"] not in LAYERS:
                raise ValueError(f"{path}:{line}: unsupported layer {row['layer']!r}")
            for key in ("log2_fc", "statistic", "p_value", "q_value"):
                if row[key]:
                    try:
                        number = float(row[key])
                    except ValueError as exc:
                        raise ValueError(f"{path}:{line}: invalid {key} {row[key]!r}") from exc
                    if not math.isfinite(number):
                        raise ValueError(f"{path}:{line}: non-finite {key}")
                    if key in ("p_value", "q_value") and not 0 <= number <= 1:
                        raise ValueError(f"{path}:{line}: {key} must be between 0 and 1")
            row["row_id"] = str(line)
            rows.append(row)
    return rows


def namespace_key(namespace: str, value: str, species: str):
    """Map a declared ID to a typed key; never match opaque feature IDs."""
    ns = namespace.strip().casefold()
    if ns == "hgnc" and species == "human":
        return "gene_symbol", value.upper()
    if ns in {"rgd", "ratgenesymbol"} and species == "rat":
        return "gene_symbol", value.upper()
    if ns == "uniprot":
        return "uniprot", value.upper()
    if ns == "refmet":
        return "refmet_id", value.upper()
    if ns == "refmetname":
        return "refmet_name", name_key(value)
    return None


def name_key(value: str) -> str:
    """Conservative metabolite-name key: case and repeated spaces only."""
    return re.sub(r"\s+", " ", value.strip()).casefold()


def identifiers(row: dict[str, str]) -> tuple[set[tuple[str, str]], set[tuple[str, str]]]:
    """Return high-confidence IDs and lower-confidence exact RefMet names."""
    high: set[tuple[str, str]] = set()
    low: set[tuple[str, str]] = set()
    declared = namespace_key(row["id_namespace"], row["feature_id"], row["species"])
    allowed = ({"gene_symbol", "uniprot"} if row["layer"] in
               {"rna", "protein", "phosphosite"} else
               {"refmet_id", "refmet_name"} if row["layer"] == "metabolite" else
               {"go"})
    if declared and declared[0] in allowed:
        (low if declared[0] == "refmet_name" else high).add(declared)
    if row["layer"] == "pathway" and row["id_namespace"].casefold() in {
            "go", "biocarta", "pathway"} and row["pathway_collection"]:
        # Collection is part of the identity. Equal pathway names in different
        # gene-set versions/collections are not automatically equivalent.
        high.add(("go", row["id_namespace"].upper() + "|" +
                  row["pathway_collection"].upper() + "|" +
                  row["feature_id"].upper()))
    if "gene_symbol" in allowed and row["gene_symbol"]:
        high.add(("gene_symbol", row["gene_symbol"].upper()))
    if "uniprot" in allowed and row["uniprot"]:
        high.add(("uniprot", row["uniprot"].upper()))
    if "refmet_id" in allowed and row["refmet_id"]:
        high.add(("refmet_id", row["refmet_id"].upper()))
    if "refmet_name" in allowed and row["refmet_name"]:
        low.add(("refmet_name", name_key(row["refmet_name"])))
    return high, low


def number(row: dict[str, str], field: str) -> float | None:
    return float(row[field]) if row.get(field, "") else None


def sign_relation(disease: dict[str, str], reference: dict[str, str]) -> str:
    d, r = number(disease, "log2_fc"), number(reference, "log2_fc")
    if (d is None or r is None) and disease["layer"] == reference["layer"] == "pathway":
        d, r = number(disease, "statistic"), number(reference, "statistic")
    if d is None or r is None or d == 0 or r == 0:
        return "zero_or_missing"
    return "same" if (d > 0) == (r > 0) else "opposite"


def sign_basis(disease: dict[str, str], reference: dict[str, str]) -> str:
    if number(disease, "log2_fc") is not None and number(reference, "log2_fc") is not None:
        return "log2_fc"
    if (disease["layer"] == reference["layer"] == "pathway" and
            number(disease, "statistic") is not None and
            number(reference, "statistic") is not None):
        return "signed_statistic"
    return ""
