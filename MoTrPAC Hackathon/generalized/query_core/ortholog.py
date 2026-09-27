"""Explicit, one-to-one cross-species gene-symbol maps."""

from __future__ import annotations

import csv
from pathlib import Path

from .schema import SPECIES, open_text

ORTHOLOG_COLUMNS = (
    "source_species", "source_gene_symbol", "target_species", "target_gene_symbol",
)


class OrthologMap:
    def __init__(self, path: Path):
        self.path = Path(path)
        self.forward: dict[tuple[str, str, str], str] = {}
        self.by_source: dict[tuple[str, str], list[tuple[str, str]]] = {}
        self.pairs: set[tuple[str, str]] = set()
        reverse: dict[tuple[str, str, str], str] = {}
        with open_text(self.path, "rt") as stream:
            reader = csv.DictReader(stream)
            if not reader.fieldnames or not set(ORTHOLOG_COLUMNS).issubset(reader.fieldnames):
                raise ValueError("Ortholog map needs " + ", ".join(ORTHOLOG_COLUMNS))
            for line, row in enumerate(reader, start=2):
                source = (row["source_species"] or "").strip().lower()
                target = (row["target_species"] or "").strip().lower()
                source_gene = (row["source_gene_symbol"] or "").strip().upper()
                target_gene = (row["target_gene_symbol"] or "").strip().upper()
                if (source not in SPECIES or target not in SPECIES or source == target or
                        not source_gene or not target_gene):
                    raise ValueError(f"{self.path}:{line}: invalid species or symbol")
                key = (source, target, source_gene)
                back = (source, target, target_gene)
                if key in self.forward and self.forward[key] != target_gene:
                    raise ValueError(f"{self.path}:{line}: one-to-many source ortholog")
                if back in reverse and reverse[back] != source_gene:
                    raise ValueError(f"{self.path}:{line}: many-to-one target ortholog")
                self.forward[key] = target_gene
                reverse[back] = source_gene
                self.pairs.add((source, target))
        for (source, target, source_gene), target_gene in self.forward.items():
            self.by_source.setdefault((source, source_gene), []).append(
                (target, target_gene))

    def targets(self, source_species: str, symbol: str):
        yield from self.by_source.get((source_species, symbol), ())

    def allows(self, source_species: str, target_species: str) -> bool:
        return (source_species, target_species) in self.pairs
