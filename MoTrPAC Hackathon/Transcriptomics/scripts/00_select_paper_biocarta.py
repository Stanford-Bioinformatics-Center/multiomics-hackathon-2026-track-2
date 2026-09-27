"""Audit the manual Table 2 transcription and select its six BioCarta rows."""

from pathlib import Path
import hashlib
import json
import pandas as pd

MODULE = Path(__file__).resolve().parents[1]
SOURCE = MODULE / "data" / "paper" / "table2_selected_pathways.csv"
OUT = MODULE / "outputs" / "paper_table2_biocarta.csv"
SETS = {
    "BIOCARTA_CCR5_PATHWAY", "BIOCARTA_IL2_PATHWAY",
    "BIOCARTA_TCR_PATHWAY", "BIOCARTA_IL12_PATHWAY",
    "BIOCARTA_TOB1_PATHWAY", "BIOCARTA_ARENRF2_PATHWAY",
}
SCORE_COLS = ["ipah_page_z", "ssc_ph_ild_page_z", "ssc_pah_page_z", "ssc_page_z"]


def main() -> None:
    manifest = json.loads((MODULE / "data" / "source_manifest.json").read_text())
    for entry in manifest.values():
        path = MODULE / "data" / entry["file"]
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != entry["sha256"]:
            raise ValueError(f"Source checksum differs from manifest: {path}")
    table = pd.read_csv(SOURCE)
    if len(table) != 8 or table["paper_label"].duplicated().any():
        raise ValueError("Expected eight distinct pathway rows transcribed from Table 2")
    if table["annotation"].value_counts().to_dict() != {
        "BioCarta": 6, "GEArray": 1, "SABiosciences": 1
    }:
        raise ValueError("Table 2 annotation counts changed")
    if table[SCORE_COLS].isna().any().any() or not (table[SCORE_COLS] < 0).all().all():
        raise ValueError("Expected four negative PAGE Z-scores on every Table 2 row")
    selected = table.loc[table["annotation"] == "BioCarta"].copy()
    if set(selected["set"]) != SETS or selected["set"].duplicated().any():
        raise ValueError("The six BioCarta names no longer match the audited set list")
    selected = selected.drop(columns="annotation")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(OUT, index=False, float_format="%.3f")
    print("Selected six BioCarta rows from eight Table 2 rows:")
    print(selected[["paper_label", "set"]].to_string(index=False))


if __name__ == "__main__":
    main()
