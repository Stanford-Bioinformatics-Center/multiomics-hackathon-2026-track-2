"""Filesystem locations. Everything outside tool/ is read-only input."""
import os
from pathlib import Path

TOOL = Path(os.environ.get("MPROBE_TOOL", Path(__file__).resolve().parents[2]))
HACK = Path(os.environ.get("MPROBE_HACKATHON", TOOL.parent))

DATA = HACK / "data"
RAW = DATA / "raw"
JOIN = DATA / "join_table_v2.csv"
SCRIPTS = HACK / "scripts"
EVERY_RAW = HACK / "everywhere" / "raw"
FIGTAB = HACK / "fig" / "tables"
MODEL_EXT = HACK / "model" / "ext"
PROTONLY_EXT = HACK / "protonly" / "ext"
DECK = HACK / "deck_extracts" / "deck_numbers.csv"

STORE = Path(os.environ.get("MPROBE_STORE", TOOL / "store"))  # override to read a staged copy (e.g. /scratch)
CONTRASTS = STORE / "contrasts.parquet"
COLUMNS = STORE / "columns.parquet"
ANNOT = STORE / "annotations.parquet"
PROVENANCE = STORE / "provenance.parquet"
LIBRARY = TOOL / "store" / "library"   # outputs stay in the repo even when the store is read from a staged copy
OUT = TOOL / "out"
SITE = TOOL / "site"
EXAMPLES = TOOL / "examples"
