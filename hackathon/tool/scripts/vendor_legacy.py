#!/usr/bin/env python3
"""Refresh src/motrpac_probe/legacy_scripts/ from hackathon/scripts/ (verbatim copies)."""
import shutil
from pathlib import Path

TOOL = Path(__file__).resolve().parents[1]
from_dir, to_dir = TOOL.parent / "scripts", TOOL / "src" / "motrpac_probe" / "legacy_scripts"
for f in ["03_join.py", "05_pah_grid.py", "06_discordance.py", "08_discordance_model.py", "10_everywhere.py"]:
    shutil.copy2(from_dir / f, to_dir / f)
    print("copied", f)
