"""Runs functions from the hackathon analysis scripts (scripts/03, 05, 06, 08, 10) without running their
pipelines, so the tool uses the same code as the narrative.

Method details: docs/METHODS.md#legacypy
"""
import ast
import hashlib
from functools import lru_cache

import numpy as np
import pandas as pd
from scipy.stats import binomtest, false_discovery_control

from pathlib import Path

from .paths import SCRIPTS

# script -> names to lift (functions or module-level assignments)
SOURCES = {
    "05_pah_grid.py": ["SIGNATURE", "MITO9", "GRID_COLS", "GRID_COLS_V2", "PATH_COLS", "WK", "sign_test", "camera_pr",
                       "col_label", "BLOOD"],
    "10_everywhere.py": ["bh", "collapse", "TP", "TP_ORDER", "GRP", "parse_human", "RAT_ORDER"],
    "03_join.py": ["H_CONTRASTS"],
    "06_discordance.py": ["rho_ci", "PAIRS"],
    "08_discordance_model.py": ["rank_rho", "boot_rho", "B", "RNG", "MIN_CLASS"],
}


VENDORED = Path(__file__).resolve().parent / "legacy_scripts"


def script_dir():
    """hackathon/scripts/ when present (the live originals), else the verbatim copies shipped in the package."""
    return SCRIPTS if all((SCRIPTS / s).exists() for s in SOURCES) else VENDORED


def _lift(script, names):
    src = (script_dir() / script).read_text(encoding="utf-8")
    tree = ast.parse(src)
    keep = []
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in names:
            keep.append(node)
        elif isinstance(node, ast.Assign):
            tg = [t.id for t in node.targets if isinstance(t, ast.Name)]
            if any(t in names for t in tg):
                keep.append(node)
    mod = ast.Module(body=keep, type_ignores=[])
    from scipy.stats import spearmanr
    ns = {"np": np, "pd": pd, "binomtest": binomtest, "false_discovery_control": false_discovery_control,
          "spearmanr": spearmanr}
    exec(compile(mod, str(script_dir() / script), "exec"), ns)
    missing = [n for n in names if n not in ns]
    if missing:
        raise ImportError(f"{script}: could not lift {missing}")
    return {n: ns[n] for n in names}, hashlib.sha256(src.encode()).hexdigest()[:16]


@lru_cache(maxsize=None)
def _all():
    out, hashes = {}, {}
    for script, names in SOURCES.items():
        got, h = _lift(script, names)
        out.update(got)
        hashes[script] = h
    return out, hashes


def get(name):
    return _all()[0][name]


def script_hashes():
    return {**dict(_all()[1]), "_source": "hackathon/scripts" if script_dir() == SCRIPTS else "vendored copies"}


def __getattr__(name):  # legacy.camera_pr etc.
    objs = _all()[0]
    if name in objs:
        return objs[name]
    raise AttributeError(name)
