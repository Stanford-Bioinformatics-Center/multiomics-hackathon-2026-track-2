"""Bridge to the analysis scripts in hackathon/scripts/.

Those scripts run their whole pipeline at module level, so they cannot be imported directly. This module
parses each script, keeps only the named top-level functions and constant assignments, and executes those
nodes in a fresh namespace. The code that runs is the scripts' own source (no copies), so a change in the
scripts propagates to the tool; `SOURCES` records which names come from where, for provenance.json.
"""
import ast
import hashlib
from functools import lru_cache

import numpy as np
import pandas as pd
from scipy.stats import binomtest, false_discovery_control

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


def _lift(script, names):
    src = (SCRIPTS / script).read_text()
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
    exec(compile(mod, str(SCRIPTS / script), "exec"), ns)
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
    return dict(_all()[1])


def __getattr__(name):  # legacy.camera_pr etc.
    objs = _all()[0]
    if name in objs:
        return objs[name]
    raise AttributeError(name)
