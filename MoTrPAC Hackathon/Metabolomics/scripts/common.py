"""Shared paths and helpers for the Metabolomics module scripts."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
WORKBENCH = DATA / "workbench" / "ST000763"
PROCESSED = DATA / "processed"
OUTPUT = ROOT / "output"
for folder in (PROCESSED, OUTPUT):
    folder.mkdir(parents=True, exist_ok=True)

GROUPS = ["Healthy", "LowRisk", "Normal Pressures", "Borderline Pressures", "PAH"]
BLOOD_TIMES = ["during_20_min", "during_40_min", "post_10_min",
               "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr"]
Q = 0.05


def bh(p):
    """Benjamini-Hochberg q values; missing P values stay missing."""
    a = np.asarray(p, dtype=float)
    q = np.full(a.shape, np.nan)
    ok = np.isfinite(a)
    if ok.any():
        q[ok] = multipletests(a[ok], method="fdr_bh")[1]
    return q


def name_key(name):
    """Case- and space-insensitive metabolite name used for exact matching."""
    return re.sub(r"\s+", " ", str(name).strip().casefold())


def load_motrpac_one_feature_per_name():
    """Script 04's MoTrPAC export with one feature per metabolite name and tissue.

    About half of MoTrPAC features have no refmet_name, but their feature_id is
    already a RefMet-style name (for example "PC 40:4"), so it is used instead
    and recorded in match_basis. When several features share a name, the
    choice ignores exercise results: targeted platform first, then the highest
    average abundance.
    """
    mo = pd.read_csv(PROCESSED / "motrpac_metabolite_contrasts.csv.gz")
    assert set(mo.source_package_version.astype(str)) == {"2.0.8"}
    mo["match_basis"] = np.where(mo.refmet_name.notna(), "refmet_name", "feature_id")
    mo["match_key"] = mo.refmet_name.fillna(mo.feature_id).map(name_key)
    choice = (mo.groupby(["tissue", "match_key", "feature_id", "platform", "match_basis"],
                         as_index=False).AveExpr.mean())
    choice["targeted"] = choice.platform.str.startswith("metab-t-")
    choice = (choice.sort_values(["tissue", "match_key", "targeted", "AveExpr", "feature_id"],
                                 ascending=[True, True, False, False, True])
                    .drop_duplicates(["tissue", "match_key"]))
    mo = mo.merge(choice[["tissue", "match_key", "feature_id"]],
                  on=["tissue", "match_key", "feature_id"])
    return mo, choice


def primary_st000763_features(features):
    """One ST000763 assay feature per RefMet name, chosen by coverage and median
    intensity, never by a test result."""
    return (features[features.refmet_name.notna()]
            .sort_values(["refmet_name", "n_measured", "median_intensity", "feature_key"],
                         ascending=[True, False, False, True])
            .drop_duplicates("refmet_name"))
