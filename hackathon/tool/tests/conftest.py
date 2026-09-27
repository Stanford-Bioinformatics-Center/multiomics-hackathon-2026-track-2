"""Shared fixtures for the motrpac_probe test suite.

* `slow` marker: tests that call out to R (limma). Skipped unless MPROBE_RUN_R=1.
* `S`: the parquet store as the array view used by every analysis (core.get_store(), cached per process).
* `deck`: golden numbers from ../deck_extracts/deck_numbers.csv as {key: value-string}.
"""
import os

import pandas as pd
import pytest

from motrpac_probe import core, paths


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: calls R / limma; skipped unless MPROBE_RUN_R=1")


def pytest_collection_modifyitems(config, items):
    if os.environ.get("MPROBE_RUN_R") == "1":
        return
    skip = pytest.mark.skip(reason="slow R cross-check; set MPROBE_RUN_R=1 to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip)


@pytest.fixture(scope="session")
def S():
    return core.get_store()


@pytest.fixture(scope="session")
def store_view(S):
    """Alias of `S` (the store fixture requested by name in some tests)."""
    return S


@pytest.fixture(scope="session")
def deck():
    d = pd.read_csv(paths.DECK, dtype=str, keep_default_na=False)
    return dict(zip(d.key, d.value))


def deck_close(value, key, deck, tol=None):
    """True if `value` equals the deck number `key` to the precision it was reported with (half a unit in the
    last printed decimal), or within `tol` if given."""
    s = deck[key]
    ref = float(s)
    if tol is None:
        dec = len(s.split(".", 1)[1]) if "." in s and "e" not in s.lower() else 0
        tol = 0.5 * 10 ** (-dec) + 1e-9
    return abs(float(value) - ref) <= tol
