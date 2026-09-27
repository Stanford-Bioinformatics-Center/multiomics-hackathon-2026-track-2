"""Discordance Option-1 adapter tests (test-first for task 2.2).

These tests define the intended contract of the NOT-YET-WRITTEN adapter
`motrpac_probe_service.discordance_casestudy` (created in task 2.2). They are
expected to fail/error at collection time (ImportError) until the adapter exists;
that is deliberate — do not stub the adapter to force them green.

Property tests use Hypothesis (the engine's PBT library) at >=100 iterations each.
Unit tests assert the exact committed per-class counts from
`demo_muscle_ee/run_summary.json` and the `parents[3]` path-coupling constraint.

Grounding source-of-truth files (read-only, committed):
  MoTrPAC Hackathon/generalized/discordance/demo_muscle_ee/catalog.csv
  MoTrPAC Hackathon/generalized/discordance/demo_muscle_ee/model_metrics.csv
  MoTrPAC Hackathon/generalized/discordance/demo_muscle_ee/run_summary.json
  MoTrPAC Hackathon/generalized/discordance/demo_ptm_parent_ee/ptm_parent_summary.csv
"""
from __future__ import annotations

import csv
import json
import math
from collections import Counter
from dataclasses import asdict, is_dataclass
from pathlib import Path
from uuid import uuid4

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

# The adapter under test. It does not exist until task 2.2 — importing it here is
# intentional so this whole module errors (collection ImportError) until then.
from motrpac_probe_service import discordance_casestudy as dc


# --------------------------------------------------------------------------- #
# Committed demo file locations (mirrors the design's planned adapter constants)
# --------------------------------------------------------------------------- #
_REPO = Path(__file__).resolve().parents[3]
_DISCORDANCE_DIR = _REPO / "MoTrPAC Hackathon" / "generalized" / "discordance"
_DEMO_DIR = _DISCORDANCE_DIR / "demo_muscle_ee"
_PTM_DIR = _DISCORDANCE_DIR / "demo_ptm_parent_ee"

_CATALOG_CSV = _DEMO_DIR / "catalog.csv"
_MODEL_METRICS_CSV = _DEMO_DIR / "model_metrics.csv"
_RUN_SUMMARY_JSON = _DEMO_DIR / "run_summary.json"

# Authoritative counts from run_summary.json.catalog_classes (verified by reading
# the committed file; see task 2.1 grounding facts).
_EXPECTED_CLASS_COUNTS = {
    "supported_concordant": 1,
    "supported_opposite": 0,
    "rna_response_protein_equivalent": 285,
    "indeterminate": 5642,
    "no_protein_measurement": 9228,
    "no_rna_measurement": 255,
}

# The fixed model grid: 2 subsets x 3 models = 6 rows.
_EXPECTED_MODEL_SUBSETS = {"all_mapped", "rna_responsive"}
_EXPECTED_MODEL_VARIANTS = {"zero", "rna_only", "temporal"}

pytestmark = pytest.mark.skipif(
    not _DEMO_DIR.is_dir(),
    reason="discordance demo outputs not present",
)


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _is_json_scalar(v) -> bool:
    return v is None or isinstance(v, (str, bool, int, float))


def _assert_json_representable(obj) -> None:
    """Every value in `obj` must be a JSON-representable type with no NaN/Infinity."""
    if isinstance(obj, dict):
        for k, val in obj.items():
            assert isinstance(k, str), f"non-string JSON key: {k!r}"
            _assert_json_representable(val)
    elif isinstance(obj, (list, tuple)):
        for val in obj:
            _assert_json_representable(val)
    else:
        assert _is_json_scalar(obj), f"non-JSON-representable value: {obj!r} ({type(obj)})"
        if isinstance(obj, float):
            assert not math.isnan(obj), "NaN is not JSON-safe"
            assert not math.isinf(obj), "Infinity is not JSON-safe"


def _to_plain(resp):
    """Normalize an adapter response (dataclass or dict) to a plain dict."""
    if is_dataclass(resp) and not isinstance(resp, type):
        return asdict(resp)
    if isinstance(resp, dict):
        return resp
    raise AssertionError(f"unexpected adapter response type: {type(resp)}")


def _catalog_class_counts_from_file() -> dict:
    with _CATALOG_CSV.open(newline="") as fh:
        counts = Counter(row["event_class"] for row in csv.DictReader(fh))
    return dict(counts)


# --------------------------------------------------------------------------- #
# Unit tests — exact committed counts + path coupling
# --------------------------------------------------------------------------- #
def test_parents3_path_coupling():
    """R8 AC7 / design: `Path(adapter).resolve().parents[3]` must equal the repo root.

    The adapter lives at apps/api/motrpac_probe_service/discordance_casestudy.py, so
    parents[3] climbs module -> motrpac_probe_service -> api -> apps -> repo root.
    This is the constraint that lets the adapter locate the committed demo files.
    """
    adapter_file = Path(dc.__file__).resolve()
    assert adapter_file.parents[3] == _REPO
    # And the constant the adapter derives with parents[3] must match the same root.
    assert Path(dc._REPO).resolve() == _REPO


# Non-zero subset of the expected counts. A Counter over catalog.csv (and
# run_summary.json.catalog_classes) only ever emits keys for classes that
# actually have >=1 row, so `supported_opposite` (0 rows in this committed demo)
# is legitimately absent. Comparing the file-derived counts against the full
# _EXPECTED_CLASS_COUNTS (which carries the zero key) is a contradiction the file
# can never satisfy. Compare against the non-zero subset to preserve R2 AC7 parity
# without the impossible key.
_EXPECTED_NONZERO_CLASS_COUNTS = {k: v for k, v in _EXPECTED_CLASS_COUNTS.items() if v > 0}


def test_committed_catalog_class_counts_exact():
    """R2 AC7 — per-class counts read from the catalog equal the committed values.

    `supported_opposite` has zero rows in this demo, so a Counter never emits it;
    parity is asserted against the non-zero expected subset."""
    counts = _catalog_class_counts_from_file()
    assert counts == _EXPECTED_NONZERO_CLASS_COUNTS


def test_run_summary_catalog_classes_exact():
    """R2 AC7 — run_summary.json.catalog_classes equals the committed values.

    Like catalog.csv, the committed run_summary omits the zero-count
    `supported_opposite`; parity is asserted against the non-zero expected subset."""
    summary = json.loads(_RUN_SUMMARY_JSON.read_text())
    assert summary["catalog_classes"] == _EXPECTED_NONZERO_CLASS_COUNTS


def test_build_catalog_summary_classification_and_coverage_counts():
    """The adapter surfaces the four classification classes (zero-filled) and the
    two coverage states separately, with values exactly as committed."""
    resp = _to_plain(dc.build_catalog_summary())
    classification = resp["classification_counts"]
    coverage = resp["coverage_counts"]

    assert classification == {
        "supported_concordant": 1,
        "supported_opposite": 0,  # valid class, zero rows in this demo
        "rna_response_protein_equivalent": 285,
        "indeterminate": 5642,
    }
    assert coverage == {
        "no_protein_measurement": 9228,
        "no_rna_measurement": 255,
    }


def test_build_model_has_exactly_six_rows_over_the_fixed_grid():
    """R2 / Property 3 — model_metrics is exactly 6 rows = 2 subsets x 3 models."""
    resp = _to_plain(dc.build_model())
    metrics = resp["metrics"]
    assert len(metrics) == 6
    seen = {(m["subset"], m["model"]) for m in metrics}
    assert seen == {
        (s, v) for s in _EXPECTED_MODEL_SUBSETS for v in _EXPECTED_MODEL_VARIANTS
    }


def test_build_ptm_parent_summary_and_candidate_count():
    """R2 — PTM-parent summary rows (one per timepoint) + committed candidate count (749)."""
    resp = _to_plain(dc.build_ptm_parent())
    assert len(resp["summary_rows"]) == 3  # three timepoints in the committed demo
    assert resp["candidate_count"] == 749


# --------------------------------------------------------------------------- #
# Property 1 — read-fidelity + JSON-safety
# Feature: finalize-discordance-mvp, Property 1: For any committed discordance
# demo output, every scalar in the adapter's response equals the value present in
# the source file (NaN/Infinity/unparseable -> null), and the entire response
# serializes to valid JSON containing only JSON-representable types.
# Validates: Requirements 2.2, 2.4, 11.1
# --------------------------------------------------------------------------- #
def _numify(cell):
    """The adapter's numeric-normalization contract: NaN/inf/unparseable -> None."""
    if cell is None or cell == "":
        return None
    try:
        f = float(cell)
    except (TypeError, ValueError):
        return None
    return None if (math.isnan(f) or math.isinf(f)) else f


@settings(max_examples=100, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(builder=st.sampled_from(["catalog", "model", "ptm"]))
def test_property1_read_fidelity_and_json_safety(builder):
    # The generated input selects which committed resource we exercise; every branch
    # is driven by the concrete committed corpus (read-fidelity is per-value).
    if builder == "catalog":
        resp = _to_plain(dc.build_catalog_summary())
    elif builder == "model":
        resp = _to_plain(dc.build_model())
    else:
        resp = _to_plain(dc.build_ptm_parent())

    # (a) JSON-safety: only JSON-representable types, no NaN/Infinity, and it must
    #     actually serialize with allow_nan=False.
    _assert_json_representable(resp)
    json.dumps(resp, allow_nan=False)

    # (b) read-fidelity: adapter scalars equal the source cells (with NaN/inf/unparseable
    #     normalized to None). Verified against the source-of-truth files.
    if builder == "model":
        with _MODEL_METRICS_CSV.open(newline="") as fh:
            src = list(csv.DictReader(fh))
        by_key = {(r["subset"], r["model"]): r for r in src}
        for m in resp["metrics"]:
            row = by_key[(m["subset"], m["model"])]
            assert int(m["n_rows"]) == int(float(row["n_rows"]))
            assert int(m["n_genes"]) == int(float(row["n_genes"]))
            for col in ("mae", "rmse", "r2"):
                assert m[col] == _numify(row[col])
    elif builder == "catalog":
        # classification + coverage counts equal the by-reading counts of the file.
        file_counts = _catalog_class_counts_from_file()
        merged = dict(resp["classification_counts"])
        merged.update(resp["coverage_counts"])
        # supported_opposite is zero-filled by the adapter but absent from the file.
        for cls, n in merged.items():
            assert n == file_counts.get(cls, 0)
    else:  # ptm
        with (_PTM_DIR / "ptm_parent_summary.csv").open(newline="") as fh:
            src = list(csv.DictReader(fh))
        assert len(resp["summary_rows"]) == len(src)


# --------------------------------------------------------------------------- #
# Property 3 — count parity + fixed 6-row grid
# Feature: finalize-discordance-mvp, Property 3: For any committed discordance demo,
# the adapter's per-class counts read from catalog.csv equal run_summary.json.
# catalog_classes, and model_metrics contains exactly six rows over
# {all_mapped, rna_responsive} x {zero, rna_only, temporal}.
# Validates: Requirements 2.7
# --------------------------------------------------------------------------- #
@settings(max_examples=100, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(_=st.integers())
def test_property3_count_parity_and_fixed_model_grid(_):
    # count parity: catalog.csv per-class counts == run_summary.json.catalog_classes
    file_counts = _catalog_class_counts_from_file()
    summary = json.loads(_RUN_SUMMARY_JSON.read_text())
    assert file_counts == summary["catalog_classes"]

    # the adapter must reflect the same parity through its own surfaces
    cat = _to_plain(dc.build_catalog_summary())
    adapter_counts = dict(cat["classification_counts"])
    adapter_counts.update(cat["coverage_counts"])
    for cls, n in summary["catalog_classes"].items():
        assert adapter_counts[cls] == n

    # fixed 6-row model grid
    model = _to_plain(dc.build_model())
    metrics = model["metrics"]
    assert len(metrics) == 6
    assert {(m["subset"], m["model"]) for m in metrics} == {
        (s, v) for s in _EXPECTED_MODEL_SUBSETS for v in _EXPECTED_MODEL_VARIANTS
    }


# --------------------------------------------------------------------------- #
# Property 4 — event-class domain gate
# Feature: finalize-discordance-mvp, Property 4: For any catalog whose required
# files exist, catalog_available() returns True exactly when every row's event_class
# is a member of KNOWN_EVENT_CLASSES; a single missing/malformed/out-of-set label
# makes catalog_available() return False and the catalog is not returned.
# Validates: Requirements 2.6, 2.8
# --------------------------------------------------------------------------- #
def _known_classes() -> frozenset:
    return frozenset(dc.KNOWN_EVENT_CLASSES)


@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(
    labels=st.lists(
        st.one_of(
            st.sampled_from(sorted(_EXPECTED_CLASS_COUNTS.keys())),  # in-set labels
            st.just(""),                                             # missing/empty
            st.text(min_size=1, max_size=12),                        # arbitrary (mostly out-of-set)
        ),
        min_size=1,
        max_size=20,
    )
)
def test_property4_event_class_domain_gate(labels, tmp_path, monkeypatch):
    """Write a synthetic catalog + run_summary into a temp demo dir, point the adapter
    at it, and assert catalog_available() is True iff every label is in the known set."""
    known = _known_classes()

    # Hypothesis reuses the function-scoped `tmp_path` across generated examples, so a fixed
    # `tmp_path / "demo"` would raise FileExistsError on the 2nd example. Give each example a
    # unique subdir.
    demo = tmp_path / f"demo_{uuid4().hex}"
    demo.mkdir()
    # minimal catalog.csv with just the event_class column the gate reads
    with (demo / "catalog.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["event_class"])
        for lab in labels:
            w.writerow([lab])
    # run_summary.json must also exist for the "required files present" precondition
    (demo / "run_summary.json").write_text(json.dumps({"catalog_classes": {}}))

    # Repoint the adapter's file constants at the synthetic demo dir.
    monkeypatch.setattr(dc, "DEMO_DIR", demo, raising=False)
    monkeypatch.setattr(dc, "CATALOG_CSV", demo / "catalog.csv", raising=False)
    monkeypatch.setattr(dc, "RUN_SUMMARY_JSON", demo / "run_summary.json", raising=False)

    expected_available = all(lab in known for lab in labels)
    assert dc.catalog_available() is expected_available

    # When the domain gate fails, the catalog must not be returned (adapter raises).
    if not expected_available:
        with pytest.raises((FileNotFoundError, ValueError)):
            dc.build_catalog_summary()


# --------------------------------------------------------------------------- #
# Endpoint tests (task 2.3) — the three GET-only discordance endpoints wired in
# app.py. Confirm each returns 200 with the expected shape when the committed
# demo is present, and that they are read-only (no write/create/update/delete).
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def client():
    from fastapi.testclient import TestClient
    from motrpac_probe_service.app import app

    return TestClient(app)


def test_endpoint_discordance_catalog(client):
    r = client.get("/api/discordance/catalog")
    assert r.status_code == 200
    body = r.json()
    assert body["analysis_type"] == "discordance_catalog"
    assert body["classification_counts"] == {
        "supported_concordant": 1,
        "supported_opposite": 0,
        "rna_response_protein_equivalent": 285,
        "indeterminate": 5642,
    }
    assert body["coverage_counts"] == {
        "no_protein_measurement": 9228,
        "no_rna_measurement": 255,
    }
    assert body["occupancy_caveat"]  # PTM caveat carried
    _assert_json_representable(body)


def test_endpoint_discordance_model(client):
    r = client.get("/api/discordance/model")
    assert r.status_code == 200
    body = r.json()
    assert body["analysis_type"] == "discordance_model"
    assert len(body["metrics"]) == 6
    assert {(m["subset"], m["model"]) for m in body["metrics"]} == {
        (s, v) for s in _EXPECTED_MODEL_SUBSETS for v in _EXPECTED_MODEL_VARIANTS
    }
    assert body["weak_prediction_note"]  # honest weak-prediction framing carried
    _assert_json_representable(body)


def test_endpoint_discordance_ptm_parent(client):
    r = client.get("/api/discordance/ptm-parent")
    assert r.status_code == 200
    body = r.json()
    assert body["analysis_type"] == "ptm_parent_audit"
    assert len(body["summary_rows"]) == 3
    assert body["candidate_count"] == 749
    assert body["occupancy_caveat"]
    _assert_json_representable(body)


def test_endpoints_are_read_only(client):
    """GET-only: the discordance routes reject write verbs (405), honoring the
    no write/create/update/delete constraint (R2 AC1)."""
    for path in (
        "/api/discordance/catalog",
        "/api/discordance/model",
        "/api/discordance/ptm-parent",
    ):
        for verb in ("post", "put", "patch", "delete"):
            resp = getattr(client, verb)(path)
            assert resp.status_code == 405, f"{verb.upper()} {path} should be 405, got {resp.status_code}"


def test_endpoint_catalog_unavailable_returns_404(client, monkeypatch):
    """When a required committed output is unavailable, the endpoint returns 404
    (FileNotFoundError -> 404), matching get_metab_casestudy — never fabricated data."""
    from motrpac_probe_service import discordance_casestudy as _dc

    monkeypatch.setattr(_dc, "catalog_available", lambda: False)
    r = client.get("/api/discordance/catalog")
    assert r.status_code == 404


# --------------------------------------------------------------------------- #
# Property 2 — read-only idempotence + byte-stability
# Feature: finalize-discordance-mvp, Property 2: For any sequence of GET requests
# to the discordance endpoints, the committed source files are never modified
# (their bytes are identical before and after), and every build of a given
# resource returns a payload equal to every other build of that resource.
# Validates: Requirements 2.1, 2.5
# --------------------------------------------------------------------------- #

# The committed source files the adapter reads. A byte-stable, read-only adapter
# must leave every one of these untouched no matter how many GETs it serves.
_COMMITTED_SOURCE_FILES = (
    _DEMO_DIR / "catalog.csv",
    _DEMO_DIR / "model_metrics.csv",
    _DEMO_DIR / "model_predictions.csv",
    _DEMO_DIR / "run_summary.json",
    _PTM_DIR / "ptm_parent_summary.csv",
    _PTM_DIR / "ptm_parent_candidates.csv",
)

# The three GET-only discordance endpoints exercised by the generated sequence.
_DISCORDANCE_ENDPOINTS = (
    "/api/discordance/catalog",
    "/api/discordance/model",
    "/api/discordance/ptm-parent",
)


def _snapshot_committed_files() -> dict:
    """Map each committed source file to its exact current bytes (byte fingerprint)."""
    return {p: p.read_bytes() for p in _COMMITTED_SOURCE_FILES if p.exists()}


# Each example reads six committed files, drives several HTTP GETs, and rebuilds every
# resource twice — legitimately I/O-bound, so the per-example timing deadline is disabled
# (byte-identity, not latency, is what this property asserts).
@settings(
    max_examples=100,
    deadline=None,
    suppress_health_check=[HealthCheck.function_scoped_fixture],
)
@given(
    # A generated sequence of GETs across the three endpoints (order + repetition vary).
    sequence=st.lists(
        st.sampled_from(_DISCORDANCE_ENDPOINTS),
        min_size=1,
        max_size=12,
    )
)
def test_property2_read_only_idempotence_and_byte_stability(client, sequence):
    # (a) byte-stability: fingerprint the committed source-file bytes BEFORE the GETs,
    #     run the generated GET sequence, then confirm every file is byte-identical.
    before = _snapshot_committed_files()
    assert before, "expected committed discordance source files to exist"

    for path in sequence:
        r = client.get(path)
        assert r.status_code == 200, f"GET {path} -> {r.status_code}"

    after = _snapshot_committed_files()
    assert after == before, "GET sequence must not create/delete/modify committed files"

    # (b) idempotence: repeated builds of the SAME resource are equal to each other.
    #     The payload is a pure function of the committed bytes, so two independent
    #     builds must be identical (no hidden state, no drift).
    assert _to_plain(dc.build_catalog_summary()) == _to_plain(dc.build_catalog_summary())
    assert _to_plain(dc.build_model()) == _to_plain(dc.build_model())
    assert _to_plain(dc.build_ptm_parent()) == _to_plain(dc.build_ptm_parent())

    # And repeated GETs of the same endpoint return byte-identical response bodies.
    for path in _DISCORDANCE_ENDPOINTS:
        first = client.get(path)
        second = client.get(path)
        assert first.status_code == 200 and second.status_code == 200
        assert first.json() == second.json(), f"repeated GET {path} must be idempotent"

    # Finally, the source files are still unchanged after the idempotence checks too.
    assert _snapshot_committed_files() == before
