"""run_analysis(): JSON-safe typed response, _path fix, deterministic run_id, frozen family."""
import json

import pytest

from motrpac_probe_service import AnalysisRequest, SignatureRow, run_analysis, SCHEMA_VERSION
from motrpac_probe_service.schema import AnalysisResponse

pytestmark = pytest.mark.skipif(
    __import__("motrpac_probe.paths", fromlist=["CONTRASTS"]).CONTRASTS.exists() is False,
    reason="no mprobe store: run `mprobe store fetch`")


def _full19():
    return AnalysisRequest(example_name="pah_muscle_malenfant2015", target_species="rat")


def _to_jsonable(resp: AnalysisResponse):
    from dataclasses import asdict
    return asdict(resp)


def test_runs_and_is_json_safe():
    resp = run_analysis(_full19())
    assert isinstance(resp, AnalysisResponse)
    # the whole response must serialize with the stdlib json encoder (no DataFrame/Store/Timer/NaN/np)
    text = json.dumps(_to_jsonable(resp), allow_nan=False)
    assert len(text) > 1000
    assert resp.schema_version == SCHEMA_VERSION
    assert resp.status == "ok"
    assert resp.n_counted_genes == 19
    assert resp.n_input_rows == 25


def test_path_keyerror_is_fixed_provenance_present():
    resp = run_analysis(_full19())
    # provenance() reads opts['_path']; if unfixed this whole call would have raised KeyError('_path')
    assert "signature" in resp.provenance and "sha256_16" in resp.provenance["signature"]
    assert resp.provenance["run_id"] == resp.run_id
    assert "store" in resp.provenance and resp.provenance["store"]["contrasts_parquet_sha256_16"]


def test_multiplicity_family_recorded():
    resp = run_analysis(_full19())
    fam = resp.multiplicity_family
    assert fam.method == "BH"
    assert fam.methods == ["BH", "Bonferroni"]
    assert fam.family_size == len(fam.column_ids) >= 40
    assert fam.n_tests >= 1
    assert fam.threshold == 0.05
    # provenance mirrors it
    assert resp.provenance["multiplicity_family"]["family_size"] == fam.family_size
    assert resp.provenance["multiplicity_family"]["methods"] == fam.methods
    assert resp.provenance["multiplicity_family"]["n_tests"] == fam.n_tests
    for col in resp.columns:
        if col.camera_p is None:
            assert col.camera_bonferroni is None
        else:
            assert col.camera_bonferroni == pytest.approx(min(1.0, col.camera_p * fam.n_tests))


def test_headline_male_rat_skmgn_protein_8wk_not_significant():
    resp = run_analysis(_full19())
    h = resp.headline["male_rat_skm_gn_protein_8wk"]
    assert h["n_opposed"] == 13 and h["n_measured"] == 19
    assert abs(h["camera_t"] - 2.2295) < 1e-3
    assert abs(h["camera_fdr_family"] - 0.0584) < 5e-4       # computed, not hardcoded
    assert h["significant_at_threshold"] is False
    assert h["interpretation"] in ("same direction", "no set-level shift")


def test_run_id_deterministic_and_excludes_presentation_options():
    a = run_analysis(_full19())
    b = run_analysis(_full19())
    assert a.run_id == b.run_id  # deterministic for identical inputs
    # presentation-only options must NOT change the run_id
    r_pres = AnalysisRequest(example_name="pah_muscle_malenfant2015", target_species="rat",
                             include_nonsignificant=False, tissue="SKM-GN", sex="male", timepoint="8w")
    assert run_analysis(r_pres).run_id == a.run_id
    # an analysis parameter (FDR threshold) MUST change the run_id
    r_fdr = AnalysisRequest(example_name="pah_muscle_malenfant2015", target_species="rat", fdr_threshold=0.10)
    assert run_analysis(r_fdr).run_id != a.run_id


def test_features_have_lineage_and_trace_to_columns():
    resp = run_analysis(_full19())
    col_ids = {c.column_id for c in resp.columns}
    assert resp.features, "must return per-feature evidence"
    for f in resp.features[:50]:
        assert f.evidence_id.startswith(resp.run_id)
        assert f.column_id in col_ids           # every feature row traces to a returned column
        if f.measured:
            assert f.opposed in (True, False)


def test_guardrails_and_contrasts_present():
    resp = run_analysis(_full19())
    assert len(resp.guardrails) == 7 and all("safe" in g and "unsafe" in g for g in resp.guardrails)
    assert "rat_train" in resp.contrasts
    c = resp.contrasts["rat_train"]
    assert c["exercise_numerator"] == "trained rat" and c["exercise_denominator"] == "sex-matched sedentary rat"
    assert "healthy" not in c["exercise_numerator"].lower()  # never "healthy gene set"


def test_disease_contrast_uses_fixture_metadata_or_generic_upload_label():
    from motrpac_probe_service.service import _disease_contrast

    pah = _disease_contrast(_full19())
    assert "PAH" in pah[0] and "control" in pah[1]
    t2d = _disease_contrast(AnalysisRequest(example_name="type2_diabetes_muscle_mootha2003"))
    assert "diabetes" in t2d[0] and "glucose" in t2d[1]
    assert "PAH" not in " ".join(t2d)
    uploaded = _disease_contrast(AnalysisRequest(signature_csv_text="gene_symbol,direction\nMYH7,-1\n"))
    assert "source comparator as supplied" == uploaded[1]


def test_signature_rows_input_path_is_content_addressed():
    # supplying rows directly must not require a server file path and must still work
    rows = [SignatureRow(gene_symbol="NDUFA9", direction="-1"),
            SignatureRow(gene_symbol="MYH7", direction="1")]
    resp = run_analysis(AnalysisRequest(signature_rows=rows, target_species="rat", signature_name="tiny"))
    assert resp.status == "ok"
    assert resp.n_counted_genes == 2


def test_feature_lineage_fields_present_and_collapse_recorded():
    """Point -> evidence -> source: every measured feature carries source_feature_id, n_collapsed,
    aggregation_method, and mapping_decision_id. MYH7 collapses 2 features and PDLIM3 collapses 3
    in the rat SKM-GN protein columns (collapse rule = max |stat|)."""
    resp = run_analysis(_full19())
    measured = [f for f in resp.features if f.measured]
    assert measured
    for f in measured:
        assert f.source_feature_id, f"{f.gene}@{f.column_id} missing source_feature_id"
        assert f.n_collapsed is not None and f.n_collapsed >= 1
        assert f.aggregation_method in ("single_feature", "max_abs_stat")
        assert f.mapping_decision_id.startswith(resp.run_id)
        # aggregation method must agree with the collapse count
        assert (f.aggregation_method == "max_abs_stat") == (f.n_collapsed > 1)
    # MYH7 / PDLIM3 collapse in a rat SKM-GN protein column
    skmgn_prot = [f for f in measured if f.dataset == "rat_train" and f.tissue == "SKM-GN" and f.layer == "PROT"]
    myh7 = [f for f in skmgn_prot if f.gene == "MYH7"]
    pdlim3 = [f for f in skmgn_prot if f.gene == "PDLIM3"]
    assert myh7 and myh7[0].n_collapsed == 2 and myh7[0].aggregation_method == "max_abs_stat"
    assert pdlim3 and pdlim3[0].n_collapsed == 3 and pdlim3[0].aggregation_method == "max_abs_stat"


def test_mapping_decision_id_ties_to_audit():
    resp = run_analysis(_full19())
    audit_genes = {r.selected_symbol for r in resp.mapping_audit if r.selected_symbol}
    for f in resp.features:
        # the gene in a mapping_decision_id must be one that appears in the mapping audit
        assert f.gene in audit_genes
