"""The generalised explorer: neutral lists in, MoTrPAC results by layer out."""
from motrpac_probe_service import explorer


def _layer(result, layer_id):
    return next(l for l in result["layers"] if l["id"] == layer_id)


def test_examples_are_neutral():
    labels = " ".join(e["label"] for e in explorer.list_examples())
    for word in ("PAH", "diabetes", "T2D", "hypertension"):
        assert word not in labels


def test_protein_list_reproduces_published_matches():
    r = explorer.analyse([{"kind": "genes", "example": "muscle_protein_9"}])
    assert r["inputs"][0]["n_mapped"] == 9
    for layer_id, n_sig_24h in (("proteomics", 0), ("transcriptomics", 7)):
        layer = _layer(r, layer_id)
        cols = layer["columns"]
        vl24 = {j for j, c in enumerate(cols) if c["tissue"] == "VL" and c["category"] == "EE-CON" and c["time"] == "24h"}
        values = [v for v in layer["values"] if v[1] in vl24]
        assert len(values) == 9
        assert sum(v[4] < 0.05 for v in values) == n_sig_24h
    # one family per layer: BH over every tested comparison, never re-run on a subset
    fam = _layer(r, "transcriptomics")["family"]
    assert fam["n_tests"] == sum(c["set_p"] is not None for c in _layer(r, "transcriptomics")["columns"])


def test_ranked_list_uses_rank_correlation():
    r = explorer.analyse([{"kind": "genes", "example": "blood_ranked"}])
    cols = [c for c in _layer(r, "transcriptomics")["columns"] if c["tissue"] == "BLOOD" and c["category"] == "EE-CON"]
    rho = {c["time"]: round(c["rho"], 3) for c in cols}
    assert rho["20m-during"] == -0.079 and rho["15-45m"] == 0.125
    assert all(c["set_t"] is None for c in cols)


def test_plain_list_without_direction_and_unmapped_ids():
    r = explorer.analyse([{"kind": "genes", "text": "PPARGC1A\nNR4A3\nNOTAGENE"}])
    inp = r["inputs"][0]
    assert inp["directed"] is False and inp["n_mapped"] == 2 and inp["unmapped"] == ["NOTAGENE"]
    assert _layer(r, "epigenomics")["available"] is False


def test_saved_example_results_match_the_live_analysis():
    """The website's saved example results must equal a fresh live run (rerun scripts_save_examples.py if not)."""
    import json
    from pathlib import Path
    saved_dir = Path(explorer.REPO) / "apps" / "web" / "public" / "examples"
    for example_id in ("blood_pathway_6", "muscle_protein_9", "plasma_metabolite_41"):
        saved = json.loads((saved_dir / f"{example_id}.json").read_text(encoding="utf-8"))
        live = explorer.analyse([{"kind": explorer.EXAMPLES[example_id]["kind"], "example": example_id}])
        saved.pop("saved")
        assert json.loads(json.dumps(live, allow_nan=False)) == saved, f"{example_id}: rerun apps/api/scripts_save_examples.py"
