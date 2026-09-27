"""Metabolite-list preview must preserve the module's matching and show misses."""

import csv
import io

import pytest

from motrpac_probe_service.metab_upload import REQUIRED, ROOT, preview_metabolite_list


HEADER = ("list_id,source,source_type,comparison,metabolite_name,refmet_name,"
          "pah_log2_effect,pah_direction,pah_q_value,is_reported_hit,caveat\n")


def test_preview_matches_committed_context_and_reports_unmapped():
    rows = (
        "demo,study,paper_reported,PAH vs control,Eicosadienoic acid,Eicosadienoic acid,1,up,0.01,True,example\n"
        "demo,study,paper_reported,PAH vs control,Unknown compound,Unknown compound,,down,,True,example\n"
    )
    result = preview_metabolite_list(HEADER + rows)
    assert result["n_input_rows"] == 2
    assert result["n_blood_matched"] == 1
    assert result["unmapped_ids"] == ["Unknown compound"]
    first = result["rows"][0]
    assert first["context_label"] == "stable in MoTrPAC blood"
    assert {r["contrast"] for r in first["observations"] if r["tissue"] == "blood"} == {
        "EE-CON", "EE-EE", "CON-CON"}
    assert first["source"].endswith("motrpac_metabolite_contrasts.csv.gz")


def test_preview_rejects_missing_columns_and_invalid_direction():
    with pytest.raises(ValueError, match="required columns"):
        preview_metabolite_list("refmet_name,pah_direction\nA,up\n")
    bad = HEADER + "demo,study,paper_reported,PAH vs control,A,A,,sideways,,True,example\n"
    with pytest.raises(ValueError, match="pah_direction"):
        preview_metabolite_list(bad)


def test_preview_context_labels_match_the_committed_module_for_its_input_list():
    with (ROOT / "output" / "06_context_table_all.csv").open(newline="", encoding="utf-8") as stream:
        expected = list(csv.DictReader(stream))
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=REQUIRED)
    writer.writeheader()
    writer.writerows({key: row.get(key, "") for key in REQUIRED} for row in expected)
    actual = preview_metabolite_list(buffer.getvalue())["rows"]
    assert len(actual) == len(expected)
    assert [row["context_label"] for row in actual] == [row["context_label"] for row in expected]
