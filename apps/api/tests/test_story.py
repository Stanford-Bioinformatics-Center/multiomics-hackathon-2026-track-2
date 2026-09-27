"""The story endpoint must carry source values and their interpretation limits."""
from motrpac_probe_service.app import get_story


def _panels():
    body = get_story().model_dump()
    return {panel["id"]: panel for panel in body["panels"]}, body["availability"]


def test_story_source_counts_and_provenance():
    panels, availability = _panels()
    protein = panels["human_muscle_protein"]["records"]
    assert len(protein) == 27
    assert all(r["q_value"] >= .05 for r in protein)
    assert all(r["source"] == "outputs/PAH protein/pah_motrpac_protein_matches.csv" for r in protein)

    rna24 = [r for r in panels["human_muscle_rna"]["records"] if r["timepoint"] == "post_24_hr"]
    assert len(rna24) == 9
    assert all(r["value"] > 0 for r in rna24)
    assert sum(r["q_value"] < .05 for r in rna24) == 7

    go = panels["blood_go"]["records"]
    assert len([r for r in go if r["record_type"] == "pathway"]) == 40
    assert sum(r["metrics"]["same_direction"] for r in go if r["record_type"] == "pathway_time_summary") == 27
    assert sum(r["metrics"]["opposite_direction"] for r in go if r["record_type"] == "pathway_time_summary") == 13
    assert availability["contexts"]
    assert all(c["contrast"] and c["panel_id"] for c in availability["contexts"])


def test_story_distinguishes_correction_families_and_null_interaction():
    panels, _ = _panels()
    rat = panels["rat_muscle_protein"]["records"]
    assert {r["correction_method"] for r in rat} == {"BH", "BY"}
    official_female8 = [r for r in rat if r["record_type"] == "official_sex_specific"
                        and r["sex"] == "female" and r["timepoint"] == "8 wk"]
    targeted_combined8 = [r for r in rat if r["record_type"] == "targeted_comparison"
                          and r["contrast"] == "eight_week_combined"]
    assert sum(r["q_value"] < .05 for r in official_female8) == 3
    assert sum(r["q_value"] < .05 for r in targeted_combined8) == 7

    exercise = panels["metabolomics_exercise"]
    assert "equivalence was not established" in exercise["interpretation"]
    rest = panels["metabolomics_rest"]
    assert "setting and/or SSc" in rest["interpretation"]
    assert all(r["source"].startswith("Metabolomics/output/") for r in rest["records"])


def test_bonferroni_only_where_the_source_bh_family_is_reproduced():
    panels, _ = _panels()
    protein = panels["human_muscle_protein"]["records"]
    # 6,211 proteins per time point: Bonferroni = min(1, p * 6,211), never below the BH q.
    assert all(r["p_value"] is not None and r["bonferroni_p"] is not None for r in protein)
    assert all(r["bonferroni_p"] >= r["q_value"] - 1e-12 for r in protein)
    rna = panels["human_muscle_rna"]["records"]
    assert all((r["bonferroni_p"] is not None) == (r["timepoint"] == "post_24_hr") for r in rna)
    ranks = panels["blood_ranks"]["records"]
    assert all(abs(r["bonferroni_p"] - min(1.0, 6 * r["p_value"])) < 1e-12 for r in ranks)
    # Families that the committed exports cannot reproduce carry no Bonferroni value.
    for panel in ("muscle_oxphos", "blood_biocarta", "metabolomics_context"):
        assert all(r["bonferroni_p"] is None for r in panels[panel]["records"])
    official = [r for r in panels["rat_muscle_protein"]["records"] if r["correction_method"] == "BY"]
    assert official and all(r["bonferroni_p"] is None for r in official)


def test_metabolite_context_carries_all_three_motrpac_contrasts():
    panels, availability = _panels()
    course = [r for r in panels["metabolomics_context"]["records"] if r["record_type"] == "metabolite_timecourse"]
    assert {r["contrast"] for r in course} == {"EE-CON", "EE-EE", "CON-CON"}
    assert len({r["molecule"] for r in course}) == 27  # matched resting hits
    assert {"EE-EE", "CON-CON"} <= set(availability["contrasts"])
