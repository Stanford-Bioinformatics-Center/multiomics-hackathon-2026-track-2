"""Read-only story panels from the committed MoTrPAC Hackathon result tables.

This adapter selects and labels existing results. It does not fit models or adjust
P values. Each record names the source table so plots can expose their provenance.
"""
from __future__ import annotations

import csv
import gzip
import json
import math
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd
from pydantic import BaseModel, Field


ROOT = Path(__file__).resolve().parents[3] / "MoTrPAC Hackathon"


class StoryRecord(BaseModel):
    id: str
    panel_id: str
    record_type: str
    species: str
    tissue: str
    layer: str
    exercise: Optional[str] = None
    contrast: str
    timepoint: Optional[str] = None
    sex: Optional[str] = None
    molecule: Optional[str] = None
    disease_contrast: Optional[str] = None
    disease_layer: Optional[str] = None
    disease_direction: Optional[str] = None
    disease_unit: Optional[str] = None
    motrpac_direction: Optional[str] = None
    disease_value: Optional[float] = None
    disease_p_value: Optional[float] = None
    disease_q_value: Optional[float] = None
    motrpac_measured: bool = True
    value: Optional[float] = None
    unit: Optional[str] = None
    n: Optional[int] = None
    p_value: Optional[float] = None
    q_value: Optional[float] = None
    bonferroni_p: Optional[float] = None
    correction_method: Optional[str] = None
    correction_family: Optional[str] = None
    label: str
    source: str
    metrics: dict[str, Any] = Field(default_factory=dict)


class StoryPanel(BaseModel):
    id: str
    title: str
    question: str
    observation: str
    interpretation: str
    hypothesis: str
    caveat: str
    records: list[StoryRecord]


class StoryResponse(BaseModel):
    panels: list[StoryPanel]
    availability: dict[str, Any]
    provenance: dict[str, str]


def _rows(source: str) -> list[dict[str, str]]:
    path = ROOT / source
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def _num(value: Any) -> Optional[float]:
    if value is None or value == "":
        return None
    result = float(value)
    return result if math.isfinite(result) else None


def _direction(value: Any) -> Optional[str]:
    number = _num(value)
    if number is None:
        return None
    return "up" if number > 0 else "down" if number < 0 else "no change"


def _table(source: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / source)


def _bh(p: np.ndarray) -> np.ndarray:
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


def _verified_bonferroni(table: pd.DataFrame, family: list[str], p: str = "p_value",
                         q: str = "adj_p_value") -> Optional[pd.Series]:
    """Bonferroni over a source family, only when BH over the same rows reproduces the source q.

    If the committed table does not reproduce the source BH values, the family is not known
    exactly and no Bonferroni value is reported.
    """
    rows = table.dropna(subset=[p])
    out = pd.Series(np.nan, index=table.index)
    groups = rows.groupby(family) if family else [(None, rows)]
    for _, group in groups:
        pv = group[p].to_numpy(float)
        if np.max(np.abs(_bh(pv) - group[q].to_numpy(float))) > 1e-8:
            return None
        out.loc[group.index] = np.minimum(pv * len(group), 1.0)
    return out


def _record(panel: str, index: int, source: str, **fields: Any) -> StoryRecord:
    return StoryRecord(id=f"{panel}:{index}", panel_id=panel, source=source, **fields)


def _panel(id: str, title: str, question: str, observation: str,
           interpretation: str, hypothesis: str, caveat: str,
           records: list[StoryRecord]) -> StoryPanel:
    return StoryPanel(id=id, title=title, question=question,
                      observation=observation, interpretation=interpretation,
                      hypothesis=hypothesis, caveat=caveat, records=records)


def _human_muscle_protein() -> StoryPanel:
    source = "outputs/PAH protein/pah_motrpac_protein_matches.csv"
    data = _rows(source)
    family = _table("data/processed/motrpac_muscle_protein_ee_con.csv.gz")
    family["bonferroni"] = _verified_bonferroni(family, ["Timepoint"])
    full = family.set_index(["feature_id", "Timepoint"])
    records = [_record(
        "human_muscle_protein", i, source,
        record_type="molecule", species="human", tissue="skeletal muscle",
        layer="proteomics", exercise="endurance", contrast="EE-CON",
        timepoint=row["motrpac_timepoint"], sex="all", molecule=row["motrpac_symbol"],
        disease_contrast="resting PAH / healthy control muscle protein",
        disease_layer="proteomics",
        disease_direction="down", disease_unit="PAH/control ratio",
        motrpac_direction=_direction(row["motrpac_logFC"]),
        disease_value=_num(row["pah_to_control_ratio"]),
        value=_num(row["motrpac_logFC"]), unit="log2 fold change",
        q_value=_num(row["motrpac_bh_adj_p_value"]),
        p_value=_num(full.loc[(row["uniprot_accession"], row["motrpac_timepoint"]), "p_value"]),
        bonferroni_p=_num(full.loc[(row["uniprot_accession"], row["motrpac_timepoint"]), "bonferroni"]),
        correction_method="BH",
        correction_family="MoTrPAC muscle protein EE-CON, 6,211 proteins per time point",
        label=f"{row['paper_symbol']} (PAH/control protein ratio {row['pah_to_control_ratio']})",
        metrics={"paper_symbol": row["paper_symbol"],
                 "uniprot": row["uniprot_accession"],
                 "paper_reported_p_value_text": row["paper_reported_p_value_text"]},
    ) for i, row in enumerate(data)]
    significant = sum(r.q_value is not None and r.q_value < 0.05 for r in records)
    return _panel(
        "human_muscle_protein", "Human muscle protein",
        "Do the nine proteins lower in resting PAH muscle change after acute exercise in healthy muscle?",
        f"{significant}/{len(records)} matched protein-time results have MoTrPAC BH q < 0.05.",
        "The measured protein layer has no statistically clear response for these nine markers at the sampled times.",
        "RNA and protein may respond on different time scales.",
        "The PAH paper and MoTrPAC use separate cohorts; nonsignificant does not mean an exact zero effect.",
        records,
    )


def _human_muscle_rna() -> StoryPanel:
    source = "outputs/PAH protein/pah_muscle_rna_timepoint_detail.csv"
    data = _rows(source)
    family = _table("data/processed/motrpac_muscle_rna_ee_con_24h.csv.gz")
    family["bonferroni"] = _verified_bonferroni(family, ["Timepoint"])
    full = family.set_index(["feature_id", "Timepoint"])

    def lookup(row: dict[str, str], column: str) -> Optional[float]:
        key = (row["motrpac_rna_feature_id"], row["motrpac_rna_timepoint"])
        return _num(full.loc[key, column]) if key in full.index else None
    records = [_record(
        "human_muscle_rna", i, source,
        record_type="molecule", species="human", tissue="skeletal muscle",
        layer="transcriptomics", exercise="endurance", contrast="EE-CON",
        timepoint=row["motrpac_rna_timepoint"], sex="all",
        molecule=row["motrpac_current_symbol"], disease_direction="down",
        disease_contrast="resting PAH / healthy control muscle protein",
        disease_layer="proteomics",
        motrpac_direction=_direction(row["motrpac_rna_logFC"]),
        value=_num(row["motrpac_rna_logFC"]), unit="log2 fold change",
        q_value=_num(row["motrpac_rna_bh_adj_p_value"]),
        p_value=lookup(row, "p_value"), bonferroni_p=lookup(row, "bonferroni"),
        correction_method="BH",
        correction_family="MoTrPAC muscle RNA EE-CON, all genes per time point",
        label=row["motrpac_current_symbol"],
        metrics={"paper_symbol": row["paper_symbol"],
                 "uniprot": row["uniprot_accession"],
                 "feature_id": row["motrpac_rna_feature_id"]},
    ) for i, row in enumerate(data)]
    summary = {r["motrpac_rna_timepoint"]: r for r in _rows(
        "outputs/PAH protein/pah_muscle_rna_timepoint_summary.csv")}
    late = summary["post_24_hr"]
    return _panel(
        "human_muscle_rna", "Human muscle RNA",
        "When do the matching nine genes change in healthy muscle RNA?",
        f"At 24 h, {late['positive_logFC_count']}/9 estimates are positive; "
        f"{late['positive_and_bh_q_lt_0_05_count']}/9 have BH q < 0.05.",
        "The 24-hour RNA pattern differs from the nonsignificant protein results at the same sampled time.",
        "A later protein response could be tested with later sampling.",
        "The 24-hour focus was chosen after seeing all three times; RNA abundance does not establish protein change.",
        records,
    )


def _muscle_oxphos() -> StoryPanel:
    source = "outputs/PAH protein/motrpac_muscle_oxphos_camera.csv"
    data = _rows(source)
    records = [_record(
        "muscle_oxphos", i, source,
        record_type="pathway", species="human", tissue="skeletal muscle",
        layer="transcriptomics", exercise="endurance", contrast=row["contrast_category"],
        timepoint=row["Timepoint"], sex="all", molecule=row["set"],
        motrpac_direction=row["direction"].lower(), value=_num(row["z.std"]),
        unit="CAMERA z", n=int(row["set_size"]), p_value=_num(row["p_value"]),
        q_value=_num(row["adj_p_value"]), correction_method="BH",
        correction_family="MoTrPAC source CAMERA family",
        label="GO oxidative phosphorylation",
        metrics={"collection": row["collection"], "database": row["database"]},
    ) for i, row in enumerate(data)]
    return _panel(
        "muscle_oxphos", "Muscle OXPHOS RNA",
        "Does the broad oxidative-phosphorylation RNA program change over time?",
        f"Down early (q {records[0].q_value:.3g}); up at 24 h (q {records[-1].q_value:.3g}).",
        "The pathway-level RNA signal changes direction over time in healthy muscle.",
        "A delayed mitochondrial transcript response might precede later protein remodeling.",
        "This is a MoTrPAC pathway test, not a PAH muscle pathway test; middle-time q is nonsignificant.",
        records,
    )


def _blood_ranks() -> StoryPanel:
    source = "outputs/PAH blood/pah_motrpac_blood_rank_association_by_time.csv"
    data = _rows(source)
    bonferroni = _verified_bonferroni(_table(source), [], p="permutation_p_two_sided",
                                      q="permutation_fdr_six_times")
    records = [_record(
        "blood_ranks", i, source,
        record_type="rank_correlation", species="human", tissue="blood",
        layer="transcriptomics", exercise="endurance", contrast="EE-CON",
        timepoint=row["Timepoint"], sex="all",
        disease_contrast="GSE33463 IPAH - healthy PBMC",
        disease_layer="transcriptomics",
        motrpac_direction="opposite ranks" if float(row["spearman_rho"]) < 0 else "same-direction ranks",
        value=_num(row["spearman_rho"]), unit="Spearman rho",
        n=int(row["shared_genes"]), p_value=_num(row["permutation_p_two_sided"]),
        q_value=_num(row["permutation_fdr_six_times"]),
        bonferroni_p=None if bonferroni is None else _num(bonferroni.iloc[i]),
        correction_method="BH of permutation P",
        correction_family="six time points; 10,000 gene-label shuffles per time",
        label=f"IPAH PBMC t vs healthy MoTrPAC whole-blood RNA z at {row['Timepoint']}",
        metrics={"permutations": int(row["permutations"]),
                 "permutation_null": row["permutation_null"]},
    ) for i, row in enumerate(data)]
    return _panel(
        "blood_ranks", "Blood RNA ranks",
        "How do IPAH PBMC gene ranks relate to healthy exercise whole-blood RNA ranks?",
        f"{records[0].n:,} matched genes/time; rho ranges from "
        f"{min(r.value for r in records):+.3f} to {max(r.value for r in records):+.3f}.",
        "The associations are small and change sign after exercise.",
        "Cell-mixture changes could contribute to the rank pattern.",
        "Separate PBMC and whole-blood cohorts; shuffling gene labels does not preserve gene dependence.",
        records,
    )


def _blood_go() -> StoryPanel:
    summary_source = "outputs/PAH blood/pah_motrpac_blood_gobp_time_summary.csv"
    match_source = "outputs/PAH blood/pah_motrpac_blood_gobp_by_time.csv.gz"
    summary = _rows(summary_source)
    records = [_record(
        "blood_go", i, summary_source,
        record_type="pathway_time_summary", species="human", tissue="blood",
        layer="transcriptomics", exercise="endurance", contrast="EE-CON",
        timepoint=row["Timepoint"], sex="all",
        value=int(row["both_bh_below_0.05"]), unit="jointly BH-significant GO labels",
        n=int(row["shared_sets"]),
        label=f"{row['both_bh_below_0.05']} of {row['shared_sets']} labels at {row['Timepoint']}",
        metrics={"same_direction": int(row["same_direction_both_bh"]),
                 "opposite_direction": int(row["opposite_direction_both_bh"])},
    ) for i, row in enumerate(summary)]
    motrpac_family = _table("outputs/PAH blood/motrpac_blood_gobp_camera.csv.gz")
    motrpac_family["bonferroni"] = _verified_bonferroni(motrpac_family, ["contrast_category", "Timepoint"])
    go_bonferroni = motrpac_family.set_index(["Timepoint", "set"])["bonferroni"]
    for row in _rows(match_source):
        if row["both_bh_below_0.05"].lower() != "true":
            continue
        records.append(_record(
            "blood_go", len(records), match_source,
            record_type="pathway", species="human", tissue="blood",
            layer="transcriptomics", exercise="endurance", contrast="EE-CON",
            timepoint=row["Timepoint"], sex="all", molecule=row["set"],
            disease_contrast="GSE33463 IPAH - healthy PBMC",
            disease_layer="transcriptomics",
            disease_direction=row["pah_direction"].lower(),
            disease_unit="PAH cameraPR z",
            motrpac_direction=row["motrpac_direction"].lower(),
            disease_value=_num(row["pah_camera_z"]),
            disease_p_value=_num(row["pah_p_value"]),
            disease_q_value=_num(row["pah_fdr"]),
            value=_num(row["motrpac_z"]), unit="MoTrPAC CAMERA z",
            p_value=_num(row["motrpac_p_value"]), q_value=_num(row["motrpac_fdr"]),
            bonferroni_p=_num(go_bonferroni.get((row["Timepoint"], row["set"]))),
            correction_method="BH in each source analysis",
            correction_family="5,183 PAH GOBP sets; original MoTrPAC GOBP family",
            label=row["set_short"],
            metrics={"direction_relation": row["direction_relation"],
                     "pah_set_genes": int(row["pah_set_genes"]),
                     "motrpac_set_genes": int(row["motrpac_set_genes"])},
        ))
    same = sum(int(r["same_direction_both_bh"]) for r in summary)
    opposite = sum(int(r["opposite_direction_both_bh"]) for r in summary)
    return _panel(
        "blood_go", "Blood GO pathways",
        "Which GO labels reach BH q < 0.05 in both separate blood analyses?",
        f"{summary[0]['shared_sets']} labels/time; {same + opposite} set-time rows: "
        f"{same} same direction, {opposite} opposite.",
        "Both-direction patterns occur; shared GO genes make labels nonindependent.",
        "Cell-composition markers are candidates for targeted follow-up.",
        "This exact-name join is descriptive, not a disease-by-exercise test.",
        records,
    )


def _blood_biocarta() -> StoryPanel:
    source = "Transcriptomics/outputs/paper_motrpac_biocarta_mapped.csv"
    data = _rows(source)
    records = [_record(
        "blood_biocarta", i, source,
        record_type="pathway", species="human", tissue="blood",
        layer="transcriptomics", exercise="endurance" if row["modality"] == "EE" else "resistance",
        contrast="EE-CON" if row["modality"] == "EE" else "RE-CON",
        timepoint=row["timepoint"], sex="all", molecule=row["set"],
        disease_contrast="Cheadle 2012 IPAH vs healthy PBMC",
        disease_layer="transcriptomics",
        disease_direction=_direction(row["ipah_page_z"]),
        disease_unit="PAGE z",
        motrpac_direction=_direction(row["z.std"]),
        disease_value=_num(row["ipah_page_z"]), value=_num(row["z.std"]),
        unit="MoTrPAC CAMERA z", n=int(row["set_size"]),
        p_value=_num(row["p_value"]), q_value=_num(row["adj_p_value"]),
        correction_method="BH", correction_family="original MoTrPAC CAMERA family",
        label=row["paper_label"],
        metrics={"paper_measure": "IPAH vs control PAGE z",
                 "motrpac_measure": "healthy whole-blood CAMERA z"},
    ) for i, row in enumerate(data)]
    positive = sum(r.value is not None and r.value > 0 for r in records)
    positive_q = sum(r.value is not None and r.value > 0 and r.q_value is not None and r.q_value < .05 for r in records)
    negative_q = sum(r.value is not None and r.value < 0 and r.q_value is not None and r.q_value < .05 for r in records)
    return _panel(
        "blood_biocarta", "Blood BioCarta",
        "How do six Cheadle 2012 IPAH-lower BioCarta labels behave in healthy blood?",
        f"{positive}/{len(records)} pathway-time z scores are positive; "
        f"{positive_q} positive and {negative_q} negative have BH q < 0.05.",
        "A positive MoTrPAC pathway z is a group tendency, not an effect in PAH patients.",
        "Different exercise modalities could produce different pathway profiles.",
        "The paper's PAGE scores and MoTrPAC CAMERA scores use different samples and may use different gene membership.",
        records,
    )


def _rat_muscle_protein() -> StoryPanel:
    official_source = "Rat Comparison PAH Proteins/output/01_official_sex_specific_weekly.csv"
    weekly_source = "Rat Comparison PAH Proteins/output/03_training_vs_sedentary_by_week.csv"
    compare_source = "Rat Comparison PAH Proteins/output/04_training_comparisons.csv"
    records: list[StoryRecord] = []
    for row in _rows(official_source):
        records.append(_record(
            "rat_muscle_protein", len(records), official_source,
            record_type="official_sex_specific", species="rat", tissue="skeletal muscle",
            layer="proteomics", exercise="endurance training", contrast="trained minus sedentary",
            timepoint=f"{row['week']} wk", sex=row["sex"], molecule=row["paper_symbol"],
            disease_contrast="resting PAH / healthy control human muscle protein",
            disease_layer="proteomics",
            disease_direction="down", motrpac_direction=_direction(row["logFC"]),
            value=_num(row["logFC"]), unit="normalized log2 protein difference",
            p_value=_num(row["p_value"]), q_value=_num(row["adj_p_value"]),
            correction_method="BY", correction_family="official MoTrPAC broad differential-abundance family",
            label=row["paper_symbol"],
            metrics={"feature_id": row["feature_ID"],
                     "selection_fdr": _num(row["selection_fdr"]),
                     "bonferroni_note": "Unavailable: official rat results use a BY-adjusted family."},
        ))
    weekly_bonferroni = _verified_bonferroni(_table(weekly_source), ["week"], q="BH_q_9_within_week")
    compare_bonferroni = _verified_bonferroni(_table(compare_source), ["comparison"], q="BH_q_within_nine")
    for i, row in enumerate(_rows(weekly_source)):
        records.append(_record(
            "rat_muscle_protein", len(records), weekly_source,
            record_type="sex_adjusted_weekly", species="rat", tissue="skeletal muscle",
            layer="proteomics", exercise="endurance training", contrast="trained minus sedentary",
            timepoint=f"{row['week']} wk", sex="both", molecule=row["paper_symbol"],
            disease_contrast="resting PAH / healthy control human muscle protein",
            disease_layer="proteomics",
            disease_direction="down", motrpac_direction=_direction(row["trained_minus_sedentary"]),
            value=_num(row["trained_minus_sedentary"]), unit="normalized log2 protein difference",
            p_value=_num(row["p_value"]), q_value=_num(row["BH_q_9_within_week"]),
            bonferroni_p=None if weekly_bonferroni is None else _num(weekly_bonferroni.iloc[i]),
            correction_method="BH", correction_family="nine selected proteins within each week",
            label=row["paper_symbol"],
            metrics={"standard_error": _num(row["standard_error"])},
        ))
    for i, row in enumerate(_rows(compare_source)):
        comparison = row["comparison"]
        records.append(_record(
            "rat_muscle_protein", len(records), compare_source,
            record_type="targeted_comparison", species="rat", tissue="skeletal muscle",
            layer="proteomics", exercise="endurance training", contrast=comparison,
            timepoint="8 wk" if comparison.startswith("eight_week") else "all sampled weeks",
            sex="female" if comparison.endswith("female") else "male" if comparison.endswith("male") else "both",
            molecule=row["paper_symbol"], disease_direction="down",
            disease_contrast="resting PAH / healthy control human muscle protein",
            disease_layer="proteomics",
            motrpac_direction=_direction(row["trained_minus_sedentary"]),
            value=_num(row["trained_minus_sedentary"]), unit="normalized log2 protein difference",
            n=int(row["n_trained"]) + int(row["n_sedentary"]),
            p_value=_num(row["p_value"]), q_value=_num(row["BH_q_within_nine"]),
            bonferroni_p=None if compare_bonferroni is None else _num(compare_bonferroni.iloc[i]),
            correction_method="BH", correction_family=f"nine selected proteins in {comparison}",
            label=row["paper_symbol"],
            metrics={"n_trained": int(row["n_trained"]),
                     "n_sedentary": int(row["n_sedentary"])},
        ))
    week8 = [r for r in records if r.record_type == "targeted_comparison" and r.contrast == "eight_week_combined"]
    week8_q = sum(r.q_value is not None and r.q_value < .05 for r in week8)
    return _panel(
        "rat_muscle_protein", "Rat muscle protein",
        "How do the same nine proteins behave over 1–8 weeks of healthy-rat training?",
        f"At 8 wk, {sum(r.value > 0 for r in week8)}/9 targeted sex-combined estimates are positive; "
        f"{week8_q}/9 have BH q < 0.05.",
        "The official sex-specific and targeted nine-protein analyses have different multiplicity families.",
        "A later protein pattern can be tested in an independent training cohort.",
        "Different rats occupy each trained week; sedentary muscle was sampled only at the 8-week endpoint.",
        records,
    )


def _metabolomics_exercise() -> StoryPanel:
    source = "Metabolomics/output/02_exercise_response_summary.csv"
    data = _rows(source)
    records = [_record(
        "metabolomics_exercise", i, source,
        record_type="interaction_summary" if row["test"].startswith("change difference") else "within_group_summary",
        species="human", tissue="blood", layer="metabolomics",
        exercise="peak exercise (modality not specified in output)",
        contrast=row["test"], timepoint="rest to peak", sex="all",
        disease_contrast=row["test"],
        disease_layer="metabolomics",
        motrpac_measured=False,
        value=int(row["n_q_lt_0_05"]), unit="features with BH q < 0.05",
        n=int(row["n_eligible"]),
        correction_method="BH", correction_family=f"eligible features in {row['test']}",
        label=row["test"],
        metrics={"raw_p_lt_0_05": int(row["n_raw_p_lt_0_05"]),
                 "minimum_q": _num(row["smallest_q"])},
    ) for i, row in enumerate(data)]
    target = next(r for r in records if r.contrast == "change difference, PAH_vs_Healthy")
    return _panel(
        "metabolomics_exercise", "PAH exercise metabolites",
        "Is the rest-to-peak metabolite change different in SSc-PAH?",
        f"PAH vs healthy: {int(target.value)}/{target.n} feature interactions have BH q < 0.05.",
        "No feature-level difference was detected under this test; equivalence was not established.",
        "A verified-pairing study could test whether responses truly match across groups.",
        "The 108 rest/peak pairs are inferred, not confirmed by participant IDs.",
        records,
    )


def _metabolomics_rest() -> StoryPanel:
    source = "Metabolomics/output/03_resting_summary.csv"
    data = _rows(source)
    records = [_record(
        "metabolomics_rest", i, source,
        record_type="resting_summary", species="human", tissue="blood",
        layer="metabolomics", contrast=row["test"] + " at rest",
        disease_contrast=row["test"] + " at rest",
        disease_layer="metabolomics",
        motrpac_measured=False,
        timepoint="rest", sex="all", value=int(row["n_q_lt_0_05"]),
        unit="features with BH q < 0.05", n=int(row["n_eligible"]),
        correction_method="BH",
        correction_family=f"eligible resting features in {row['test']}",
        label=row["test"],
        metrics={"unique_refmet_hits": int(row["n_unique_refmet"]),
                 "minimum_q": _num(row["smallest_q"])},
    ) for i, row in enumerate(data)]
    healthy = next(r for r in records if r.label == "PAH_vs_Healthy")
    normal = next(r for r in records if r.label == "PAH_vs_Normal_Pressures")
    return _panel(
        "metabolomics_rest", "Resting PAH metabolites",
        "Do resting PAH plasma levels differ from samples collected in the same setting?",
        f"PAH vs healthy: {int(healthy.value)}/{healthy.n} features; "
        f"PAH vs normal-pressure SSc: {int(normal.value)}/{normal.n} at BH q < 0.05.",
        "The PAH–healthy difference tracks sampling setting and/or SSc, which this design cannot separate.",
        "Sampling procedures may alter measured metabolite levels.",
        "All PAH samples are Cath and all healthy samples Non-invasive; a zero-hit comparison does not prove equality.",
        records,
    )


def _metabolomics_context() -> StoryPanel:
    source = "Metabolomics/output/06_context_table_hits.csv"
    summary_source = "Metabolomics/output/06_context_summary.json"
    with (ROOT / summary_source).open(encoding="utf-8") as stream:
        summary = json.load(stream)["lists"]["st000763_rest_pah_vs_healthy"]
    data = _rows(source)
    context_all = _table("Metabolomics/output/06_context_table_all.csv")
    hits = context_all[context_all.is_reported_hit.astype(str).str.lower() == "true"].set_index("refmet_name")
    records = [_record(
        "metabolomics_context", i, source,
        record_type="metabolite_context", species="human", tissue="blood",
        layer="metabolomics", exercise="endurance", contrast="EE-CON",
        timepoint="during_20_min", sex="all", molecule=row["refmet_name"],
        disease_contrast="ST000763 SSc-PAH - healthy at rest (Cath vs Non-invasive)",
        disease_layer="metabolomics",
        motrpac_measured=row["context_label"] != "no MoTrPAC blood match",
        disease_direction=row["pah_direction"],
        disease_unit="log2 plasma level difference",
        motrpac_direction=_direction(row["blood_ee_con_logFC_during_20_min"]),
        disease_value=_num(row["pah_log2_effect"]),
        disease_q_value=_num(row["pah_q_value"]),
        value=_num(row["blood_ee_con_logFC_during_20_min"]),
        q_value=_num(hits.loc[row["refmet_name"], "blood_ee_con_q_during_20_min"]),
        unit="log2 fold change", correction_method="BH plus |log2 FC| >= 0.5 for context label",
        correction_family="original source families; context label is a size-and-q rule",
        label=row["context_label"],
        metrics={"motrpac_blood_match_basis": row["motrpac_blood_match_basis"],
                 "blood_peak_time": row["blood_peak_time"],
                 "blood_peak_logFC": _num(row["blood_peak_logFC"]),
                 "blood_max_abs_control_logFC": _num(row["blood_max_abs_control_logFC"]),
                 "pah_effect_over_control_drift": _num(row["pah_effect_over_control_drift"]),
                 "exercise_vs_pah_direction": row["exercise_vs_pah_direction"],
                 "bonferroni_note": "Unavailable: the exported rows do not reproduce MoTrPAC's metabolomics BH family."},
    ) for i, row in enumerate(data)]
    course_source = "Metabolomics/data/processed/motrpac_metabolite_contrasts.csv.gz"
    course = _table(course_source)
    course = course[course.tissue == "blood"]
    matched = hits[hits.motrpac_blood_feature_id.notna()]
    exercise = {"EE-CON": "endurance", "EE-EE": "endurance", "CON-CON": "control"}
    for name, hit in matched.iterrows():
        rows = course[course.feature_id == hit.motrpac_blood_feature_id]
        for row in rows.itertuples():
            if row.contrast_category == "EE-CON" and row.Timepoint == "during_20_min":
                continue  # already shown as the metabolite's context row
            records.append(_record(
                "metabolomics_context", len(records), course_source,
                record_type="metabolite_timecourse", species="human", tissue="blood",
                layer="metabolomics", exercise=exercise[row.contrast_category],
                contrast=row.contrast_category, timepoint=row.Timepoint, sex="all", molecule=name,
                disease_contrast="ST000763 SSc-PAH - healthy at rest (Cath vs Non-invasive)",
                disease_layer="metabolomics", disease_direction=hit.pah_direction,
                disease_unit="log2 plasma level difference",
                motrpac_direction=_direction(row.logFC), disease_value=_num(hit.pah_log2_effect),
                disease_q_value=_num(hit.pah_q_value), value=_num(row.logFC), unit="log2 fold change",
                p_value=_num(row.p_value), q_value=_num(row.adj_p_value), correction_method="BH",
                correction_family="MoTrPAC source metabolomics family",
                label=f"{name}: {row.contrast_category} {row.Timepoint}",
                metrics={"feature_id": row.feature_id, "platform": row.platform,
                         "bonferroni_note": "Unavailable: the exported rows do not reproduce MoTrPAC's metabolomics BH family."},
            ))
    return _panel(
        "metabolomics_context", "MoTrPAC metabolite context",
        "Are the resting PAH–healthy hits exercise-sensitive or prone to control drift in healthy blood?",
        f"{summary['hit_labels']['stable in MoTrPAC blood']}/{summary['n_hits_matched_blood']} matched hits are "
        f"stable; median PAH difference is {summary['median_pah_effect_over_control_drift_hits']:.1f}× "
        "largest MoTrPAC control drift.",
        "Most matched hits do not meet the combined q-and-effect-size rule for exercise or resting drift.",
        "Matched sampling conditions could sharpen future PAH biomarker studies.",
        "These are separate cohorts; 14 resting hits have no MoTrPAC blood match, and time-of-day effects cannot be ruled out.",
        records,
    )


def build_story() -> StoryResponse:
    panels = [
        _human_muscle_protein(), _human_muscle_rna(), _muscle_oxphos(),
        _blood_ranks(), _blood_go(), _blood_biocarta(), _rat_muscle_protein(),
        _metabolomics_exercise(), _metabolomics_rest(), _metabolomics_context(),
    ]
    keys = ("panel_id", "species", "tissue", "layer", "exercise", "contrast", "timepoint", "sex")
    contexts = [dict(zip(keys, values)) for values in sorted(
        {tuple(getattr(r, key) for key in keys) for panel in panels for r in panel.records},
        key=lambda values: tuple("" if value is None else str(value) for value in values),
    )]
    availability = {
        "species": sorted({r.species for p in panels for r in p.records}),
        "tissues": sorted({r.tissue for p in panels for r in p.records}),
        "layers": sorted({r.layer for p in panels for r in p.records}),
        "exercises": sorted({r.exercise for p in panels for r in p.records if r.exercise}),
        "contrasts": sorted({r.contrast for p in panels for r in p.records}),
        "timepoints": sorted({r.timepoint for p in panels for r in p.records if r.timepoint}),
        "sexes": sorted({r.sex for p in panels for r in p.records if r.sex}),
        "contexts": contexts,
    }
    return StoryResponse(
        panels=panels, availability=availability,
        provenance={"root": "MoTrPAC Hackathon", "human_package_version": "2.0.8",
                    "human_collection": "c2.0", "method": "read-only source-output adapter"},
    )
