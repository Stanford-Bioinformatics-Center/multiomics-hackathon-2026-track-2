"""Match Table 2 BioCarta labels to MoTrPAC's precomputed blood RNA tests."""

from pathlib import Path
import pandas as pd

MODULE = Path(__file__).resolve().parents[1]
PAPER = MODULE / "outputs" / "paper_table2_biocarta.csv"
CAMERA = MODULE / "outputs" / "motrpac_biocarta_camera.csv"
OUT = MODULE / "outputs"


def main() -> None:
    paper = pd.read_csv(PAPER)
    camera = pd.read_csv(CAMERA)
    expected_sets = {
        "BIOCARTA_CCR5_PATHWAY", "BIOCARTA_IL2_PATHWAY",
        "BIOCARTA_TCR_PATHWAY", "BIOCARTA_IL12_PATHWAY",
        "BIOCARTA_TOB1_PATHWAY", "BIOCARTA_ARENRF2_PATHWAY",
    }
    if set(paper["set"]) != expected_sets or paper["set"].duplicated().any():
        raise ValueError("Curated Table 2 BioCarta set list changed unexpectedly")
    if len(paper) != 6:
        raise ValueError("Expected all six BioCarta pathways from Table 2")
    if not (paper.filter(like="_page_z") < 0).all().all():
        raise ValueError("Expected negative published PAGE scores")
    keys = ["set", "modality", "timepoint"]
    if camera.duplicated(keys).any() or len(camera) != 60:
        raise ValueError("Expected one CAMERA result per set, modality and timepoint")
    mapped = camera.merge(paper, on="set", how="left", validate="many_to_one",
                          indicator=True)
    if (mapped["_merge"] != "both").any():
        raise ValueError("Some MoTrPAC sets did not match Table 2")
    mapped = mapped.drop(columns="_merge")
    mapped["ipah_vs_ee_re_direction"] = mapped["z.std"].map(
        lambda z: "opposite" if z > 0 else ("same" if z < 0 else "zero")
    )
    # The paper score and CAMERA score come from different cohorts and methods.
    # Direction comparison is descriptive, not a joint significance test.
    mapped.to_csv(OUT / "paper_motrpac_biocarta_mapped.csv", index=False)
    qc = (mapped.groupby(["set", "paper_label"],
                         as_index=False)
          .agg(contrasts=("z.std", "size"),
               ee_contrasts=("modality", lambda v: int((v == "EE").sum())),
               re_contrasts=("modality", lambda v: int((v == "RE").sum())),
               bh_below_0_05=("adj_p_value", lambda v: int((v < 0.05).sum())),
               min_camera_z=("z.std", "min"),
               max_camera_z=("z.std", "max")))
    if not ((qc["contrasts"] == 10) & (qc["ee_contrasts"] == 6) &
            (qc["re_contrasts"] == 4)).all():
        raise ValueError("Incomplete set/timepoint coverage")
    qc.to_csv(OUT / "pathway_mapping_qc.csv", index=False)
    print(qc.to_string(index=False))


if __name__ == "__main__":
    main()
