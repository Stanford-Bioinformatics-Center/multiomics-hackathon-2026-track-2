"""Render query_core outputs from a small synthetic disease/reference pair."""

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

STAGING = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STAGING))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from plot_results import bubble_size, feature_label, numeric, render, time_key  # noqa: E402
from query_core.engine import run_query  # noqa: E402


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


class PlotResultsTests(unittest.TestCase):
    def test_time_order_and_missing_q_size(self) -> None:
        self.assertLess(time_key("during_20_min", "", 1),
                        time_key("post_24_hr", "", 0))
        self.assertLess(time_key("post_15_30_45_min", "", 1),
                        time_key("post_3.5_4_hr", "", 0))
        self.assertLess(time_key("post_24_hr", "2", 0),
                        time_key("during_20_min", "3", 1))
        self.assertEqual(bubble_size(None), 30)
        self.assertLess(bubble_size(None), bubble_size(0.05))
        self.assertEqual(feature_label("protein", "Q16795", "NDUFA9", False),
                         "NDUFA9 (Q16795)")
        self.assertEqual(feature_label("rna", "NDUFA9", "NDUFA9", True),
                         "rna: NDUFA9")
        self.assertTrue(pd.isna(numeric(pd.Series(["", "1.2"]))[0]))

    def test_render_from_query_core_with_ambiguous_and_missing_q(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            disease = []
            for gene, fc, stat in (("A", "1", "4"), ("B", "-1", "-3"),
                                   ("C", "0.5", "2"), ("D", "0.2", "1")):
                disease.append({
                    "study_id": "StudyX", "species": "human", "tissue": "muscle",
                    "layer": "rna", "feature_id": gene, "id_namespace": "HGNC",
                    "log2_fc": fc, "statistic": stat, "p_value": "0.01",
                    "q_value": "0.04", "contrast": "case minus control",
                })
            reference = []
            for layer in ("rna", "protein"):
                for time in ("post_24_hr", "during_20_min"):
                    for gene, fc in (("A", "0.8"), ("B", "-0.4"),
                                     ("C", "0.2"), ("D", "-0.3")):
                        reference.append({
                            "study_id": "MoTrPAC", "species": "human", "tissue": "muscle",
                            "layer": layer, "feature_id": gene, "id_namespace": "HGNC",
                            "log2_fc": fc, "statistic": fc, "p_value": "0.01",
                            "q_value": "" if gene == "C" else "0.03",
                            "contrast": f"exercise vs control {time}", "timepoint": time,
                            "contrast_category": "EE-CON",
                        })
            # D has two protein candidates at 20 minutes, so it must not
            # appear as a bubble for that layer/time cell.
            reference.append({**reference[-1], "timepoint": "during_20_min",
                              "contrast": "exercise vs control during_20_min",
                              "feature_id": "D_extra", "gene_symbol": "D"})
            write_csv(root / "disease.csv", disease)
            write_csv(root / "reference.csv", reference)
            query_dir = root / "query"
            run_query(root / "disease.csv", root / "reference.csv", query_dir)
            manifest = render(query_dir, root / "plots", top=3)
            self.assertEqual(len(manifest["groups"]), 1)
            group = manifest["groups"][0]
            folder = root / "plots" / group["directory"]
            for name in ("bubble_heatmap.png", "bubble_heatmap.svg",
                         "bubble_plot_data.csv", "bubble_selection.csv.gz",
                         "rank_timecourse.png", "rank_plot_data.csv",
                         "coverage.png", "coverage_plot_data.csv",
                         "plot_metadata.json"):
                self.assertTrue((folder / name).exists(), name)
            bubbles = pd.read_csv(folder / "bubble_plot_data.csv")
            self.assertTrue((bubbles.mapping_status == "unique").all())
            self.assertNotIn("D", set(bubbles.disease_feature_id))  # top 3 only
            self.assertEqual(set(bubbles.disease_feature_id), {"A", "B", "C"})
            self.assertIn("rna\nduring_20_min", group["bubble"]["column_order"])
            self.assertLess(group["bubble"]["column_order"].index("rna\nduring_20_min"),
                            group["bubble"]["column_order"].index("rna\npost_24_hr"))
            self.assertEqual(group["rank"]["status"], "rendered")
            self.assertEqual(group["coverage"]["status"], "rendered")

    def test_pathway_plot_uses_statistic_and_preserves_methods(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            disease = [{
                "study_id": "DiseaseY", "species": "human", "tissue": "blood",
                "layer": "pathway", "feature_id": "BIOCARTA_TOB1_PATHWAY",
                "id_namespace": "BIOCARTA", "pathway_collection": "C2/BIOCARTA",
                "log2_fc": "", "statistic": "-2.4", "statistic_type": "PAGE_Z",
                "effect_scale": "PAGE pathway Z", "p_value": "0.01",
                "q_value": "0.03", "contrast": "case minus control",
            }]
            reference = []
            for time, score in (("during_20_min", "1.7"), ("post_24_hr", "-1.1")):
                reference.append({
                    "study_id": "MoTrPAC", "species": "human", "tissue": "blood",
                    "layer": "pathway", "feature_id": "BIOCARTA_TOB1_PATHWAY",
                    "id_namespace": "BIOCARTA", "pathway_collection": "C2/BIOCARTA",
                    "log2_fc": "", "statistic": score,
                    "statistic_type": "signed CAMERA z.std",
                    "effect_scale": "CAMERA pathway Z", "p_value": "0.02",
                    "q_value": "0.04", "contrast": f"exercise minus control {time}",
                    "timepoint": time, "contrast_category": "EE-CON",
                    "reference_scope": "partial BioCarta selection",
                })
            write_csv(root / "disease.csv", disease)
            write_csv(root / "reference.csv", reference)
            query_dir = root / "query"
            run_query(root / "disease.csv", root / "reference.csv", query_dir)
            manifest = render(query_dir, root / "plots")
            group = manifest["groups"][0]
            folder = root / "plots" / group["directory"]
            self.assertEqual(group["bubble"]["status"], "no_matches")
            self.assertEqual(group["pathway_bubble"]["status"], "rendered")
            self.assertEqual(group["pathway_bubble"]["color_field"],
                             "reference_statistic")
            self.assertEqual(group["pathway_bubble"]["color_label"],
                             "MoTrPAC CAMERA Z")
            self.assertEqual(group["pathway_bubble"]["disease_statistic_methods"],
                             ["PAGE_Z"])
            self.assertTrue((folder / "pathway_bubble_heatmap.png").exists())
            plotted = pd.read_csv(folder / "pathway_bubble_plot_data.csv")
            self.assertEqual(set(plotted.reference_plot_value), {1.7, -1.1})
            self.assertTrue(plotted.reference_log2_fc.isna().all())


if __name__ == "__main__":
    unittest.main()
