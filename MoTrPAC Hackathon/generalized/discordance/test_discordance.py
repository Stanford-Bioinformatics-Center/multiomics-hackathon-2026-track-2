"""Synthetic end-to-end checks; no patient data or MoTrPAC package required."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from build_discordance import read_effects


SCRIPT = Path(__file__).with_name("build_discordance.py")
TRANSFER = Path(__file__).with_name("transfer_model.py")
PTM_AUDIT = Path(__file__).with_name("audit_ptm_parent.py")
TIMES = ("post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")


def make_synthetic() -> pd.DataFrame:
    rng = np.random.default_rng(17)
    records = []
    for i in range(45):
        gene = f"GENE{i:03d}"
        accession = f"P{i:05d}"
        early_rna = rng.normal(0, 0.5)
        mid_rna = early_rna + rng.normal(0, 0.15)
        early_ptm = rng.normal(0, 0.35)
        target_protein = 0.55 * mid_rna + 0.3 * early_ptm + rng.normal(0, 0.12)
        target_rna = rng.normal(0, 0.35)
        if i == 0:
            target_rna, target_protein = 0.8, -0.6
        if i == 1:
            target_rna, target_protein = 0.8, 0.02
        if i == 2:
            target_rna, target_protein = 0.8, 0.01
        for j, time in enumerate(TIMES):
            rna_fc = (early_rna, mid_rna, target_rna)[j]
            protein_fc = (0.01, 0.1 * target_protein, target_protein)[j]
            for layer, feature, effect, ave in (
                ("rna", f"rna_{gene}", rna_fc, 10.0),
                ("protein", f"protein_{gene}", protein_fc, 9.0),
                ("phosphosite", f"site_{gene}", protein_fc + early_ptm, 5.0),
            ):
                radius = 0.06 if not (i == 2 and j == 2 and layer == "protein") else 0.5
                q = 0.001 if abs(effect) > 0.2 and radius < 0.2 else 0.8
                records.append({
                    "species": "human", "tissue": "muscle", "layer": layer,
                    "contrast_type": "exercise_with_controls", "contrast_category": "EE-CON",
                    "contrast_short": f"Endur.{time} - Control.{time}", "timepoint": time,
                    "feature_id": feature, "gene_symbol": gene, "uniprot": accession if layer != "rna" else "",
                    "log2_fc": effect, "ci_lower": effect - radius, "ci_upper": effect + radius,
                    "q_value": q, "ave_expr": ave,
                })
        if i == 3:
            for time in TIMES:
                records.append({
                    "species": "human", "tissue": "muscle", "layer": "protein",
                    "contrast_type": "exercise_with_controls", "contrast_category": "EE-CON",
                    "contrast_short": f"Endur.{time} - Control.{time}", "timepoint": time,
                    "feature_id": "protein_GENE003_lower_abundance", "gene_symbol": gene,
                    "uniprot": accession, "log2_fc": 3.0, "ci_lower": 2.9,
                    "ci_upper": 3.1, "q_value": 0.001, "ave_expr": 1.0,
                })
        if i == 4:
            records.append({
                "species": "human", "tissue": "muscle", "layer": "protein",
                "contrast_type": "exercise_with_controls", "contrast_category": "EE-CON",
                "contrast_short": f"Endur.{TIMES[-1]} - Control.{TIMES[-1]}", "timepoint": TIMES[-1],
                "feature_id": "protein_GENE004_second_accession", "gene_symbol": gene,
                "uniprot": "P99999", "log2_fc": 0.3, "ci_lower": 0.2,
                "ci_upper": 0.4, "q_value": 0.001, "ave_expr": 5.0,
            })
    return pd.DataFrame(records)


class DiscordanceTests(unittest.TestCase):
    def test_native_aliases(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "native.csv"
            pd.DataFrame([{
                "tissue": "muscle", "assay": "prot-pr", "contrast_category": "EE-CON",
                "Timepoint": TIMES[-1], "feature_id": "P1", "gene_symbol": "G1",
                "uniprot": "P1", "logFC": 0.2, "CI.L_calculated": 0.1, "CI.R_calculated": 0.3,
                "adj_p_value": 0.01,
            }]).to_csv(path, index=False)
            frame = read_effects(path)
            self.assertEqual(frame.loc[0, "layer"], "protein")
            self.assertEqual(frame.loc[0, "timepoint"], TIMES[-1])
            self.assertEqual(frame.loc[0, "log2_fc"], 0.2)
            self.assertEqual(frame.loc[0, "ci_lower"], 0.1)

    def test_catalog_and_grouped_model(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "effects.csv"
            output = root / "out"
            make_synthetic().to_csv(source, index=False)
            command = [
                sys.executable, str(SCRIPT), "--input", str(source), "--out-dir", str(output),
                "--tissue", "muscle", "--contrast-category", "EE-CON",
                "--target-time", TIMES[-1], "--earlier-times", TIMES[0], TIMES[1],
            ]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            catalog = pd.read_csv(output / "catalog.csv")
            predictions = pd.read_csv(output / "model_predictions.csv")
            metrics = pd.read_csv(output / "model_metrics.csv")
            summary = json.loads((output / "run_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(catalog.loc[catalog.gene_symbol == "GENE000", "event_class"].iloc[0], "supported_opposite")
            self.assertEqual(catalog.loc[catalog.gene_symbol == "GENE001", "event_class"].iloc[0], "rna_response_protein_equivalent")
            self.assertEqual(catalog.loc[catalog.gene_symbol == "GENE002", "event_class"].iloc[0], "indeterminate")
            selected = catalog.loc[catalog.gene_symbol == "GENE003"].iloc[0]
            self.assertEqual(selected["protein_feature_id"], "protein_GENE003")
            self.assertEqual(selected["protein_feature_count"], 2)
            self.assertEqual(len(catalog.loc[catalog.gene_symbol == "GENE004"]), 2)
            self.assertFalse(predictions.gene_symbol.eq("GENE004").any())
            self.assertEqual(summary["model"]["model_status"], "fitted")
            self.assertEqual(set(metrics.model), {"zero", "rna_only", "temporal"})
            self.assertTrue(predictions["prediction_temporal"].notna().all())
            np.testing.assert_allclose(
                predictions["observed_discordance_residual"],
                predictions["rna_log2_fc"] - predictions["protein_log2_fc"],
            )
            self.assertTrue((output / "discordance_overview.png").is_file())

            ptm_output = root / "ptm"
            result = subprocess.run([
                sys.executable, str(PTM_AUDIT), "--input", str(source),
                "--out-dir", str(ptm_output), "--tissue", "muscle",
                "--contrast-category", "EE-CON", "--times", *TIMES,
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            ptm_summary = pd.read_csv(ptm_output / "ptm_parent_summary.csv")
            self.assertEqual(len(ptm_summary), 3)
            self.assertGreater(ptm_summary.loc[0, "supported_phosphosite_and_parent_equivalent_rows"], 0)
            self.assertTrue((ptm_output / "ptm_parent_overview.png").is_file())

            other = make_synthetic()
            other["contrast_category"] = "RE-CON"
            other["contrast_short"] = other["contrast_short"].str.replace("Endur.", "Resist.", regex=False)
            other_source = root / "other_effects.csv"
            other_output = root / "other_out"
            other.to_csv(other_source, index=False)
            other_command = [
                sys.executable, str(SCRIPT), "--input", str(other_source), "--out-dir", str(other_output),
                "--tissue", "muscle", "--contrast-category", "RE-CON",
                "--target-time", TIMES[-1], "--earlier-times", TIMES[0], TIMES[1],
            ]
            result = subprocess.run(other_command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            transfer_output = root / "transfer"
            result = subprocess.run([
                sys.executable, str(TRANSFER), "--train-dir", str(output),
                "--test-dir", str(other_output), "--out-dir", str(transfer_output),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            transfer_metrics = pd.read_csv(transfer_output / "transfer_metrics.csv")
            self.assertEqual(set(transfer_metrics.model), {"zero", "rna_only", "temporal"})


if __name__ == "__main__":
    unittest.main()
