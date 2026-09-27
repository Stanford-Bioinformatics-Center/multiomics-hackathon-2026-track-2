"""Synthetic Workbench-shaped inputs exercise parsing and study-design gates."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from workbench_adapter import run, split_output  # noqa: E402


def write_json(path: Path, value: object) -> Path:
    path.write_text(json.dumps(value), encoding="utf-8")
    return path


class WorkbenchAdapterTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.sample_ids = [f"S{i}" for i in range(1, 7)]
        self.factors = write_json(self.root / "factors.json", {
            str(i): {"study_id": "ST123456", "local_sample_id": sid,
                     "factors": f"Group: {'Case' if i <= 3 else 'Control'} | Time: Rest"}
            for i, sid in enumerate(self.sample_ids, 1)
        })
        self.data = write_json(self.root / "data.json", {
            "1": {"study_id": "ST123456", "analysis_id": "AN1", "metabolite_id": "ME1",
                  "metabolite_name": "Glucose", "refmet_name": "Glucose", "units": "peak area",
                  "DATA": dict(zip(self.sample_ids, [8, 16, 32, 2, 4, 8]))},
            "2": {"study_id": "ST123456", "analysis_id": "AN1", "metabolite_id": "ME2",
                  "metabolite_name": "unknown feature", "refmet_name": "", "units": "peak area",
                  "DATA": dict(zip(self.sample_ids, [2, 3, 4, 2, 3, 4]))},
        })
        self.design = {
            "study_id": "ST123456", "species": "human", "tissue": "plasma",
            "contrast": "Case minus Control at rest",
            "design": {"type": "unpaired", "group_factor": "Group",
                       "case_values": ["Case"], "control_values": ["Control"],
                       "time_factor": "Time", "time_value": "Rest"},
            "normalization": {"input_scale": "linear", "pseudocount": 0,
                              "sample_center": "none"},
            "min_samples_per_arm": 3, "min_fraction_per_arm": 1,
        }

    def tearDown(self) -> None:
        self.temp.cleanup()

    def test_unpaired_signed_log2_effect_and_unmapped_feature(self) -> None:
        design_path = write_json(self.root / "design.json", self.design)
        result = run(self.factors, self.data, design_path)
        glucose = result.set_index("source_feature_id").loc["AN1:ME1"]
        self.assertEqual(glucose.feature_id, "Glucose")
        self.assertEqual(glucose.id_namespace, "RefMetName")
        self.assertAlmostEqual(glucose.log2_fc, 2.0)
        self.assertGreater(glucose.statistic, 0)
        self.assertTrue(0 <= glucose.q_value <= 1)
        unknown = result.set_index("source_feature_id").loc["AN1:ME2"]
        self.assertTrue(unknown.status == "tested")
        self.assertTrue(unknown.id_namespace is None)

    def test_paired_requires_real_subject_ids(self) -> None:
        design = dict(self.design)
        design["design"] = {"type": "paired_change", "group_factor": "Group",
                            "case_values": ["Case"], "time_factor": "Time",
                            "baseline_value": "Rest", "endpoint_value": "Peak"}
        design_path = write_json(self.root / "paired.json", design)
        with self.assertRaisesRegex(ValueError, "subject_id_factor or subject_map_csv"):
            run(self.factors, self.data, design_path)

    def test_paired_uses_subject_factor_not_sample_order(self) -> None:
        pairs = [
            ("A_peak", "Case", "Peak", "A", 8),
            ("B_rest", "Case", "Rest", "B", 3),
            ("A_rest", "Case", "Rest", "A", 4),
            ("C_peak", "Case", "Peak", "C", 20),
            ("B_peak", "Case", "Peak", "B", 12),
            ("C_rest", "Case", "Rest", "C", 5),
            ("D_rest", "Case", "Rest", "D", 7),
        ]
        factors = write_json(self.root / "paired_factors.json", {
            str(i): {"study_id": "ST123456", "local_sample_id": sid,
                     "factors": f"Group: {group} | Time: {time} | Subject ID: {person}"}
            for i, (sid, group, time, person, _) in enumerate(pairs)
        })
        data = write_json(self.root / "paired_data.json", {
            "1": {"study_id": "ST123456", "analysis_id": "AN1", "metabolite_id": "ME1",
                  "metabolite_name": "Glucose", "refmet_name": "Glucose", "units": "log2 intensity",
                  "DATA": {sid: value for sid, _, _, _, value in pairs}},
        })
        design = dict(self.design)
        design["design"] = {"type": "paired_change", "group_factor": "Group",
                            "case_values": ["Case"], "time_factor": "Time",
                            "baseline_value": "Rest", "endpoint_value": "Peak",
                            "subject_id_factor": "Subject ID"}
        design["normalization"] = {"input_scale": "log2"}
        design_path = write_json(self.root / "paired.json", design)
        result = run(factors, data, design_path)
        self.assertAlmostEqual(result.iloc[0].log2_fc, (4 + 9 + 15) / 3)
        self.assertEqual(result.iloc[0].n_case, 3)
        self.assertEqual(result.iloc[0].n_design_case, 3)

    def test_unsupported_design_is_rejected(self) -> None:
        design = dict(self.design)
        design["design"] = {"type": "guess_pairs", "group_factor": "Group",
                            "case_values": ["Case"]}
        design_path = write_json(self.root / "bad.json", design)
        with self.assertRaisesRegex(ValueError, "Unsupported design"):
            run(self.factors, self.data, design_path)

    def test_difference_in_changes_uses_both_identified_groups(self) -> None:
        samples = []
        values = {}
        for arm, changes in (("Case", [2, 3, 4]), ("Control", [0, 1, 2])):
            for i, delta in enumerate(changes, 1):
                person = f"{arm}_{i}"
                for time, value in (("Rest", 5 + i), ("Peak", 5 + i + delta)):
                    sid = f"{person}_{time}"
                    samples.append({"study_id": "ST123456", "local_sample_id": sid,
                                    "factors": f"Group: {arm} | Time: {time} | Person: {person}"})
                    values[sid] = value
        factors = write_json(self.root / "did_factors.json", samples[::-1])
        data = write_json(self.root / "did_data.json", [
            {"study_id": "ST123456", "analysis_id": "AN1", "metabolite_id": "ME1",
             "metabolite_name": "Glucose", "refmet_name": "Glucose", "refmet_id": "RM1",
             "units": "log2 intensity", "DATA": values}
        ])
        design = dict(self.design)
        design["design"] = {"type": "difference_in_changes", "group_factor": "Group",
                            "case_values": ["Case"], "control_values": ["Control"],
                            "time_factor": "Time", "baseline_value": "Rest",
                            "endpoint_value": "Peak", "subject_id_factor": "Person"}
        design["normalization"] = {"input_scale": "log2"}
        design_path = write_json(self.root / "did.json", design)
        result = run(factors, data, design_path)
        self.assertAlmostEqual(result.iloc[0].log2_fc, 2)
        self.assertEqual(result.iloc[0].feature_id, "RM1")
        self.assertEqual(result.iloc[0].id_namespace, "RefMet")

    def test_duplicate_policy_uses_coverage_not_p_value(self) -> None:
        result = pd.DataFrame({
            "id_namespace": ["RefMetName", "RefMetName"],
            "feature_id": ["Glucose", "Glucose"],
            "source_feature_id": ["AN2:ME2", "AN1:ME1"],
            "n_case": [3, 2], "n_control": [3, 2],
            "p_value": [0.9, 0.001], "duplicate_canonical_id": [True, True],
            "status": ["tested", "tested"],
        })
        with self.assertRaisesRegex(ValueError, "duplicate canonical IDs"):
            split_output(result, "error")
        chosen, _, audit = split_output(result, "highest_coverage")
        self.assertEqual(chosen.iloc[0].source_feature_id, "AN2:ME2")
        self.assertEqual(int(audit.chosen_for_signature.sum()), 1)


if __name__ == "__main__":
    unittest.main()
