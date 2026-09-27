"""App-runner smoke tests with a new, non-PAH study ID."""

import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from run_config import load_jobs, run_config


FIELDS = ["study_id", "species", "tissue", "layer", "feature_id",
          "id_namespace", "log2_fc", "statistic", "p_value", "q_value",
          "contrast", "timepoint", "contrast_category", "gene_symbol"]


def write_rows(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


class RunConfigTest(unittest.TestCase):
    def test_config_relative_paths_and_query_outputs(self):
        with TemporaryDirectory() as root:
            base = Path(root)
            disease = base / "study.csv"
            reference = base / "reference.csv"
            row = dict(study_id="UNRELATED_STUDY", species="human", tissue="muscle",
                       layer="rna", feature_id="STAT3", id_namespace="HGNC",
                       log2_fc="0.5", statistic="2.0", p_value="0.01",
                       q_value="0.02", contrast="disease - comparator",
                       timepoint="", contrast_category="", gene_symbol="STAT3")
            write_rows(disease, [row])
            write_rows(reference, [{**row, "study_id": "REFERENCE",
                                    "contrast": "exercise - control",
                                    "timepoint": "post_24_hr",
                                    "contrast_category": "EE-CON"}])
            config = base / "jobs.json"
            config.write_text(json.dumps({"schema_version": 1, "jobs": [{
                "id": "unrelated_study", "disease": "study.csv",
                "reference": "reference.csv", "plots": False,
                "filters": {"tissue": "muscle", "reference_contrast_category": "EE-CON"},
            }]}), encoding="utf-8")
            manifest = run_config(config, base / "results")
            self.assertEqual(manifest["jobs"][0]["status"], "complete")
            self.assertEqual(manifest["jobs"][0]["counts"]["matched_candidate_rows"], 1)
            self.assertTrue((base / "results" / "unrelated_study" / "matches.csv.gz").is_file())
            self.assertTrue((base / "results" / "manifest.json").is_file())

    def test_rejects_unsafe_or_duplicate_job_ids(self):
        with TemporaryDirectory() as root:
            base = Path(root)
            config = base / "jobs.json"
            config.write_text(json.dumps({"schema_version": 1, "jobs": [
                {"id": "../escape", "disease": "absent.csv", "reference": "absent.csv"}
            ]}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "safe"):
                load_jobs(config)

    def test_rejects_boolean_plot_top_and_bad_threshold(self):
        with TemporaryDirectory() as root:
            base = Path(root)
            config = base / "jobs.json"
            job = {"id": "valid_id", "disease": "unused.csv",
                   "reference": "unused.csv", "plots": {"top": True}}
            config.write_text(json.dumps({"schema_version": 1, "jobs": [job]}),
                              encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "positive integer top"):
                load_jobs(config)
            job["plots"] = False
            job["filters"] = {"q_threshold": True}
            config.write_text(json.dumps({"schema_version": 1, "jobs": [job]}),
                              encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "q_threshold"):
                load_jobs(config)

    def test_failed_job_does_not_hide_completed_job(self):
        with TemporaryDirectory() as root:
            base = Path(root)
            disease, reference = base / "study.csv", base / "reference.csv"
            row = dict(study_id="UNRELATED_STUDY", species="human", tissue="muscle",
                       layer="rna", feature_id="STAT3", id_namespace="HGNC",
                       log2_fc="0.5", statistic="2", p_value="0.01",
                       q_value="0.02", contrast="disease - control",
                       timepoint="", contrast_category="", gene_symbol="STAT3")
            write_rows(disease, [row])
            write_rows(reference, [{**row, "timepoint": "post_24_hr",
                                    "contrast_category": "EE-CON"}])
            bad = base / "bad.csv"
            bad.write_text("wrong,header\n1,2\n", encoding="utf-8")
            config = base / "jobs.json"
            config.write_text(json.dumps({"schema_version": 1, "jobs": [
                {"id": "good", "disease": "study.csv", "reference": "reference.csv"},
                {"id": "bad", "disease": "study.csv", "reference": "bad.csv"},
                {"id": "missing", "disease": "study.csv", "reference": "not_uploaded.csv"},
            ]}), encoding="utf-8")
            manifest = run_config(config, base / "results")
            self.assertEqual([job["status"] for job in manifest["jobs"]],
                             ["complete", "error", "error"])
            self.assertTrue((base / "results/good/matches.csv.gz").is_file())
            self.assertIn("missing required columns", manifest["jobs"][1]["error"])
            self.assertIn("FileNotFoundError", manifest["jobs"][2]["error"])


if __name__ == "__main__":
    unittest.main()
