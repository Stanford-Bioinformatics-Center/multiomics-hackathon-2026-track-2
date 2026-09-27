"""Behavioral checks for rank-sum enrichment and explicit universe gates."""

from __future__ import annotations

import csv
import gzip
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ranked_enrichment import run  # noqa: E402


def write_input(root: Path, *, layer: str, namespace: str, names: list[str]) -> tuple[Path, Path]:
    signature = root / "signature.csv"
    fields = ("study_id", "species", "tissue", "layer", "feature_id", "id_namespace",
              "log2_fc", "statistic", "p_value", "q_value", "contrast")
    with signature.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for i, name in enumerate(names):
            writer.writerow({
                "study_id": "DiseaseX", "species": "human", "tissue": "blood",
                "layer": layer, "feature_id": name, "id_namespace": namespace,
                "log2_fc": "", "statistic": i - 20,
                "p_value": "", "q_value": "", "contrast": "case minus control",
            })
    provenance = root / "universe.json"
    provenance.write_text(json.dumps({
        "source": "Fictional full assay table", "complete_measured_universe": True,
        "expected_feature_count": len(names), "inclusion_rule": "All assayed, QC-passing features",
        "feature_id_policy": "One identifier per measured feature",
        "gmt_source": "Fictional collection", "gmt_version": "test-1",
    }), encoding="utf-8")
    return signature, provenance


class RankedEnrichmentTests(unittest.TestCase):
    def test_rna_signed_test_and_bh_over_eligible_sets(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            names = [f"GENE{i:02d}" for i in range(40)]
            signature, provenance = write_input(root, layer="rna", namespace="HGNC",
                                                names=names)
            gmt = root / "sets.gmt"
            gmt.write_text(
                "HIGH\tTop genes\t" + "\t".join(names[30:]) + "\n" +
                "LOW\tBottom genes\t" + "\t".join(names[:10]) + "\n" +
                "TOO_SMALL\tNot tested\tGENE01\tGENE02\n",
                encoding="utf-8")
            manifest = run(signature, gmt, provenance, root / "out",
                           gmt_namespace="HGNC", pathway_namespace="Pathway",
                           collection="fictional/test")
            self.assertEqual(manifest["n_eligible_sets"], 2)
            with (root / "out" / "pathway_tests.csv").open(newline="") as stream:
                rows = {row["feature_id"]: row for row in csv.DictReader(stream)}
            self.assertGreater(float(rows["HIGH"]["statistic"]), 0)
            self.assertLess(float(rows["LOW"]["statistic"]), 0)
            self.assertLess(float(rows["HIGH"]["q_value"]), 0.05)
            self.assertEqual(rows["HIGH"]["log2_fc"], "")
            self.assertEqual(rows["TOO_SMALL"]["test_status"], "below_min_size")
            self.assertEqual(rows["TOO_SMALL"]["q_value"], "")
            with (root / "out" / "pathway_signature.csv").open(newline="") as stream:
                self.assertEqual(len(list(csv.DictReader(stream))), 2)

    def test_metabolite_refmet_name_gmt_and_universe_validation(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            names = [f"Metab {i:02d}" for i in range(40)]
            signature, provenance = write_input(root, layer="metabolite",
                                                namespace="RefMetName", names=names)
            gmt = root / "sets.gmt.gz"
            with gzip.open(gmt, "wt", encoding="utf-8") as stream:
                stream.write("METAB_SET\tTest names\t" + "\t".join(
                    name.upper() for name in names[30:]) + "\n")
            manifest = run(signature, gmt, provenance, root / "out",
                           gmt_namespace="RefMetName", pathway_namespace="Pathway",
                           collection="RefMet/test")
            self.assertEqual(manifest["n_eligible_sets"], 1)
            with (root / "out" / "pathway_signature.csv").open(newline="") as stream:
                row = next(csv.DictReader(stream))
            self.assertEqual(row["n_overlap"], "10")
            self.assertGreater(float(row["statistic"]), 0)
            with self.assertRaisesRegex(ValueError, "GMT member namespace"):
                run(signature, gmt, provenance, root / "invalid",
                    gmt_namespace="RefMet", pathway_namespace="Pathway",
                    collection="RefMet/test")
            data = json.loads(provenance.read_text(encoding="utf-8"))
            data["complete_measured_universe"] = False
            provenance.write_text(json.dumps(data), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "complete_measured_universe"):
                run(signature, gmt, provenance, root / "invalid",
                    gmt_namespace="RefMetName", pathway_namespace="Pathway",
                    collection="RefMet/test")

    def test_optional_r_camera_pr_smoke(self) -> None:
        rscript = shutil.which("Rscript")
        if rscript is None and sys.platform == "win32":
            candidates = sorted(Path("C:/Program Files/R").glob("R-*/bin/Rscript.exe"))
            rscript = str(candidates[-1]) if candidates else None
        if rscript is None:
            self.skipTest("Rscript is not installed")
        available = subprocess.run(
            [rscript, "-e", "quit(status=if(requireNamespace('limma',quietly=TRUE) && requireNamespace('jsonlite',quietly=TRUE) && requireNamespace('digest',quietly=TRUE)) 0 else 1)"],
            capture_output=True, text=True)
        if available.returncode:
            self.skipTest("Optional R dependencies are not installed")
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            names = [f"GENE{i:02d}" for i in range(40)]
            signature, provenance = write_input(root, layer="rna", namespace="HGNC",
                                                names=names)
            gmt = root / "sets.gmt"
            gmt.write_text(
                "HIGH\tTop genes\t" + "\t".join(names[30:]) + "\n" +
                "LOW\tBottom genes\t" + "\t".join(names[:10]) + "\n",
                encoding="utf-8")
            command = [
                rscript, str(Path(__file__).resolve().parents[1] / "camera_pr.R"),
                "--signature", str(signature), "--gmt", str(gmt),
                "--universe-provenance", str(provenance), "--gmt-namespace", "HGNC",
                "--pathway-namespace", "Pathway", "--collection", "fictional/test",
                "--out-dir", str(root / "out"),
            ]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with (root / "out" / "pathway_signature.csv").open(newline="") as stream:
                rows = {row["feature_id"]: row for row in csv.DictReader(stream)}
            self.assertGreater(float(rows["HIGH"]["statistic"]), 0)
            self.assertLess(float(rows["LOW"]["statistic"]), 0)
            self.assertEqual(rows["HIGH"]["log2_fc"], "")
            manifest = json.loads((root / "out" / "manifest.json").read_text())
            self.assertEqual(manifest["n_eligible_sets"], 2)
            self.assertEqual(manifest["inter_gene_correlation"], 0.01)


if __name__ == "__main__":
    unittest.main()
