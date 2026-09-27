"""Behavioral tests for matching, uncertainty, and legacy normalization."""

from __future__ import annotations

import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from query_core.adapter_legacy import convert
from query_core.engine import run_query
from query_core.merge_reference import merge
from query_core.schema import DISEASE_COLUMNS, REFERENCE_COLUMNS


def write_rows(path: Path, rows: list[dict]) -> None:
    fields = list(dict.fromkeys(key for row in rows for key in row))
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_rows(path: Path) -> list[dict]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def disease(feature: str, layer="rna", **extra):
    base = dict.fromkeys(DISEASE_COLUMNS, "")
    base.update(study_id="Disease1", species="human", tissue="blood",
                layer=layer, feature_id=feature, id_namespace="HGNC",
                log2_fc="1", statistic="2", p_value="0.01",
                q_value="0.02", contrast="case - control")
    base.update(extra)
    return base


def reference(feature: str, layer="protein", **extra):
    base = dict.fromkeys(REFERENCE_COLUMNS, "")
    base.update(study_id="MoTrPAC", species="human", tissue="blood",
                layer=layer, feature_id=feature, id_namespace="HGNC",
                log2_fc="1", statistic="2", p_value="0.01",
                q_value="0.02", contrast="exercise - control",
                timepoint="20m", contrast_category="EE-CON")
    base.update(extra)
    return base


class QueryCoreTests(unittest.TestCase):
    def test_matches_audit_rank_and_pathway_scale(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            ds = [disease("A", gene_symbol="A", biospecimen="PBMC"),
                  disease("B", gene_symbol="B", log2_fc="2"),
                  disease("C", gene_symbol="C", log2_fc="-1"),
                  disease("X", gene_symbol="X"),
                  disease("RM001", layer="metabolite", id_namespace="RefMet",
                          refmet_id="RM001", refmet_name="citrate"),
                  disease("BIOCARTA_TOB1_PATHWAY", layer="pathway",
                          id_namespace="BIOCARTA", pathway_collection="C2/BIOCARTA",
                          effect_scale="PAGE Z", statistic_type="signed PAGE Z",
                          log2_fc="", statistic="-2"),
                  disease("BIOCARTA_OTHER_PATHWAY", layer="pathway",
                          id_namespace="BIOCARTA", pathway_collection="C2/OTHER",
                          log2_fc="", statistic="-2")]
            rs = [reference("A", gene_symbol="A", biospecimen="whole blood"),
                  reference("B", gene_symbol="B", log2_fc="-1"),
                  reference("C", gene_symbol="C", log2_fc="2"),
                  reference("A:site1", layer="phosphosite", id_namespace="MoTrPAC",
                            gene_symbol="A"),
                  reference("A:site2", layer="phosphosite", id_namespace="MoTrPAC",
                            gene_symbol="A"),
                  reference("RM001", layer="metabolite", id_namespace="RefMet",
                            refmet_id="RM001", refmet_name="citrate"),
                  reference("RM999", layer="metabolite", id_namespace="RefMet",
                            refmet_id="RM999", refmet_name="citrate"),
                  reference("BIOCARTA_TOB1_PATHWAY", layer="pathway",
                            id_namespace="BIOCARTA", pathway_collection="C2/BIOCARTA",
                            effect_scale="CAMERA Z", statistic_type="signed CAMERA Z",
                            log2_fc="", statistic="3"),
                  reference("BIOCARTA_TOB1_PATHWAY", layer="pathway",
                            id_namespace="BIOCARTA", pathway_collection="C2/OTHER",
                            log2_fc="", statistic="4")]
            disease_csv, reference_csv = root / "d.csv", root / "r.csv"
            write_rows(disease_csv, ds)
            write_rows(reference_csv, rs)
            out = root / "out"
            summary = run_query(disease_csv, reference_csv, out,
                                rank_field="log2_fc")
            matches = read_rows(out / "matches.csv.gz")
            a_protein = [m for m in matches if m["disease_feature_id"] == "A" and
                         m["reference_layer"] == "protein"]
            self.assertEqual(len(a_protein), 1)
            self.assertEqual(a_protein[0]["biospecimen_compatibility"], "different")
            self.assertEqual(a_protein[0]["supported_direction"], "concordant")
            self.assertEqual(a_protein[0]["match_relation"], "gene_linked_cross_layer")
            sites = [m for m in matches if m["reference_layer"] == "phosphosite"]
            self.assertEqual(len(sites), 2)
            self.assertTrue(all(m["mapping_status"] == "ambiguous" for m in sites))
            metab = [m for m in matches if m["disease_feature_id"] == "RM001"]
            self.assertEqual(len(metab), 1)
            self.assertEqual(metab[0]["matched_on"], "refmet_id")
            pathway = [m for m in matches if m["disease_feature_id"] ==
                       "BIOCARTA_TOB1_PATHWAY"]
            self.assertEqual(len(pathway), 1)
            self.assertEqual(pathway[0]["sign_basis"], "signed_statistic")
            self.assertEqual(pathway[0]["supported_direction"], "discordant")
            self.assertEqual(pathway[0]["match_relation"], "pathway_name_overlap")
            self.assertEqual(pathway[0]["disease_log2_fc"], "")
            self.assertEqual({m["feature_id"] for m in read_rows(out / "unmatched.csv")},
                             {"X", "BIOCARTA_OTHER_PATHWAY"})
            rank = [r for r in read_rows(out / "rank_correlation.csv")
                    if r["reference_layer"] == "protein"][0]
            self.assertEqual(rank["n_used"], "3")
            self.assertAlmostEqual(float(rank["spearman_rho"]), -1.0)
            self.assertEqual(summary["counts"]["different_biospecimen_candidate_rows"], 1)

    def test_missing_q_is_not_called_no_change(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            d = disease("A", q_value="", source_q_value="0.01",
                        source_q_scope="probe-level BH")
            write_rows(root / "d.csv", [d])
            write_rows(root / "r.csv", [reference("A")])
            run_query(root / "d.csv", root / "r.csv", root / "out")
            match = read_rows(root / "out/matches.csv.gz")[0]
            self.assertEqual(match["evidence_status"], "missing_q")
            self.assertEqual(match["supported_direction"], "")
            self.assertEqual(match["source_q_value"] if "source_q_value" in match
                             else match["disease_source_q_value"], "0.01")

    def test_cross_species_requires_explicit_one_to_one_map(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            write_rows(root / "d.csv", [disease("CCN1", gene_symbol="CCN1")])
            write_rows(root / "r.csv", [reference("Ccn1", species="rat",
                                                      id_namespace="RatGeneSymbol",
                                                      gene_symbol="Ccn1")])
            run_query(root / "d.csv", root / "r.csv", root / "without")
            self.assertEqual(len(read_rows(root / "without/matches.csv.gz")), 0)
            map_path = root / "ortholog.csv"
            write_rows(map_path, [dict(source_species="human", source_gene_symbol="CCN1",
                                       target_species="rat", target_gene_symbol="Ccn1")])
            mapped_summary = run_query(root / "d.csv", root / "r.csv", root / "with",
                                       ortholog_map=map_path)
            match = read_rows(root / "with/matches.csv.gz")[0]
            self.assertEqual((match["species"], match["reference_species"]),
                             ("human", "rat"))
            self.assertEqual(match["match_confidence"], "ortholog")
            self.assertIn("supplied explicit one-to-one ortholog map",
                          mapped_summary["interpretation"])

    def test_protein_gene_match_is_not_called_same_analyte(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            write_rows(root / "d.csv", [disease("P11111", layer="protein",
                                                  id_namespace="UniProt", gene_symbol="GENE",
                                                  uniprot="P11111")])
            write_rows(root / "r.csv", [reference("P22222", layer="protein",
                                                   id_namespace="UniProt", gene_symbol="GENE",
                                                   uniprot="P22222")])
            run_query(root / "d.csv", root / "r.csv", root / "out")
            match = read_rows(root / "out/matches.csv.gz")[0]
            self.assertEqual(match["matched_on"], "gene_symbol")
            self.assertEqual(match["match_relation"], "gene_linked_same_layer")

    def test_mixed_known_reference_release_is_rejected_per_study_species(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            write_rows(root / "d.csv", [disease("A")])
            human = reference("A", source_package_version="2.0.8",
                              source_collection="c2.0")
            changed = reference("B", source_package_version="2.0.9",
                                source_collection="c2.0")
            write_rows(root / "mixed.csv", [human, changed])
            with self.assertRaisesRegex(ValueError, "Mixed known reference releases"):
                run_query(root / "d.csv", root / "mixed.csv", root / "rejected")
            rat = reference("C", species="rat", id_namespace="RatGeneSymbol",
                            source_package_version="2.0.0", source_collection="rat")
            unknown = reference("D", source_package_version="", source_collection="")
            write_rows(root / "separate.csv", [human, rat, unknown])
            summary = run_query(root / "d.csv", root / "separate.csv", root / "accepted")
            self.assertEqual(len(summary["reference_releases"]), 2)
            human_audit = next(x for x in summary["reference_releases"]
                               if x["species"] == "human")
            self.assertEqual(human_audit["package_versions"], ["2.0.8"])
            self.assertEqual(human_audit["rows_without_package_version"], 1)

    def test_legacy_adapter_preserves_raw_feature_and_scale(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "data/processed/motrpac_muscle_rna_ee_con_24h.csv.gz"
            source.parent.mkdir(parents=True)
            write_rows(source, [dict(tissue="muscle", assay="transcript-rna-seq",
                                     contrast_category="EE-CON", contrast_short="EE - CON",
                                     Timepoint="24h", feature_id="ENSG1", gene_symbol="GENE1",
                                     logFC="0.5", **{"z.std": "2.3"}, p_value="0.02",
                                     adj_p_value="0.03", source_package_version="2.0.8",
                                     source_collection="c2.0")])
            output = root / "reference.csv.gz"
            counts = convert(root, output, skip_missing=True)
            self.assertEqual(counts["muscle_rna"], 1)
            row = read_rows(output)[0]
            self.assertEqual((row["feature_id"], row["id_namespace"],
                              row["source_feature_id"]), ("GENE1", "HGNC", "ENSG1"))
            self.assertEqual(row["statistic"], "2.3")

    def test_opt_in_metabolite_feature_name_fallback_is_audited(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "Metabolomics/data/processed/motrpac_metabolite_contrasts.csv.gz"
            source.parent.mkdir(parents=True)
            write_rows(source, [dict(tissue="blood", platform="metab-u", feature_id="PC 40:4",
                                     refmet_name="NA", refmet_id="NA", contrast_category="EE-CON",
                                     Timepoint="20m", logFC="0.4", AveExpr="3",
                                     p_value="0.02", adj_p_value="0.03",
                                     source_package_version="2.0.8", source_collection="c2.0")])
            without, with_fallback = root / "without.csv.gz", root / "with.csv.gz"
            convert(root, without, skip_missing=True)
            convert(root, with_fallback, skip_missing=True,
                    metabolite_feature_name_fallback=True)
            self.assertEqual(read_rows(without)[0]["id_namespace"], "MoTrPAC")
            row = read_rows(with_fallback)[0]
            self.assertEqual(row["id_namespace"], "RefMetName")
            self.assertEqual(row["match_basis"], "feature_id_name_fallback")

    def test_merge_replaces_duplicates_and_rejects_mixed_release(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            base = reference("A", layer="rna", tissue="muscle",
                             source_feature_id="ENSG1", assay="transcript-rna-seq",
                             source_package_version="2.0.8", source_collection="c2.0")
            newer = {**base, "log2_fc": "2"}
            write_rows(root / "base.csv", [base])
            write_rows(root / "full.csv", [newer])
            result = merge([root / "base.csv", root / "full.csv"],
                           root / "merged.csv", root / "audit.csv",
                           expected_package_version="2.0.8",
                           expected_collection="c2.0")
            self.assertEqual(result["duplicates_replaced"], 1)
            self.assertEqual(result["duplicates_with_different_numeric_values"], 1)
            self.assertEqual(read_rows(root / "merged.csv")[0]["log2_fc"], "2")
            write_rows(root / "wrong.csv", [{**newer, "source_collection": "c1.3"}])
            with self.assertRaisesRegex(ValueError, "does not match"):
                merge([root / "base.csv", root / "wrong.csv"],
                      root / "bad.csv", root / "bad_audit.csv",
                      expected_package_version="2.0.8",
                      expected_collection="c2.0")


if __name__ == "__main__":
    unittest.main()
