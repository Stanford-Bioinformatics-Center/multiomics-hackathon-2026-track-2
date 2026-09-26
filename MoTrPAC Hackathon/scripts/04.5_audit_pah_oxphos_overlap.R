#!/usr/bin/env Rscript

# Identify which genes corresponding to the PAH paper's nine lower muscle
# proteins belong to the GO biological-process OXPHOS set. This audit comes
# after the nine-gene lookup (script 04) and before the pathway result (05).
# Run from any directory with:
#   Rscript scripts/04.5_audit_pah_oxphos_overlap.R
#
# Inputs: outputs/PAH protein/pah_protein_motrpac_rna_matches.csv (script 04),
# and the published MotrpacHumanPreSuspensionAnalysis 0.2.4
# MOLECULAR_SIGNATURES object. Script 05 is not an input.
# This is an exact membership lookup. It does not perform enrichment testing
# and does not estimate an exercise effect in PAH patients.

script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_flag) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_flag), mustWork = TRUE)
project_root <- dirname(dirname(script_path))

package_name <- "MotrpacHumanPreSuspensionAnalysis"
if (!requireNamespace(package_name, quietly = TRUE)) {
  stop("Install MotrpacHumanPreSuspensionAnalysis before running this script.")
}
package_version <- as.character(utils::packageVersion(package_name))
if (package_version != "0.2.4") {
  stop("Expected MotrpacHumanPreSuspensionAnalysis 0.2.4; found ", package_version)
}
suppressPackageStartupMessages(library(MotrpacHumanPreSuspensionAnalysis))

matches_path <- file.path(project_root, "outputs", "PAH protein",
                          "pah_protein_motrpac_rna_matches.csv")
matches <- read.csv(matches_path, stringsAsFactors = FALSE)

if (nrow(matches) != 9L || anyDuplicated(matches$uniprot_accession) ||
    anyNA(matches$motrpac_current_symbol) ||
    any(matches$motrpac_rna_timepoint != "post_24_hr") ||
    any(matches$motrpac_package_version != package_version)) {
  stop("Expected nine distinct 24-hour RNA matches from package version 0.2.4.")
}

set_name <- "GOBP_OXIDATIVE_PHOSPHORYLATION"
# Membership in the source GO gene set is the question here. The pathway test
# exported by script 05 is a separate result for the measured genes in that
# set; its q value does not apply to any one member or to the nine PAH genes.
go_oxphos_genes <- MOLECULAR_SIGNATURES[["GOBP"]][[set_name]]
if (is.null(go_oxphos_genes) || anyDuplicated(go_oxphos_genes)) {
  stop("GO OXPHOS gene set is missing or has duplicate symbols.")
}
member <- matches$motrpac_current_symbol %in% go_oxphos_genes

output <- data.frame(
  paper_symbol = matches$paper_symbol,
  uniprot_accession = matches$uniprot_accession,
  pah_to_control_protein_ratio = matches$pah_to_control_protein_ratio,
  motrpac_current_symbol = matches$motrpac_current_symbol,
  go_oxphos_member = ifelse(member, "yes", "no"),
  motrpac_24h_rna_logFC = matches$motrpac_rna_logFC,
  motrpac_24h_rna_bh_q = matches$motrpac_rna_bh_adj_p_value,
  go_set_name = set_name,
  source_doi = matches$source_doi,
  motrpac_package_version = package_version,
  stringsAsFactors = FALSE
)

output_path <- file.path(project_root, "outputs", "PAH protein",
                         "pah_proteins_go_oxphos_membership.csv")
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.csv(output, output_path, row.names = FALSE, na = "")

cat("GO OXPHOS members among nine PAH-linked genes:", sum(member), "/9\n")
cat("Members:", paste(output$motrpac_current_symbol[member], collapse = ", "), "\n")
cat("GO set symbols:", length(go_oxphos_genes), "\n")
cat("Wrote", nrow(output), "membership rows to", output_path, "\n")
