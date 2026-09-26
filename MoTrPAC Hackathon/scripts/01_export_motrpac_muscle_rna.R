#!/usr/bin/env Rscript

# Export the published MoTrPAC summary statistics used for the 24-hour
# skeletal-muscle RNA comparison. Run with:
#   Rscript scripts/01_export_motrpac_muscle_rna.R
# The script locates the project root from its own location.

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

# Load published model summaries with the package's documented interface.
# combine_with_featgene adds gene symbols and Ensembl IDs from the package's
# feature-to-gene mapping. No model is fitted here.
results <- load_differential_analysis(
  selected_tissues = "muscle",
  selected_omes = "transcript-rna-seq",
  single_matrix = TRUE,
  combine_with_featgene = TRUE,
  verbose = FALSE
)
results <- as.data.frame(results)
keep <- as.character(results$contrast_category) == "EE-CON" &
  as.character(results$Timepoint) == "post_24_hr"
results <- results[keep, , drop = FALSE]
if (nrow(results) != 16509L) stop("Unexpected number of 24-hour RNA features.")
if (anyDuplicated(results$feature_id)) stop("Duplicate RNA feature IDs.")
results$source_package_version <- package_version

columns <- c(
  "tissue", "assay", "contrast_type", "contrast_category", "contrast", "contrast_short",
  "Timepoint", "feature_id", "gene_symbol", "ensembl_gene",
  "logFC", "z.std", "p_value", "adj_p_value", "full_model",
  "source_package_version"
)
output <- results[, columns, drop = FALSE]
output_dir <- file.path(project_root, "data", "processed")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(output_dir, "motrpac_muscle_rna_ee_con_24h.csv.gz")
connection <- gzfile(output_path, open = "wt")
write.csv(output, connection, row.names = FALSE, na = "")
close(connection)

cat("Wrote", nrow(output), "MoTrPAC muscle RNA rows to", output_path, "\n")
