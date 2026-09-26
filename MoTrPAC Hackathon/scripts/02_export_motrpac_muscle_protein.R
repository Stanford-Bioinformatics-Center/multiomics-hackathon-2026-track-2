#!/usr/bin/env Rscript

# Export published MoTrPAC skeletal-muscle protein summary statistics at
# the three post-exercise sampling times. Run with:
#   Rscript scripts/02_export_motrpac_muscle_protein.R

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

# These are precomputed model results, not individual protein measurements.
# The documented loader also adds the package's UniProt and gene annotations.
results <- load_differential_analysis(
  selected_tissues = "muscle",
  selected_omes = "prot-pr",
  single_matrix = TRUE,
  combine_with_featgene = TRUE,
  verbose = FALSE
)
results <- as.data.frame(results)
times <- c("post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")
keep <- as.character(results$contrast_category) == "EE-CON" &
  as.character(results$Timepoint) %in% times
results <- results[keep, , drop = FALSE]
if (nrow(results) != 18591L) stop("Unexpected number of protein result rows.")
if (anyDuplicated(paste(results$Timepoint, results$feature_id))) {
  stop("Duplicate protein feature IDs within a time point.")
}
if (anyNA(results$uniprot)) stop("At least one protein lacks a UniProt accession.")
results$source_package_version <- package_version

columns <- c(
  "tissue", "assay", "contrast_type", "contrast_category", "contrast", "contrast_short",
  "Timepoint", "feature_id", "uniprot", "gene_symbol",
  "logFC", "z.std", "p_value", "adj_p_value", "full_model",
  "source_package_version"
)
output <- results[, columns, drop = FALSE]
output_dir <- file.path(project_root, "data", "processed")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(output_dir, "motrpac_muscle_protein_ee_con.csv.gz")
connection <- gzfile(output_path, open = "wt")
write.csv(output, connection, row.names = FALSE, na = "")
close(connection)

cat("Wrote", nrow(output), "MoTrPAC muscle protein rows to", output_path, "\n")
