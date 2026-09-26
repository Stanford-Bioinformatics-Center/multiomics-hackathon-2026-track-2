#!/usr/bin/env Rscript

# Export MoTrPAC's precomputed muscle RNA CAMERA-PR result for the GO
# biological-process oxidative-phosphorylation set at three endurance times.
# This is an exploratory pathway follow-up to the PAH paper's lower muscle
# proteins. The genes corresponding to five of those nine proteins belong to
# this GO OXPHOS set (NDUFA9,
# UQCRC2, UQCRC1, ATP5MG, ATP5F1B); the other four do not. Script 00
# separately counts RNA results for all nine PAH-linked genes at all three
# times. This script tests neither those nine genes as a set nor PAH patients.
# OXPHOS and the 24-hour emphasis were not prospectively selected here.
# Run with: Rscript scripts/05_export_motrpac_oxphos_pathway.R

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

# CAMERA_RESULTS contains pathway tests MoTrPAC already ran. Filtering after
# the test preserves their original BH-adjusted p-values and test universe.
results <- as.data.frame(CAMERA_RESULTS)
keep <- as.character(results$tissue) == "muscle" &
  as.character(results$assay) == "transcript-rna-seq" &
  as.character(results$contrast_type) == "exercise_with_controls" &
  as.character(results$collection) == "C5" &
  as.character(results$database) == "GOBP" &
  as.character(results$set) == "GOBP_OXIDATIVE_PHOSPHORYLATION" &
  startsWith(as.character(results$contrast_short), "Endur.")
results <- results[keep, , drop = FALSE]

# Use the package's contrast dictionary to recover the time point and
# verify that these really are the EE-CON comparisons.
dictionary <- as.data.frame(CONTRAST_CONVERTER)
index <- match(as.character(results$contrast), as.character(dictionary$contrast))
if (anyNA(index)) stop("A CAMERA contrast is absent from CONTRAST_CONVERTER.")
results$Timepoint <- as.character(dictionary$Timepoint[index])
results$contrast_category <- as.character(dictionary$contrast_category[index])
times <- c("post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")
if (nrow(results) != 3L || !setequal(results$Timepoint, times) ||
    any(results$contrast_category != "EE-CON")) {
  stop("Expected one endurance-versus-control OXPHOS result at each muscle time.")
}
results <- results[match(times, results$Timepoint), , drop = FALSE]
results$source_package_version <- package_version

columns <- c(
  "tissue", "assay", "contrast_category", "contrast", "contrast_short",
  "Timepoint", "collection", "database", "set_id", "set", "set_size",
  "direction", "z.std", "p_value", "adj_p_value", "source_package_version"
)
output <- results[, columns, drop = FALSE]
output_dir <- file.path(project_root, "outputs", "PAH protein")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(output_dir, "motrpac_muscle_oxphos_camera.csv")
write.csv(output, output_path, row.names = FALSE, na = "")

cat("Wrote", nrow(output), "precomputed pathway rows to", output_path, "\n")
