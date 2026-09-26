#!/usr/bin/env Rscript

# Export MoTrPAC's precomputed GO Biological Process (GOBP) pathway results
# for healthy-volunteer blood RNA after endurance exercise, each time versus
# its time-matched control. This creates the exercise-side pathway table; it
# does not compare the table with PAH or estimate exercise effects in PAH.
# Run from any directory with:
#   Rscript scripts/10_export_motrpac_blood_gobp.R
#
# Source: MotrpacHumanPreSuspensionAnalysis 0.2.4 CAMERA_RESULTS
# Outputs: outputs/PAH blood/motrpac_blood_gobp_camera.csv.gz and QC

script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_flag) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_flag), mustWork = TRUE)
project_root <- dirname(dirname(script_path))
output_dir <- file.path(project_root, "outputs", "PAH blood")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

package_name <- "MotrpacHumanPreSuspensionAnalysis"
if (!requireNamespace(package_name, quietly = TRUE)) {
  stop("Install MotrpacHumanPreSuspensionAnalysis 0.2.4 first.")
}
package_version <- as.character(utils::packageVersion(package_name))
if (package_version != "0.2.4") {
  stop("Expected MotrpacHumanPreSuspensionAnalysis 0.2.4; found ",
       package_version)
}
suppressPackageStartupMessages(library(MotrpacHumanPreSuspensionAnalysis))

# Filter the package's existing CAMERA tests after they were calculated.
# This keeps their original P values and adjusted P values intact.
results <- as.data.frame(CAMERA_RESULTS)
required <- c("tissue", "assay", "contrast_type", "collection", "database",
              "contrast", "contrast_short", "set", "set_short", "set_size",
              "direction", "z.std", "p_value", "adj_p_value")
if (!all(required %in% names(results))) {
  stop("CAMERA_RESULTS lacks expected pathway columns.")
}
keep <- as.character(results$tissue) == "blood" &
  as.character(results$assay) == "transcript-rna-seq" &
  as.character(results$contrast_type) == "exercise_with_controls" &
  as.character(results$collection) == "C5" &
  as.character(results$database) == "GOBP" &
  startsWith(as.character(results$contrast_short), "Endur.")
results <- results[keep, , drop = FALSE]

# The contrast dictionary establishes the six EE-CON time points explicitly.
dictionary <- as.data.frame(CONTRAST_CONVERTER)
index <- match(as.character(results$contrast),
               as.character(dictionary$contrast))
if (anyNA(index)) stop("A CAMERA contrast is absent from CONTRAST_CONVERTER.")
results$Timepoint <- as.character(dictionary$Timepoint[index])
results$contrast_category <- as.character(
  dictionary$contrast_category[index])
times <- c("during_20_min", "during_40_min", "post_10_min",
           "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")
per_time <- table(factor(results$Timepoint, levels = times))
if (nrow(results) != 28368L ||
    !setequal(unique(results$Timepoint), times) ||
    !all(as.integer(per_time) == 4728L) ||
    any(results$contrast_category != "EE-CON") ||
    anyDuplicated(results[, c("Timepoint", "set")])) {
  stop("Unexpected MoTrPAC blood RNA GOBP set or time-point count.")
}

results <- results[order(match(results$Timepoint, times), results$set), ]
results$source_package_version <- package_version
columns <- c("tissue", "assay", "contrast_category", "contrast_short",
             "Timepoint", "collection", "database", "set", "set_short",
             "set_size", "direction", "z.std", "p_value", "adj_p_value",
             "source_package_version")
write.csv(results[, columns],
          gzfile(file.path(output_dir, "motrpac_blood_gobp_camera.csv.gz")),
          row.names = FALSE, na = "")

qc <- data.frame(
  Timepoint = times,
  pathway_rows = as.integer(per_time),
  pathways_with_bh_below_0.05 = vapply(times, function(time) {
    sum(results$Timepoint == time & results$adj_p_value < 0.05)
  }, integer(1)),
  source_package_version = package_version
)
write.csv(qc, file.path(output_dir, "motrpac_blood_gobp_camera_qc.csv"),
          row.names = FALSE)

print(qc)
cat("Exported precomputed MoTrPAC blood RNA pathway results.\n")
