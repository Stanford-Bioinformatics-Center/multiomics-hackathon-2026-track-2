# Export MoTrPAC blood and muscle metabolomics summary results.
#
# Three endurance-exercise contrasts, all precomputed by MoTrPAC:
#   EE-CON  (exercise post - pre) - (control post - pre): exercise effect
#   EE-EE   exercise post - pre only: comparable to patient peak - rest
#   CON-CON control post - pre only: drift without exercise
# No model is fitted here.

library(MotrpacHumanPreSuspensionAnalysis)

package_version <- as.character(packageVersion("MotrpacHumanPreSuspensionAnalysis"))
if (package_version != "2.0.8") {
  stop("Expected MotrpacHumanPreSuspensionAnalysis 2.0.8; found ", package_version)
}

# Find the module folder from this script's own location.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
module_dir <- if (length(script_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
} else {
  getwd()
}
processed_dir <- file.path(module_dir, "data", "processed")
output_dir <- file.path(module_dir, "output")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

results <- as.data.frame(load_differential_analysis(
  selected_tissues = c("blood", "muscle"),
  selected_omes = "all",
  single_matrix = TRUE,
  epigen = FALSE,
  combine_with_featgene = TRUE,
  verbose = FALSE
))

metab <- subset(
  results,
  assay == "metab" &
    contrast_category %in% c("EE-CON", "EE-EE", "CON-CON") &
    Timepoint != "pre_exercise"
)
metab <- metab[, c("tissue", "platform", "feature_id", "refmet_name", "refmet_id",
                   "contrast_category", "Timepoint", "logFC", "AveExpr",
                   "p_value", "adj_p_value")]
metab$source_package_version <- package_version
metab$source_collection <- "c2.0"

connection <- gzfile(file.path(processed_dir, "motrpac_metabolite_contrasts.csv.gz"), open = "wt")
write.csv(metab, connection, row.names = FALSE)
close(connection)

counts <- aggregate(
  feature_id ~ tissue + contrast_category + Timepoint,
  data = metab,
  FUN = length
)
names(counts)[names(counts) == "feature_id"] <- "n_features"
write.csv(counts, file.path(output_dir, "04_motrpac_export_counts.csv"), row.names = FALSE)

cat("Rows written:", nrow(metab), "\n")
print(counts)
