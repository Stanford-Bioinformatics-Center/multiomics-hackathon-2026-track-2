# Export MoTrPAC's precomputed CAMERA pathway results (MotrpacHumanPreSuspensionAnalysis 2.0.8,
# copied in MoTrPAC Hackathon/Transcriptomics/data/motrpac/CAMERA_RESULTS.rda) for the explorer's
# pathway-name input. Keeps exercise-vs-control, exercisers-only and control-only contrasts.
# Run from the repository root:
#   Rscript apps/api/scripts_build_pathways.R
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else "."
e <- new.env()
load(file.path(root, "MoTrPAC Hackathon/Transcriptomics/data/motrpac/CAMERA_RESULTS.rda"), envir = e)
x <- as.data.frame(e$CAMERA_RESULTS)
x <- x[x$contrast_type %in% c("exercise_with_controls", "exercise_no_controls", "control_only"),
       c("tissue", "assay", "contrast_type", "contrast_short", "collection", "database", "set", "set_short",
         "set_size", "direction", "z.std", "p_value", "adj_p_value")]
out <- file.path(root, "apps/api/data/motrpac_camera_pathways.csv.gz")
con <- gzfile(out, "wt"); write.csv(x, con, row.names = FALSE); close(con)
cat("rows:", nrow(x), "->", out, "\n")
