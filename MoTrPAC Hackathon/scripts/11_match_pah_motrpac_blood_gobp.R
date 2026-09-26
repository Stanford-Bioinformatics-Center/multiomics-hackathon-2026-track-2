#!/usr/bin/env Rscript

# Match PAH PBMC and healthy MoTrPAC blood RNA GOBP pathway results by exact
# GO set name, separately at each exercise time.
# Run from any directory with:
#   Rscript scripts/11_match_pah_motrpac_blood_gobp.R
#
# Inputs: outputs/PAH blood/ pathway tables from scripts 09 and 10
# Outputs: outputs/PAH blood/ (matched pathways and summaries)
#
# Each FDR belongs to its own cohort and test universe. Matching two
# significant results is descriptive; it is not a joint interaction test or
# evidence that exercise changes these pathways in PAH patients.

script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_flag) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_flag), mustWork = TRUE)
project_root <- dirname(dirname(script_path))
output_dir <- file.path(project_root, "outputs", "PAH blood")

pah <- read.csv(file.path(output_dir, "pah_gobp_camera.csv"),
                check.names = FALSE)
motrpac <- read.csv(gzfile(file.path(output_dir,
                                     "motrpac_blood_gobp_camera.csv.gz")),
                     check.names = FALSE)
if (!all(c("set", "NGenes", "Direction", "PValue", "FDR",
           "pah_camera_z") %in% names(pah)) ||
    !all(c("set", "Timepoint", "set_short", "set_size", "direction",
           "z.std", "p_value", "adj_p_value") %in% names(motrpac))) {
  stop("Required columns are missing from a pathway input table.")
}
times <- c("during_20_min", "during_40_min", "post_10_min",
           "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")
if (anyDuplicated(pah$set) ||
    anyDuplicated(motrpac[, c("Timepoint", "set")]) ||
    !setequal(unique(motrpac$Timepoint), times)) {
  stop("Duplicate pathway names or unexpected MoTrPAC time points.")
}

pah <- pah[, c("set", "NGenes", "Direction", "PValue", "FDR",
               "pah_camera_z")]
names(pah) <- c("set", "pah_set_genes", "pah_direction", "pah_p_value",
                "pah_fdr", "pah_camera_z")
motrpac <- motrpac[, c("set", "Timepoint", "set_short", "set_size",
                       "direction", "z.std", "p_value", "adj_p_value")]
names(motrpac) <- c("set", "Timepoint", "set_short", "motrpac_set_genes",
                    "motrpac_direction", "motrpac_z", "motrpac_p_value",
                    "motrpac_fdr")

matched <- merge(pah, motrpac, by = "set", sort = FALSE)
if (nrow(matched) != 20148L ||
    anyDuplicated(matched[, c("Timepoint", "set")]) ||
    !all(as.integer(table(factor(matched$Timepoint, levels = times))) ==
         3358L) ||
    !all(matched$pah_direction %in% c("Up", "Down")) ||
    !all(matched$motrpac_direction %in% c("Up", "Down")) ||
    !all(is.finite(matched$pah_fdr)) ||
    !all(is.finite(matched$motrpac_fdr))) {
  stop("Unexpected shared pathway universe or missing pathway statistics.")
}

matched$direction_relation <- ifelse(
  matched$pah_direction == matched$motrpac_direction, "same", "opposite")
matched$both_bh_below_0.05 <-
  matched$pah_fdr < 0.05 & matched$motrpac_fdr < 0.05
matched <- matched[order(match(matched$Timepoint, times),
                         matched$motrpac_fdr, matched$set), ]
write.csv(matched,
          gzfile(file.path(output_dir,
                           "pah_motrpac_blood_gobp_by_time.csv.gz")),
          row.names = FALSE)

# Summarize both kinds of direction relationship, rather than treating an
# opposite direction as the goal. Related GO labels often share many genes.
time_summary <- do.call(rbind, lapply(times, function(time) {
  x <- matched[matched$Timepoint == time, ]
  data.frame(
    Timepoint = time,
    shared_sets = nrow(x),
    both_bh_below_0.05 = sum(x$both_bh_below_0.05),
    same_direction_both_bh = sum(x$both_bh_below_0.05 &
                                   x$direction_relation == "same"),
    opposite_direction_both_bh = sum(x$both_bh_below_0.05 &
                                       x$direction_relation == "opposite")
  )
}))
write.csv(time_summary,
          file.path(output_dir, "pah_motrpac_blood_gobp_time_summary.csv"),
          row.names = FALSE)

set_summary <- aggregate(
  cbind(
    same_direction_both_bh = as.integer(
      matched$both_bh_below_0.05 & matched$direction_relation == "same"),
    opposite_direction_both_bh = as.integer(
      matched$both_bh_below_0.05 & matched$direction_relation == "opposite")
  ),
  by = list(set = matched$set), FUN = sum
)
names(set_summary)[2:3] <- c("same_direction_bh_timepoints",
                             "opposite_direction_bh_timepoints")
set_summary$both_bh_timepoints <-
  set_summary$same_direction_bh_timepoints +
  set_summary$opposite_direction_bh_timepoints
pah_details <- unique(matched[, c("set", "pah_direction", "pah_fdr")])
if (anyDuplicated(pah_details$set)) stop("PAH pathway details vary by time.")
set_summary <- merge(set_summary, pah_details, by = "set", sort = FALSE)
set_summary <- set_summary[order(-set_summary$both_bh_timepoints,
                                 set_summary$pah_fdr, set_summary$set), ]
write.csv(set_summary,
          file.path(output_dir, "pah_motrpac_blood_gobp_set_summary.csv"),
          row.names = FALSE)

print(time_summary)
cat("GO labels with >=2 opposite-direction BH time points:",
    sum(set_summary$opposite_direction_bh_timepoints >= 2L), "\n")
cat("GO labels with >=2 same-direction BH time points:",
    sum(set_summary$same_direction_bh_timepoints >= 2L), "\n")
