#!/usr/bin/env Rscript

# Check whether the PAH-versus-MoTrPAC blood RNA rank correlation changes
# when two blood-related GO sets are excluded or a different probe is chosen.
# Run from any directory with:
#   Rscript scripts/12_check_pah_motrpac_blood_rank_sensitivity.R
#
# Inputs: outputs/PAH blood/ probe and gene rankings (script 07)
#         data/processed/blood_gene_ranked.csv.gz (script 07.5 output)
#         MotrpacHumanPreSuspensionAnalysis 0.2.4 GOBP memberships
# Output: outputs/PAH blood/pah_motrpac_blood_rank_sensitivity.csv
#
# These checks test two specific analysis choices. They do not adjust for
# blood-cell composition, batch, or other differences between cohorts.

script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_flag) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_flag), mustWork = TRUE)
project_root <- dirname(dirname(script_path))
output_dir <- file.path(project_root, "outputs", "PAH blood")

package_name <- "MotrpacHumanPreSuspensionAnalysis"
if (!requireNamespace(package_name, quietly = TRUE)) {
  stop("Install MotrpacHumanPreSuspensionAnalysis 0.2.4 first.")
}
package_version <- as.character(utils::packageVersion(package_name))
if (package_version != "0.2.4") {
  stop("Expected MotrpacHumanPreSuspensionAnalysis 0.2.4; found ",
       package_version)
}

primary <- read.csv(gzfile(file.path(output_dir,
                                      "pah_blood_gene_ranked.csv.gz")),
                    check.names = FALSE)
probes <- read.csv(gzfile(file.path(output_dir,
                                     "pah_blood_probe_ranked.csv.gz")),
                   check.names = FALSE)
exercise <- read.csv(gzfile(file.path(project_root, "data", "processed",
                                         "blood_gene_ranked.csv.gz")),
                     check.names = FALSE)
if (!all(c("gene_symbol", "t") %in% names(primary)) ||
    !all(c("gene_symbol", "probe_id", "t", "AveExpr") %in% names(probes)) ||
    !all(c("tissue", "assay", "Timepoint", "gene_symbol", "z.std") %in%
         names(exercise))) {
  stop("A required input column is missing.")
}

# Primary script 07 retained the largest |moderated t| probe for each gene.
# This alternate rule uses the largest overall average measured expression,
# then probe ID to break ties. It does not select directly on the group t.
probes <- probes[!is.na(probes$gene_symbol) & is.finite(probes$t) &
                   is.finite(probes$AveExpr), ]
probes <- probes[order(probes$gene_symbol, -probes$AveExpr,
                       probes$probe_id), ]
alternate <- probes[!duplicated(probes$gene_symbol), ]
if (anyDuplicated(primary$gene_symbol) ||
    anyDuplicated(alternate$gene_symbol) ||
    !setequal(primary$gene_symbol, alternate$gene_symbol)) {
  stop("Primary and alternate probe rules have different gene universes.")
}

signatures <- suppressPackageStartupMessages(getExportedValue(
  package_name, "MOLECULAR_SIGNATURES"))
gobp <- signatures[["GOBP"]]
sets_to_exclude <- c("GOBP_OXYGEN_TRANSPORT",
                     "GOBP_HYDROGEN_PEROXIDE_CATABOLIC_PROCESS")
if (!all(sets_to_exclude %in% names(gobp))) {
  stop("A sensitivity-check GO gene set is missing.")
}
excluded_genes <- unique(unlist(gobp[sets_to_exclude], use.names = FALSE))

exercise <- exercise[exercise$tissue == "blood" &
                       exercise$assay == "transcript-rna-seq", ]
times <- c("during_20_min", "during_40_min", "post_10_min",
           "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")
if (anyDuplicated(exercise[, c("Timepoint", "gene_symbol")]) ||
    !setequal(unique(exercise$Timepoint), times)) {
  stop("Duplicate exercise genes or unexpected MoTrPAC time points.")
}

rows <- vector("list", length(times) * 3L)
k <- 0L
for (time in times) {
  ex <- exercise[exercise$Timepoint == time,
                 c("gene_symbol", "z.std")]
  for (analysis in c("primary_max_abs_t",
                     "exclude_oxygen_transport_and_peroxide_catabolism",
                     "highest_mean_expression_probe")) {
    disease <- if (analysis == "highest_mean_expression_probe")
      alternate else primary
    x <- merge(disease[, c("gene_symbol", "t")], ex,
               by = "gene_symbol", sort = FALSE)
    x <- x[is.finite(x$t) & is.finite(x$z.std), ]
    if (analysis == "exclude_oxygen_transport_and_peroxide_catabolism") {
      x <- x[!(x$gene_symbol %in% excluded_genes), ]
    }
    if (anyDuplicated(x$gene_symbol) || nrow(x) < 13000L) {
      stop("Unexpected shared-gene universe at ", time, " / ", analysis)
    }
    k <- k + 1L
    rows[[k]] <- data.frame(
      Timepoint = time,
      analysis = analysis,
      shared_genes = nrow(x),
      spearman_rho = stats::cor(x$t, x$z.std, method = "spearman")
    )
  }
}
result <- do.call(rbind, rows)
write.csv(result,
          file.path(output_dir, "pah_motrpac_blood_rank_sensitivity.csv"),
          row.names = FALSE)

print(result)
cat("These two sensitivity checks do not remove all cell-composition effects.\n")
