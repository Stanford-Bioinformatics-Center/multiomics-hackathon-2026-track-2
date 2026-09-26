#!/usr/bin/env Rscript

# Test GO Biological Process (GOBP) sets against the full GSE33463 PBMC
# IPAH-minus-healthy gene ranking. This uses PAH data only; the exercise
# pathway comparison will be a separate module.
# Run from any directory with:
#   Rscript scripts/09_test_pah_gobp_pathways.R
#
# Input:  outputs/PAH blood/pah_blood_gene_ranked.csv.gz (script 07)
# Source: MotrpacHumanPreSuspensionAnalysis 0.2.4 GOBP gene sets
# Outputs: outputs/PAH blood/pah_gobp_camera.csv and QC
#
# A competitive pathway result describes a pattern across measured genes.
# It does not show that every member changed or that a pathway is causal.

script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_flag) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_flag), mustWork = TRUE)
project_root <- dirname(dirname(script_path))
output_dir <- file.path(project_root, "outputs", "PAH blood")

if (!requireNamespace("limma", quietly = TRUE) ||
    !requireNamespace("MotrpacHumanPreSuspensionAnalysis", quietly = TRUE)) {
  stop("Install limma and MotrpacHumanPreSuspensionAnalysis 0.2.4 first.")
}
package_version <- as.character(utils::packageVersion(
  "MotrpacHumanPreSuspensionAnalysis"))
if (package_version != "0.2.4") {
  stop("Expected MotrpacHumanPreSuspensionAnalysis 0.2.4; found ",
       package_version)
}

gene <- read.csv(gzfile(file.path(output_dir,
                                   "pah_blood_gene_ranked.csv.gz")),
                 check.names = FALSE)
if (!all(c("gene_symbol", "t") %in% names(gene)) ||
    nrow(gene) != 19593L || anyNA(gene$gene_symbol) ||
    anyDuplicated(gene$gene_symbol) || !all(is.finite(gene$t))) {
  stop("Unexpected GSE33463 gene ranking from script 07.")
}

# Use the same GOBP set definitions as the original GMT-based analysis.
# Package 0.2.4 was checked locally against its source GMT: all 7,647 GOBP
# names, order, and gene memberships matched. Use numerical indices so the
# gene statistics and set membership follow the same row order.
signatures <- suppressPackageStartupMessages(getExportedValue(
  "MotrpacHumanPreSuspensionAnalysis", "MOLECULAR_SIGNATURES"))
gene_sets <- lapply(signatures[["GOBP"]], unique)
if (length(gene_sets) != 7647L ||
    anyDuplicated(names(gene_sets))) {
  stop("Unexpected GOBP collection in package 0.2.4.")
}
index <- lapply(gene_sets, function(s) which(gene$gene_symbol %in% s))
# Use tested-gene membership, not the unfiltered size of each GO set.
index <- index[lengths(index) >= 10L & lengths(index) <= 500L]
if (length(index) <= 1000L) stop("Too few GOBP sets passed size filters.")

# cameraPR tests whether genes in a set tend toward more extreme signed t
# statistics than genes outside it. FDR is BH across the tested GOBP sets.
pah_camera <- limma::cameraPR(gene$t, index = index)
pah_camera$set <- rownames(pah_camera)
pah_camera$pah_camera_z <-
  ifelse(pah_camera$Direction == "Up", 1, -1) *
  stats::qnorm(pmax(pah_camera$PValue, .Machine$double.xmin) / 2,
               lower.tail = FALSE)
write.csv(pah_camera, file.path(output_dir, "pah_gobp_camera.csv"),
          row.names = FALSE)

qc <- data.frame(
  metric = c("pah_genes_ranked", "gobp_sets_available",
             "gobp_sets_tested_10_to_500_genes", "gobp_sets_bh_below_0.05",
             "motrpac_package_version"),
  value = c(nrow(gene), length(gene_sets), length(index),
            sum(pah_camera$FDR < 0.05), package_version)
)
write.csv(qc, file.path(output_dir, "pah_gobp_camera_qc.csv"),
          row.names = FALSE)

print(qc)
cat("GOBP results use the IPAH-minus-healthy PBMC gene ranking.\n")
