#!/usr/bin/env Rscript

# Rebuild the MoTrPAC blood gene ranking used by scripts 08 and 12 from
# published package summary statistics. This follows the gene-collapse rule
# in the original MoTrPAC project's 02_build_exercise_signature.R: within
# each assay, time, and symbol, retain the feature with the largest |z.std|.
# It exports RNA and Olink rows; the PAH comparison uses RNA rows only.
# No exercise model is refitted here.
# Run: Rscript scripts/07.5_export_motrpac_blood_gene_ranks.R

script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_flag) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_flag), mustWork = TRUE)
project_root <- dirname(dirname(script_path))

package_name <- "MotrpacHumanPreSuspensionAnalysis"
if (!requireNamespace(package_name, quietly = TRUE) ||
    !requireNamespace("dplyr", quietly = TRUE)) {
  stop("Install MotrpacHumanPreSuspensionAnalysis 0.2.4 and dplyr first.")
}
package_version <- as.character(utils::packageVersion(package_name))
if (package_version != "0.2.4") {
  stop("Expected MotrpacHumanPreSuspensionAnalysis 0.2.4; found ",
       package_version)
}
suppressPackageStartupMessages(library(MotrpacHumanPreSuspensionAnalysis))
suppressPackageStartupMessages(library(dplyr))

da <- as.data.frame(load_differential_analysis(
    selected_tissues = "blood",
    selected_omes = c("transcript-rna-seq", "prot-ol"),
    single_matrix = TRUE,
    combine_with_featgene = TRUE,
    verbose = FALSE
  ))
da <- da[as.character(da$contrast_type) == "exercise_with_controls" &
           as.character(da$contrast_category) == "EE-CON" &
           as.character(da$tissue) == "blood" &
           as.character(da$assay) %in% c("transcript-rna-seq", "prot-ol"), ]
# The original all-assay export has a platform column from metabolomics.
# It is empty for these two assays; the focused loader omits it entirely.
if (!"platform" %in% names(da)) da$platform <- NA_character_

clean_identifier <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}
for (column in c("tissue", "assay", "platform", "Timepoint",
                 "feature_id", "gene_symbol")) {
  da[[column]] <- clean_identifier(da[[column]])
}
times <- c("during_20_min", "during_40_min", "post_10_min",
           "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")
if (!setequal(unique(da$Timepoint), times) || anyNA(da$z.std)) {
  stop("Unexpected blood exercise time points or missing z.std values.")
}

ranked <- da %>%
  filter(!is.na(gene_symbol)) %>%
  mutate(abs_z = abs(z.std)) %>%
  arrange(tissue, assay, Timepoint, gene_symbol, desc(abs_z),
          adj_p_value, p_value, feature_id) %>%
  group_by(tissue, assay, Timepoint, gene_symbol) %>%
  mutate(source_feature_count = n()) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  group_by(tissue, assay, Timepoint) %>%
  mutate(rank_signed = min_rank(desc(z.std)),
         rank_absolute = min_rank(desc(abs_z))) %>%
  ungroup() %>%
  arrange(tissue, assay, Timepoint, desc(z.std)) %>%
  select(tissue, assay, platform, Timepoint, gene_symbol, entrez_gene,
         ensembl_gene, uniprot, feature_id, source_feature_count,
         logFC, z.std, abs_z, p_value, adj_p_value, rank_signed,
         rank_absolute, contrast_short)

per_assay_time <- table(ranked$assay, ranked$Timepoint)
if (nrow(ranked) != 110274L ||
    any(as.integer(per_assay_time["transcript-rna-seq", ]) != 16974L) ||
    any(as.integer(per_assay_time["prot-ol", ]) != 1405L) ||
    anyDuplicated(ranked[, c("assay", "Timepoint", "gene_symbol")])) {
  stop("The published blood gene ranking has an unexpected size or duplicate.")
}

output_dir <- file.path(project_root, "data", "processed")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(output_dir, "blood_gene_ranked.csv.gz")
connection <- gzfile(output_path, open = "wt")
write.csv(ranked, connection, row.names = FALSE, na = "")
close(connection)
print(per_assay_time)
cat("Wrote", nrow(ranked), "published-result gene ranks to",
    output_path, "\n")
