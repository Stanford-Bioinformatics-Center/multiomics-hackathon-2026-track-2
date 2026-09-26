#!/usr/bin/env Rscript

# Reproduce the exploratory time-point check for the nine proteins reported
# lower in PAH muscle by Malenfant et al. (2015), Table 2.
# Run from any directory with:
#   Rscript scripts/00_audit_pah_muscle_rna_timepoints.R
#
# Inputs: data/from paper/pah_lower_proteins_malenfant2015.csv and the published
# MotrpacHumanPreSuspensionAnalysis 0.2.4 differential-analysis summaries.
# No statistical model is fitted by this script. Each MoTrPAC BH q value was
# calculated by the source package within its full tissue/assay/contrast set.

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

pah_path <- file.path(project_root, "data", "from paper",
                      "pah_lower_proteins_malenfant2015.csv")
pah <- read.csv(pah_path, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(pah) != 9L || anyDuplicated(pah$uniprot_accession) ||
    anyNA(pah$uniprot_accession) || anyNA(pah$paper_symbol)) {
  stop("Expected nine distinct PAH paper proteins with UniProt accessions.")
}

times <- c("post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")

# Use exact UniProt matches in the package's protein results to obtain current
# gene symbols. This resolves the paper's ATP5L and ATP5B names to ATP5MG
# and ATP5F1B rather than silently assuming the names still match.
protein <- as.data.frame(load_differential_analysis(
  selected_tissues = "muscle",
  selected_omes = "prot-pr",
  single_matrix = TRUE,
  combine_with_featgene = TRUE,
  verbose = FALSE
))
protein <- protein[
  as.character(protein$contrast_category) == "EE-CON" &
    as.character(protein$Timepoint) %in% times &
    as.character(protein$uniprot) %in% pah$uniprot_accession,
  , drop = FALSE
]
current_symbol <- character(nrow(pah))
for (i in seq_len(nrow(pah))) {
  matches <- protein[as.character(protein$uniprot) == pah$uniprot_accession[i],
                     , drop = FALSE]
  if (nrow(matches) != 3L ||
      !setequal(as.character(matches$Timepoint), times) ||
      anyNA(matches$gene_symbol) ||
      length(unique(as.character(matches$gene_symbol))) != 1L) {
    stop("Expected one consistent protein-to-gene mapping at each time for ",
         pah$uniprot_accession[i])
  }
  current_symbol[i] <- as.character(matches$gene_symbol[1])
}
rm(protein)

# This loader returns already computed, control-adjusted EE-CON model results.
# Keep all three sampled post-exercise muscle RNA times for this audit.
rna <- as.data.frame(load_differential_analysis(
  selected_tissues = "muscle",
  selected_omes = "transcript-rna-seq",
  single_matrix = TRUE,
  combine_with_featgene = TRUE,
  verbose = FALSE
))
rna <- rna[
  as.character(rna$contrast_category) == "EE-CON" &
    as.character(rna$Timepoint) %in% times &
    as.character(rna$gene_symbol) %in% current_symbol,
  , drop = FALSE
]

detail <- vector("list", nrow(pah) * length(times))
k <- 0L
for (time in times) {
  for (i in seq_len(nrow(pah))) {
    matched <- rna[
      as.character(rna$Timepoint) == time &
        as.character(rna$gene_symbol) == current_symbol[i],
      , drop = FALSE
    ]
    if (nrow(matched) != 1L || is.na(matched$logFC) ||
        is.na(matched$adj_p_value)) {
      stop("Expected exactly one complete RNA result for ",
           current_symbol[i], " at ", time)
    }
    k <- k + 1L
    positive <- as.numeric(matched$logFC) > 0
    positive_q <- positive && as.numeric(matched$adj_p_value) < 0.05
    detail[[k]] <- data.frame(
      paper_symbol = pah$paper_symbol[i],
      uniprot_accession = pah$uniprot_accession[i],
      motrpac_current_symbol = current_symbol[i],
      motrpac_rna_timepoint = time,
      motrpac_rna_feature_id = as.character(matched$feature_id),
      motrpac_rna_logFC = as.numeric(matched$logFC),
      motrpac_rna_bh_adj_p_value = as.numeric(matched$adj_p_value),
      positive_logFC = positive,
      positive_and_bh_q_lt_0_05 = positive_q,
      source_doi = pah$source_doi[i],
      motrpac_package_version = package_version,
      stringsAsFactors = FALSE
    )
  }
}
detail <- do.call(rbind, detail)

summary <- do.call(rbind, lapply(times, function(time) {
  at_time <- detail[detail$motrpac_rna_timepoint == time, , drop = FALSE]
  data.frame(
    motrpac_rna_timepoint = time,
    matching_genes = nrow(at_time),
    positive_logFC_count = sum(at_time$positive_logFC),
    positive_and_bh_q_lt_0_05_count = sum(at_time$positive_and_bh_q_lt_0_05),
    stringsAsFactors = FALSE
  )
}))
if (!identical(summary$matching_genes, rep(9L, length(times)))) {
  stop("Expected nine gene results at each time.")
}

output_dir <- file.path(project_root, "outputs", "PAH protein")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
detail_path <- file.path(output_dir, "pah_muscle_rna_timepoint_detail.csv")
summary_path <- file.path(output_dir, "pah_muscle_rna_timepoint_summary.csv")
write.csv(detail, detail_path, row.names = FALSE, na = "")
write.csv(summary, summary_path, row.names = FALSE, na = "")

print(summary, row.names = FALSE)
cat("Wrote", nrow(detail), "gene-by-time rows to", detail_path, "\n")
cat("Wrote the three-timepoint summary to", summary_path, "\n")
