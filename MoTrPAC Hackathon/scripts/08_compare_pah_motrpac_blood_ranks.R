#!/usr/bin/env Rscript

# Compare the GSE33463 IPAH-minus-healthy PBMC ranking with MoTrPAC blood RNA
# endurance-minus-time-matched-control rankings at six exercise times.
# Run from any directory with:
#   Rscript scripts/08_compare_pah_motrpac_blood_ranks.R
#
# Inputs: outputs/PAH blood/pah_blood_gene_ranked.csv.gz (script 07)
#         data/processed/blood_gene_ranked.csv.gz
#         (MoTrPAC control-adjusted export; scripts 01-02 in the source project)
# Outputs: outputs/PAH blood/ (matched genes and rank association)
#
# This is a cross-cohort rank association, not an exercise effect in PAH.
# Gene-label permutations disrupt gene-gene dependence, so their P values
# are exploratory; the size and sign of Spearman rho are the main results.

script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_flag) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_flag), mustWork = TRUE)
project_root <- dirname(dirname(script_path))
output_dir <- file.path(project_root, "outputs", "PAH blood")

pah_path <- file.path(output_dir, "pah_blood_gene_ranked.csv.gz")
exercise_path <- file.path(project_root, "data", "processed",
                           "blood_gene_ranked.csv.gz")
pah <- read.csv(gzfile(pah_path), check.names = FALSE)
exercise <- read.csv(gzfile(exercise_path), check.names = FALSE)

if (!all(c("gene_symbol", "logFC", "t", "p_value", "adj_p_value") %in%
         names(pah)) ||
    !all(c("tissue", "assay", "Timepoint", "gene_symbol", "logFC",
           "z.std", "p_value", "adj_p_value") %in% names(exercise))) {
  stop("A required PAH or MoTrPAC ranking column is missing.")
}
exercise <- exercise[exercise$tissue == "blood" &
                       exercise$assay == "transcript-rna-seq", ]
times <- c("during_20_min", "during_40_min", "post_10_min",
           "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")
if (anyDuplicated(pah$gene_symbol) ||
    anyDuplicated(exercise[, c("Timepoint", "gene_symbol")]) ||
    !setequal(unique(exercise$Timepoint), times)) {
  stop("Duplicate gene symbols or unexpected MoTrPAC RNA time points.")
}

set.seed(33463)
joined <- vector("list", length(times))
summary_rows <- vector("list", length(times))
for (i in seq_along(times)) {
  time <- times[i]
  ex <- exercise[exercise$Timepoint == time,
                 c("gene_symbol", "logFC", "z.std", "p_value", "adj_p_value")]
  names(ex)[-1L] <- paste0("exercise_", names(ex)[-1L])
  x <- merge(pah[, c("gene_symbol", "logFC", "t", "p_value",
                      "adj_p_value")], ex, by = "gene_symbol", sort = FALSE)
  x <- x[is.finite(x$t) & is.finite(x$exercise_z.std), ]
  if (anyDuplicated(x$gene_symbol) || nrow(x) != 13055L) {
    stop("Unexpected shared-gene universe at ", time, ": ", nrow(x))
  }
  names(x)[2:5] <- paste0("pah_", names(x)[2:5])
  x$Timepoint <- time
  joined[[i]] <- x

  # Spearman rho is Pearson correlation of ranks within the shared genes.
  # Negative: opposing ranks; positive: same-direction ranks. Neither sign
  # establishes a mechanism or a treatment response in patients.
  a <- as.numeric(scale(rank(x$pah_t)))
  b <- as.numeric(scale(rank(x$exercise_z.std)))
  rho <- sum(a * b) / (length(a) - 1L)

  # Retain the original analysis's 10,000 label shuffles for reproducibility.
  # They do not preserve coexpression, so these P values are descriptive.
  permuted <- numeric(10000L)
  for (j in seq_along(permuted)) {
    permuted[j] <- sum(a * b[sample.int(length(b))]) / (length(a) - 1L)
  }
  summary_rows[[i]] <- data.frame(
    Timepoint = time,
    shared_genes = nrow(x),
    spearman_rho = rho,
    permutation_p_two_sided = (1 + sum(abs(permuted) >= abs(rho))) /
      (length(permuted) + 1),
    permutations = length(permuted),
    gene_universe = "genes with finite PAH t and MoTrPAC blood RNA z.std",
    permutation_null = "shuffled gene labels; gene-gene dependence not preserved"
  )
  cat(time, "shared genes", nrow(x), "rho", round(rho, 4), "\n")
}

summary <- do.call(rbind, summary_rows)
summary$permutation_fdr_six_times <-
  p.adjust(summary$permutation_p_two_sided, method = "BH")
write.csv(summary,
          file.path(output_dir,
                    "pah_motrpac_blood_rank_association_by_time.csv"),
          row.names = FALSE)
write.csv(do.call(rbind, joined),
          gzfile(file.path(output_dir,
                           "pah_motrpac_blood_matched_genes.csv.gz")),
          row.names = FALSE)
cat("Saved six time points; negative rho indicates opposing ranks.\n")
