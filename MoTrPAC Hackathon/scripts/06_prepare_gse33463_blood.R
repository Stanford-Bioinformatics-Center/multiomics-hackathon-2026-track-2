#!/usr/bin/env Rscript

# Prepare GEO GSE33463 PBMC expression for an IPAH-versus-healthy comparison.
# GEO source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33463
# Run from any directory with:
#   Rscript scripts/06_prepare_gse33463_blood.R
#
# Input: data/processed/GSE33463/GSE33463_series_matrix.txt.gz
# The submitter already median-scaled and log2-transformed these expression
# values. This script selects samples and checks the matrix; it does not
# normalize again or fit a differential-expression model.

script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_flag) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_flag), mustWork = TRUE)
project_root <- dirname(dirname(script_path))

input_path <- file.path(project_root, "data", "processed", "GSE33463",
                        "GSE33463_series_matrix.txt.gz")
output_dir <- file.path(project_root, "outputs", "PAH blood")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# GEO puts sample metadata before !series_matrix_table_begin. Read those lines
# once, then use the same count to skip directly to the expression table.
connection <- gzfile(input_path, open = "rt")
metadata <- character()
repeat {
  line <- readLines(connection, n = 1L, warn = FALSE)
  if (!length(line)) stop("Series Matrix table marker was not found.")
  metadata <- c(metadata, line)
  if (identical(line, "!series_matrix_table_begin")) break
}
close(connection)

field <- function(key) {
  rows <- metadata[startsWith(metadata, paste0(key, "\t"))]
  if (length(rows) != 1L) stop("Expected one metadata field: ", key)
  value <- strsplit(rows, "\t", fixed = TRUE)[[1L]][-1L]
  sub('"$', "", sub('^"', "", value))
}

samples <- data.frame(
  sample_id = field("!Sample_geo_accession"),
  title = field("!Sample_title"),
  source = field("!Sample_source_name_ch1"),
  disease_status = sub("^disease status: ", "",
                       field("!Sample_characteristics_ch1")),
  description = field("!Sample_description"),
  stringsAsFactors = FALSE
)
if (nrow(samples) != 140L || anyDuplicated(samples$sample_id) ||
    any(samples$source != "PBMC")) {
  stop("Expected 140 unique PBMC samples.")
}

# The other 69 samples include systemic-sclerosis and other PH groups; they
# are excluded from this specific idiopathic-PAH-versus-healthy comparison.
samples$group <- ifelse(
  samples$disease_status == "control", "Healthy",
  ifelse(samples$disease_status == "Idiopathic Pulmonary Arterial Hypertension",
         "IPAH", "Other")
)
if (sum(samples$group == "Healthy") != 41L ||
    sum(samples$group == "IPAH") != 30L ||
    sum(samples$group == "Other") != 69L ||
    any(samples$description[samples$group == "Healthy"] != "control") ||
    any(samples$description[samples$group == "IPAH"] !=
        "Idiopathic Pulmonary Arterial Hypertension")) {
  stop("Unexpected GSE33463 group labels or counts.")
}

all_samples_path <- file.path(output_dir, "gse33463_samples_all.csv")
selected_samples_path <- file.path(output_dir,
                                   "gse33463_samples_ipah_healthy.csv")
write.csv(samples, all_samples_path, row.names = FALSE)
selected <- samples$group != "Other"
selected_samples <- samples[selected, , drop = FALSE]
write.csv(selected_samples, selected_samples_path, row.names = FALSE)

matrix_table <- read.delim(gzfile(input_path), skip = length(metadata),
                           nrows = 48803L, check.names = FALSE)
if (nrow(matrix_table) != 48803L ||
    !identical(names(matrix_table)[-1L], samples$sample_id) ||
    anyDuplicated(matrix_table$ID_REF)) {
  stop("Unexpected probe count, sample order, or duplicate probe IDs.")
}
expression <- as.matrix(matrix_table[, -1L])
storage.mode(expression) <- "double"
rownames(expression) <- matrix_table$ID_REF
if (anyNA(expression) || !all(is.finite(expression))) {
  stop("Expression matrix contains missing or non-finite values.")
}
expression <- expression[, selected, drop = FALSE]
if (!identical(colnames(expression), selected_samples$sample_id) ||
    ncol(expression) != 71L) {
  stop("Selected expression columns do not match the sample sheet.")
}

quantiles <- t(apply(expression, 2L, quantile,
                     probs = c(0, 0.25, 0.5, 0.75, 1)))
qc <- data.frame(
  sample_id = colnames(expression),
  group = selected_samples$group,
  min = quantiles[, 1L],
  q1 = quantiles[, 2L],
  median = quantiles[, 3L],
  q3 = quantiles[, 4L],
  max = quantiles[, 5L],
  missing_fraction = colMeans(is.na(expression)),
  stringsAsFactors = FALSE
)
qc_path <- file.path(output_dir, "gse33463_sample_qc.csv")
write.csv(qc, qc_path, row.names = FALSE)

plot_path <- file.path(output_dir, "gse33463_expression_boxplot.png")
png(plot_path, width = 1800, height = 900, res = 150)
boxplot(as.data.frame(expression), outline = FALSE, xaxt = "n",
        col = ifelse(selected_samples$group == "IPAH", "#bd6c58", "#6e9fc2"),
        border = "#59636d", main = "GSE33463 processed PBMC expression by sample",
        ylab = "log2 expression", xlab = "41 healthy controls (blue), 30 IPAH (orange)")
abline(h = 8, col = "#555555", lty = 2)
dev.off()

prepared_path <- file.path(output_dir, "gse33463_prepared.rds")
saveRDS(list(expression = expression, samples = selected_samples),
        prepared_path, compress = TRUE)

cat("Prepared", nrow(expression), "probes from", ncol(expression), "samples\n")
print(table(selected_samples$group))
cat("Median expression across samples:", range(qc$median), "\n")
cat("Saved selected matrix to", prepared_path, "\n")
