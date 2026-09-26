#!/usr/bin/env Rscript

# Fit the GSE33463 IPAH-minus-healthy PBMC comparison, one test per probe.
# Run from any directory with:
#   Rscript scripts/07_fit_gse33463_blood_limma.R
#
# Inputs: outputs/PAH blood/gse33463_prepared.rds (script 06)
#         data/processed/GPL6947.annot.gz (GEO platform annotation)
# Outputs: outputs/PAH blood/ (design, probe/gene rankings, mapping QC)
#
# The gene table retains the probe with the largest absolute moderated t per
# gene. Its adjusted P value belongs to that selected probe; it is not a
# formally controlled gene-level FDR.

script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_flag) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_flag), mustWork = TRUE)
project_root <- dirname(dirname(script_path))

if (!requireNamespace("limma", quietly = TRUE)) {
  stop("The Bioconductor package 'limma' is required.")
}
suppressPackageStartupMessages(library(limma))

prepared_path <- file.path(project_root, "outputs", "PAH blood",
                           "gse33463_prepared.rds")
annotation_path <- file.path(project_root, "data", "processed",
                             "GPL6947.annot.gz")
output_dir <- file.path(project_root, "outputs", "PAH blood")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input <- readRDS(prepared_path)
expression <- input$expression
samples <- input$samples
if (!is.matrix(expression) || !is.numeric(expression) ||
    nrow(expression) != 48803L || ncol(expression) != 71L ||
    anyNA(expression) || !all(is.finite(expression)) ||
    anyDuplicated(rownames(expression)) ||
    !identical(colnames(expression), samples$sample_id) ||
    sum(samples$group == "Healthy") != 41L ||
    sum(samples$group == "IPAH") != 30L ||
    !identical(sort(unique(samples$group)), c("Healthy", "IPAH"))) {
  stop("Prepared matrix or IPAH/healthy sample labels differ from script 06.")
}

# With Healthy as the reference, the groupIPAH coefficient is the difference
# in mean log2 signal: IPAH minus Healthy. There is no paired sample or time
# variable in this cross-sectional model.
group <- factor(samples$group, levels = c("Healthy", "IPAH"))
design <- model.matrix(~ group)
rownames(design) <- samples$sample_id
if (!identical(colnames(design), c("(Intercept)", "groupIPAH"))) {
  stop("Unexpected model-matrix columns.")
}
design_export <- data.frame(
  sample_id = samples$sample_id,
  group = as.character(group),
  intercept = unname(design[, "(Intercept)"]),
  IPAH_indicator = unname(design[, "groupIPAH"])
)
write.csv(design_export, file.path(output_dir, "gse33463_limma_design.csv"),
          row.names = FALSE)

# limma fits a linear model to every probe. eBayes stabilizes the probe-wise
# residual variance estimates; it does not change the direction of logFC.
# BH correction here is across the 48,803 tested probes.
fit <- eBayes(lmFit(expression, design))
probe <- topTable(fit, coef = "groupIPAH", number = Inf,
                  sort.by = "none", adjust.method = "BH")
probe$probe_id <- rownames(probe)

# GPL6947 is the array's probe annotation. Multiple probes may map to one
# symbol, while symbols containing delimiters do not identify one clear gene.
annotation <- read.delim(gzfile(annotation_path), skip = 28L,
                         nrows = 48803L, check.names = FALSE, quote = "")
if (nrow(annotation) != 48803L ||
    !all(c("ID", "Gene symbol", "Gene ID") %in% names(annotation))) {
  stop("Unexpected GPL6947 annotation format.")
}
annotation <- annotation[, c("ID", "Gene symbol", "Gene ID")]
names(annotation) <- c("probe_id", "gene_symbol", "entrez_gene")
if (anyDuplicated(annotation$probe_id) ||
    !setequal(annotation$probe_id, rownames(expression))) {
  stop("GPL6947 probe IDs do not match the prepared expression matrix.")
}

probe <- merge(probe, annotation, by = "probe_id", all.x = TRUE,
               sort = FALSE)
names(probe)[names(probe) == "P.Value"] <- "p_value"
names(probe)[names(probe) == "adj.P.Val"] <- "adj_p_value"
probe$gene_symbol <- trimws(probe$gene_symbol)
probe$gene_symbol[probe$gene_symbol == "" |
                  grepl("///|//|;|,", probe$gene_symbol)] <- NA
probe <- probe[order(-probe$t, probe$probe_id), ]
probe <- probe[, c("probe_id", "gene_symbol", "entrez_gene", "logFC", "t",
                   "p_value", "adj_p_value", "AveExpr", "B")]
write.csv(probe, gzfile(file.path(output_dir,
                                  "pah_blood_probe_ranked.csv.gz")),
          row.names = FALSE)

# Collapse only after fitting and adjusting the probe tests. This makes a
# convenient one-row-per-symbol ranking, but the probe selected by |t| is
# chosen using the same outcome being tested.
gene <- probe[!is.na(probe$gene_symbol) & is.finite(probe$t), ]
gene <- gene[order(gene$gene_symbol, -abs(gene$t), gene$probe_id), ]
gene$source_probe_count <- as.integer(ave(gene$probe_id, gene$gene_symbol,
                                          FUN = length))
gene <- gene[!duplicated(gene$gene_symbol), ]
gene <- gene[order(-gene$t, gene$gene_symbol), ]
gene <- gene[, c("gene_symbol", "probe_id", "source_probe_count",
                 "entrez_gene", "logFC", "t", "p_value", "adj_p_value",
                 "AveExpr", "B")]
write.csv(gene, gzfile(file.path(output_dir,
                                 "pah_blood_gene_ranked.csv.gz")),
          row.names = FALSE)

qc <- data.frame(
  metric = c("tested_probes", "annotated_probes", "unique_genes",
             "genes_with_multiple_probes", "collapsed_extra_probes",
             "selected_genes_with_probe_bh_below_0.05"),
  value = c(nrow(probe), sum(!is.na(probe$gene_symbol)), nrow(gene),
            sum(gene$source_probe_count > 1L),
            sum(gene$source_probe_count - 1L),
            sum(gene$adj_p_value < 0.05))
)
write.csv(qc, file.path(output_dir, "gse33463_probe_mapping_qc.csv"),
          row.names = FALSE)

print(qc)
cat("Model coefficient: groupIPAH = IPAH minus Healthy\n")
cat("Probe-level BH values are not gene-level FDR values.\n")
