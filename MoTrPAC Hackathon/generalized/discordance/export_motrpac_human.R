#!/usr/bin/env Rscript

# Export package 2.0.8's human summary-level effects to the canonical long CSV
# read by the query core and build_discordance.py.
# Usage:
# Rscript export_motrpac_human.R --tissue muscle --contrast-category EE-CON \
#   --out reference_effects.csv.gz

argv <- commandArgs(trailingOnly = TRUE)
option <- function(name, default = NULL) {
  where <- match(name, argv)
  if (is.na(where)) return(default)
  if (where == length(argv) || startsWith(argv[where + 1L], "--")) {
    stop("Missing value for ", name)
  }
  argv[where + 1L]
}

tissue_arg <- option("--tissue", "muscle")
category_arg <- option("--contrast-category", "EE-CON")
type_arg <- option("--contrast-type", "exercise_with_controls")
assay_arg <- option("--assays", "rna,protein,phosphosite")
output_path <- option("--out")
if (is.null(output_path)) {
  stop("Usage: Rscript export_motrpac_human.R --tissue muscle|all --contrast-category EE-CON|all --assays rna,protein,phosphosite|all --out reference_effects.csv.gz")
}
parse_choice <- function(value, choices, name) {
  if (value == "all") return(choices)
  values <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  if (length(values) == 0L || any(!values %in% choices)) stop("Invalid ", name, ": ", value)
  unique(values)
}
tissues <- parse_choice(tissue_arg, c("muscle", "adipose", "blood"), "--tissue")
categories <- parse_choice(category_arg, c("EE-CON", "RE-CON", "EE-EE", "RE-RE", "EE-RE", "CON-CON"), "--contrast-category")
types <- parse_choice(type_arg, c("exercise_with_controls", "exercise_no_controls", "Endur_vs_Resist", "baseline", "control_only"), "--contrast-type")
if (assay_arg == "all") {
  selected_omes <- "all"
} else {
  assays <- parse_choice(assay_arg, c("rna", "protein", "phosphosite", "olink", "metabolite"), "--assays")
  ome_map <- list(
    rna = "transcript-rna-seq",
    protein = "prot-pr",
    phosphosite = "prot-ph",
    olink = "prot-ol",
    metabolite = "metab"
  )
  selected_omes <- unique(unlist(ome_map[assays], use.names = FALSE))
}

package_name <- "MotrpacHumanPreSuspensionAnalysis"
if (!requireNamespace(package_name, quietly = TRUE)) {
  stop("Install MotrpacHumanPreSuspensionAnalysis 2.0.8 first")
}
package_version <- as.character(utils::packageVersion(package_name))
if (package_version != "2.0.8") {
  stop("Expected MotrpacHumanPreSuspensionAnalysis 2.0.8; found ", package_version)
}
suppressPackageStartupMessages(library(MotrpacHumanPreSuspensionAnalysis))

results <- as.data.frame(load_differential_analysis(
  selected_tissues = tissues,
  selected_omes = selected_omes,
  single_matrix = TRUE,
  combine_with_featgene = TRUE,
  verbose = FALSE
))
results <- results[
  as.character(results$contrast_type) %in% types &
    as.character(results$contrast_category) %in% categories,
  , drop = FALSE
]
if (nrow(results) == 0L) stop("No rows found for the requested selection")

needed <- c(
  "tissue", "assay", "contrast_type", "contrast_category", "contrast_short",
  "Timepoint", "feature_id", "gene_symbol", "uniprot", "refmet_id",
  "refmet_name", "platform", "logFC",
  "p_value", "adj_p_value", "AveExpr", "z.std"
)
for (column in setdiff(needed, colnames(results))) results[[column]] <- NA
lower_column <- if ("CI.L_calculated" %in% colnames(results)) "CI.L_calculated" else "CI.L"
upper_column <- if ("CI.R_calculated" %in% colnames(results)) "CI.R_calculated" else "CI.R"
if (!lower_column %in% colnames(results) || !upper_column %in% colnames(results)) {
  warning("This package result has no confidence limits; equivalence categories will be unavailable")
  results$CI.L <- NA_real_
  results$CI.R <- NA_real_
  lower_column <- "CI.L"
  upper_column <- "CI.R"
}

layer <- c(
  "transcript-rna-seq" = "rna",
  "prot-pr" = "protein",
  "prot-ph" = "phosphosite",
  "prot-ol" = "protein",
  "metab" = "metabolite"
)[as.character(results$assay)]
source_gene_symbol <- as.character(results$gene_symbol)
source_uniprot <- as.character(results$uniprot)
valid_gene <- !is.na(source_gene_symbol) & nzchar(source_gene_symbol) &
  !grepl("[,;|]", source_gene_symbol)
valid_uniprot <- !is.na(source_uniprot) & nzchar(source_uniprot) &
  !grepl("[,;|]", source_uniprot)
gene_symbol <- ifelse(valid_gene, source_gene_symbol, "")
uniprot <- ifelse(valid_uniprot, source_uniprot, "")
source_refmet_id <- as.character(results$refmet_id)
source_refmet_name <- as.character(results$refmet_name)
valid_refmet_id <- !is.na(source_refmet_id) & nzchar(source_refmet_id) &
  !grepl("[,;|]", source_refmet_id)
valid_refmet_name <- !is.na(source_refmet_name) & nzchar(source_refmet_name) &
  !grepl("[;|]", source_refmet_name)
refmet_id <- ifelse(valid_refmet_id, source_refmet_id, "")
refmet_name <- ifelse(valid_refmet_name, source_refmet_name, "")
namespace <- ifelse(layer == "rna" & valid_gene, "HGNC",
                    ifelse(layer == "protein" & valid_uniprot, "UniProt",
                           ifelse(layer == "metabolite" & valid_refmet_id, "RefMet",
                                  ifelse(layer == "metabolite" & valid_refmet_name,
                                         "RefMetName", "MoTrPAC"))))
entity_id <- ifelse(layer == "rna" & valid_gene, gene_symbol,
                    ifelse(layer == "protein" & valid_uniprot, uniprot,
                           ifelse(layer == "metabolite" & valid_refmet_id, refmet_id,
                                  ifelse(layer == "metabolite" & valid_refmet_name,
                                         refmet_name, as.character(results$feature_id)))))
biospecimen <- c(
  "muscle" = "skeletal muscle",
  "adipose" = "subcutaneous adipose tissue"
)[as.character(results$tissue)]
biospecimen[as.character(results$tissue) == "blood" & layer == "rna"] <- "whole blood"
biospecimen[as.character(results$tissue) == "blood" & layer != "rna"] <- "plasma"

output <- data.frame(
  study_id = "MoTrPAC",
  species = "human",
  tissue = as.character(results$tissue),
  biospecimen = unname(biospecimen),
  layer = unname(layer),
  assay = as.character(results$assay),
  platform = as.character(results$platform),
  id_namespace = unname(namespace),
  entity_id = unname(entity_id),
  contrast_type = as.character(results$contrast_type),
  contrast_category = as.character(results$contrast_category),
  contrast_short = as.character(results$contrast_short),
  contrast = as.character(results$contrast_short),
  timepoint = as.character(results$Timepoint),
  feature_id = unname(entity_id),
  source_feature_id = as.character(results$feature_id),
  gene_symbol = unname(gene_symbol),
  source_gene_symbol = source_gene_symbol,
  uniprot = unname(uniprot),
  source_uniprot = source_uniprot,
  refmet_id = unname(refmet_id),
  refmet_name = unname(refmet_name),
  log2_fc = as.numeric(results$logFC),
  effect_scale = "log2 fold change",
  statistic = as.numeric(results$z.std),
  statistic_type = "z.std",
  ci_lower = as.numeric(results[[lower_column]]),
  ci_upper = as.numeric(results[[upper_column]]),
  p_value = as.numeric(results$p_value),
  q_value = as.numeric(results$adj_p_value),
  ave_expr = as.numeric(results$AveExpr),
  source_package_version = package_version,
  source_collection = "c2.0",
  stringsAsFactors = FALSE
)
if (anyNA(output$layer)) stop("Unknown assay in MoTrPAC loader output")

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
connection <- if (grepl("\\.gz$", output_path, ignore.case = TRUE)) {
  gzfile(output_path, open = "wt")
} else {
  file(output_path, open = "wt")
}
utils::write.csv(output, connection, row.names = FALSE, na = "")
close(connection)
cat("Wrote", nrow(output), "rows to", normalizePath(output_path, winslash = "/", mustWork = FALSE), "\n")
print(table(output$layer, output$timepoint))
