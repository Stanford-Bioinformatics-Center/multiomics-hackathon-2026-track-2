#!/usr/bin/env Rscript

# Export all measured MoTrPAC rat gastrocnemius protein contrasts. This
# exporter uses the bundled official objects, not the nine PAH-selected rows.
args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(2L, 3L))) {
  stop("Usage: Rscript export_rat_training_reference.R DATA_RAT_DIR OUTPUT.csv.gz [ORTHOLOG.csv]")
}
rat_dir <- normalizePath(args[1], mustWork = TRUE)
output_path <- args[2]

read_object <- function(name) {
  environment <- new.env(parent = emptyenv())
  loaded <- load(file.path(rat_dir, paste0(name, ".rda")), envir = environment)
  if (!identical(loaded, name)) stop("Unexpected object in ", name, ".rda")
  as.data.frame(environment[[name]])
}

unique_or_blank <- function(values) {
  values <- unique(as.character(values[!is.na(values) & nzchar(values)]))
  if (length(values) == 1L) values else ""
}

da <- read_object("PROT_SKMGN_DA")
feature_map <- read_object("FEATURE_TO_GENE_FILT")
orthologs <- read_object("RAT_TO_HUMAN_GENE")
needed <- c("feature_ID", "sex", "comparison_group", "logFC", "tscore",
            "p_value", "adj_p_value", "assay_code", "tissue")
if (!all(needed %in% names(da))) stop("Rat DA schema changed")
da <- da[as.character(da$assay_code) == "prot-pr" &
           as.character(da$tissue) == "SKM-GN", , drop = FALSE]
if (!nrow(da)) stop("No rat gastrocnemius protein contrasts")

feature_map$feature_ID <- as.character(feature_map$feature_ID)
feature_map$gene_symbol <- as.character(feature_map$gene_symbol)
symbol_by_feature <- vapply(split(feature_map$gene_symbol, feature_map$feature_ID),
                            unique_or_blank, character(1))
orthologs$RAT_SYMBOL <- as.character(orthologs$RAT_SYMBOL)
orthologs$HUMAN_ORTHOLOG_SYMBOL <- as.character(orthologs$HUMAN_ORTHOLOG_SYMBOL)
human_by_rat <- vapply(split(orthologs$HUMAN_ORTHOLOG_SYMBOL,
                             orthologs$RAT_SYMBOL), unique_or_blank, character(1))

original_feature <- as.character(da$feature_ID)
rat_symbol <- unname(symbol_by_feature[original_feature])
rat_symbol[is.na(rat_symbol)] <- ""
human_symbol <- unname(human_by_rat[rat_symbol])
human_symbol[is.na(human_symbol)] <- ""
sex <- as.character(da$sex)
week <- as.character(da$comparison_group)
if (!all(sex %in% c("female", "male")) ||
    !all(grepl("^[0-9]+w$", week))) {
  stop("Unexpected sex or training-week labels")
}

out <- data.frame(
  study_id = "MoTrPAC-rat-training",
  species = "rat",
  tissue = "muscle",
  biospecimen = "gastrocnemius",
  layer = "protein",
  feature_id = ifelse(nzchar(rat_symbol), rat_symbol, original_feature),
  id_namespace = ifelse(nzchar(rat_symbol), "RatGeneSymbol", "MoTrPAC"),
  log2_fc = da$logFC,
  statistic = da$tscore,
  p_value = da$p_value,
  q_value = da$adj_p_value,
  contrast = paste("trained minus sedentary", sex, week),
  contrast_category = paste0("TRAIN-SED-", sex),
  timepoint = paste0("week_", sub("w$", "", week)),
  gene_symbol = rat_symbol,
  human_ortholog_symbol = human_symbol,
  uniprot = "",
  refmet_id = "",
  refmet_name = "",
  source_feature_id = original_feature,
  source = "MotrpacRatTraining6moData 2.0.0; bundled PROT_SKMGN_DA",
  evidence_note = ifelse(nzchar(rat_symbol),
                         "Single rat symbol for this protein feature.",
                         "Ambiguous or missing feature-to-gene mapping."),
  stringsAsFactors = FALSE
)
if (anyDuplicated(out[, c("source_feature_id", "contrast_category", "timepoint")])) {
  stop("Duplicate rat feature within sex and training week")
}
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
connection <- gzfile(output_path, open = "wt")
write.csv(out, connection, row.names = FALSE, na = "")
close(connection)
cat("Exported", nrow(out), "rat protein rows;",
    sum(!nzchar(rat_symbol)), "rows lack unique rat-gene mapping.\n")

if (length(args) == 3L) {
  pair <- unique(orthologs[, c("HUMAN_ORTHOLOG_SYMBOL", "RAT_SYMBOL")])
  pair <- pair[!is.na(pair$HUMAN_ORTHOLOG_SYMBOL) &
                 !is.na(pair$RAT_SYMBOL) &
                 nzchar(pair$HUMAN_ORTHOLOG_SYMBOL) &
                 nzchar(pair$RAT_SYMBOL), , drop = FALSE]
  names(pair) <- c("source_gene_symbol", "target_gene_symbol")
  human_n <- table(pair$source_gene_symbol)
  rat_n <- table(pair$target_gene_symbol)
  pair <- pair[human_n[pair$source_gene_symbol] == 1L &
                 rat_n[pair$target_gene_symbol] == 1L, , drop = FALSE]
  pair$source_species <- "human"
  pair$target_species <- "rat"
  pair$mapping_source <- "MotrpacRatTraining6moData 2.0.0 RAT_TO_HUMAN_GENE"
  pair <- pair[, c("source_species", "source_gene_symbol", "target_species",
                   "target_gene_symbol", "mapping_source")]
  pair <- pair[order(pair$source_gene_symbol), , drop = FALSE]
  write.csv(pair, args[3], row.names = FALSE)
  cat("Exported", nrow(pair), "one-to-one human-to-rat ortholog rows.\n")
}
