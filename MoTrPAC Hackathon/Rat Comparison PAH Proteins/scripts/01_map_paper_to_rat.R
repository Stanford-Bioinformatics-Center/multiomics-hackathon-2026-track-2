# Run from the Rat Comparison PAH Proteins folder with Rscript.
# Table 2 paper identifiers -> current human symbols -> rat ortholog -> one
# gastrocnemius protein feature. Abort rather than silently choose a duplicate.
dir.create('output', showWarnings = FALSE)
read_object <- function(name) {
  e <- new.env(parent = emptyenv())
  load(file.path('data', 'rat', paste0(name, '.rda')), envir = e)
  e[[name]]
}
pah <- read.csv(file.path('data', 'paper', 'pah_lower_proteins_malenfant2015.csv'),
                stringsAsFactors = FALSE)
stopifnot(nrow(pah) == 9L, all(pah$pah_to_control_ratio < 1),
          !anyDuplicated(pah$paper_symbol), !anyDuplicated(pah$uniprot_accession))
# Two paper-era gene names have changed. Keep this crosswalk explicit.
current <- c(NDUFA9 = 'NDUFA9', UQCRC2 = 'UQCRC2', UQCRC1 = 'UQCRC1',
             ATP5L = 'ATP5MG', ATP5B = 'ATP5F1B', IDH2 = 'IDH2',
             OGDH = 'OGDH', SLC25A4 = 'SLC25A4', ECH1 = 'ECH1')
stopifnot(identical(sort(pah$paper_symbol), sort(names(current))))
orthologs <- read_object('RAT_TO_HUMAN_GENE')
feature_map <- read_object('FEATURE_TO_GENE_FILT')
abundance <- read_object('PROT_SKMGN_NORM_DATA')
da <- read_object('PROT_SKMGN_DA')

mapping <- do.call(rbind, lapply(seq_len(nrow(pah)), function(i) {
  human <- unname(current[pah$paper_symbol[i]])
  rat <- unique(na.omit(orthologs$RAT_SYMBOL[
    orthologs$HUMAN_ORTHOLOG_SYMBOL == human]))
  if (length(rat) != 1L) stop('Expected one rat ortholog for ', human)
  feature <- unique(na.omit(feature_map$feature_ID[
    feature_map$gene_symbol == rat &
      feature_map$feature_ID %in% abundance$feature_ID &
      feature_map$feature_ID %in% da$feature_ID]))
  if (length(feature) != 1L) stop('Expected one measured rat protein feature for ', human)
  data.frame(uniprot_accession = pah$uniprot_accession[i],
             paper_symbol = pah$paper_symbol[i],
             current_human_symbol = human, rat_symbol = rat,
             rat_feature_ID = feature,
             pah_to_control_ratio = pah$pah_to_control_ratio[i],
             reported_p_value_text = pah$reported_p_value_text[i],
             source_doi = pah$source_doi[i], stringsAsFactors = FALSE)
}))
write.csv(mapping, file.path('output', '01_paper_to_rat_mapping.csv'), row.names = FALSE)

# These are MoTrPAC's published, sex-specific trained-minus-sedentary results.
official <- do.call(rbind, lapply(seq_len(nrow(mapping)), function(i) {
  x <- da[da$feature_ID == mapping$rat_feature_ID[i],
          c('feature_ID', 'sex', 'comparison_group', 'logFC', 'logFC_se',
            'p_value', 'adj_p_value', 'selection_fdr')]
  if (nrow(x) != 8L) stop('Expected four weeks x two sexes for ', mapping$paper_symbol[i])
  x$paper_symbol <- mapping$paper_symbol[i]
  x
}))
official$week <- as.integer(sub('w$', '', official$comparison_group))
official <- official[order(match(official$paper_symbol, mapping$paper_symbol),
                           official$sex, official$week), ]
rownames(official) <- NULL
write.csv(official, file.path('output', '01_official_sex_specific_weekly.csv'),
          row.names = FALSE)
cat('Mapped', nrow(mapping), 'paper proteins to unique rat features; exported',
    nrow(official), 'official sex-specific weekly contrasts.\n')
