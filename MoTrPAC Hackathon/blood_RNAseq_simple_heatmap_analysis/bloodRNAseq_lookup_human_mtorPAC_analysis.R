library(MotrpacHumanPreSuspensionAnalysis)
library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)

# given a set of genes or proteins, query into Motrpac
geneset_fname <- 'BIOCARTA_TCR_PATHWAY.v2026.txt'
# BIOCARTA_CCR5_PATHWAY.v2026.1.txt
# BIOCARTA_IL2_PATHWAY.v2026.txt
# BIOCARTA_TCR_PATHWAY.v2026.txt


set <- readLines(geneset_fname)
tissue_type <- 'blood'


# get corresponding RNA transcripts that are present in MoTrPAC
rna_results <- load_differential_analysis(
  selected_tissues = tissue_type,
  selected_omes = 'transcript-rna-seq',
  single_matrix = TRUE,
  combine_with_featgene = TRUE,
  verbose = FALSE
)

rna_results <- rna_results %>% filter(contrast_category %in% c('RE-CON', 'EE-CON')) # we only look at exercise v/s control
rna_results <- rna_results %>% filter(gene_symbol %in% set)
rna_results <- rna_results %>% select(c('gene_symbol', 'contrast_category', 'Timepoint', 'logFC', 'p_value', 'adj_p_value' ))

# Further filter columns to transform the data from longform to wide table
value_to_plot <- 'adj_p_value' # can be logFC, p_value, or adj_p_value
rna_results <- rna_results %>% select(c('gene_symbol', 'contrast_category', 'Timepoint', value_to_plot ))

rna_results_wide <- pivot_wider(rna_results, names_from = c('contrast_category', 'Timepoint'), values_from = value_to_plot)
rna_results_wide <- rna_results_wide %>% relocate(`EE-CON_pre_exercise`, .after = gene_symbol)
rna_results_wide <- rna_results_wide %>% relocate(`RE-CON_pre_exercise`, .after = `EE-CON_post_24_hr`)

rownames(rna_results_wide) <- NULL
rna_results_wide <- column_to_rownames(rna_results_wide, var='gene_symbol')

# Plot heatmap of differential expression results (padj, logFC, or nominal p-value)
plot_title <- strsplit(geneset_fname, '[.]')[[1]][1]
Heatmap(as.matrix(rna_results_wide),
        name = value_to_plot,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_title = plot_title
        )

