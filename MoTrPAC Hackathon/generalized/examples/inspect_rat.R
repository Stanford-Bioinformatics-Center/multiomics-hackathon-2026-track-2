root <- commandArgs(trailingOnly = TRUE)[1]
for (name in c("PROT_SKMGN_DA", "FEATURE_TO_GENE_FILT", "RAT_TO_HUMAN_GENE")) {
  env <- new.env(parent = emptyenv())
  load(file.path(root, paste0(name, ".rda")), envir = env)
  x <- env[[name]]
  cat("\n", name, nrow(x), ncol(x), "\n", paste(names(x), collapse = ", "), "\n")
  print(utils::head(x[, seq_len(min(8L, ncol(x))), drop = FALSE], 2L))
}
