#!/usr/bin/env Rscript
# Cache the MSigDB C2:CGP gene sets used by make_disease_examples.py.
# The sets come from the msigdbr package that ships with the motrpac conda env
# (offline; no download). Output: tool/examples/sources/msigdbr_cgp_subset.csv
#
# Run:  ~/miniconda3/envs/motrpac/bin/Rscript tool/scripts/fetch_msigdbr_cgp.R
#
# MOOTHA_VOXPHOS is the signature used for the type 2 diabetes example.
# THUM_SYSTOLIC_HEART_FAILURE_* and KAYO_AGING_MUSCLE_* are cached only as
# reference sets for the sanity checks in make_disease_examples.py.
suppressMessages(library(msigdbr))

sets <- c("MOOTHA_VOXPHOS",
          "THUM_SYSTOLIC_HEART_FAILURE_UP", "THUM_SYSTOLIC_HEART_FAILURE_DN",
          "KAYO_AGING_MUSCLE_UP", "KAYO_AGING_MUSCLE_DN")

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
out <- file.path(dirname(here), "examples", "sources", "msigdbr_cgp_subset.csv")

x <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP")
x <- as.data.frame(x[x$gs_name %in% sets,
                     c("gs_name", "gene_symbol", "ncbi_gene", "gs_pmid", "gs_description")])
x <- x[!duplicated(x[, c("gs_name", "gene_symbol")]), ]
x$msigdbr_version <- as.character(packageVersion("msigdbr"))
x$db_version <- tryCatch(unique(msigdbr::msigdbr_collections()$db_version)[1],
                         error = function(e) NA)
missing <- setdiff(sets, unique(x$gs_name))
if (length(missing)) stop("sets not found in msigdbr: ", paste(missing, collapse = ", "))
write.csv(x, out, row.names = FALSE)
cat("wrote", out, "\n")
print(table(x$gs_name))
