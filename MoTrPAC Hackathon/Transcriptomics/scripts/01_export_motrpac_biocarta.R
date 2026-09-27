#!/usr/bin/env Rscript

# Export the official, precomputed c2.0 CAMERA results. No gene-set test is
# recomputed here. The paper's PAGE values are a different statistic.
args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(args) != 1L) stop("Run this file with Rscript.")
script <- normalizePath(sub("^--file=", "", args), mustWork = TRUE)
module <- dirname(dirname(script))
rda <- file.path(module, "data", "motrpac", "CAMERA_RESULTS.rda")
if (!file.exists(rda)) stop("CAMERA_RESULTS.rda not found: ", rda)
description <- file.path(module, "data", "motrpac", "DESCRIPTION")
if (!file.exists(description)) stop("2.0.8 DESCRIPTION not found: ", description)
version <- read.dcf(description, fields = "Version")[1, 1]
if (version != "2.0.8") stop("Expected source version 2.0.8; got ", version)

paper <- read.csv(file.path(module, "outputs", "paper_table2_biocarta.csv"),
                  check.names = FALSE)
if (nrow(paper) != 6L || anyDuplicated(paper$set)) {
  stop("Expected six distinct Table 2 BioCarta pathways.")
}
env <- new.env(parent = emptyenv())
loaded <- load(rda, envir = env)
if (!identical(loaded, "CAMERA_RESULTS")) {
  stop("RDA must contain only CAMERA_RESULTS.")
}
camera <- as.data.frame(env$CAMERA_RESULTS)
needed <- c("tissue", "assay", "contrast_type", "contrast_short",
            "collection", "database", "set", "set_size", "direction",
            "z.std", "p_value", "adj_p_value")
if (!all(needed %in% names(camera))) stop("CAMERA_RESULTS schema changed.")

keep <- as.character(camera$tissue) == "blood" &
  as.character(camera$assay) == "transcript-rna-seq" &
  as.character(camera$contrast_type) == "exercise_with_controls" &
  as.character(camera$collection) == "C2" &
  as.character(camera$database) == "BIOCARTA" &
  as.character(camera$set) %in% paper$set
x <- camera[keep, needed, drop = FALSE]
x[] <- lapply(x, function(col) if (is.factor(col)) as.character(col) else col)
x$modality <- ifelse(startsWith(x$contrast_short, "Endur."), "EE",
                     ifelse(startsWith(x$contrast_short, "Resist."), "RE", NA))
x$timepoint <- sub("^(Endur|Resist)\\.([^ ]+).*", "\\2", x$contrast_short)
ee <- c("during_20_min", "during_40_min", "post_10_min",
        "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")
re <- ee[3:6]
if (nrow(x) != 60L || anyNA(x$modality) ||
    anyDuplicated(x[c("set", "modality", "timepoint")]) ||
    !setequal(x$set, paper$set) ||
    !setequal(unique(x$timepoint[x$modality == "EE"]), ee) ||
    !setequal(unique(x$timepoint[x$modality == "RE"]), re) ||
    any(!is.finite(x$z.std)) || any(!is.finite(x$adj_p_value)) ||
    any(x$adj_p_value < 0 | x$adj_p_value > 1)) {
  stop("Unexpected BioCarta pathway rows or statistics.")
}
x <- x[order(match(x$set, paper$set), match(x$modality, c("EE", "RE")),
             match(x$timepoint, ee)), ]
x$source_package_version <- version
dir.create(file.path(module, "outputs"), showWarnings = FALSE, recursive = TRUE)
write.csv(x, file.path(module, "outputs", "motrpac_biocarta_camera.csv"),
          row.names = FALSE)
cat("Exported", nrow(x), "rows: 6 sets x (6 EE + 4 RE) time points.\n")
