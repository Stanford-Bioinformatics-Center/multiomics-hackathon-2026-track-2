#!/usr/bin/env Rscript

# Generic pre-ranked RNA pathway test. This intentionally uses a preset
# inter-gene correlation: sample-level covariance cannot be estimated here.

parse_args <- function(args) {
  if (length(args) == 0L || "--help" %in% args) {
    cat(paste0(
      "Usage: Rscript camera_pr.R --signature FULL_RANK.csv[.gz] ",
      "--gmt SETS.gmt[.gz] --universe-provenance universe.json ",
      "--gmt-namespace HGNC --pathway-namespace GO --collection C5/GOBP ",
      "--out-dir DIR [--score-field statistic|log2_fc] ",
      "[--min-size 10] [--max-size 500] [--inter-gene-cor 0.01]\n"))
    if (length(args) == 0L) stop("Required arguments are missing.")
    quit(save = "no", status = 0L)
  }
  if (length(args) %% 2L != 0L || any(!startsWith(args[seq(1L, length(args), 2L)], "--")))
    stop("Arguments must be --name value pairs.")
  keys <- sub("^--", "", args[seq(1L, length(args), 2L)])
  values <- args[seq(2L, length(args), 2L)]
  if (anyDuplicated(keys)) stop("Repeated CLI option.")
  options <- as.list(stats::setNames(values, keys))
  allowed <- c("signature", "gmt", "universe-provenance", "gmt-namespace",
               "pathway-namespace", "collection", "out-dir", "score-field",
               "min-size", "max-size", "inter-gene-cor")
  unknown <- setdiff(keys, allowed)
  if (length(unknown)) stop("Unknown CLI option(s): ", paste(unknown, collapse = ", "))
  required <- allowed[seq_len(7L)]
  missing <- setdiff(required, keys)
  if (length(missing)) stop("Missing CLI option(s): ", paste(missing, collapse = ", "))
  options
}

read_table <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  utils::read.csv(con, check.names = FALSE, stringsAsFactors = FALSE,
                  na.strings = c("", "NA", "NaN", "null"))
}

read_gmt <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  if (!length(lines)) stop("GMT contains no pathways.")
  fields <- strsplit(lines, "\t", fixed = TRUE)
  if (any(lengths(fields) < 3L)) stop("GMT line must contain name, description, and member(s).")
  names <- trimws(vapply(fields, `[[`, character(1L), 1L))
  if (any(!nzchar(names)) || anyDuplicated(toupper(names))) stop("GMT pathway names must be unique and nonempty.")
  members <- lapply(fields, function(parts) {
    ids <- toupper(trimws(parts[-c(1L, 2L)]))
    unique(ids[nzchar(ids)])
  })
  if (any(lengths(members) == 0L)) stop("GMT has an empty pathway.")
  descriptions <- vapply(fields, `[[`, character(1L), 2L)
  names(members) <- names
  list(members = members, description = descriptions)
}

one_value <- function(frame, field) {
  values <- unique(trimws(as.character(frame[[field]])))
  if (length(values) != 1L || is.na(values) || !nzchar(values))
    stop("Signature requires one nonempty ", field, " per run.")
  values
}

main <- function() {
  options <- parse_args(commandArgs(trailingOnly = TRUE))
  for (package in c("limma", "jsonlite", "digest"))
    if (!requireNamespace(package, quietly = TRUE)) stop("Install R package: ", package)
  score_field <- if (is.null(options[["score-field"]])) "statistic" else options[["score-field"]]
  if (!(score_field %in% c("statistic", "log2_fc"))) stop("Invalid --score-field.")
  min_size <- if (is.null(options[["min-size"]])) 10L else as.integer(options[["min-size"]])
  max_size <- if (is.null(options[["max-size"]])) 500L else as.integer(options[["max-size"]])
  correlation <- if (is.null(options[["inter-gene-cor"]])) 0.01 else as.numeric(options[["inter-gene-cor"]])
  if (is.na(min_size) || is.na(max_size) || min_size < 2L || max_size < min_size)
    stop("Require 2 <= min-size <= max-size.")
  if (!is.finite(correlation) || correlation < 0 || correlation >= 1)
    stop("inter-gene-cor must be in [0, 1).")
  if (!(options[["pathway-namespace"]] %in% c("GO", "BIOCARTA", "Pathway")) ||
      !nzchar(trimws(options[["collection"]])))
    stop("Declare GO/BIOCARTA/Pathway and a nonempty collection.")
  signature_path <- options[["signature"]]
  gmt_path <- options[["gmt"]]
  provenance_path <- options[["universe-provenance"]]
  signature <- read_table(signature_path)
  required <- c("study_id", "species", "tissue", "layer", "feature_id", "id_namespace",
                "log2_fc", "statistic", "p_value", "q_value", "contrast")
  if (length(setdiff(required, names(signature))))
    stop("Signature lacks canonical columns: ", paste(setdiff(required, names(signature)), collapse = ", "))
  if (!nrow(signature)) stop("Signature is empty.")
  context <- lapply(c("study_id", "species", "tissue", "layer", "id_namespace", "contrast"),
                    function(field) one_value(signature, field))
  names(context) <- c("study_id", "species", "tissue", "layer", "id_namespace", "contrast")
  if (context$layer != "rna") stop("cameraPR adapter requires full-rank RNA.")
  if (!(tolower(context$id_namespace) %in% c("hgnc", "ratgenesymbol", "rgd")))
    stop("Unsupported RNA id_namespace.")
  if (tolower(options[["gmt-namespace"]]) != tolower(context$id_namespace))
    stop("GMT namespace must equal signature id_namespace; map upstream.")
  biospecimen <- if ("biospecimen" %in% names(signature)) one_value(signature, "biospecimen") else ""
  identifiers <- toupper(trimws(as.character(signature$feature_id)))
  if (anyNA(identifiers) || any(!nzchar(identifiers)) || anyDuplicated(identifiers))
    stop("RNA feature_id must be nonempty and unique (case-insensitive).")
  scores <- suppressWarnings(as.numeric(signature[[score_field]]))
  if (any(!is.finite(scores))) stop("Every measured gene needs a finite signed ranking score.")
  if (length(scores) < min_size + 2L) stop("Universe is too small for requested set size.")
  names(scores) <- identifiers
  provenance <- jsonlite::fromJSON(provenance_path, simplifyVector = TRUE)
  fields <- c("source", "inclusion_rule", "feature_id_policy", "gmt_source", "gmt_version")
  if (!is.list(provenance) ||
      any(vapply(fields, function(x) is.null(provenance[[x]]) ||
                   is.na(provenance[[x]]) || !nzchar(trimws(as.character(provenance[[x]]))),
                 logical(1L)))) stop("Universe provenance lacks required source/inclusion/GMT fields.")
  if (!identical(provenance$complete_measured_universe, TRUE))
    stop("Universe provenance must declare complete_measured_universe=true.")
  if (!is.null(provenance$expected_feature_count) &&
      as.integer(provenance$expected_feature_count) != length(scores))
    stop("expected_feature_count does not match signature rows.")
  gmt <- read_gmt(gmt_path)
  index <- lapply(gmt$members, function(members) which(identifiers %in% members))
  overlap <- lengths(index)
  outside <- length(scores) - overlap
  status <- ifelse(overlap < min_size, "below_min_size",
                   ifelse(overlap > max_size, "above_max_size",
                          ifelse(outside < 2L, "outside_universe_too_small", "eligible")))
  eligible <- which(status == "eligible")
  empty <- rep(NA_real_, length(index))
  result <- data.frame(
    study_id = context$study_id, species = context$species, tissue = context$tissue,
    layer = "pathway", feature_id = names(gmt$members),
    id_namespace = options[["pathway-namespace"]], log2_fc = empty,
    statistic = empty, p_value = empty, q_value = empty,
    contrast = context$contrast, biospecimen = biospecimen,
    pathway_membership_version = as.character(provenance$gmt_version),
    pathway_collection = options[["collection"]],
    effect_scale = "cameraPR competitive signed Z equivalent",
    statistic_type = "cameraPR signed two-sided P-to-Z",
    score_field = score_field, n_set_members = lengths(gmt$members),
    n_overlap = overlap, n_universe = length(scores),
    gmt_description = gmt$description,
    camera_inter_gene_cor = correlation, camera_direction = "",
    q_value_scope = "BH over eligible sets in this GMT and contrast",
    test_status = status, stringsAsFactors = FALSE
  )
  if (length(eligible)) {
    camera <- limma::cameraPR(scores, index = index[eligible], use.ranks = FALSE,
                              inter.gene.cor = correlation, sort = FALSE)
    if (!identical(rownames(camera), names(index)[eligible]))
      stop("cameraPR returned unexpected pathway order.")
    p <- pmin(1, pmax(0, as.numeric(camera$PValue)))
    q <- stats::p.adjust(p, method = "BH")
    direction <- as.character(camera$Direction)
    z <- ifelse(direction == "Up", 1, -1) *
      stats::qnorm(pmax(p, .Machine$double.xmin) / 2, lower.tail = FALSE)
    result$statistic[eligible] <- z
    result$p_value[eligible] <- p
    result$q_value[eligible] <- q
    result$camera_direction[eligible] <- direction
    if (!identical(as.integer(camera$NGenes), as.integer(overlap[eligible])))
      stop("cameraPR set sizes disagree with overlap audit.")
  }
  out_dir <- options[["out-dir"]]
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(result, file.path(out_dir, "pathway_tests.csv"), row.names = FALSE, na = "")
  utils::write.csv(result[eligible, , drop = FALSE], file.path(out_dir, "pathway_signature.csv"),
                   row.names = FALSE, na = "")
  source_hash <- function(path) digest::digest(file = path, algo = "sha256")
  manifest <- list(
    signature = normalizePath(signature_path), signature_sha256 = source_hash(signature_path),
    gmt = normalizePath(gmt_path), gmt_sha256 = source_hash(gmt_path),
    universe_provenance = normalizePath(provenance_path),
    universe_provenance_sha256 = source_hash(provenance_path),
    user_universe_declaration = provenance,
    universe_verification = "declared by user; completeness cannot be proved from a rank table",
    source_context = context, gmt_member_namespace = options[["gmt-namespace"]],
    pathway_namespace = options[["pathway-namespace"]],
    pathway_collection = options[["collection"]], score_field = score_field,
    n_measured_genes = length(scores), n_gmt_sets = length(index),
    n_eligible_sets = length(eligible), min_size = min_size, max_size = max_size,
    test = "limma::cameraPR; two-sided; use.ranks=FALSE; signed P-to-Z output",
    limma_version = as.character(utils::packageVersion("limma")),
    inter_gene_correlation = correlation,
    inter_gene_correlation_basis = "preset sensitivity parameter; not estimated from sample data",
    q_value_scope = "BH over all eligible sets in this GMT and contrast",
    limits = c(
      "The preset correlation may not equal the true correlation in any gene set.",
      "Pre-ranked results cannot estimate set-specific correlation from sample-level expression.",
      "Signed P-to-Z is not a fold change and is not directly comparable to PAGE or other set scores."
    )
  )
  jsonlite::write_json(manifest, file.path(out_dir, "manifest.json"),
                       pretty = TRUE, auto_unbox = TRUE, na = "null")
  cat("cameraPR tested", length(eligible), "of", length(index),
      "sets against", length(scores), "measured genes\n")
}

main()
