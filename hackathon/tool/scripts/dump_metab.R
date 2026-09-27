# dump_metab.R -- dump MoTrPAC metabolomics differential-analysis (DA) summary statistics + feature maps
# to CSV for the mprobe METAB layer.
# Run from hackathon/tool/:  ~/miniconda3/envs/motrpac/bin/Rscript scripts/dump_metab.R
# (or via Slurm: sbatch -p cpu -w danilogin --mem=16G --time=00:30:00 -o /mnt/nfs/shared/%x-%j.out ...)
#
# Idempotent / write-once: a file that already exists is never rewritten (bump the _v suffix instead).
# Outputs (store/raw_metab/):
#   rat_metab_DA_v1.csv.gz             METAB_<TISSUE>_DA (MotrpacRatTraining6moData), all tissues, per-platform rows
#   rat_metab_DA_metareg_v1.csv.gz     METAB_<TISSUE>_DA_METAREG (repeated features merged by meta-regression)
#   rat_metab_feature_map_v1.csv.gz    METAB_FEATURE_ID_MAP + RefMet classes/IDs
#   human_metab_DA_v1.csv.gz           <TISSUE>_METAB_DA (MotrpacHumanPreSuspensionAnalysis), muscle/blood/adipose
#   human_metab_clinical_DA_v1.csv.gz  BLOOD_METAB_T_CLINICAL_DA (clinical chemistry panel; kept separate)
#   human_metab_feature_map_v1.csv.gz  DA feature_id -> refmet_name (+ METABOLOMICS_CVS, HUMAN_FEATURE_TO_GENE, RefMet)
#   refmet_classes_v1.csv.gz           RefMet bulk download (Metabolomics Workbench) subset to names used above

suppressMessages({
  library(data.table)
  library(MotrpacHumanPreSuspensionAnalysis)
})
out_dir <- "store/raw_metab"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
f_out <- function(x) file.path(out_dir, x)
rpkg <- "MotrpacRatTraining6moData"
hpkg <- "MotrpacHumanPreSuspensionAnalysis"

write_once <- function(x, f) {
  if (file.exists(f)) { message("EXISTS, not overwriting: ", f); return(invisible(FALSE)) }
  fwrite(x, f)
  Sys.chmod(f, "0444")
  message(sprintf("wrote %s  (%d rows x %d cols)", f, nrow(x), ncol(x)))
  invisible(TRUE)
}
get_obj <- function(name, pkg) {
  e <- new.env()
  data(list = name, package = pkg, envir = e)
  x <- as.data.table(e[[name]])
  for (v in names(x)) if (is.factor(x[[v]])) set(x, j = v, value = as.character(x[[v]]))
  x
}
cat("MotrpacRatTraining6moData", as.character(packageVersion(rpkg)),
    "| MotrpacHumanPreSuspensionAnalysis", as.character(packageVersion(hpkg)),
    "| MotrpacBicQC", as.character(packageVersion("MotrpacBicQC")),
    "| data.table", as.character(packageVersion("data.table")), "|", R.version.string, "\n")

# ------------------------------------------------------------------ rat DA
all_obj <- data(package = rpkg)$results[, "Item"]
da_objs <- grep("^METAB_.*_DA$", all_obj, value = TRUE)
mr_objs <- grep("^METAB_.*_DA_METAREG$", all_obj, value = TRUE)
rat_cols <- c("tissue", "assay", "dataset", "site", "is_targeted", "feature_ID", "feature", "metabolite_refmet",
              "metabolite", "sex", "comparison_group", "logFC", "logFC_se", "tscore", "zscore", "p_value",
              "adj_p_value", "selection_fdr", "numNAs", "comparison_average_intensity",
              "reference_average_intensity", "meta_reg_het_p", "meta_reg_pvalue")
load_rat <- function(objs) rbindlist(lapply(objs, function(o) {
  d <- get_obj(o, rpkg)
  cat(sprintf("%-26s %7d rows %6d feature_IDs  datasets=%d\n", o, nrow(d), uniqueN(d$feature_ID),
              uniqueN(d$dataset)))
  d <- d[, intersect(rat_cols, names(d)), with = FALSE]
  # stat_type: which test statistic the row carries (tscore = limma per platform, zscore = meta-regression)
  d[, stat_type := fifelse(dataset == "meta-reg", "zscore", "tscore")]
  d[, source_object := o]
  d
}), fill = TRUE)
rat_da <- load_rat(da_objs)
rat_mr <- load_rat(mr_objs)
cat(sprintf("## rat METAB_*_DA: %d tissues, %d rows | METAREG: %d rows\n", uniqueN(rat_da$tissue), nrow(rat_da),
            nrow(rat_mr)))

# ------------------------------------------------------------------ human DA
h_objs <- c("MUSCLE_METAB_DA", "BLOOD_METAB_DA", "ADIPOSE_METAB_DA")
hcols <- c("tissue", "assay", "platform", "contrast_short", "contrast_category", "contrast_type",
           "randomGroupCode", "Timepoint", "feature_id", "logFC", "CI.L_calculated", "CI.R_calculated", "t",
           "z.std", "degrees_of_freedom", "AveExpr", "p_value", "adj_p_value")
load_human <- function(objs) rbindlist(lapply(objs, function(o) {
  d <- get_obj(o, hpkg)
  cat(sprintf("%-26s %7d rows %6d features %3d contrasts\n", o, nrow(d), uniqueN(d$feature_id),
              uniqueN(d$contrast_short)))
  d <- d[, intersect(hcols, names(d)), with = FALSE]
  d[, source_object := o]
  d
}), fill = TRUE)
hum_da <- load_human(h_objs)
hum_clin <- load_human("BLOOD_METAB_T_CLINICAL_DA")

# ------------------------------------------------------------------ feature annotation inputs
rat_map <- get_obj("METAB_FEATURE_ID_MAP", rpkg)
cvs <- get_obj("METABOLOMICS_CVS", hpkg)
hf2g <- get_obj("HUMAN_FEATURE_TO_GENE", hpkg)[assay == "metab"]
hf2g <- hf2g[, which(vapply(hf2g, function(v) any(!is.na(v)), logical(1))), with = FALSE]

# ------------------------------------------------------------------ RefMet (Metabolomics Workbench bulk download)
refmet_file <- f_out("refmet_classes_v1.csv.gz")
names_needed <- sort(unique(na.omit(c(rat_map$metabolite_refmet, rat_da$metabolite_refmet, rat_mr$metabolite_refmet,
                                      hum_da$feature_id, as.character(cvs$refmet_name), hf2g$refmet_name))))
names_needed <- names_needed[names_needed != ""]
if (file.exists(refmet_file)) {
  refmet <- fread(refmet_file, colClasses = "character")
} else {
  url <- "https://www.metabolomicsworkbench.org/databases/refmet/refmet_download.php"
  tf <- tempfile(fileext = ".csv")
  ok <- tryCatch({ download.file(url, tf, quiet = TRUE, mode = "wb"); TRUE }, error = function(e) {
    message("RefMet download FAILED: ", conditionMessage(e)); FALSE })
  if (ok) {
    full <- fread(tf, colClasses = "character", encoding = "UTF-8")
    for (v in names(full)) set(full, which(full[[v]] == ""), v, NA_character_)
    unlink(tf)
    cat(sprintf("RefMet bulk download: %d rows\n", nrow(full)))
    # match by exact refmet_name; then case-insensitive; then (human) by refmet_id from HUMAN_FEATURE_TO_GENE
    q <- data.table(query_name = names_needed)
    q[full, on = .(query_name = refmet_name), refmet_id := i.refmet_id]
    q[!is.na(refmet_id), match_type := "exact_name"]
    full[, lc := tolower(refmet_name)]
    ci <- full[!duplicated(lc)]
    q[is.na(refmet_id), refmet_id := ci[match(tolower(query_name), lc), refmet_id]]
    q[!is.na(refmet_id) & is.na(match_type), match_type := "case_insensitive_name"]
    hid <- unique(hf2g[!is.na(refmet_id), .(refmet_name, refmet_id)])
    q[is.na(refmet_id), refmet_id := hid[match(query_name, refmet_name), refmet_id]]
    q[!is.na(refmet_id) & is.na(match_type), match_type := "refmet_id_from_HUMAN_FEATURE_TO_GENE"]
    q[!is.na(refmet_id) & !(refmet_id %in% full$refmet_id), `:=`(refmet_id = NA, match_type = NA)]
    full[, lc := NULL]
    refmet <- merge(q, full, by = "refmet_id", all.x = TRUE, sort = FALSE)
    setcolorder(refmet, c("query_name", "match_type", "refmet_id", "refmet_name"))
    refmet[, retrieved := format(Sys.Date())]
    refmet[, source_url := url]
    setorder(refmet, query_name)
    cat(sprintf("RefMet: %d query names, %d matched (%s)\n", nrow(refmet), sum(!is.na(refmet$refmet_id)),
                paste(names(table(refmet$match_type)), table(refmet$match_type), sep = "=", collapse = ", ")))
    write_once(refmet, refmet_file)
  } else {
    refmet <- NULL
  }
}

# BicQC local RefMet-derived dictionary (MoTrPAC named metabolites, revised 2021-05-05) as a fallback
bic_f <- system.file("extdata", "motrpac-metabolomics-named-revised-20210505.csv", package = "MotrpacBicQC")
bic <- unique(fread(bic_f, colClasses = "character")[, .(refmet_name, super_class, main_class, sub_class, formula,
                                                         pubchem_cid, kegg_id, hmdb_id, chebi_id, inchi_key)],
              by = "refmet_name")
ann_cols <- c("refmet_id", "super_class", "main_class", "sub_class", "formula", "pubchem_cid", "kegg_id", "hmdb_id",
              "chebi_id", "lipidmaps_id", "inchi_key")
annotate <- function(dt, name_col) {
  # adds RefMet columns (prefixed refmet_ann_ where needed) keyed on dt[[name_col]]
  a <- data.table(key_name = unique(na.omit(dt[[name_col]])))
  if (!is.null(refmet)) {
    r <- refmet[!is.na(refmet_id), c("query_name", "refmet_name", ann_cols), with = FALSE]
    setnames(r, "refmet_name", "refmet_name_current")
    a <- merge(a, r, by.x = "key_name", by.y = "query_name", all.x = TRUE)
    a[!is.na(refmet_id), class_source := "RefMet_bulk_download"]
  } else {
    for (v in c("refmet_name_current", ann_cols)) a[, (v) := NA_character_]
  }
  b <- bic[match(a$key_name, bic$refmet_name)]
  fill <- is.na(a$super_class) & !is.na(b$super_class)
  for (v in intersect(names(b), ann_cols)) a[fill, (v) := b[[v]][fill]]
  a[fill, class_source := "MotrpacBicQC_named_20210505"]
  clash <- intersect(setdiff(names(a), "key_name"), names(dt))
  if (length(clash)) setnames(a, clash, paste0("refmet_", clash))
  merge(dt, a, by.x = name_col, by.y = "key_name", all.x = TRUE, sort = FALSE)
}

# ------------------------------------------------------------------ write rat
write_once(rat_da, f_out("rat_metab_DA_v1.csv.gz"))
write_once(rat_mr, f_out("rat_metab_DA_metareg_v1.csv.gz"))
rm_f <- annotate(copy(rat_map), "metabolite_refmet")
setcolorder(rm_f, c("tissue", "dataset", "feature_ID_da", "metabolite_refmet"))
setorder(rm_f, tissue, dataset, feature_ID_da)
cat(sprintf("rat feature map: %d rows; DA (tissue,dataset,feature_ID) found in map: %.4f\n", nrow(rm_f),
            mean(paste(rat_da$tissue, rat_da$dataset, rat_da$feature_ID) %in%
                   paste(rm_f$tissue, rm_f$dataset, rm_f$feature_ID_da))))
write_once(rm_f, f_out("rat_metab_feature_map_v1.csv.gz"))

# ------------------------------------------------------------------ write human
write_once(hum_da, f_out("human_metab_DA_v1.csv.gz"))
write_once(hum_clin, f_out("human_metab_clinical_DA_v1.csv.gz"))
# DA feature_id is the RefMet name of the lowest-CV platform feature (see load_differential_analysis message);
# METABOLOMICS_CVS rows with lowest_CV == "yes" give the originating platform feature name.
hm <- unique(hum_da[, .(tissue, platform, feature_id)])
low <- cvs[lowest_CV == "yes", .(tissue, platform = assay, refmet_name = as.character(refmet_name),
                                 platform_feature_id = feature_id, site, feature_cv)]
low <- low[, .SD[1], by = .(tissue, platform, refmet_name)]
hm <- merge(hm, low, by.x = c("tissue", "platform", "feature_id"), by.y = c("tissue", "platform", "refmet_name"),
            all.x = TRUE, sort = FALSE)
hm[, refmet_name := fifelse(feature_id %in% as.character(cvs$refmet_name), feature_id, NA_character_)]
h2 <- unique(hf2g[, .(refmet_name, f2g_refmet_id = refmet_id, f2g_kegg_id = kegg_id)], by = "refmet_name")
hm <- merge(hm, h2, by = "refmet_name", all.x = TRUE, sort = FALSE)
hm <- annotate(hm, "refmet_name")
setcolorder(hm, c("tissue", "platform", "feature_id", "refmet_name"))
setorder(hm, tissue, platform, feature_id)
cat(sprintf("human feature map: %d rows, %d unique feature_id, refmet_name non-NA %.4f\n", nrow(hm),
            uniqueN(hm$feature_id), mean(!is.na(hm$refmet_name))))
write_once(hm, f_out("human_metab_feature_map_v1.csv.gz"))
cat("done\n")
