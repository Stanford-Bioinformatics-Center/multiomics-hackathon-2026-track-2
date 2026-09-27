# Match each normalized gastrocnemius protein value to one rat and training group.
dir.create('output', showWarnings = FALSE)
read_object <- function(name) {
  e <- new.env(parent = emptyenv())
  load(file.path('data', 'rat', paste0(name, '.rda')), envir = e)
  e[[name]]
}
mapping <- read.csv(file.path('output', '01_paper_to_rat_mapping.csv'),
                    stringsAsFactors = FALSE)
abundance <- read_object('PROT_SKMGN_NORM_DATA')
pheno <- read_object('PHENO')
sample_ids <- names(abundance)[-(1:4)]
meta <- pheno[match(sample_ids, pheno$viallabel),
              c('viallabel', 'pid', 'registration.sex', 'key.intervention',
                'key.sacrificetime')]
if (anyNA(meta$viallabel) || anyDuplicated(meta$pid)) stop('Sample/rat match failed.')
meta$group <- ifelse(meta$key.intervention == 'control', 'sedentary',
                     sub('^([1248]) weeks? of training.*$', '\\1w',
                         meta$key.sacrificetime))
if (!all(meta$group %in% c('sedentary', '1w', '2w', '4w', '8w')) ||
    !all(grepl('^8 weeks', meta$key.sacrificetime[meta$group == 'sedentary'])))
  stop('Training group or sedentary endpoint check failed.')
measured_week <- rep(8L, nrow(meta))
trained <- meta$group != 'sedentary'
measured_week[trained] <- as.integer(sub('w$', '', meta$group[trained]))

samples <- do.call(rbind, lapply(seq_len(nrow(mapping)), function(i) {
  row <- abundance[abundance$feature_ID == mapping$rat_feature_ID[i], ]
  if (nrow(row) != 1L) stop('Missing or duplicate abundance feature.')
  data.frame(paper_symbol = mapping$paper_symbol[i],
             rat_feature_ID = mapping$rat_feature_ID[i],
             sample_id = sample_ids, rat_id = meta$pid,
             sex = meta$registration.sex, group = meta$group,
             measured_week = measured_week,
             normalized_log2_protein = as.numeric(row[1, sample_ids]),
             stringsAsFactors = FALSE)
}))
stopifnot(nrow(samples) == 9L * 57L)
write.csv(samples, file.path('output', '02_rat_protein_samples.csv'),
          row.names = FALSE, na = '')
counts <- as.data.frame(with(meta, table(sex = registration.sex, group = group)))
names(counts)[3] <- 'n_rats'
write.csv(counts, file.path('output', '02_sample_counts.csv'), row.names = FALSE)
print(counts, row.names = FALSE)
