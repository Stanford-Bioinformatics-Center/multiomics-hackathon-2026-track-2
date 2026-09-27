# Exploratory targeted tests: all training weeks pooled and a matched 8-week
# endpoint sensitivity, with female, male, and sex-adjusted combined models.
dir.create('output', showWarnings = FALSE)
samples <- read.csv(file.path('output', '02_rat_protein_samples.csv'),
                    stringsAsFactors = FALSE)
mapping <- read.csv(file.path('output', '01_paper_to_rat_mapping.csv'),
                    stringsAsFactors = FALSE)
results <- do.call(rbind, lapply(mapping$paper_symbol, function(symbol) {
  x <- samples[samples$paper_symbol == symbol &
                 !is.na(samples$normalized_log2_protein), ]
  x$sex <- factor(x$sex, levels = c('female', 'male'))
  x$trained <- factor(ifelse(x$group == 'sedentary', 'untrained', 'trained'),
                      levels = c('untrained', 'trained'))
  make_row <- function(d, label, adjust_sex) {
    fit <- if (adjust_sex) lm(normalized_log2_protein ~ sex + trained, data = d)
           else lm(normalized_log2_protein ~ trained, data = d)
    z <- summary(fit)$coefficients['trainedtrained', ]
    data.frame(paper_symbol = symbol, comparison = label,
               n_trained = sum(d$trained == 'trained'),
               n_sedentary = sum(d$trained == 'untrained'),
               trained_minus_sedentary = unname(z['Estimate']),
               p_value = unname(z['Pr(>|t|)']))
  }
  eight <- x[x$group %in% c('sedentary', '8w'), ]
  rbind(make_row(x, 'all_weeks_combined', TRUE),
        make_row(eight, 'eight_week_combined', TRUE),
        make_row(eight[eight$sex == 'female', ], 'eight_week_female', FALSE),
        make_row(eight[eight$sex == 'male', ], 'eight_week_male', FALSE))
}))
results$BH_q_within_nine <- ave(results$p_value, results$comparison,
                                FUN = function(p) p.adjust(p, 'BH'))
rownames(results) <- NULL
write.csv(results, file.path('output', '04_training_comparisons.csv'),
          row.names = FALSE)
print(aggregate(BH_q_within_nine < 0.05 ~ comparison, results, sum), row.names = FALSE)
