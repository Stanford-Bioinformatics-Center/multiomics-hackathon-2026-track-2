# Fit normalized log2 protein ~ sex + group for each of the nine proteins.
# Equal female/male prediction weights give sex-adjusted group means.
dir.create('output', showWarnings = FALSE)
samples <- read.csv(file.path('output', '02_rat_protein_samples.csv'),
                    stringsAsFactors = FALSE)
mapping <- read.csv(file.path('output', '01_paper_to_rat_mapping.csv'),
                    stringsAsFactors = FALSE)
groups <- c('sedentary', '1w', '2w', '4w', '8w')
all_means <- list()
all_differences <- list()
for (i in seq_len(nrow(mapping))) {
  symbol <- mapping$paper_symbol[i]
  x <- samples[samples$paper_symbol == symbol, ]
  x$sex <- factor(x$sex, levels = c('female', 'male'))
  x$group <- factor(x$group, levels = groups)
  fit <- lm(normalized_log2_protein ~ sex + group, data = x)
  critical <- qt(0.975, df.residual(fit))
  all_means[[i]] <- do.call(rbind, lapply(groups, function(g) {
    new <- data.frame(sex = factor(c('female', 'male'), levels = levels(x$sex)),
                      group = factor(rep(g, 2), levels = groups))
    design <- model.matrix(delete.response(terms(fit)), new,
                           contrasts.arg = fit$contrasts)
    w <- colMeans(design)
    estimate <- sum(w * coef(fit))
    se <- sqrt(as.numeric(t(w) %*% vcov(fit) %*% w))
    data.frame(paper_symbol = symbol, group = g,
               measured_week = if (g == 'sedentary') 8L else
                 as.integer(sub('w$', '', g)),
               n_rats = sum(!is.na(x$normalized_log2_protein) & x$group == g),
               sex_adjusted_mean = estimate,
               ci95_low = estimate - critical * se,
               ci95_high = estimate + critical * se,
               sedentary_endpoint_reference = g == 'sedentary')
  }))
  coefficients <- summary(fit)$coefficients
  all_differences[[i]] <- do.call(rbind, lapply(groups[-1], function(g) {
    z <- coefficients[paste0('group', g), ]
    data.frame(paper_symbol = symbol, week = as.integer(sub('w$', '', g)),
               trained_minus_sedentary = unname(z['Estimate']),
               standard_error = unname(z['Std. Error']),
               p_value = unname(z['Pr(>|t|)']))
  }))
}
means <- do.call(rbind, all_means)
differences <- do.call(rbind, all_differences)
differences$BH_q_9_within_week <- ave(differences$p_value, differences$week,
                                     FUN = function(p) p.adjust(p, 'BH'))
patterns <- do.call(rbind, lapply(mapping$paper_symbol, function(symbol) {
  z <- means[means$paper_symbol == symbol & means$group != 'sedentary', ]
  z <- z[order(z$measured_week), ]
  x <- samples[samples$paper_symbol == symbol & samples$group != 'sedentary', ]
  fit <- lm(normalized_log2_protein ~ sex + measured_week,
            data = transform(x, sex = factor(sex, levels = c('female', 'male'))))
  slope <- summary(fit)$coefficients['measured_week', ]
  data.frame(paper_symbol = symbol,
             strict_increase_1_to_8w = all(diff(z$sex_adjusted_mean) > 0),
             strict_decrease_1_to_8w = all(diff(z$sex_adjusted_mean) < 0),
             mean_change_8w_minus_1w = z$sex_adjusted_mean[4] - z$sex_adjusted_mean[1],
             trained_only_slope_per_week = unname(slope['Estimate']),
             linear_trend_p = unname(slope['Pr(>|t|)']))
}))
patterns$linear_trend_BH_q_9 <- p.adjust(patterns$linear_trend_p, 'BH')
rownames(means) <- NULL
rownames(differences) <- NULL
rownames(patterns) <- NULL
stopifnot(nrow(means) == 45L, nrow(differences) == 36L)
write.csv(means, file.path('output', '03_sex_adjusted_group_means.csv'),
          row.names = FALSE)
write.csv(differences, file.path('output', '03_training_vs_sedentary_by_week.csv'),
          row.names = FALSE)
write.csv(patterns, file.path('output', '03_training_duration_patterns.csv'),
          row.names = FALSE)
cat('Wrote 45 group means, 36 weekly contrasts, and nine duration patterns.\n')
