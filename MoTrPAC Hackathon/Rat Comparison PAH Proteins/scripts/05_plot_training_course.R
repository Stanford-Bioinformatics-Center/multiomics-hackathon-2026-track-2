# The magenta line joins trained groups of DIFFERENT rats at weeks 1,2,4,8.
# Blue is the ONE measured sedentary 8-week endpoint, repeated as a benchmark.
dir.create('output', showWarnings = FALSE)
means <- read.csv(file.path('output', '03_sex_adjusted_group_means.csv'),
                  stringsAsFactors = FALSE)
mapping <- read.csv(file.path('output', '01_paper_to_rat_mapping.csv'),
                    stringsAsFactors = FALSE)
stopifnot(nrow(means) == 45L)
pink <- '#B33676'
blue <- '#2674A6'
png(file.path('output', '05_rat_protein_training_course.png'),
    width = 2400, height = 1800, res = 220)
par(mfrow = c(3, 3), mar = c(3.5, 4.0, 2.2, 1.0), oma = c(5.7, 2, 3, 0))
for (symbol in mapping$paper_symbol) {
  z <- means[means$paper_symbol == symbol, ]
  trained <- z[z$group != 'sedentary', ]
  trained <- trained[order(trained$measured_week), ]
  sedentary <- z[z$group == 'sedentary', ]
  stopifnot(nrow(trained) == 4L, nrow(sedentary) == 1L)
  yrange <- range(c(z$ci95_low, z$ci95_high))
  pad <- max(0.045, diff(yrange) * 0.12)
  plot(NA, xlim = c(1, 8), ylim = yrange + c(-pad, pad), xaxt = 'n',
       xlab = '', ylab = '', main = symbol, cex.main = 1.1)
  axis(1, at = c(1, 2, 4, 8))
  rect(1, sedentary$ci95_low, 8, sedentary$ci95_high,
       col = adjustcolor(blue, alpha.f = 0.10), border = NA)
  abline(h = sedentary$sex_adjusted_mean, col = blue, lwd = 2, lty = 2)
  points(8, sedentary$sex_adjusted_mean, pch = 19, col = blue, cex = 1.1)
  segments(trained$measured_week, trained$ci95_low,
           trained$measured_week, trained$ci95_high,
           col = adjustcolor(pink, alpha.f = 0.65), lwd = 1.7)
  lines(trained$measured_week, trained$sex_adjusted_mean,
        type = 'b', pch = 19, col = pink, lwd = 2)
  mtext('Normalized log2 protein', side = 2, line = 2.7, cex = 0.72)
}
mtext('Training duration (weeks)', side = 1, outer = TRUE, line = 0.4)
mtext('Trained: magenta line and 95% CIs     Sedentary: blue 8-week reference and CI band',
      side = 1, outer = TRUE, line = 2.0, cex = 0.82)
mtext('The blue reference was measured only at the 8-week endpoint; its flat line is repeated for comparison.',
      side = 1, outer = TRUE, line = 3.3, cex = 0.76)
mtext('Rat gastrocnemius proteins: sex-combined training course versus sedentary reference',
      side = 3, outer = TRUE, line = 0.8, cex = 1.2)
dev.off()
cat('Wrote output/05_rat_protein_training_course.png\n')
