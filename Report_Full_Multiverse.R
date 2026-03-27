#### Full Multiverse Report — Comprehensive Analysis & Plots ####
#
# Reads multiverse_full_raw.csv and generates:
#   - One or more plots per variable (dataset + ReReReRe)
#   - Interaction plots for key variable pairs
#   - A text report with summary statistics
#
# Output: archive/multiverse_full/report/ (all plots + report.txt)
# ===

setwd("C:/Users/vitto/Desktop/ReReReRe")

# --- Load data ---
raw <- read.csv("archive/multiverse_full/multiverse_full_raw.csv")
best <- read.csv("archive/multiverse_full/multiverse_full_best.csv")

outdir <- "archive/multiverse_full/report"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("Loaded %d raw rows, %d best rows\n", nrow(raw), nrow(best)))

# --- Color palettes ---
pal_nf <- colorRampPalette(c("#fee8c8", "#fdbb84", "#e34a33", "#67000d"))(7)
pal_ipf <- c("3" = "#e34a33", "6" = "steelblue", "10" = "#31a354")
pal_n <- c("100" = "#fc8d59", "300" = "#91bfdb", "500" = "#4575b4")
pal_pct <- c("0.05" = "#1b9e77", "0.1" = "#d95f02", "0.25" = "#7570b3")
pal_cp <- c("0.03" = "#e7298a", "0.05" = "#66a61e", "0.1" = "#e6ab02")

NF_LEVELS <- c(4, 8, 12, 16, 20, 25, 30)
IPF_LEVELS <- c(3, 6, 10)
N_LEVELS <- c(100, 300, 500)
PCT_LEVELS <- c(0.05, 0.10, 0.25)
CP_LEVELS <- c(0.03, 0.05, 0.10)
Z_RANGE <- seq(0.1, 3.0, by = 0.2)

# --- Helper: save plot ---
save_plot <- function(name, width = 1400, height = 900, res = 150) {
  png(file.path(outdir, paste0(name, ".png")), width = width, height = height, res = res)
}

# --- Text report sink ---
report_file <- file.path(outdir, "report.txt")

report <- file(report_file, open = "w")
writeLines(paste(rep("=", 80), collapse = ""), report)
writeLines("FULL MULTIVERSE REPORT — ReReReRe Parameter Analysis", report)
writeLines(paste("Generated:", Sys.time()), report)
writeLines(paste(rep("=", 80), collapse = ""), report)
writeLines("", report)
writeLines(sprintf("Total raw rows: %d", nrow(raw)), report)
writeLines(sprintf("Total best rows: %d", nrow(best)), report)
writeLines(sprintf("Replications per cell: %d", max(raw$rep_id)), report)
writeLines(sprintf("Dataset conditions: %d", nrow(best) / max(raw$rep_id) / length(CP_LEVELS)), report)
writeLines("", report)

# ============================================================
# PLOT 1: MCC by nFactors (main effect, oracle best)
# ============================================================
agg_nf <- aggregate(mcc ~ nFactors, data = best, FUN = function(x) c(mean = mean(x), sd = sd(x), median = median(x)))
agg_nf <- data.frame(nFactors = agg_nf$nFactors, mean = agg_nf$mcc[,1], sd = agg_nf$mcc[,2], median = agg_nf$mcc[,3])

save_plot("01_mcc_by_nFactors")
par(mar = c(5, 5, 4, 2))
plot(agg_nf$nFactors, agg_nf$mean, type = "b", pch = 16, col = "#e34a33", lwd = 3, cex = 1.5,
     xlab = "Number of Factors (nF)", ylab = "Mean MCC (oracle best z)",
     main = "D1: Effect of nFactors on Detection Performance",
     ylim = c(0, max(agg_nf$mean + agg_nf$sd) * 1.1))
arrows(agg_nf$nFactors, agg_nf$mean - agg_nf$sd, agg_nf$nFactors, agg_nf$mean + agg_nf$sd,
       angle = 90, code = 3, length = 0.05, col = "#e34a33", lwd = 1.5)
lines(agg_nf$nFactors, agg_nf$median, type = "b", pch = 1, col = "gray50", lwd = 1.5, lty = 2)
legend("topleft", legend = c("Mean +/- SD", "Median"),
       pch = c(16, 1), col = c("#e34a33", "gray50"), lwd = c(3, 1.5), lty = c(1, 2), bg = "white")
abline(h = seq(0, 0.6, 0.1), col = "gray90", lty = 3)
mtext("Averaged over all ipf, n, pct, corProp | R=10 reps", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()

writeLines("--- D1: nFactors ---", report)
writeLines(sprintf("  nF=%2d: mean=%.3f, sd=%.3f, median=%.3f", agg_nf$nFactors, agg_nf$mean, agg_nf$sd, agg_nf$median), report)
writeLines("", report)

# ============================================================
# PLOT 2: MCC by items_per_factor (main effect)
# ============================================================
agg_ipf <- aggregate(mcc ~ items_per_factor, data = best, FUN = function(x) c(mean = mean(x), sd = sd(x)))
agg_ipf <- data.frame(ipf = agg_ipf$items_per_factor, mean = agg_ipf$mcc[,1], sd = agg_ipf$mcc[,2])

save_plot("02_mcc_by_items_per_factor")
par(mar = c(5, 5, 4, 2))
bp <- barplot(agg_ipf$mean, names.arg = agg_ipf$ipf, col = unname(pal_ipf),
              xlab = "Items per Factor", ylab = "Mean MCC (oracle best z)",
              main = "D2: Effect of Items per Factor on Detection",
              ylim = c(0, max(agg_ipf$mean + agg_ipf$sd) * 1.15), border = NA)
arrows(bp, agg_ipf$mean - agg_ipf$sd, bp, agg_ipf$mean + agg_ipf$sd,
       angle = 90, code = 3, length = 0.05, lwd = 1.5)
text(bp, agg_ipf$mean + agg_ipf$sd + 0.01, sprintf("%.3f", agg_ipf$mean), cex = 0.9, font = 2)
mtext("Averaged over all nF, n, pct, corProp | R=10 reps", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()

writeLines("--- D2: items_per_factor ---", report)
writeLines(sprintf("  ipf=%2d: mean=%.3f, sd=%.3f", agg_ipf$ipf, agg_ipf$mean, agg_ipf$sd), report)
writeLines("", report)

# ============================================================
# PLOT 3: MCC by n_respondents (main effect)
# ============================================================
agg_n <- aggregate(mcc ~ n_respondents, data = best, FUN = function(x) c(mean = mean(x), sd = sd(x)))
agg_n <- data.frame(n = agg_n$n_respondents, mean = agg_n$mcc[,1], sd = agg_n$mcc[,2])

save_plot("03_mcc_by_n_respondents")
par(mar = c(5, 5, 4, 2))
plot(agg_n$n, agg_n$mean, type = "b", pch = 16, col = "#4575b4", lwd = 3, cex = 1.5,
     xlab = "Sample Size (n)", ylab = "Mean MCC (oracle best z)",
     main = "D3: Effect of Sample Size on Detection Performance",
     ylim = c(0, max(agg_n$mean + agg_n$sd) * 1.1))
arrows(agg_n$n, agg_n$mean - agg_n$sd, agg_n$n, agg_n$mean + agg_n$sd,
       angle = 90, code = 3, length = 0.05, col = "#4575b4", lwd = 1.5)
text(agg_n$n, agg_n$mean + agg_n$sd + 0.01, sprintf("%.3f", agg_n$mean), cex = 0.85, font = 2)
abline(h = seq(0, 0.6, 0.1), col = "gray90", lty = 3)
mtext("Averaged over all nF, ipf, pct, corProp | R=10 reps", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()

writeLines("--- D3: n_respondents ---", report)
writeLines(sprintf("  n=%4d: mean=%.3f, sd=%.3f", agg_n$n, agg_n$mean, agg_n$sd), report)
writeLines("", report)

# ============================================================
# PLOT 4: MCC by pct_careless (main effect)
# ============================================================
agg_pct <- aggregate(mcc ~ pct_careless, data = best, FUN = function(x) c(mean = mean(x), sd = sd(x)))
agg_pct <- data.frame(pct = agg_pct$pct_careless, mean = agg_pct$mcc[,1], sd = agg_pct$mcc[,2])

save_plot("04_mcc_by_pct_careless")
par(mar = c(5, 5, 4, 2))
bp <- barplot(agg_pct$mean, names.arg = paste0(agg_pct$pct * 100, "%"),
              col = unname(pal_pct), border = NA,
              xlab = "% Careless Respondents", ylab = "Mean MCC (oracle best z)",
              main = "D4: Effect of Careless Base Rate on Detection",
              ylim = c(0, max(agg_pct$mean + agg_pct$sd) * 1.15))
arrows(bp, agg_pct$mean - agg_pct$sd, bp, agg_pct$mean + agg_pct$sd,
       angle = 90, code = 3, length = 0.05, lwd = 1.5)
text(bp, agg_pct$mean + agg_pct$sd + 0.01, sprintf("%.3f", agg_pct$mean), cex = 0.9, font = 2)
mtext("Averaged over all nF, ipf, n, corProp | R=10 reps", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()

writeLines("--- D4: pct_careless ---", report)
writeLines(sprintf("  pct=%.0f%%: mean=%.3f, sd=%.3f", agg_pct$pct * 100, agg_pct$mean, agg_pct$sd), report)
writeLines("", report)

# ============================================================
# PLOT 5: MCC by corProp (main effect)
# ============================================================
agg_cp <- aggregate(mcc ~ corProp, data = best, FUN = function(x) c(mean = mean(x), sd = sd(x)))
agg_cp <- data.frame(cp = agg_cp$corProp, mean = agg_cp$mcc[,1], sd = agg_cp$mcc[,2])

save_plot("05_mcc_by_corProp")
par(mar = c(5, 5, 4, 2))
bp <- barplot(agg_cp$mean, names.arg = agg_cp$cp,
              col = unname(pal_cp), border = NA,
              xlab = "corProp", ylab = "Mean MCC (oracle best z)",
              main = "R1: Effect of corProp on Detection Performance",
              ylim = c(0, max(agg_cp$mean + agg_cp$sd) * 1.15))
arrows(bp, agg_cp$mean - agg_cp$sd, bp, agg_cp$mean + agg_cp$sd,
       angle = 90, code = 3, length = 0.05, lwd = 1.5)
text(bp, agg_cp$mean + agg_cp$sd + 0.01, sprintf("%.3f", agg_cp$mean), cex = 0.9, font = 2)
mtext("Averaged over all nF, ipf, n, pct | R=10 reps | Oracle best z", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()

writeLines("--- R1: corProp ---", report)
writeLines(sprintf("  corProp=%.2f: mean=%.3f, sd=%.3f", agg_cp$cp, agg_cp$mean, agg_cp$sd), report)
writeLines("", report)

# ============================================================
# PLOT 6: MCC by z_threshold (main effect, across all conditions)
# ============================================================
agg_z <- aggregate(mcc ~ z_threshold, data = raw, FUN = function(x) c(mean = mean(x), sd = sd(x)))
agg_z <- data.frame(z = agg_z$z_threshold, mean = agg_z$mcc[,1], sd = agg_z$mcc[,2])

save_plot("06_mcc_by_z_threshold")
par(mar = c(5, 5, 4, 2))
plot(agg_z$z, agg_z$mean, type = "b", pch = 16, col = "#7570b3", lwd = 3, cex = 1.3,
     xlab = "z_threshold", ylab = "Mean MCC",
     main = "R2: Effect of z_threshold on Detection (all conditions)",
     ylim = c(0, max(agg_z$mean) * 1.3))
abline(v = 1.5, col = "red", lty = 2, lwd = 1.5)
text(1.5, max(agg_z$mean) * 1.15, "z=1.5\n(recommended)", col = "red", cex = 0.8)
abline(h = seq(0, 0.3, 0.05), col = "gray90", lty = 3)
mtext("Averaged over all nF, ipf, n, pct, corProp | R=10 reps", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()

writeLines("--- R2: z_threshold ---", report)
writeLines(sprintf("  z=%.1f: mean=%.3f", agg_z$z, agg_z$mean), report)
writeLines("", report)

# ============================================================
# PLOT 7: z_threshold curves by nFactors (interaction D1 x R2)
# ============================================================
save_plot("07_z_curves_by_nFactors", width = 1600, height = 900)
par(mar = c(5, 5, 4, 8), xpd = TRUE)

agg_nf_z <- aggregate(mcc ~ nFactors + z_threshold, data = raw[raw$corProp == 0.05, ], FUN = mean)

plot(NULL, xlim = range(Z_RANGE), ylim = c(0, 0.55),
     xlab = "z_threshold", ylab = "Mean MCC",
     main = "D1 x R2: z_threshold Sensitivity by nFactors (corProp=0.05)")
abline(h = seq(0, 0.6, 0.05), col = "gray92", lty = 3)
abline(v = 1.5, col = "gray70", lty = 2)

for (i in seq_along(NF_LEVELS)) {
  nf <- NF_LEVELS[i]
  sub <- agg_nf_z[agg_nf_z$nFactors == nf, ]
  lines(sub$z_threshold, sub$mcc, type = "l", col = pal_nf[i], lwd = 2.5)
  points(sub$z_threshold[which.max(sub$mcc)], max(sub$mcc), pch = 16, col = pal_nf[i], cex = 1.5)
}
legend("topright", inset = c(-0.12, 0), legend = paste("nF =", NF_LEVELS),
       col = pal_nf, lwd = 2.5, cex = 0.8, bg = "white", title = "nFactors")
mtext("Dots mark optimal z per nF | corProp=0.05 | R=10 reps", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()

# ============================================================
# PLOT 8: Heatmap nFactors x items_per_factor (interaction D1 x D2)
# ============================================================
save_plot("08_heatmap_nF_x_ipf")
par(mar = c(5, 5, 4, 6))

heat_data <- aggregate(mcc ~ nFactors + items_per_factor, data = best, FUN = mean)
mat <- matrix(NA, nrow = length(NF_LEVELS), ncol = length(IPF_LEVELS))
for (i in seq_along(NF_LEVELS)) {
  for (j in seq_along(IPF_LEVELS)) {
    val <- heat_data$mcc[heat_data$nFactors == NF_LEVELS[i] & heat_data$items_per_factor == IPF_LEVELS[j]]
    if (length(val) == 1) mat[i, j] <- val
  }
}

n_colors <- 100
colors <- colorRampPalette(c("#f7fbff", "#6baed6", "#2171b5", "#08306b"))(n_colors)
val_range <- range(mat, na.rm = TRUE)

plot(NULL, xlim = c(0.5, length(IPF_LEVELS) + 0.5), ylim = c(0.5, length(NF_LEVELS) + 0.5),
     xlab = "Items per Factor", ylab = "Number of Factors (nF)",
     main = "D1 x D2: Oracle MCC — nFactors x Items/Factor",
     xaxt = "n", yaxt = "n")
axis(1, at = seq_along(IPF_LEVELS), labels = IPF_LEVELS)
axis(2, at = seq_along(NF_LEVELS), labels = NF_LEVELS, las = 1)

for (i in seq_along(NF_LEVELS)) {
  for (j in seq_along(IPF_LEVELS)) {
    if (!is.na(mat[i, j])) {
      idx <- round((mat[i, j] - val_range[1]) / (val_range[2] - val_range[1]) * (n_colors - 1)) + 1
      rect(j - 0.45, i - 0.45, j + 0.45, i + 0.45,
           col = colors[max(1, min(n_colors, idx))], border = "white", lwd = 2)
      text(j, i, sprintf("%.3f", mat[i, j]), cex = 0.85, font = 2,
           col = if (mat[i, j] > 0.35) "white" else "black")
    }
  }
}
mtext("Averaged over n, pct, corProp | R=10 reps | Oracle best z", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()

# ============================================================
# PLOT 9: MCC by total_items (the unifying variable)
# ============================================================
best$total_items <- best$nFactors * best$items_per_factor

save_plot("09_mcc_by_total_items")
par(mar = c(5, 5, 4, 2))

agg_ti <- aggregate(mcc ~ total_items, data = best, FUN = function(x) c(mean = mean(x), sd = sd(x)))
agg_ti <- data.frame(total = agg_ti$total_items, mean = agg_ti$mcc[,1], sd = agg_ti$mcc[,2])
agg_ti <- agg_ti[order(agg_ti$total), ]

plot(agg_ti$total, agg_ti$mean, type = "b", pch = 16, col = "#2171b5", lwd = 2, cex = 1.2,
     xlab = "Total Items (nF x ipf)", ylab = "Mean MCC (oracle best z)",
     main = "Combined Effect: MCC by Total Questionnaire Length",
     ylim = c(0, max(agg_ti$mean) * 1.15))
arrows(agg_ti$total, agg_ti$mean - agg_ti$sd, agg_ti$total, agg_ti$mean + agg_ti$sd,
       angle = 90, code = 3, length = 0.03, col = "#2171b580", lwd = 1)

# Add LOESS smooth
lo <- loess(mean ~ total, data = agg_ti, span = 0.6)
pred_x <- seq(min(agg_ti$total), max(agg_ti$total), length.out = 200)
pred_y <- predict(lo, newdata = data.frame(total = pred_x))
lines(pred_x, pred_y, col = "#e34a33", lwd = 3, lty = 1)

abline(h = 0.3, col = "gray60", lty = 2)
text(max(agg_ti$total) * 0.5, 0.32, "MCC = 0.3 threshold", col = "gray50", cex = 0.7)
abline(h = seq(0, 0.6, 0.1), col = "gray92", lty = 3)
mtext("Blue: raw means | Red: LOESS smooth | R=10 reps", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()

writeLines("--- Total Items Effect ---", report)
writeLines(sprintf("  %3d items: mean=%.3f, sd=%.3f", agg_ti$total, agg_ti$mean, agg_ti$sd), report)
writeLines("", report)

# ============================================================
# PLOT 10: nFactors x n_respondents interaction (D1 x D3)
# ============================================================
save_plot("10_nF_x_n_interaction", width = 1400, height = 900)
par(mar = c(5, 5, 4, 8), xpd = TRUE)

agg_nf_n <- aggregate(mcc ~ nFactors + n_respondents, data = best, FUN = mean)

plot(NULL, xlim = range(NF_LEVELS), ylim = c(0, max(agg_nf_n$mcc) * 1.1),
     xlab = "Number of Factors (nF)", ylab = "Mean MCC (oracle best z)",
     main = "D1 x D3: nFactors x Sample Size Interaction")
abline(h = seq(0, 0.6, 0.05), col = "gray92", lty = 3)

for (i in seq_along(N_LEVELS)) {
  n_val <- N_LEVELS[i]
  sub <- agg_nf_n[agg_nf_n$n_respondents == n_val, ]
  lines(sub$nFactors, sub$mcc, type = "b", pch = 16, col = pal_n[as.character(n_val)], lwd = 2.5, cex = 1.3)
}
legend("topright", inset = c(-0.12, 0),
       legend = paste("n =", N_LEVELS), col = unname(pal_n), pch = 16, lwd = 2.5,
       cex = 0.9, bg = "white", title = "Sample Size")
mtext("Averaged over ipf, pct, corProp | R=10 reps", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()

# ============================================================
# PLOT 11: nFactors x pct_careless interaction (D1 x D4)
# ============================================================
save_plot("11_nF_x_pct_interaction", width = 1400, height = 900)
par(mar = c(5, 5, 4, 8), xpd = TRUE)

agg_nf_pct <- aggregate(mcc ~ nFactors + pct_careless, data = best, FUN = mean)

plot(NULL, xlim = range(NF_LEVELS), ylim = c(0, max(agg_nf_pct$mcc) * 1.1),
     xlab = "Number of Factors (nF)", ylab = "Mean MCC (oracle best z)",
     main = "D1 x D4: nFactors x Careless Base Rate Interaction")
abline(h = seq(0, 0.6, 0.05), col = "gray92", lty = 3)

for (i in seq_along(PCT_LEVELS)) {
  pct_val <- PCT_LEVELS[i]
  sub <- agg_nf_pct[agg_nf_pct$pct_careless == pct_val, ]
  lines(sub$nFactors, sub$mcc, type = "b", pch = 16,
        col = pal_pct[as.character(pct_val)], lwd = 2.5, cex = 1.3)
}
legend("topright", inset = c(-0.12, 0),
       legend = paste0(PCT_LEVELS * 100, "%"), col = unname(pal_pct), pch = 16, lwd = 2.5,
       cex = 0.9, bg = "white", title = "% Careless")
mtext("Averaged over ipf, n, corProp | R=10 reps", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()

# ============================================================
# PLOT 12: nFactors x corProp interaction (D1 x R1)
# ============================================================
save_plot("12_nF_x_corProp_interaction", width = 1400, height = 900)
par(mar = c(5, 5, 4, 8), xpd = TRUE)

agg_nf_cp <- aggregate(mcc ~ nFactors + corProp, data = best, FUN = mean)

plot(NULL, xlim = range(NF_LEVELS), ylim = c(0, max(agg_nf_cp$mcc) * 1.1),
     xlab = "Number of Factors (nF)", ylab = "Mean MCC (oracle best z)",
     main = "D1 x R1: nFactors x corProp Interaction")
abline(h = seq(0, 0.6, 0.05), col = "gray92", lty = 3)

for (i in seq_along(CP_LEVELS)) {
  cp_val <- CP_LEVELS[i]
  sub <- agg_nf_cp[agg_nf_cp$corProp == cp_val, ]
  lines(sub$nFactors, sub$mcc, type = "b", pch = 16,
        col = pal_cp[as.character(cp_val)], lwd = 2.5, cex = 1.3)
}
legend("topright", inset = c(-0.12, 0),
       legend = paste("corProp =", CP_LEVELS), col = unname(pal_cp), pch = 16, lwd = 2.5,
       cex = 0.9, bg = "white", title = "corProp")
mtext("Averaged over ipf, n, pct | R=10 reps", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()

# ============================================================
# PLOT 13: Auto-z vs Fixed z vs Oracle (comprehensive comparison)
# ============================================================
save_plot("13_auto_z_vs_fixed_vs_oracle", width = 1600, height = 900)
par(mar = c(5, 5, 4, 2))

# Get auto-z data
auto_raw <- read.csv("archive/multiverse_full/multiverse_full_auto_z.csv")
auto_sub <- auto_raw[auto_raw$corProp == 0.05, ]
auto_agg <- aggregate(mcc ~ nFactors, data = auto_sub, FUN = mean)

# Fixed z levels
z_vals <- c(0.9, 1.5, 1.9)  # closest to 1.0, 1.5, 2.0 in our grid
z_labels <- c("z=1.0", "z=1.5", "z=2.0")
z_colors <- c("#fc8d59", "#91bfdb", "#4575b4")

fixed_aggs <- list()
for (i in seq_along(z_vals)) {
  sub <- raw[raw$corProp == 0.05 & abs(raw$z_threshold - z_vals[i]) < 0.05, ]
  fixed_aggs[[i]] <- aggregate(mcc ~ nFactors, data = sub, FUN = mean)
}

# Oracle
oracle_sub <- best[best$corProp == 0.05, ]
oracle_agg <- aggregate(mcc ~ nFactors, data = oracle_sub, FUN = mean)

plot(NULL, xlim = range(NF_LEVELS), ylim = c(0, max(oracle_agg$mcc) * 1.15),
     xlab = "Number of Factors (nF)", ylab = "Mean MCC",
     main = "R3: Auto-z vs Fixed z vs Oracle — The Threshold Comparison")
abline(h = seq(0, 0.6, 0.05), col = "gray92", lty = 3)

# Oracle (dashed)
lines(oracle_agg$nFactors, oracle_agg$mcc, type = "b", pch = 4, col = "black",
      lwd = 2, cex = 1.2, lty = 2)

# Fixed z
for (i in seq_along(z_vals)) {
  lines(fixed_aggs[[i]]$nFactors, fixed_aggs[[i]]$mcc, type = "b", pch = 15,
        col = z_colors[i], lwd = 2, cex = 1.1)
}

# Auto-z (bold)
lines(auto_agg$nFactors, auto_agg$mcc, type = "b", pch = 17, col = "#e34a33",
      lwd = 3, cex = 1.5)

legend("topleft",
       legend = c("Oracle best z", "Auto-z (LOESS)", z_labels),
       pch = c(4, 17, 15, 15, 15),
       col = c("black", "#e34a33", z_colors),
       lwd = c(2, 3, 2, 2, 2), lty = c(2, 1, 1, 1, 1),
       cex = 0.85, bg = "white")
mtext("corProp=0.05 | Averaged over ipf, n, pct | R=10 reps", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()

# ============================================================
# PLOT 14: Sensitivity vs Specificity trade-off by z_threshold
# ============================================================
save_plot("14_sens_spec_tradeoff", width = 1400, height = 900)
par(mar = c(5, 5, 4, 2))

agg_ss <- aggregate(cbind(sensitivity, specificity) ~ z_threshold, data = raw[raw$corProp == 0.05, ], FUN = mean)

plot(agg_ss$z_threshold, agg_ss$sensitivity, type = "b", pch = 16, col = "#e34a33", lwd = 3, cex = 1.3,
     xlab = "z_threshold", ylab = "Rate",
     main = "R2: Sensitivity-Specificity Trade-off by z_threshold",
     ylim = c(0, 1))
lines(agg_ss$z_threshold, agg_ss$specificity, type = "b", pch = 17, col = "#4575b4", lwd = 3, cex = 1.3)

# Crossover point
cross_idx <- which.min(abs(agg_ss$sensitivity - agg_ss$specificity))
abline(v = agg_ss$z_threshold[cross_idx], col = "gray50", lty = 2)
text(agg_ss$z_threshold[cross_idx] + 0.15, 0.5,
     sprintf("Crossover\nz=%.1f", agg_ss$z_threshold[cross_idx]), col = "gray40", cex = 0.8)

legend("right", legend = c("Sensitivity", "Specificity"),
       pch = c(16, 17), col = c("#e34a33", "#4575b4"), lwd = 3, cex = 0.9, bg = "white")
abline(h = seq(0, 1, 0.1), col = "gray92", lty = 3)
mtext("corProp=0.05 | Averaged over all conditions | R=10 reps", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()

# ============================================================
# PLOT 15: Boxplot distribution of MCC by nFactors (variability)
# ============================================================
save_plot("15_boxplot_mcc_by_nF", width = 1400, height = 900)
par(mar = c(5, 5, 4, 2))

boxplot(mcc ~ nFactors, data = best, col = pal_nf, border = "gray30",
        xlab = "Number of Factors (nF)", ylab = "MCC (oracle best z)",
        main = "D1: Distribution of MCC by nFactors (all conditions)",
        outline = FALSE)
# Add means
means <- tapply(best$mcc, best$nFactors, mean)
points(seq_along(NF_LEVELS), means, pch = 18, col = "white", cex = 1.8)
points(seq_along(NF_LEVELS), means, pch = 18, col = "#e34a33", cex = 1.5)
legend("topleft", legend = "Mean", pch = 18, col = "#e34a33", cex = 0.9, bg = "white")
mtext("Outliers hidden | All corProp, ipf, n, pct | R=10 reps", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()

# ============================================================
# PLOT 16: Heatmap nFactors x z_threshold (best corProp=0.05)
# ============================================================
save_plot("16_heatmap_nF_x_z", width = 1800, height = 900)
par(mar = c(5, 5, 4, 6))

heat_nf_z <- aggregate(mcc ~ nFactors + z_threshold, data = raw[raw$corProp == 0.05, ], FUN = mean)
mat_nfz <- matrix(NA, nrow = length(NF_LEVELS), ncol = length(Z_RANGE))
for (i in seq_along(NF_LEVELS)) {
  for (j in seq_along(Z_RANGE)) {
    val <- heat_nf_z$mcc[heat_nf_z$nFactors == NF_LEVELS[i] & abs(heat_nf_z$z_threshold - Z_RANGE[j]) < 0.05]
    if (length(val) == 1) mat_nfz[i, j] <- val
  }
}

colors_heat <- colorRampPalette(c("#fff7ec", "#fee8c8", "#fdbb84", "#fc8d59",
                                   "#ef6548", "#d7301f", "#990000"))(n_colors)
val_range2 <- range(mat_nfz, na.rm = TRUE)

plot(NULL, xlim = c(0.5, length(Z_RANGE) + 0.5), ylim = c(0.5, length(NF_LEVELS) + 0.5),
     xlab = "z_threshold", ylab = "Number of Factors (nF)",
     main = "D1 x R2: Mean MCC Heatmap — nFactors x z_threshold (corProp=0.05)",
     xaxt = "n", yaxt = "n")
axis(1, at = seq_along(Z_RANGE), labels = sprintf("%.1f", Z_RANGE), cex.axis = 0.7)
axis(2, at = seq_along(NF_LEVELS), labels = NF_LEVELS, las = 1)

for (i in seq_along(NF_LEVELS)) {
  for (j in seq_along(Z_RANGE)) {
    if (!is.na(mat_nfz[i, j])) {
      idx <- round((mat_nfz[i, j] - val_range2[1]) / (val_range2[2] - val_range2[1]) * (n_colors - 1)) + 1
      rect(j - 0.48, i - 0.48, j + 0.48, i + 0.48,
           col = colors_heat[max(1, min(n_colors, idx))], border = NA)
      text(j, i, sprintf("%.2f", mat_nfz[i, j]), cex = 0.5,
           col = if (mat_nfz[i, j] > 0.3) "white" else "black")
    }
  }
}

# Mark optimal z per nF
for (i in seq_along(NF_LEVELS)) {
  best_j <- which.max(mat_nfz[i, ])
  rect(best_j - 0.48, i - 0.48, best_j + 0.48, i + 0.48, border = "cyan", lwd = 2.5)
}
mtext("Cyan borders = optimal z per nF | corProp=0.05 | R=10 reps", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()

# ============================================================
# PLOT 17: items_per_factor effect within each nF level
# ============================================================
save_plot("17_ipf_within_nF", width = 1600, height = 900)
par(mar = c(5, 5, 4, 2))

agg_nf_ipf <- aggregate(mcc ~ nFactors + items_per_factor, data = best, FUN = mean)

# Grouped bar chart
nf_labels <- paste("nF =", NF_LEVELS)
bar_data <- matrix(NA, nrow = length(IPF_LEVELS), ncol = length(NF_LEVELS))
for (i in seq_along(NF_LEVELS)) {
  for (j in seq_along(IPF_LEVELS)) {
    val <- agg_nf_ipf$mcc[agg_nf_ipf$nFactors == NF_LEVELS[i] & agg_nf_ipf$items_per_factor == IPF_LEVELS[j]]
    if (length(val) == 1) bar_data[j, i] <- val
  }
}
rownames(bar_data) <- paste("ipf =", IPF_LEVELS)
colnames(bar_data) <- nf_labels

barplot(bar_data, beside = TRUE, col = unname(pal_ipf), border = NA,
        xlab = "", ylab = "Mean MCC (oracle best z)",
        main = "D1 x D2: Items per Factor Effect Within Each nF Level",
        ylim = c(0, max(bar_data, na.rm = TRUE) * 1.15),
        cex.names = 0.8, las = 1)
legend("topleft", legend = paste("ipf =", IPF_LEVELS),
       fill = unname(pal_ipf), border = NA, cex = 0.9, bg = "white")
mtext("Averaged over n, pct, corProp | R=10 reps", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()

# ============================================================
# PLOT 18: Mean z-scores (good vs careless) by nFactors
# ============================================================
save_plot("18_z_separation_by_nF")
par(mar = c(5, 5, 4, 2))

# From best data (one RR call per cell)
agg_zsep <- aggregate(cbind(mean_z_good, mean_z_careless) ~ nFactors,
                       data = best[best$corProp == 0.05, ], FUN = mean)

plot(agg_zsep$nFactors, agg_zsep$mean_z_good, type = "b", pch = 16, col = "#4575b4",
     lwd = 3, cex = 1.3,
     xlab = "Number of Factors (nF)", ylab = "Mean z-score",
     main = "Signal Separation: z-scores of Good vs Careless Respondents",
     ylim = range(c(agg_zsep$mean_z_good, agg_zsep$mean_z_careless)) * c(0.9, 1.1))
lines(agg_zsep$nFactors, agg_zsep$mean_z_careless, type = "b", pch = 17, col = "#e34a33",
      lwd = 3, cex = 1.3)

# Shade the gap
polygon(c(agg_zsep$nFactors, rev(agg_zsep$nFactors)),
        c(agg_zsep$mean_z_good, rev(agg_zsep$mean_z_careless)),
        col = "#fddbc720", border = NA)

legend("topleft", legend = c("Good respondents", "Careless respondents"),
       pch = c(16, 17), col = c("#4575b4", "#e34a33"), lwd = 3, cex = 0.9, bg = "white")
abline(h = 1.5, col = "gray50", lty = 2)
text(min(NF_LEVELS) + 2, 1.6, "z=1.5 threshold", col = "gray40", cex = 0.7)
mtext("corProp=0.05 | Averaged over ipf, n, pct | R=10 reps", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()

writeLines("--- Z-score Separation ---", report)
writeLines(sprintf("  nF=%2d: good=%.2f, careless=%.2f, gap=%.2f",
                   agg_zsep$nFactors, agg_zsep$mean_z_good, agg_zsep$mean_z_careless,
                   agg_zsep$mean_z_good - agg_zsep$mean_z_careless), report)
writeLines("", report)

# ============================================================
# FINAL REPORT SUMMARY
# ============================================================
writeLines(paste(rep("=", 80), collapse = ""), report)
writeLines("OVERALL SUMMARY", report)
writeLines(paste(rep("=", 80), collapse = ""), report)
writeLines("", report)

writeLines(sprintf("Overall mean MCC (oracle best z): %.3f", mean(best$mcc)), report)
writeLines(sprintf("Overall median MCC: %.3f", median(best$mcc)), report)
writeLines(sprintf("Max MCC observed: %.3f", max(best$mcc)), report)
writeLines(sprintf("Cells with MCC >= 0.3: %.1f%%", mean(best$mcc >= 0.3) * 100), report)
writeLines(sprintf("Cells with MCC >= 0.5: %.1f%%", mean(best$mcc >= 0.5) * 100), report)
writeLines("", report)

writeLines("VARIABLE IMPORTANCE (eta-squared from one-way ANOVAs on oracle MCC):", report)
vars <- c("nFactors", "items_per_factor", "n_respondents", "pct_careless", "corProp")
for (v in vars) {
  fit <- aov(as.formula(paste("mcc ~", v)), data = best)
  ss <- summary(fit)[[1]]
  eta2 <- ss$`Sum Sq`[1] / sum(ss$`Sum Sq`)
  writeLines(sprintf("  %20s: eta2 = %.3f (%.1f%% variance explained)", v, eta2, eta2 * 100), report)
}
writeLines("", report)

writeLines("RECOMMENDED DEFAULTS:", report)
writeLines("  corProp = 0.03 (best overall)", report)
writeLines("  z_threshold = 1.5 (robust universal default)", report)
writeLines("  auto_z = optional (no performance gain, adds convenience)", report)
writeLines("  align_signs = TRUE (mandatory)", report)
writeLines("  min_pairs = 15", report)
writeLines("  iterations = 100", report)
writeLines("", report)
writeLines("PRACTICAL RECOMMENDATION:", report)
writeLines("  Use ReReReRe when total_items >= 60 (e.g., nF>=10 with ipf>=6)", report)
writeLines("  Sweet spot: total_items >= 120, n >= 300", report)
writeLines("  Below 60 items: detection is weak (MCC < 0.20)", report)

close(report)
cat(sprintf("\nReport saved to %s\n", report_file))
cat(sprintf("Total plots generated: 18\n"))
cat(sprintf("All output in: %s/\n", outdir))
cat("\n=== REPORT COMPLETE ===\n")
