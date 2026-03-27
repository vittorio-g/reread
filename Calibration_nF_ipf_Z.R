#### Calibration: nF x items_per_factor → optimal z_threshold ####
#
# Creates a 2D lookup surface for auto-calibrated z_threshold
# based on BOTH nFactors and items_per_factor.
#
# Design:
#   nFactors:        4, 6, 8, 10, 12, 15, 18, 20, 25, 30  (10)
#   items_per_factor: 3, 4, 6, 8, 10, 12                    (6)
#   z_threshold:     0.1-3.0 step 0.1                       (30, post-hoc)
#   Reps:            30 per cell
#
# Fixed: n=300, pct=10%, corProp=0.03, iterations=100, min_pairs=15
#
# Total: 60 cells x 30 reps = 1,800 RR calls
# ===

setwd("C:/Users/vitto/Desktop/ReReReRe")
source("Synthetic_Good_Responses_2.R")
source("Careless_machine_2.R")
source("ReReReRe.R")

# --- Design ---
NF_LEVELS  <- c(4, 6, 8, 10, 12, 15, 18, 20, 25, 30)
IPF_LEVELS <- c(3, 4, 6, 8, 10, 12)
Z_RANGE    <- seq(0.1, 3.0, by = 0.1)
REPS       <- 30

# Fixed
N_RESP     <- 300
PCT        <- 0.10
CORPROP    <- 0.03
ITERATIONS <- 100
MIN_PAIRS  <- 15
SEED_BASE  <- 20260327
CARELESS_LEVELS <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5)

# --- MCC ---
compute_mcc <- function(tp, tn, fp, fn) {
  num <- as.double(tp)*as.double(tn) - as.double(fp)*as.double(fn)
  den <- sqrt(as.double(tp+fp)*as.double(tp+fn)*as.double(tn+fp)*as.double(tn+fn))
  if (den == 0) return(0)
  num / den
}

# --- Grid ---
grid <- expand.grid(nF = NF_LEVELS, ipf = IPF_LEVELS, stringsAsFactors = FALSE)
grid$total_items <- grid$nF * grid$ipf
n_cells <- nrow(grid)
n_total <- n_cells * REPS

cat(sprintf("=== Calibration nF x ipf → z ===\n"))
cat(sprintf("  Cells: %d (nF x ipf)\n", n_cells))
cat(sprintf("  Reps per cell: %d\n", REPS))
cat(sprintf("  Total RR calls: %d\n", n_total))
cat(sprintf("  z levels: %d (post-hoc, step 0.1)\n", length(Z_RANGE)))
cat(sprintf("  Total result rows: %d\n", n_total * length(Z_RANGE)))
cat("\n")

# --- Run ---
all_results <- list()
result_idx <- 0
call_count <- 0
t_start <- Sys.time()

for (ci in seq_len(n_cells)) {
  nf  <- grid$nF[ci]
  ipf <- grid$ipf[ci]
  ti  <- grid$total_items[ci]

  for (rep in seq_len(REPS)) {
    seed <- SEED_BASE + (ci - 1) * REPS + rep

    # Generate data
    dat <- tryCatch({
      simulated_good_responses(nf, rep(ipf, nf), n = N_RESP, seed = seed)
    }, error = function(e) NULL)
    if (is.null(dat)) next

    # Inject careless
    inj <- inject_careless(dat, PCT, careless_levels = CARELESS_LEVELS, seed = seed + 500000)
    labels <- inj$labels

    # Run ReReReRe
    rr <- tryCatch({
      ReReReRe(inj$data_corrupted, corProp = CORPROP, iterations = ITERATIONS,
               min_pairs = MIN_PAIRS, align_signs = TRUE, progress = FALSE)
    }, error = function(e) NULL)
    if (is.null(rr)) next

    call_count <- call_count + 1

    # Evaluate all z thresholds
    for (z_thr in Z_RANGE) {
      flagged <- as.integer(rr$z_score <= z_thr)
      tp <- sum(flagged == 1 & labels == 1)
      tn <- sum(flagged == 0 & labels == 0)
      fp <- sum(flagged == 1 & labels == 0)
      fn <- sum(flagged == 0 & labels == 1)
      mcc <- compute_mcc(tp, tn, fp, fn)
      sens <- if (tp+fn > 0) tp/(tp+fn) else 0
      spec <- if (tn+fp > 0) tn/(tn+fp) else 0

      result_idx <- result_idx + 1
      all_results[[result_idx]] <- data.frame(
        nFactors = nf, items_per_factor = ipf, total_items = ti,
        z_threshold = z_thr, rep_id = rep,
        mcc = mcc, sensitivity = sens, specificity = spec,
        tp = tp, tn = tn, fp = fp, fn = fn,
        mean_z_good = mean(rr$z_score[labels == 0], na.rm = TRUE),
        mean_z_careless = mean(rr$z_score[labels == 1], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }

  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  pct_done <- ci / n_cells * 100
  eta <- if (ci > 0) elapsed / ci * (n_cells - ci) else NA
  cat(sprintf("[%5.1f%%] nF=%2d ipf=%2d (%3d items) done | %d calls | %.1f min | ETA %.1f min\n",
              pct_done, nf, ipf, ti, call_count, elapsed, eta))
}

elapsed_total <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat(sprintf("\n=== DONE in %.1f min | %d calls | %d rows ===\n",
            elapsed_total, call_count, result_idx))

# --- Save raw ---
raw_df <- do.call(rbind, all_results)
outdir <- "archive/calibration_nF_ipf_z"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
write.csv(raw_df, file.path(outdir, "calibration_raw.csv"), row.names = FALSE)
cat(sprintf("Saved calibration_raw.csv (%d rows)\n", nrow(raw_df)))

# --- Best z per rep per cell ---
best_df <- do.call(rbind, lapply(
  split(raw_df, paste(raw_df$nFactors, raw_df$items_per_factor, raw_df$rep_id)),
  function(chunk) chunk[which.max(chunk$mcc), ]
))
best_df <- best_df[order(best_df$nFactors, best_df$items_per_factor, best_df$rep_id), ]
write.csv(best_df, file.path(outdir, "calibration_best.csv"), row.names = FALSE)
cat(sprintf("Saved calibration_best.csv (%d rows)\n", nrow(best_df)))

# --- Lookup table: mean optimal z per nF x ipf ---
lookup <- aggregate(cbind(z_threshold, mcc) ~ nFactors + items_per_factor + total_items,
                     data = best_df, FUN = function(x) c(mean = mean(x), median = median(x), sd = sd(x)))
lookup_flat <- data.frame(
  nFactors = lookup$nFactors,
  items_per_factor = lookup$items_per_factor,
  total_items = lookup$total_items,
  mean_best_z = lookup$z_threshold[, "mean"],
  median_best_z = lookup$z_threshold[, "median"],
  sd_best_z = lookup$z_threshold[, "sd"],
  mean_mcc = lookup$mcc[, "mean"],
  median_mcc = lookup$mcc[, "median"],
  sd_mcc = lookup$mcc[, "sd"]
)
lookup_flat <- lookup_flat[order(lookup_flat$nFactors, lookup_flat$items_per_factor), ]
write.csv(lookup_flat, file.path(outdir, "calibration_lookup.csv"), row.names = FALSE)
cat(sprintf("Saved calibration_lookup.csv (%d rows)\n", nrow(lookup_flat)))

# --- Print lookup table ---
cat("\n=== LOOKUP TABLE: mean optimal z by nF x ipf ===\n")
cat(sprintf("%4s", "nF"))
for (ipf in IPF_LEVELS) cat(sprintf(" | ipf=%2d (z/MCC)", ipf))
cat("\n")
for (nf in NF_LEVELS) {
  cat(sprintf("%4d", nf))
  for (ipf in IPF_LEVELS) {
    row <- lookup_flat[lookup_flat$nFactors == nf & lookup_flat$items_per_factor == ipf, ]
    if (nrow(row) == 1) {
      cat(sprintf(" | %4.1f / %.3f  ", row$mean_best_z, row$mean_mcc))
    } else {
      cat(" |      ---      ")
    }
  }
  cat("\n")
}

# ============================================================
# PLOTS
# ============================================================
cat("\nGenerating plots...\n")

# --- PLOT 1: Heatmap of optimal z (nF x ipf) ---
png(file.path(outdir, "plot_01_heatmap_optimal_z.png"), width = 1400, height = 1000, res = 150)
par(mar = c(5, 5, 4, 7))

mat_z <- matrix(NA, nrow = length(NF_LEVELS), ncol = length(IPF_LEVELS))
for (i in seq_along(NF_LEVELS)) {
  for (j in seq_along(IPF_LEVELS)) {
    val <- lookup_flat$mean_best_z[lookup_flat$nFactors == NF_LEVELS[i] &
                                     lookup_flat$items_per_factor == IPF_LEVELS[j]]
    if (length(val) == 1) mat_z[i, j] <- val
  }
}

n_col <- 100
cols_z <- colorRampPalette(c("#f7fcf5", "#74c476", "#238b45", "#00441b"))(n_col)
vr <- range(mat_z, na.rm = TRUE)

plot(NULL, xlim = c(0.5, length(IPF_LEVELS) + 0.5), ylim = c(0.5, length(NF_LEVELS) + 0.5),
     xlab = "Items per Factor", ylab = "Number of Factors (nF)",
     main = "Optimal z_threshold by nFactors x Items/Factor\n(mean across 30 reps)",
     xaxt = "n", yaxt = "n")
axis(1, at = seq_along(IPF_LEVELS), labels = IPF_LEVELS)
axis(2, at = seq_along(NF_LEVELS), labels = NF_LEVELS, las = 1)

for (i in seq_along(NF_LEVELS)) {
  for (j in seq_along(IPF_LEVELS)) {
    if (!is.na(mat_z[i, j])) {
      idx <- round((mat_z[i, j] - vr[1]) / (vr[2] - vr[1]) * (n_col - 1)) + 1
      rect(j-0.45, i-0.45, j+0.45, i+0.45,
           col = cols_z[max(1, min(n_col, idx))], border = "white", lwd = 2)
      text(j, i, sprintf("%.1f", mat_z[i, j]), cex = 0.85, font = 2,
           col = if (mat_z[i, j] > 1.5) "white" else "black")
    }
  }
}
mtext("n=300, pct=10%, corProp=0.03 | R=30 reps", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()
cat("Saved plot_01_heatmap_optimal_z.png\n")

# --- PLOT 2: Heatmap of MCC at optimal z ---
png(file.path(outdir, "plot_02_heatmap_mcc.png"), width = 1400, height = 1000, res = 150)
par(mar = c(5, 5, 4, 7))

mat_mcc <- matrix(NA, nrow = length(NF_LEVELS), ncol = length(IPF_LEVELS))
for (i in seq_along(NF_LEVELS)) {
  for (j in seq_along(IPF_LEVELS)) {
    val <- lookup_flat$mean_mcc[lookup_flat$nFactors == NF_LEVELS[i] &
                                  lookup_flat$items_per_factor == IPF_LEVELS[j]]
    if (length(val) == 1) mat_mcc[i, j] <- val
  }
}

cols_mcc <- colorRampPalette(c("#f7fbff", "#6baed6", "#2171b5", "#08306b"))(n_col)
vr_m <- range(mat_mcc, na.rm = TRUE)

plot(NULL, xlim = c(0.5, length(IPF_LEVELS) + 0.5), ylim = c(0.5, length(NF_LEVELS) + 0.5),
     xlab = "Items per Factor", ylab = "Number of Factors (nF)",
     main = "Mean MCC (at optimal z) by nFactors x Items/Factor\n(mean across 30 reps)",
     xaxt = "n", yaxt = "n")
axis(1, at = seq_along(IPF_LEVELS), labels = IPF_LEVELS)
axis(2, at = seq_along(NF_LEVELS), labels = NF_LEVELS, las = 1)

for (i in seq_along(NF_LEVELS)) {
  for (j in seq_along(IPF_LEVELS)) {
    if (!is.na(mat_mcc[i, j])) {
      idx <- round((mat_mcc[i, j] - vr_m[1]) / (vr_m[2] - vr_m[1]) * (n_col - 1)) + 1
      rect(j-0.45, i-0.45, j+0.45, i+0.45,
           col = cols_mcc[max(1, min(n_col, idx))], border = "white", lwd = 2)
      text(j, i, sprintf("%.2f", mat_mcc[i, j]), cex = 0.8, font = 2,
           col = if (mat_mcc[i, j] > 0.35) "white" else "black")
    }
  }
}
mtext("n=300, pct=10%, corProp=0.03 | R=30 reps | Oracle best z", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()
cat("Saved plot_02_heatmap_mcc.png\n")

# --- PLOT 3: Scatter cloud — optimal z by total_items (all reps) ---
png(file.path(outdir, "plot_03_scatter_z_by_total_items.png"), width = 1600, height = 1000, res = 150)
par(mar = c(5, 5, 4, 2))

set.seed(42)
jitter_x <- runif(nrow(best_df), -0.02, 0.02) * best_df$total_items
jitter_y <- runif(nrow(best_df), -0.03, 0.03)

# Point size proportional to MCC
cex_vals <- 0.3 + best_df$mcc * 3

plot(best_df$total_items + jitter_x, best_df$z_threshold + jitter_y,
     pch = 16, col = rgb(0.2, 0.4, 0.8, 0.3), cex = cex_vals,
     xlab = "Total Items (nF x ipf)", ylab = "Optimal z_threshold",
     main = "Optimal z_threshold by Total Questionnaire Length\n(each dot = 1 replication, size = MCC)")

# Mean per total_items
mean_by_ti <- aggregate(cbind(z_threshold, mcc) ~ total_items, data = best_df, FUN = mean)
mean_by_ti <- mean_by_ti[order(mean_by_ti$total_items), ]
points(mean_by_ti$total_items, mean_by_ti$z_threshold,
       pch = 16, col = "red", cex = 0.5 + mean_by_ti$mcc * 5)

# LOESS smooth
lo <- loess(z_threshold ~ total_items, data = mean_by_ti, span = 0.5)
pred_x <- seq(min(mean_by_ti$total_items), max(mean_by_ti$total_items), length.out = 200)
pred_y <- predict(lo, newdata = data.frame(total_items = pred_x))
lines(pred_x, pred_y, col = "red", lwd = 3, lty = 2)

abline(h = 1.5, col = "gray50", lty = 3)
text(max(best_df$total_items) * 0.8, 1.6, "z=1.5 (fixed default)", col = "gray40", cex = 0.7)

legend("topright",
       legend = c("Individual reps (size=MCC)", "Mean per total_items (size=MCC)", "LOESS smooth"),
       pch = c(16, 16, NA), lty = c(NA, NA, 2), lwd = c(NA, NA, 3),
       col = c(rgb(0.2, 0.4, 0.8, 0.5), "red", "red"),
       pt.cex = c(1.5, 2, NA), cex = 0.8, bg = "white")
mtext("n=300, pct=10%, corProp=0.03 | R=30 reps per cell", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()
cat("Saved plot_03_scatter_z_by_total_items.png\n")

# --- PLOT 4: Scatter cloud — optimal z by nF, colored by ipf ---
png(file.path(outdir, "plot_04_scatter_z_by_nF_colored.png"), width = 1600, height = 1000, res = 150)
par(mar = c(5, 5, 4, 2))

ipf_colors <- c("3" = "#e34a33", "4" = "#fc8d59", "6" = "#4575b4",
                 "8" = "#91bfdb", "10" = "#31a354", "12" = "#006d2c")

set.seed(42)
jx <- runif(nrow(best_df), -0.3, 0.3)
jy <- runif(nrow(best_df), -0.03, 0.03)

plot(best_df$nFactors + jx, best_df$z_threshold + jy,
     pch = 16, cex = 0.3 + best_df$mcc * 3,
     col = adjustcolor(ipf_colors[as.character(best_df$items_per_factor)], alpha = 0.4),
     xlab = "Number of Factors (nF)", ylab = "Optimal z_threshold",
     main = "Optimal z by nFactors, colored by items/factor\n(size = MCC)")

# Mean overlay per nF x ipf
for (ipf_val in IPF_LEVELS) {
  sub <- best_df[best_df$items_per_factor == ipf_val, ]
  agg <- aggregate(cbind(z_threshold, mcc) ~ nFactors, data = sub, FUN = mean)
  lines(agg$nFactors, agg$z_threshold, col = ipf_colors[as.character(ipf_val)], lwd = 2.5)
  points(agg$nFactors, agg$z_threshold, pch = 16, cex = 0.5 + agg$mcc * 4,
         col = ipf_colors[as.character(ipf_val)])
}

legend("topright",
       legend = paste("ipf =", IPF_LEVELS),
       col = unname(ipf_colors), pch = 16, lwd = 2.5,
       cex = 0.85, bg = "white", title = "Items/Factor")
mtext("Dots: individual reps | Lines: means | Size = MCC", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()
cat("Saved plot_04_scatter_z_by_nF_colored.png\n")

# --- PLOT 5: z threshold curves (MCC vs z) for selected nF x ipf combos ---
png(file.path(outdir, "plot_05_z_curves_selected.png"), width = 1600, height = 1000, res = 150)
par(mar = c(5, 5, 4, 8), xpd = TRUE)

# Select representative combos
combos <- data.frame(
  nf = c(4, 8, 12, 20, 30, 4, 8, 12, 20, 30),
  ipf = c(3, 3, 3, 3, 3, 10, 10, 10, 10, 10)
)
combos$label <- paste0("nF=", combos$nf, " ipf=", combos$ipf)
combos$total <- combos$nf * combos$ipf

combo_cols <- colorRampPalette(c("#fc8d59", "#b30000"))(5)
combo_cols <- c(combo_cols, colorRampPalette(c("#91bfdb", "#023858"))(5))

plot(NULL, xlim = range(Z_RANGE), ylim = c(0, 0.65),
     xlab = "z_threshold", ylab = "Mean MCC",
     main = "MCC vs z_threshold for Selected nF x ipf Combinations")
abline(h = seq(0, 0.6, 0.05), col = "gray92", lty = 3)

for (k in 1:nrow(combos)) {
  sub <- raw_df[raw_df$nFactors == combos$nf[k] & raw_df$items_per_factor == combos$ipf[k], ]
  agg <- aggregate(mcc ~ z_threshold, data = sub, FUN = mean)
  lty_val <- if (combos$ipf[k] == 3) 2 else 1
  lines(agg$z_threshold, agg$mcc, col = combo_cols[k], lwd = 2, lty = lty_val)
  # Mark peak
  best_idx <- which.max(agg$mcc)
  points(agg$z_threshold[best_idx], agg$mcc[best_idx], pch = 16, col = combo_cols[k], cex = 1.5)
}

legend("topright", inset = c(-0.15, 0),
       legend = paste0(combos$label, " (", combos$total, ")"),
       col = combo_cols, lwd = 2,
       lty = c(rep(2, 5), rep(1, 5)),
       cex = 0.65, bg = "white", title = "Combo (total items)")
dev.off()
cat("Saved plot_05_z_curves_selected.png\n")

# --- PLOT 6: SD of optimal z (uncertainty heatmap) ---
png(file.path(outdir, "plot_06_heatmap_z_variability.png"), width = 1400, height = 1000, res = 150)
par(mar = c(5, 5, 4, 7))

mat_sd <- matrix(NA, nrow = length(NF_LEVELS), ncol = length(IPF_LEVELS))
for (i in seq_along(NF_LEVELS)) {
  for (j in seq_along(IPF_LEVELS)) {
    val <- lookup_flat$sd_best_z[lookup_flat$nFactors == NF_LEVELS[i] &
                                    lookup_flat$items_per_factor == IPF_LEVELS[j]]
    if (length(val) == 1) mat_sd[i, j] <- val
  }
}

cols_sd <- colorRampPalette(c("#f7f7f7", "#fdae61", "#d73027", "#67001f"))(n_col)
vr_sd <- range(mat_sd, na.rm = TRUE)

plot(NULL, xlim = c(0.5, length(IPF_LEVELS) + 0.5), ylim = c(0.5, length(NF_LEVELS) + 0.5),
     xlab = "Items per Factor", ylab = "Number of Factors (nF)",
     main = "Variability (SD) of Optimal z_threshold\n(higher = less stable calibration)",
     xaxt = "n", yaxt = "n")
axis(1, at = seq_along(IPF_LEVELS), labels = IPF_LEVELS)
axis(2, at = seq_along(NF_LEVELS), labels = NF_LEVELS, las = 1)

for (i in seq_along(NF_LEVELS)) {
  for (j in seq_along(IPF_LEVELS)) {
    if (!is.na(mat_sd[i, j])) {
      idx <- round((mat_sd[i, j] - vr_sd[1]) / (vr_sd[2] - vr_sd[1]) * (n_col - 1)) + 1
      rect(j-0.45, i-0.45, j+0.45, i+0.45,
           col = cols_sd[max(1, min(n_col, idx))], border = "white", lwd = 2)
      text(j, i, sprintf("%.2f", mat_sd[i, j]), cex = 0.8, font = 2,
           col = if (mat_sd[i, j] > 0.6) "white" else "black")
    }
  }
}
mtext("High SD = optimal z varies a lot across reps (noisy calibration)", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()
cat("Saved plot_06_heatmap_z_variability.png\n")

# --- PLOT 7: Optimal z as function of total_items (the key relationship) ---
png(file.path(outdir, "plot_07_z_vs_total_items_loess.png"), width = 1400, height = 900, res = 150)
par(mar = c(5, 5, 4, 2))

# Use lookup means
plot(lookup_flat$total_items, lookup_flat$mean_best_z,
     pch = 16, cex = 0.5 + lookup_flat$mean_mcc * 4,
     col = ipf_colors[as.character(lookup_flat$items_per_factor)],
     xlab = "Total Items (nF x ipf)", ylab = "Mean Optimal z_threshold",
     main = "Optimal z as Function of Total Items\n(size = MCC, color = ipf)")

# LOESS on means
lo2 <- loess(mean_best_z ~ total_items, data = lookup_flat, span = 0.5)
pred_x2 <- seq(min(lookup_flat$total_items), max(lookup_flat$total_items), length.out = 300)
pred_y2 <- predict(lo2, newdata = data.frame(total_items = pred_x2))
lines(pred_x2, pred_y2, col = "black", lwd = 3)

abline(h = 1.5, col = "gray50", lty = 2)
legend("topright",
       legend = c(paste("ipf =", IPF_LEVELS), "LOESS", "z=1.5 default"),
       pch = c(rep(16, 6), NA, NA),
       lty = c(rep(NA, 6), 1, 2),
       lwd = c(rep(NA, 6), 3, 1),
       col = c(unname(ipf_colors), "black", "gray50"),
       pt.cex = c(rep(1.5, 6), NA, NA),
       cex = 0.75, bg = "white")
mtext("Each point = mean of 30 reps | LOESS = smooth calibration curve", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()
cat("Saved plot_07_z_vs_total_items_loess.png\n")

# --- PLOT 8: MCC at fixed z=1.5 vs auto-z (from lookup) ---
png(file.path(outdir, "plot_08_fixed_vs_auto_z_mcc.png"), width = 1400, height = 900, res = 150)
par(mar = c(5, 5, 4, 2))

# MCC at fixed z=1.5
fixed_z <- raw_df[abs(raw_df$z_threshold - 1.5) < 0.05, ]
fixed_agg <- aggregate(mcc ~ nFactors + items_per_factor + total_items, data = fixed_z, FUN = mean)

# MCC at auto-z (use nearest z in grid to the lookup mean)
auto_mcc_list <- list()
for (k in 1:nrow(lookup_flat)) {
  nf <- lookup_flat$nFactors[k]
  ipf <- lookup_flat$items_per_factor[k]
  auto_z_val <- lookup_flat$mean_best_z[k]
  nearest_z <- Z_RANGE[which.min(abs(Z_RANGE - auto_z_val))]
  sub <- raw_df[raw_df$nFactors == nf & raw_df$items_per_factor == ipf &
                  abs(raw_df$z_threshold - nearest_z) < 0.05, ]
  auto_mcc_list[[k]] <- data.frame(
    nFactors = nf, items_per_factor = ipf,
    total_items = nf * ipf,
    mcc_auto = mean(sub$mcc),
    stringsAsFactors = FALSE
  )
}
auto_agg <- do.call(rbind, auto_mcc_list)

comp <- merge(fixed_agg, auto_agg, by = c("nFactors", "items_per_factor", "total_items"))
comp <- comp[order(comp$total_items), ]

plot(comp$total_items, comp$mcc, type = "b", pch = 15, col = "#4575b4", lwd = 2, cex = 1.2,
     xlab = "Total Items", ylab = "Mean MCC",
     main = "Fixed z=1.5 vs Auto-calibrated z by Total Items",
     ylim = c(0, max(comp$mcc_auto, comp$mcc) * 1.1))
lines(comp$total_items, comp$mcc_auto, type = "b", pch = 17, col = "#e34a33", lwd = 2, cex = 1.2)

# Oracle
oracle_agg <- aggregate(mcc ~ total_items, data = best_df, FUN = mean)
oracle_agg <- oracle_agg[order(oracle_agg$total_items), ]
lines(oracle_agg$total_items, oracle_agg$mcc, type = "b", pch = 4, col = "gray40", lwd = 1.5, lty = 2)

legend("topleft",
       legend = c("Fixed z=1.5", "Auto-z (2D lookup)", "Oracle best z"),
       pch = c(15, 17, 4), col = c("#4575b4", "#e34a33", "gray40"),
       lwd = 2, lty = c(1, 1, 2), cex = 0.85, bg = "white")
mtext("Auto-z uses mean optimal z from 2D nF x ipf lookup table", side = 1, line = 3.5, cex = 0.7, col = "gray40")
dev.off()
cat("Saved plot_08_fixed_vs_auto_z_mcc.png\n")

# --- Text report ---
report_file <- file.path(outdir, "calibration_report.txt")
rpt <- file(report_file, open = "w")
writeLines(paste(rep("=", 80), collapse = ""), rpt)
writeLines("CALIBRATION REPORT: nF x ipf -> optimal z_threshold", rpt)
writeLines(paste("Generated:", Sys.time()), rpt)
writeLines(paste(rep("=", 80), collapse = ""), rpt)
writeLines("", rpt)
writeLines(sprintf("Total RR calls: %d", call_count), rpt)
writeLines(sprintf("Runtime: %.1f minutes", elapsed_total), rpt)
writeLines(sprintf("Reps per cell: %d", REPS), rpt)
writeLines(sprintf("Fixed: n=%d, pct=%.0f%%, corProp=%.2f", N_RESP, PCT*100, CORPROP), rpt)
writeLines("", rpt)

writeLines("--- LOOKUP TABLE ---", rpt)
writeLines(sprintf("%4s %4s %5s | %6s %6s %5s | %6s %6s",
                   "nF", "ipf", "total", "z_mean", "z_med", "z_sd", "MCC_m", "MCC_sd"), rpt)
for (k in 1:nrow(lookup_flat)) {
  writeLines(sprintf("%4d %4d %5d | %6.2f %6.1f %5.2f | %6.3f %6.3f",
                     lookup_flat$nFactors[k], lookup_flat$items_per_factor[k], lookup_flat$total_items[k],
                     lookup_flat$mean_best_z[k], lookup_flat$median_best_z[k], lookup_flat$sd_best_z[k],
                     lookup_flat$mean_mcc[k], lookup_flat$sd_mcc[k]), rpt)
}

writeLines("", rpt)
writeLines("--- KEY FINDINGS ---", rpt)
writeLines("", rpt)

# Does total_items alone predict optimal z well?
cor_ti_z <- cor(lookup_flat$total_items, lookup_flat$mean_best_z)
writeLines(sprintf("Correlation(total_items, optimal_z) = %.3f", cor_ti_z), rpt)

# Do nF and ipf have independent effects on z?
if (nrow(lookup_flat) > 10) {
  fit <- lm(mean_best_z ~ nFactors + items_per_factor, data = lookup_flat)
  writeLines(sprintf("Linear model: z = %.3f + %.4f*nF + %.4f*ipf (R2=%.3f)",
                     coef(fit)[1], coef(fit)[2], coef(fit)[3], summary(fit)$r.squared), rpt)
  fit2 <- lm(mean_best_z ~ total_items, data = lookup_flat)
  writeLines(sprintf("Simple model: z = %.3f + %.4f*total_items (R2=%.3f)",
                     coef(fit2)[1], coef(fit2)[2], summary(fit2)$r.squared), rpt)
}

# Fixed z=1.5 vs auto-z comparison
writeLines("", rpt)
writeLines("--- FIXED z=1.5 vs AUTO-z COMPARISON ---", rpt)
writeLines(sprintf("%5s | %8s | %8s | %8s | %5s",
                   "total", "z=1.5", "auto-z", "oracle", "gain"), rpt)
for (k in 1:nrow(comp)) {
  gain <- comp$mcc_auto[k] - comp$mcc[k]
  oracle_val <- oracle_agg$mcc[oracle_agg$total_items == comp$total_items[k]]
  if (length(oracle_val) == 0) oracle_val <- NA
  writeLines(sprintf("%5d | %8.3f | %8.3f | %8.3f | %+.3f",
                     comp$total_items[k], comp$mcc[k], comp$mcc_auto[k],
                     oracle_val, gain), rpt)
}

overall_fixed <- mean(comp$mcc)
overall_auto <- mean(comp$mcc_auto)
overall_oracle <- mean(oracle_agg$mcc)
writeLines(sprintf("\nOverall: fixed=%.3f, auto=%.3f, oracle=%.3f",
                   overall_fixed, overall_auto, overall_oracle), rpt)
writeLines(sprintf("Auto-z gain over fixed: %+.3f MCC", overall_auto - overall_fixed), rpt)
writeLines(sprintf("Oracle ceiling gap: %.3f MCC", overall_oracle - overall_auto), rpt)

close(rpt)
cat(sprintf("\nReport saved to %s\n", report_file))
cat(sprintf("Total plots: 8\n"))
cat("\n=== ALL DONE ===\n")
