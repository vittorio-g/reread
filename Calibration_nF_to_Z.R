#### Calibration: nF → optimal z_threshold ####
#
# Purpose: Build a lookup table mapping number of factors to optimal z_threshold.
# Design: items/factor FIXED at 6 to isolate the pure nF effect.
#
# Parameters:
#   nF: 2 to 40 (step 1) = 39 levels
#   z_threshold: 0.1 to 3.0 (step 0.2) = 15 levels (post-hoc, free)
#   R: 30 replications per nF
#   n_respondents: 300 (fixed, in "recommended" zone)
#   pct_careless: 0.10 (realistic base rate)
#   corProp: 0.05 (recommended default)
#   items_per_factor: 6 (fixed)
#   iterations: 100 (permutation iterations)
#
# Output:
#   calibration_nF_z_raw.csv      — all reps × nF × z_threshold
#   calibration_nF_z_best.csv     — optimal z per rep per nF (for scatter plot)
#   calibration_nF_z_lookup.csv   — mean optimal z per nF (the actual lookup table)
#   plot_calibration_scatter.png   — scatter cloud + mean overlay
# ===

source("Synthetic_Good_Responses_2.R")
source("Careless_machine_2.R")
source("ReReReRe.R")

# --- Fixed parameters ---
ITEMS_PER_FACTOR <- 6
N_RESPONDENTS    <- 300
PCT_CARELESS     <- 0.10
COR_PROP         <- 0.05
ITERATIONS       <- 100
REPS             <- 30
SEED_BASE        <- 2026
CARELESS_LEVELS  <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5)

NF_RANGE         <- 2:40
Z_RANGE          <- seq(0.1, 3.0, by = 0.2)

# --- MCC helper ---
compute_mcc <- function(tp, tn, fp, fn) {
  num <- as.double(tp) * as.double(tn) - as.double(fp) * as.double(fn)
  den <- sqrt(as.double(tp + fp) * as.double(tp + fn) *
              as.double(tn + fp) * as.double(tn + fn))
  if (den == 0) return(0)
  num / den
}

# --- Total work estimate ---
total_calls <- length(NF_RANGE) * REPS
cat(sprintf("=== Calibration nF -> z_threshold ===\n"))
cat(sprintf("  nF range: %d to %d (%d levels)\n", min(NF_RANGE), max(NF_RANGE), length(NF_RANGE)))
cat(sprintf("  z range: %.1f to %.1f (step 0.2, %d levels, post-hoc)\n",
            min(Z_RANGE), max(Z_RANGE), length(Z_RANGE)))
cat(sprintf("  Replications: %d\n", REPS))
cat(sprintf("  Total ReReReRe calls: %d\n", total_calls))
cat(sprintf("  Fixed: items/factor=%d, n=%d, pct_careless=%.0f%%, corProp=%.2f\n",
            ITEMS_PER_FACTOR, N_RESPONDENTS, PCT_CARELESS * 100, COR_PROP))
cat(sprintf("  Estimated time: ~100 minutes\n\n"))

# --- Results storage ---
all_results <- list()
result_idx <- 0

t_start <- Sys.time()

for (nf in NF_RANGE) {
  nItems_vec <- rep(ITEMS_PER_FACTOR, nf)
  total_items <- sum(nItems_vec)

  for (rep in seq_len(REPS)) {
    seed <- SEED_BASE + (nf - min(NF_RANGE)) * REPS + rep

    # 1. Generate clean data
    dat <- simulated_good_responses(
      nConstructs = nf,
      nItems = nItems_vec,
      n = N_RESPONDENTS,
      seed = seed
    )

    # 2. Inject careless
    inj <- inject_careless(
      dat,
      pct_careless = PCT_CARELESS,
      careless_levels = CARELESS_LEVELS,
      seed = seed + 100000
    )

    labels <- inj$labels  # 1 = careless, 0 = clean

    # 3. Run ReReReRe (once per dataset)
    rr <- ReReReRe(
      inj$data_corrupted,
      corProp = COR_PROP,
      iterations = ITERATIONS,
      align_signs = TRUE
    )

    # 4. Evaluate across all z_thresholds (post-hoc, free)
    for (z_thr in Z_RANGE) {
      flagged <- as.integer(rr$z_score <= z_thr)

      tp <- sum(flagged == 1 & labels == 1)
      tn <- sum(flagged == 0 & labels == 0)
      fp <- sum(flagged == 1 & labels == 0)
      fn <- sum(flagged == 0 & labels == 1)

      mcc <- compute_mcc(tp, tn, fp, fn)
      sens <- if (tp + fn > 0) tp / (tp + fn) else 0
      spec <- if (tn + fp > 0) tn / (tn + fp) else 0

      result_idx <- result_idx + 1
      all_results[[result_idx]] <- data.frame(
        nFactors = nf,
        total_items = total_items,
        rep_id = rep,
        z_threshold = z_thr,
        mcc = mcc,
        sensitivity = sens,
        specificity = spec,
        tp = tp, tn = tn, fp = fp, fn = fn,
        mean_z_good = mean(rr$z_score[labels == 0]),
        mean_z_careless = mean(rr$z_score[labels == 1]),
        stringsAsFactors = FALSE
      )
    }
  }

  # Progress
  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  done_nf <- which(NF_RANGE == nf)
  pct <- done_nf / length(NF_RANGE) * 100
  eta <- elapsed / done_nf * (length(NF_RANGE) - done_nf)
  cat(sprintf("[%5.1f%%] nF=%2d done | %.1f min elapsed | ETA %.1f min\n",
              pct, nf, elapsed, eta))
}

elapsed_total <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat(sprintf("\n=== DONE in %.1f minutes ===\n", elapsed_total))

# --- Combine results ---
raw_df <- do.call(rbind, all_results)
write.csv(raw_df, "calibration_nF_z_raw.csv", row.names = FALSE)
cat(sprintf("Saved calibration_nF_z_raw.csv (%d rows)\n", nrow(raw_df)))

# --- Best z per rep per nF ---
best_df <- do.call(rbind, lapply(split(raw_df, paste(raw_df$nFactors, raw_df$rep_id)), function(chunk) {
  best_idx <- which.max(chunk$mcc)
  chunk[best_idx, ]
}))
best_df <- best_df[order(best_df$nFactors, best_df$rep_id), ]
write.csv(best_df, "calibration_nF_z_best.csv", row.names = FALSE)
cat(sprintf("Saved calibration_nF_z_best.csv (%d rows)\n", nrow(best_df)))

# --- Lookup table: mean optimal z per nF ---
lookup_df <- do.call(rbind, lapply(split(best_df, best_df$nFactors), function(chunk) {
  data.frame(
    nFactors = chunk$nFactors[1],
    total_items = chunk$total_items[1],
    mean_best_z = mean(chunk$z_threshold),
    median_best_z = median(chunk$z_threshold),
    sd_best_z = sd(chunk$z_threshold),
    mean_best_mcc = mean(chunk$mcc),
    median_best_mcc = median(chunk$mcc),
    sd_mcc = sd(chunk$mcc),
    mean_z_good = mean(chunk$mean_z_good),
    mean_z_careless = mean(chunk$mean_z_careless),
    stringsAsFactors = FALSE
  )
}))
lookup_df <- lookup_df[order(lookup_df$nFactors), ]
write.csv(lookup_df, "calibration_nF_z_lookup.csv", row.names = FALSE)
cat(sprintf("Saved calibration_nF_z_lookup.csv (%d rows)\n", nrow(lookup_df)))

# --- Print lookup table ---
cat("\n=== LOOKUP TABLE: nF -> optimal z_threshold ===\n")
cat(sprintf("%4s %5s %8s %8s %8s %8s\n",
            "nF", "items", "mean_z", "med_z", "mean_MCC", "sd_MCC"))
for (i in seq_len(nrow(lookup_df))) {
  cat(sprintf("%4d %5d %8.2f %8.2f %8.3f %8.3f\n",
              lookup_df$nFactors[i],
              lookup_df$total_items[i],
              lookup_df$mean_best_z[i],
              lookup_df$median_best_z[i],
              lookup_df$mean_best_mcc[i],
              lookup_df$sd_mcc[i]))
}

# --- Scatter cloud plot ---
cat("\nGenerating plot...\n")
png("plot_calibration_scatter.png", width = 1400, height = 800, res = 150)

par(mar = c(5, 5, 3, 1))

# Jitter for visibility
set.seed(42)
jx <- best_df$nFactors + runif(nrow(best_df), -0.3, 0.3)
jy <- best_df$z_threshold + runif(nrow(best_df), -0.08, 0.08)

# Scale point size by MCC (range: 0.5 to 3)
mcc_range <- range(best_df$mcc, na.rm = TRUE)
if (mcc_range[2] > mcc_range[1]) {
  pt_size <- 0.5 + 2.5 * (best_df$mcc - mcc_range[1]) / (mcc_range[2] - mcc_range[1])
} else {
  pt_size <- rep(1.5, nrow(best_df))
}

# Cloud points (semi-transparent)
plot(jx, jy,
     pch = 16, cex = pt_size,
     col = adjustcolor("steelblue", alpha.f = 0.35),
     xlab = "Number of Factors (nF)",
     ylab = "Optimal z threshold",
     main = "Calibration: Optimal z_threshold by Number of Factors",
     xlim = c(1, 41), ylim = c(0, 3.2),
     xaxt = "n", las = 1)

# X-axis with all nF values
axis(1, at = seq(2, 40, by = 2), labels = seq(2, 40, by = 2), cex.axis = 0.7)

# Mean points (highlighted)
mcc_mean_range <- range(lookup_df$mean_best_mcc, na.rm = TRUE)
if (mcc_mean_range[2] > mcc_mean_range[1]) {
  mean_pt_size <- 1 + 3 * (lookup_df$mean_best_mcc - mcc_mean_range[1]) /
    (mcc_mean_range[2] - mcc_mean_range[1])
} else {
  mean_pt_size <- rep(2, nrow(lookup_df))
}

points(lookup_df$nFactors, lookup_df$mean_best_z,
       pch = 21, cex = mean_pt_size,
       bg = "red", col = "darkred", lwd = 1.5)

# Connect means with a line
lines(lookup_df$nFactors, lookup_df$mean_best_z,
      col = "darkred", lwd = 2, lty = 2)

# Legend
legend("topright",
       legend = c("Individual rep (size = MCC)", "Mean across reps (size = mean MCC)"),
       pch = c(16, 21),
       col = c(adjustcolor("steelblue", 0.6), "darkred"),
       pt.bg = c(NA, "red"),
       pt.cex = c(1.5, 2),
       cex = 0.8,
       bg = "white")

# Subtitle
mtext(sprintf("R=%d reps | n=%d | items/factor=%d | corProp=%.2f | pct_careless=%.0f%%",
              REPS, N_RESPONDENTS, ITEMS_PER_FACTOR, COR_PROP, PCT_CARELESS * 100),
      side = 1, line = 3.8, cex = 0.75, col = "gray40")

dev.off()
cat("Saved plot_calibration_scatter.png\n")
cat("\n=== ALL DONE ===\n")
