#### Full Multiverse: optimal parameters for auto-calibrated ReReReRe ####
#
# Varies both dataset properties and ReReReRe parameters.
# z_threshold and auto_z are evaluated post-hoc (free).
#
# Design:
#   D1: nFactors        = 4, 8, 12, 16, 20, 25, 30  (7)
#   D2: items_per_factor = 3, 6, 10                   (3)
#   D3: n_respondents   = 100, 300, 500               (3)
#   D4: pct_careless    = 0.05, 0.10, 0.25            (3)
#   R1: corProp         = 0.03, 0.05, 0.10            (3)
#   R2: z_threshold     = 0.1-3.0 step 0.2            (15, post-hoc)
#   R3: auto_z          = evaluated post-hoc from LOESS
#   Reps: 10
#
# Fixed: align_signs=TRUE, iterations=100, min_pairs=15,
#        careless types 1/3 each, corruption 50-100%
#
# Total RR calls: 7*3*3*3*3*10 = 5,670
# Estimated runtime: ~85 minutes
# ===

setwd("C:/Users/vitto/Desktop/ReReReRe")
source("Synthetic_Good_Responses_2.R")
source("Careless_machine_2.R")
source("ReReReRe.R")

# --- Design parameters ---
NF_LEVELS     <- c(4, 8, 12, 16, 20, 25, 30)
IPF_LEVELS    <- c(3, 6, 10)
N_LEVELS      <- c(100, 300, 500)
PCT_LEVELS    <- c(0.05, 0.10, 0.25)
CORPROP_LEVELS <- c(0.03, 0.05, 0.10)

Z_RANGE       <- seq(0.1, 3.0, by = 0.2)
REPS          <- 10
ITERATIONS    <- 100
MIN_PAIRS     <- 15
SEED_BASE     <- 20260326
CARELESS_LEVELS <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5)

# --- MCC helper ---
compute_mcc <- function(tp, tn, fp, fn) {
  num <- as.double(tp) * as.double(tn) - as.double(fp) * as.double(fn)
  den <- sqrt(as.double(tp + fp) * as.double(tp + fn) *
              as.double(tn + fp) * as.double(tn + fn))
  if (den == 0) return(0)
  num / den
}

# --- Build condition grid (dataset level) ---
data_grid <- expand.grid(
  nFactors = NF_LEVELS,
  items_per_factor = IPF_LEVELS,
  n_respondents = N_LEVELS,
  pct_careless = PCT_LEVELS,
  stringsAsFactors = FALSE
)
n_data_conds <- nrow(data_grid)
n_rr_calls <- n_data_conds * length(CORPROP_LEVELS) * REPS
cat(sprintf("=== Full Multiverse ===\n"))
cat(sprintf("  Dataset conditions: %d\n", n_data_conds))
cat(sprintf("  corProp levels: %d\n", length(CORPROP_LEVELS)))
cat(sprintf("  Reps: %d\n", REPS))
cat(sprintf("  Total RR calls: %d\n", n_rr_calls))
cat(sprintf("  z_threshold levels: %d (post-hoc)\n", length(Z_RANGE)))
cat(sprintf("  Estimated time: ~90 minutes\n\n"))

# --- Results storage ---
all_results <- list()
result_idx <- 0
call_count <- 0

t_start <- Sys.time()

for (rep in seq_len(REPS)) {
  for (dc in seq_len(n_data_conds)) {
    nf  <- data_grid$nFactors[dc]
    ipf <- data_grid$items_per_factor[dc]
    n   <- data_grid$n_respondents[dc]
    pct <- data_grid$pct_careless[dc]
    total_items <- nf * ipf

    seed <- SEED_BASE + (rep - 1) * n_data_conds + dc

    # --- Generate data ---
    dat <- tryCatch({
      simulated_good_responses(
        nConstructs = nf,
        nItems = rep(ipf, nf),
        n = n,
        seed = seed
      )
    }, error = function(e) {
      cat(sprintf("  SKIP: nF=%d ipf=%d n=%d rep=%d — data generation failed: %s\n",
                  nf, ipf, n, rep, e$message))
      NULL
    })

    if (is.null(dat)) next

    # --- Inject careless ---
    inj <- inject_careless(
      dat,
      pct_careless = pct,
      careless_levels = CARELESS_LEVELS,
      seed = seed + 1000000
    )
    labels <- inj$labels

    # --- Get auto_z threshold (from true nF via LOESS, simulating perfect PA) ---
    auto_z_threshold <- get_calibrated_z(nf)

    # --- Run ReReReRe for each corProp ---
    for (cp in CORPROP_LEVELS) {
      rr <- tryCatch({
        ReReReRe(
          inj$data_corrupted,
          corProp = cp,
          iterations = ITERATIONS,
          min_pairs = MIN_PAIRS,
          align_signs = TRUE,
          progress = FALSE
        )
      }, error = function(e) {
        cat(sprintf("  SKIP RR: nF=%d ipf=%d n=%d pct=%.2f cp=%.2f rep=%d — %s\n",
                    nf, ipf, n, pct, cp, rep, e$message))
        NULL
      })

      if (is.null(rr)) next
      call_count <- call_count + 1

      # --- Evaluate at all z_thresholds (post-hoc) ---
      for (z_thr in Z_RANGE) {
        flagged <- as.integer(rr$z_score <= z_thr)

        tp <- sum(flagged == 1 & labels == 1)
        tn <- sum(flagged == 0 & labels == 0)
        fp <- sum(flagged == 1 & labels == 0)
        fn <- sum(flagged == 0 & labels == 1)

        mcc  <- compute_mcc(tp, tn, fp, fn)
        sens <- if (tp + fn > 0) tp / (tp + fn) else 0
        spec <- if (tn + fp > 0) tn / (tn + fp) else 0

        result_idx <- result_idx + 1
        all_results[[result_idx]] <- data.frame(
          nFactors = nf,
          items_per_factor = ipf,
          total_items = total_items,
          n_respondents = n,
          pct_careless = pct,
          corProp = cp,
          z_threshold = z_thr,
          is_auto_z = (abs(z_thr - auto_z_threshold) < 0.05),
          auto_z_value = auto_z_threshold,
          rep_id = rep,
          mcc = mcc,
          sensitivity = sens,
          specificity = spec,
          tp = tp, tn = tn, fp = fp, fn = fn,
          mean_z_good = mean(rr$z_score[labels == 0], na.rm = TRUE),
          mean_z_careless = mean(rr$z_score[labels == 1], na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # --- Progress per rep ---
  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  pct_done <- rep / REPS * 100
  eta <- elapsed / rep * (REPS - rep)
  cat(sprintf("[Rep %2d/%d | %5.1f%%] %d RR calls | %.1f min elapsed | ETA %.1f min\n",
              rep, REPS, pct_done, call_count, elapsed, eta))
}

elapsed_total <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat(sprintf("\n=== DONE in %.1f minutes | %d RR calls | %d result rows ===\n",
            elapsed_total, call_count, result_idx))

# --- Save raw results ---
raw_df <- do.call(rbind, all_results)
write.csv(raw_df, "multiverse_full_raw.csv", row.names = FALSE)
cat(sprintf("Saved multiverse_full_raw.csv (%d rows)\n", nrow(raw_df)))

# --- Best z per condition (oracle) ---
best_df <- do.call(rbind, lapply(
  split(raw_df, paste(raw_df$nFactors, raw_df$items_per_factor, raw_df$n_respondents,
                       raw_df$pct_careless, raw_df$corProp, raw_df$rep_id)),
  function(chunk) {
    best_idx <- which.max(chunk$mcc)
    chunk[best_idx, ]
  }
))
best_df <- best_df[order(best_df$nFactors, best_df$items_per_factor,
                          best_df$n_respondents, best_df$pct_careless,
                          best_df$corProp, best_df$rep_id), ]
write.csv(best_df, "multiverse_full_best.csv", row.names = FALSE)
cat(sprintf("Saved multiverse_full_best.csv (%d rows)\n", nrow(best_df)))

# --- Auto-z performance (extract rows where z_threshold matches auto_z) ---
auto_df <- raw_df[raw_df$is_auto_z == TRUE, ]
write.csv(auto_df, "multiverse_full_auto_z.csv", row.names = FALSE)
cat(sprintf("Saved multiverse_full_auto_z.csv (%d rows)\n", nrow(auto_df)))

# --- Summary tables ---
cat("\n=== SUMMARY: Mean MCC by nFactors x items_per_factor (oracle best z, corProp=0.05) ===\n")
sub <- best_df[best_df$corProp == 0.05, ]
agg <- aggregate(mcc ~ nFactors + items_per_factor, data = sub, FUN = mean)
agg_wide <- reshape(agg, idvar = "nFactors", timevar = "items_per_factor",
                     direction = "wide")
names(agg_wide) <- gsub("mcc\\.", "ipf_", names(agg_wide))
print(agg_wide, row.names = FALSE)

cat("\n=== SUMMARY: Mean MCC by corProp (all conditions) ===\n")
agg_cp <- aggregate(mcc ~ corProp, data = best_df, FUN = mean)
print(agg_cp, row.names = FALSE)

cat("\n=== SUMMARY: Auto-z vs oracle best z (corProp=0.05) ===\n")
# Auto-z MCC per nF
auto_sub <- auto_df[auto_df$corProp == 0.05, ]
auto_agg <- aggregate(mcc ~ nFactors, data = auto_sub, FUN = mean)
# Oracle MCC per nF
oracle_sub <- best_df[best_df$corProp == 0.05, ]
oracle_agg <- aggregate(mcc ~ nFactors, data = oracle_sub, FUN = mean)
comparison <- merge(auto_agg, oracle_agg, by = "nFactors", suffixes = c("_auto", "_oracle"))
comparison$gap <- comparison$mcc_oracle - comparison$mcc_auto
print(comparison, row.names = FALSE)

cat("\n=== SUMMARY: Auto-z MCC by nFactors x items_per_factor (corProp=0.05) ===\n")
auto_agg2 <- aggregate(mcc ~ nFactors + items_per_factor, data = auto_sub, FUN = mean)
auto_wide <- reshape(auto_agg2, idvar = "nFactors", timevar = "items_per_factor",
                      direction = "wide")
names(auto_wide) <- gsub("mcc\\.", "ipf_", names(auto_wide))
print(auto_wide, row.names = FALSE)

# --- Heatmap plot ---
cat("\nGenerating plots...\n")

# Plot 1: Oracle best MCC heatmap (nF x ipf, corProp=0.05, aggregated over n and pct)
png("plot_multiverse_oracle_heatmap.png", width = 1200, height = 800, res = 150)
par(mar = c(5, 5, 4, 6))

oracle_heat <- aggregate(mcc ~ nFactors + items_per_factor, data = oracle_sub, FUN = mean)
mat_heat <- matrix(NA, nrow = length(NF_LEVELS), ncol = length(IPF_LEVELS))
for (i in seq_along(NF_LEVELS)) {
  for (j in seq_along(IPF_LEVELS)) {
    val <- oracle_heat$mcc[oracle_heat$nFactors == NF_LEVELS[i] &
                            oracle_heat$items_per_factor == IPF_LEVELS[j]]
    if (length(val) == 1) mat_heat[i, j] <- val
  }
}
rownames(mat_heat) <- NF_LEVELS
colnames(mat_heat) <- IPF_LEVELS

# Color scale
n_colors <- 100
colors <- colorRampPalette(c("#fee8c8", "#fdbb84", "#e34a33", "#b30000"))(n_colors)
val_range <- range(mat_heat, na.rm = TRUE)
get_color <- function(v) {
  idx <- round((v - val_range[1]) / (val_range[2] - val_range[1]) * (n_colors - 1)) + 1
  colors[max(1, min(n_colors, idx))]
}

plot(NULL, xlim = c(0.5, length(IPF_LEVELS) + 0.5),
     ylim = c(0.5, length(NF_LEVELS) + 0.5),
     xlab = "Items per factor", ylab = "Number of factors (nF)",
     main = "Oracle best MCC by nFactors x items/factor\n(corProp=0.05, averaged over n and pct_careless)",
     xaxt = "n", yaxt = "n", asp = NA)

axis(1, at = seq_along(IPF_LEVELS), labels = IPF_LEVELS)
axis(2, at = seq_along(NF_LEVELS), labels = NF_LEVELS, las = 1)

for (i in seq_along(NF_LEVELS)) {
  for (j in seq_along(IPF_LEVELS)) {
    if (!is.na(mat_heat[i, j])) {
      rect(j - 0.45, i - 0.45, j + 0.45, i + 0.45,
           col = get_color(mat_heat[i, j]), border = "white", lwd = 1)
      text(j, i, sprintf("%.2f", mat_heat[i, j]), cex = 0.7,
           col = if (mat_heat[i, j] > 0.35) "white" else "black")
    }
  }
}

mtext(sprintf("R=%d reps | Oracle best z per cell", REPS),
      side = 1, line = 3.5, cex = 0.75, col = "gray40")
dev.off()
cat("Saved plot_multiverse_oracle_heatmap.png\n")

# Plot 2: Auto-z vs Oracle comparison by nF
png("plot_auto_z_vs_oracle.png", width = 1200, height = 700, res = 150)
par(mar = c(5, 5, 3, 1))

plot(comparison$nFactors, comparison$mcc_oracle, type = "b",
     pch = 16, col = "darkred", lwd = 2, cex = 1.3,
     xlab = "Number of Factors (nF)", ylab = "Mean MCC",
     main = "Auto-calibrated z vs Oracle best z (corProp=0.05)",
     ylim = c(0, max(comparison$mcc_oracle) * 1.1))
lines(comparison$nFactors, comparison$mcc_auto, type = "b",
      pch = 17, col = "steelblue", lwd = 2, cex = 1.3)
legend("topleft",
       legend = c("Oracle best z", "Auto-calibrated z"),
       pch = c(16, 17), col = c("darkred", "steelblue"),
       lwd = 2, cex = 0.9, bg = "white")

# Add gap labels
for (i in seq_len(nrow(comparison))) {
  if (comparison$gap[i] > 0.01) {
    text(comparison$nFactors[i], comparison$mcc_auto[i] - 0.015,
         sprintf("-%.2f", comparison$gap[i]), cex = 0.6, col = "gray50")
  }
}

mtext(sprintf("R=%d reps | Averaged over n, pct_careless, items/factor", REPS),
      side = 1, line = 3.5, cex = 0.75, col = "gray40")
dev.off()
cat("Saved plot_auto_z_vs_oracle.png\n")

# Plot 3: Effect of items_per_factor on auto-z
png("plot_ipf_effect_auto_z.png", width = 1200, height = 700, res = 150)
par(mar = c(5, 5, 3, 1))

ipf_cols <- c("3" = "#e34a33", "6" = "steelblue", "10" = "#31a354")
auto_by_ipf <- aggregate(mcc ~ nFactors + items_per_factor,
                          data = auto_sub, FUN = mean)

plot(NULL, xlim = range(NF_LEVELS), ylim = c(0, 0.6),
     xlab = "Number of Factors (nF)", ylab = "Mean MCC (auto-z)",
     main = "Auto-z performance by items/factor")

for (ipf_val in IPF_LEVELS) {
  sub_ipf <- auto_by_ipf[auto_by_ipf$items_per_factor == ipf_val, ]
  lines(sub_ipf$nFactors, sub_ipf$mcc, type = "b",
        pch = 16, col = ipf_cols[as.character(ipf_val)], lwd = 2, cex = 1.2)
}

legend("topleft",
       legend = paste("ipf =", IPF_LEVELS),
       col = unname(ipf_cols), pch = 16, lwd = 2, cex = 0.9, bg = "white")

mtext(sprintf("R=%d reps | corProp=0.05 | Averaged over n, pct_careless", REPS),
      side = 1, line = 3.5, cex = 0.75, col = "gray40")
dev.off()
cat("Saved plot_ipf_effect_auto_z.png\n")

cat("\n=== ALL DONE ===\n")
