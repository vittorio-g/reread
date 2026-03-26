#### Multiverse: nFactors x Z-threshold ####
# Varia solo:
#   - nFactors: 2 a 40 (step 1)
#   - z_threshold: 0.1 a 3 (step 0.2)
# Tutti gli altri parametri fissi.
# Output: per ogni nFactors, il z_threshold che massimizza MCC.
# Plot: heatmap + line plot del best z per nFactors.

rm(list = ls())

library(dplyr)
library(magrittr)
library(lavaan)
library(psych)
library(ggplot2)

# setwd — use script location when running from RStudio, or current dir otherwise
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  # When running via Rscript, use the script's own directory
  args <- commandArgs(trailingOnly = FALSE)
  script_arg <- grep("--file=", args, value = TRUE)
  if (length(script_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("--file=", "", script_arg)))
    setwd(script_dir)
  }
}

source("Synthetic_Good_Responses_2.R")
source("Careless_machine_2.R")
source("ReReReRe.R")

# ============================================================
# PARAMETERS
# ============================================================

nFactors_grid   <- 2:40                          # 39 levels
z_threshold_grid <- seq(0.1, 3, by = 0.2)        # 15 levels

# Fixed parameters
N_RESPONDENTS   <- 300
PCT_CARELESS    <- 0.20
COR_PROP        <- 0.05
ITERATIONS      <- 100
MIN_PAIRS       <- 15
CARELESS_TYPES  <- c(random = 1/3, longstring = 1/3, mixed = 1/3)
CARELESS_LEVELS <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5)
ITEMS_CYCLE     <- c(10, 6, 3)
EVAL_THRESHOLD  <- 0.01
R_REPLICATIONS  <- 10   # replications per nFactors for stability
SEED_BASE       <- 123

make_nItems <- function(nFactors) rep_len(ITEMS_CYCLE, nFactors)

# MCC helper
calc_mcc <- function(tp, tn, fp, fn) {
  denom <- sqrt(as.numeric(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  if (denom == 0) return(0)
  (tp * tn - fp * fn) / denom
}

# ============================================================
# DESIGN
# ============================================================

total_cells <- length(nFactors_grid) * R_REPLICATIONS
cat(sprintf("nFactors: %d-%d (%d levels)\n", min(nFactors_grid), max(nFactors_grid), length(nFactors_grid)))
cat(sprintf("z_threshold: %.1f-%.1f (%d levels)\n", min(z_threshold_grid), max(z_threshold_grid), length(z_threshold_grid)))
cat(sprintf("Replications: %d\n", R_REPLICATIONS))
cat(sprintf("Total ReReReRe calls: %d\n", total_cells))
cat(sprintf("Fixed: n=%d, pct_careless=%.2f, corProp=%.2f\n\n", N_RESPONDENTS, PCT_CARELESS, COR_PROP))

# ============================================================
# MAIN LOOP
# ============================================================

all_results <- list()
row_idx <- 0
time_start <- Sys.time()

for (nF in nFactors_grid) {

  nItems_vec <- make_nItems(nF)
  total_items <- sum(nItems_vec)

  for (rep_id in 1:R_REPLICATIONS) {

    row_idx <- row_idx + 1
    elapsed <- as.numeric(difftime(Sys.time(), time_start, units = "mins"))

    if (row_idx > 1) {
      eta <- elapsed / (row_idx - 1) * (total_cells - row_idx + 1)
      eta_str <- if (eta > 60) sprintf("%.1f h left", eta/60) else sprintf("%.0f min left", eta)
    } else {
      eta_str <- "estimating..."
    }

    cat(sprintf("[%4.1f%%] nF=%2d rep=%2d | items=%3d | %.1f min | %s\n",
                100 * (row_idx - 1) / total_cells, nF, rep_id, total_items, elapsed, eta_str))

    cell_seed <- SEED_BASE + (nF - min(nFactors_grid)) * R_REPLICATIONS + rep_id
    set.seed(cell_seed)

    # Generate data
    clean_data <- tryCatch(
      simulated_good_responses(nConstructs = nF, nItems = nItems_vec, n = N_RESPONDENTS),
      error = function(e) { cat("  WARN gen:", e$message, "\n"); NULL }
    )
    if (is.null(clean_data)) next

    # Inject careless
    injection <- tryCatch(
      inject_careless(data = clean_data, pct_careless = PCT_CARELESS,
                      pct_types = CARELESS_TYPES, careless_levels = CARELESS_LEVELS),
      error = function(e) { cat("  WARN inject:", e$message, "\n"); NULL }
    )
    if (is.null(injection)) next

    corrupted_data <- injection$data_corrupted
    labels <- injection$labels
    is_careless <- labels$careless_pct >= EVAL_THRESHOLD

    # Run ReReReRe
    rr <- tryCatch(
      ReReReRe(data = corrupted_data, corProp = COR_PROP,
               cutOff = 0, iterations = ITERATIONS,
               min_pairs = MIN_PAIRS, align_signs = TRUE, progress = FALSE),
      error = function(e) { cat("  WARN RR:", e$message, "\n"); NULL }
    )
    if (is.null(rr)) next

    # Evaluate all z_thresholds
    for (zt in z_threshold_grid) {
      flagged <- rr$z_score <= zt
      tp <- sum(flagged & is_careless)
      tn <- sum(!flagged & !is_careless)
      fp <- sum(flagged & !is_careless)
      fn <- sum(!flagged & is_careless)

      mcc <- calc_mcc(tp, tn, fp, fn)
      sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA
      spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA

      all_results[[length(all_results) + 1]] <- data.frame(
        nFactors = nF,
        total_items = total_items,
        rep_id = rep_id,
        z_threshold = zt,
        mcc = mcc,
        sensitivity = sens,
        specificity = spec,
        tp = tp, tn = tn, fp = fp, fn = fn,
        stringsAsFactors = FALSE
      )
    }
  }
}

elapsed_total <- difftime(Sys.time(), time_start, units = "mins")
cat(sprintf("\nDone in %.1f minutes.\n", as.numeric(elapsed_total)))

# ============================================================
# AGGREGATE RESULTS
# ============================================================

results_df <- do.call(rbind, all_results)
write.csv(results_df, "multiverse_nF_x_Z_raw.csv", row.names = FALSE)

# Average MCC across replications for each (nFactors, z_threshold)
agg <- results_df %>%
  group_by(nFactors, total_items, z_threshold) %>%
  summarise(
    mean_mcc = mean(mcc, na.rm = TRUE),
    sd_mcc = sd(mcc, na.rm = TRUE),
    mean_sens = mean(sensitivity, na.rm = TRUE),
    mean_spec = mean(specificity, na.rm = TRUE),
    n_reps = n(),
    .groups = "drop"
  )

write.csv(agg, "multiverse_nF_x_Z_aggregated.csv", row.names = FALSE)

# Best z_threshold per nFactors (highest mean MCC)
best_z <- agg %>%
  group_by(nFactors, total_items) %>%
  slice_max(mean_mcc, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(nFactors)

write.csv(best_z, "multiverse_nF_x_Z_best.csv", row.names = FALSE)

cat("\nBest z_threshold per nFactors:\n")
print(as.data.frame(best_z[, c("nFactors", "total_items", "z_threshold", "mean_mcc")]), row.names = FALSE)

# ============================================================
# PLOTS
# ============================================================

# 1. Line plot: best z_threshold per nFactors
p1 <- ggplot(best_z, aes(x = nFactors, y = z_threshold)) +
  geom_line(color = "#2563EB", linewidth = 1.2) +
  geom_point(aes(size = mean_mcc), color = "#2563EB", alpha = 0.8) +
  scale_size_continuous(name = "MCC", range = c(1.5, 5)) +
  scale_x_continuous(breaks = seq(2, 40, by = 2)) +
  scale_y_continuous(breaks = z_threshold_grid) +
  labs(
    title = "Optimal Z-threshold by Number of Factors",
    subtitle = sprintf("n=%d, pct_careless=%.0f%%, corProp=%.2f, %d reps per cell",
                        N_RESPONDENTS, PCT_CARELESS*100, COR_PROP, R_REPLICATIONS),
    x = "Number of Factors",
    y = "Best Z-threshold (maximizes MCC)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave("plot_best_z_by_nFactors.png", p1, width = 12, height = 6, dpi = 200)
cat("\nSaved: plot_best_z_by_nFactors.png\n")

# 2. Heatmap: MCC for all (nFactors, z_threshold) combinations
p2 <- ggplot(agg, aes(x = nFactors, y = z_threshold, fill = mean_mcc)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean MCC", option = "inferno", limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(2, 40, by = 2)) +
  scale_y_continuous(breaks = z_threshold_grid) +
  # Overlay best z per nFactors
  geom_point(data = best_z, aes(x = nFactors, y = z_threshold),
             color = "white", shape = 4, size = 2.5, stroke = 1.2, inherit.aes = FALSE) +
  labs(
    title = "MCC Heatmap: Number of Factors x Z-threshold",
    subtitle = sprintf("White X = optimal z per nFactors | n=%d, %d reps",
                        N_RESPONDENTS, R_REPLICATIONS),
    x = "Number of Factors",
    y = "Z-threshold"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave("plot_heatmap_nF_x_Z.png", p2, width = 14, height = 7, dpi = 200)
cat("Saved: plot_heatmap_nF_x_Z.png\n")

cat("\nAll done!\n")
