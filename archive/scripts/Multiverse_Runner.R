#### Multiverse Runner for ReReReRe ####
# Runs the full simulation design and saves results incrementally.
#
# Pipeline per cell:
#   simulated_good_responses() -> inject_careless() -> ReReReRe() -> metrics
#
# Optimization: ReReReRe is run once per (dataset x corProp).
# z_threshold values are applied post-hoc (free), saving compute.
#
# Evaluation design:
#   - Careless respondents have corruption at moderate-to-high levels (50-100%)
#   - EVAL_THRESHOLD = 0.01: all injected careless respondents count as positives
#   - This matches the original multiverse setup for validation purposes
#
# === REVISION 2026-03-10f: SWITCH TO Z-SCORE METRIC ===
# After fixing the longstring bug (2026-03-10e), the metric comparison reversed:
# - Z-score/indCors now beat percentile for nF>=10 (Cohen's d up to 3.1 vs 1.7)
# - The previous percentile "advantage" was an artifact of the contiguous-block bug
# - Z-score = (coupled_cor - mean_random) / sd_random: how many SDs above
#   the respondent's own random baseline is their coupled correlation
# - Flagging: flag if z_score <= z_threshold (low z = likely careless)
#
# Previous percentile-based multiverse archived in archive/pre-revision-2026-03-10f/
# ===

rm(list = ls())

library(dplyr)
library(magrittr)
library(lavaan)
library(psych)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# --- Source current functions ---
source("Synthetic_Good_Responses_2.R")
source("Careless_machine_2.R")
source("ReReReRe.R")

# ============================================================
# MULTIVERSE PARAMETERS
# ============================================================

params <- list(
  nFactors       = c(4, 8, 10, 15, 20, 25),    # 29–158 items (cycling 10,6,3)
  n_respondents  = c(50, 100, 300, 1000),
  pct_careless   = c(0.05, 0.10, 0.25, 0.50),
  corProp        = c(0.05, 0.10, 0.15, 0.20, 0.30),
  z_threshold    = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5)  # applied post-hoc (free)
)
MIN_PAIRS       <- 15  # passed to ReReReRe()

# Fixed parameters
CARELESS_TYPES  <- c(random = 1/3, longstring = 1/3, mixed = 1/3)
CARELESS_LEVELS <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5)  # 50-100% corruption
ITERATIONS      <- 100
ITEMS_CYCLE     <- c(10, 6, 3)  # items per factor, cycling in order
EVAL_THRESHOLD  <- 0.01  # effectively all injected careless = positive
                         # (minimum careless_levels is 0.50, so 0.01 catches them all)

# Helper: generate nItems vector for a given nFactors
make_nItems <- function(nFactors) {
  rep_len(ITEMS_CYCLE, nFactors)
}

# ============================================================
# HELPER: compute classification metrics using z_score + fixed threshold
# ============================================================

compute_metrics <- function(z_scores, labels_df, z_threshold_values, eval_threshold) {
  # z_scores:           numeric vector of ReReReRe z-scores per respondent
  # labels_df:          data.frame with columns: careless (logical), careless_pct, pattern
  # z_threshold_values: vector of z-score thresholds to try
  # eval_threshold:     minimum corruption % to count as "truly careless"
  #
  # Fixed threshold flagging: flag if z_score <= z_threshold
  # Z-score = (coupled_cor - mean_random) / sd_random. Lower = more likely careless.

  actual <- labels_df$careless_pct >= eval_threshold  # TRUE = truly careless

  lapply(z_threshold_values, function(zt) {
    flagged <- z_scores <= zt  # TRUE = predicted careless (low z-score)

    TP <- sum(flagged & actual)
    TN <- sum(!flagged & !actual)
    FP <- sum(flagged & !actual)
    FN <- sum(!flagged & actual)
    total <- length(flagged)

    sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
    specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
    precision   <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
    f1 <- if (!is.na(precision) && !is.na(sensitivity) && (precision + sensitivity) > 0)
            2 * precision * sensitivity / (precision + sensitivity) else NA_real_
    youden <- if (!is.na(sensitivity) && !is.na(specificity))
                sensitivity + specificity - 1 else NA_real_

    # Matthews Correlation Coefficient
    mcc_denom <- sqrt(as.numeric(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    mcc <- if (mcc_denom > 0) (TP * TN - FP * FN) / mcc_denom else NA_real_

    data.frame(
      z_threshold = zt,
      TP = TP, TN = TN, FP = FP, FN = FN,
      n_truly_careless = sum(actual),
      accuracy = (TP + TN) / total,
      sensitivity = sensitivity,
      specificity = specificity,
      precision = precision,
      f1 = f1,
      youden = youden,
      mcc = mcc,
      false_positive_rate = if ((FP + TN) > 0) FP / (FP + TN) else NA_real_,
      false_negative_rate = if ((TP + FN) > 0) FN / (TP + FN) else NA_real_,
      n_flagged = sum(flagged),
      mean_z_good = mean(z_scores[!actual], na.rm = TRUE),
      mean_z_careless = mean(z_scores[actual], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }) %>% do.call(rbind, .)
}

# Also compute detection rate by careless pattern type
compute_metrics_by_type <- function(z_scores, labels_df, z_threshold_values, eval_threshold) {
  pattern <- labels_df$pattern
  truly_careless <- labels_df$careless_pct >= eval_threshold

  lapply(z_threshold_values, function(zt) {
    flagged <- z_scores <= zt

    # Only look at truly careless respondents (above eval threshold)
    if (sum(truly_careless) == 0) return(NULL)

    types <- unique(pattern[truly_careless])
    types <- types[types != "clean"]  # exclude clean label
    lapply(types, function(typ) {
      mask <- truly_careless & (pattern == typ)
      if (sum(mask) == 0) return(NULL)
      data.frame(
        z_threshold = zt,
        pattern = typ,
        n_cases = sum(mask),
        n_detected = sum(flagged[mask]),
        detection_rate = mean(flagged[mask]),
        mean_z = mean(z_scores[mask], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}

# ============================================================
# BUILD THE GRID
# ============================================================

# Data-generation conditions (each needs a unique dataset + ReReReRe run per corProp)
data_grid <- expand.grid(
  nFactors = params$nFactors,
  n_respondents = params$n_respondents,
  pct_careless = params$pct_careless,
  stringsAsFactors = FALSE
)

# Pre-compute total items and nItems vectors for each nFactors level
data_grid$total_items <- sapply(data_grid$nFactors, function(nf) sum(make_nItems(nf)))
data_grid$n_lt_p <- data_grid$n_respondents < data_grid$total_items

cat("=== MULTIVERSE DESIGN (z-score + fixed threshold) ===\n")
cat("Data-generation conditions:", nrow(data_grid), "\n")
cat("corProp levels:", length(params$corProp), "\n")
cat("z_threshold levels:", length(params$z_threshold), "(applied post-hoc)\n")
cat("Total ReReReRe runs:", nrow(data_grid) * length(params$corProp), "\n")
cat("Total result rows:", nrow(data_grid) * length(params$corProp) * length(params$z_threshold), "\n")
cat("Careless levels:", paste(CARELESS_LEVELS, collapse=", "), "\n")
cat("Eval threshold:", EVAL_THRESHOLD, "(only corruption >=", EVAL_THRESHOLD*100, "% = truly careless)\n\n")

# ============================================================
# OUTPUT FILES
# ============================================================

results_file    <- "multiverse_results.csv"
by_type_file    <- "multiverse_results_by_type.csv"
checkpoint_file <- "multiverse_checkpoint.rds"

# Check for existing checkpoint to resume from
start_row <- 1
if (file.exists(checkpoint_file)) {
  cp <- readRDS(checkpoint_file)
  start_row <- cp$next_row
  cat("Resuming from checkpoint: row", start_row, "of", nrow(data_grid), "\n\n")
}

# Initialize output files if starting fresh
if (start_row == 1) {
  header_main <- paste(
    "nFactors", "n_respondents", "total_items", "n_lt_p",
    "pct_careless", "corProp", "z_threshold",
    "TP", "TN", "FP", "FN", "n_truly_careless",
    "accuracy", "sensitivity", "specificity",
    "precision", "f1", "youden", "mcc",
    "false_positive_rate", "false_negative_rate", "n_flagged",
    "mean_z_good", "mean_z_careless",
    sep = ","
  )
  writeLines(header_main, results_file)

  header_type <- paste(
    "nFactors", "n_respondents", "pct_careless", "corProp", "z_threshold",
    "pattern", "n_cases", "n_detected", "detection_rate", "mean_z",
    sep = ","
  )
  writeLines(header_type, by_type_file)
}

# ============================================================
# MAIN LOOP
# ============================================================

total_runs <- nrow(data_grid) * length(params$corProp)
run_counter <- (start_row - 1) * length(params$corProp)
time_start <- Sys.time()

for (dg_idx in start_row:nrow(data_grid)) {

  dg <- data_grid[dg_idx, ]

  # --- Progress report ---
  pct_done <- round(100 * (run_counter / total_runs), 1)
  elapsed <- as.numeric(difftime(Sys.time(), time_start, units = "mins"))
  if (run_counter > 0) {
    rate <- elapsed / run_counter  # minutes per run
    remaining <- rate * (total_runs - run_counter)
    eta_str <- sprintf("%.0f min remaining", remaining)
  } else {
    eta_str <- "estimating..."
  }

  nItems_vec <- make_nItems(dg$nFactors)
  total_items <- dg$total_items
  n_lt_p <- dg$n_lt_p

  cat(sprintf("[%5.1f%%] Row %d/%d | nF=%d items=%s n=%d pct=%.2f | elapsed=%.1fmin | %s\n",
              pct_done, dg_idx, nrow(data_grid),
              dg$nFactors, paste(nItems_vec, collapse=","), dg$n_respondents, dg$pct_careless,
              elapsed, eta_str))

  # --- Generate clean dataset ---
  clean_data <- tryCatch(
    simulated_good_responses(
      nConstructs = dg$nFactors,
      nItems = nItems_vec,
      n = dg$n_respondents
    ),
    error = function(e) {
      cat("  WARNING: data generation failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (is.null(clean_data)) {
    run_counter <- run_counter + length(params$corProp)
    next
  }

  # --- Inject careless responses ---
  injection <- tryCatch(
    inject_careless(
      data = clean_data,
      pct_careless = dg$pct_careless,
      pct_types = CARELESS_TYPES,
      careless_levels = CARELESS_LEVELS
    ),
    error = function(e) {
      cat("  WARNING: careless injection failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (is.null(injection)) {
    run_counter <- run_counter + length(params$corProp)
    next
  }

  corrupted_data <- injection$data_corrupted
  labels <- injection$labels

  # --- Run ReReReRe for each corProp (z_threshold applied post-hoc) ---
  for (cp in params$corProp) {

    rr_result <- tryCatch(
      ReReReRe(
        data = corrupted_data,
        corProp = cp,
        cutOff = 0,  # dummy, we apply z_threshold post-hoc
        iterations = ITERATIONS,
        min_pairs = MIN_PAIRS,
        progress = FALSE
      ),
      error = function(e) {
        cat("  WARNING: ReReReRe failed (corProp=", cp, "):", conditionMessage(e), "\n")
        return(NULL)
      }
    )

    if (is.null(rr_result)) {
      run_counter <- run_counter + 1
      next
    }

    # Apply all z_threshold values post-hoc using z-scores
    metrics <- compute_metrics(rr_result$z_score, labels, params$z_threshold, EVAL_THRESHOLD)
    metrics$nFactors <- dg$nFactors
    metrics$n_respondents <- dg$n_respondents
    metrics$total_items <- total_items
    metrics$n_lt_p <- n_lt_p
    metrics$pct_careless <- dg$pct_careless
    metrics$corProp <- cp

    # Reorder columns
    metrics <- metrics[, c("nFactors", "n_respondents",
                           "total_items", "n_lt_p",
                           "pct_careless", "corProp", "z_threshold",
                           "TP", "TN", "FP", "FN", "n_truly_careless",
                           "accuracy", "sensitivity", "specificity",
                           "precision", "f1", "youden", "mcc",
                           "false_positive_rate", "false_negative_rate",
                           "n_flagged",
                           "mean_z_good", "mean_z_careless")]

    # Append to CSV
    write.table(metrics, results_file, append = TRUE, sep = ",",
                row.names = FALSE, col.names = FALSE)

    # By-type metrics
    type_metrics <- compute_metrics_by_type(rr_result$z_score, labels, params$z_threshold, EVAL_THRESHOLD)
    if (!is.null(type_metrics) && nrow(type_metrics) > 0) {
      type_metrics$nFactors <- dg$nFactors
      type_metrics$n_respondents <- dg$n_respondents
      type_metrics$pct_careless <- dg$pct_careless
      type_metrics$corProp <- cp

      type_metrics <- type_metrics[, c("nFactors", "n_respondents",
                                       "pct_careless", "corProp", "z_threshold",
                                       "pattern", "n_cases", "n_detected",
                                       "detection_rate", "mean_z")]

      write.table(type_metrics, by_type_file, append = TRUE, sep = ",",
                  row.names = FALSE, col.names = FALSE)
    }

    run_counter <- run_counter + 1
  }

  # --- Save checkpoint after each data-generation row ---
  saveRDS(list(next_row = dg_idx + 1, timestamp = Sys.time()), checkpoint_file)
}

# ============================================================
# DONE
# ============================================================

elapsed_total <- difftime(Sys.time(), time_start, units = "mins")

cat("\n=== MULTIVERSE COMPLETE ===\n")
cat("Total time:", round(as.numeric(elapsed_total), 1), "minutes\n")
cat("Results saved to:", results_file, "\n")
cat("By-type results saved to:", by_type_file, "\n")

# Clean up checkpoint
if (file.exists(checkpoint_file)) file.remove(checkpoint_file)
