#### Final Multiverse Runner for ReReReRe ####
# Publication-quality simulation with replications and finer parameter grid.
#
# Key differences from pilot (Multiverse_Runner.R):
#   - R replications per cell (different random datasets per replication)
#   - Finer nFactors grid (9 levels) — pins down crossover with Mahalanobis
#   - Finer n_respondents grid (6 levels)
#   - Finer corProp grid (6 levels) — resolution around optimal 0.05
#   - More z_threshold levels (11, post-hoc)
#   - rep_id column in output for replication-level analysis
#   - Per-replication AUC for ReReReRe (threshold-free metric)
#   - Reproducible seeds per cell × replication
#   - Better ETA estimation (weighted by completion)
#   - Benchmark methods (Mahalanobis, PTC, LongString, IRV) on same datasets
#
# Pipeline per cell × replication:
#   simulated_good_responses() → inject_careless() →
#     benchmarks (once per dataset) + ReReReRe (× corProp) → metrics
#
# Output files:
#   1. final_multiverse_results.csv       — ReReReRe results (per corProp × z_threshold)
#   2. final_multiverse_results_by_type.csv — detection by careless pattern type
#   3. final_benchmark_results.csv        — benchmark methods (once per dataset, oracle threshold)
#
# === FINAL ANALYSES (2026-03-10h) ===

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
# BENCHMARK METHOD FUNCTIONS (computed on same datasets)
# ============================================================
# These run once per dataset (independent of corProp/z_threshold).
# Allows direct comparison with ReReReRe on identical data.

# 1. Mahalanobis distance — HIGH = outlier = likely careless
compute_mahad <- function(data) {
  x <- as.matrix(data)
  center <- colMeans(x)
  cv <- cov(x)
  tryCatch(
    mahalanobis(x, center, cv),
    error = function(e) {
      # Regularize for singular covariance (n < p)
      tryCatch(
        mahalanobis(x, center, cv + diag(1e-3, ncol(x))),
        error = function(e2) rep(NA_real_, nrow(x))
      )
    }
  )
}

# 2. Person-total correlation — LOW = likely careless
compute_ptc <- function(data) {
  x <- as.matrix(data)
  col_means <- colMeans(x)
  apply(x, 1, function(row) cor(row, col_means))
}

# 3. LongString — max consecutive identical responses — HIGH = likely careless
compute_longstring <- function(data) {
  x <- as.matrix(data)
  apply(x, 1, function(row) {
    r <- rle(as.integer(row))
    max(r$lengths)
  })
}

# 4. IRV — intra-individual response variability (row SD) — LOW = likely careless
compute_irv <- function(data) {
  apply(as.matrix(data), 1, sd)
}

# ============================================================
# BENCHMARK EVALUATION HELPERS
# ============================================================

calc_mcc <- function(tp, tn, fp, fn) {
  denom <- sqrt(as.numeric(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  if (denom == 0) return(0)
  (tp * tn - fp * fn) / denom
}

calc_auc_generic <- function(scores, is_careless, higher_is_careless) {
  if (any(is.na(scores))) return(NA_real_)
  if (!higher_is_careless) scores <- -scores
  n_pos <- sum(is_careless)
  n_neg <- sum(!is_careless)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  ranks <- rank(scores, ties.method = "average")
  (sum(ranks[is_careless]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

find_best_threshold <- function(scores, is_careless, higher_is_careless,
                                n_candidates = 200) {
  if (any(is.na(scores))) {
    return(list(mcc = NA_real_, threshold = NA_real_,
                sensitivity = NA_real_, specificity = NA_real_))
  }
  candidates <- quantile(scores, probs = seq(0.005, 0.995, length.out = n_candidates),
                          na.rm = TRUE)
  candidates <- unique(candidates)
  best_mcc <- -Inf; best_thresh <- NA_real_; best_sens <- NA_real_; best_spec <- NA_real_
  for (thresh in candidates) {
    flagged <- if (higher_is_careless) scores >= thresh else scores <= thresh
    tp <- sum(flagged & is_careless); tn <- sum(!flagged & !is_careless)
    fp <- sum(flagged & !is_careless); fn <- sum(!flagged & is_careless)
    mcc <- calc_mcc(tp, tn, fp, fn)
    if (mcc > best_mcc) {
      best_mcc <- mcc; best_thresh <- thresh
      n_pos <- tp + fn
      best_sens <- if (n_pos > 0) tp / n_pos else NA_real_
      best_spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    }
  }
  list(mcc = best_mcc, threshold = best_thresh,
       sensitivity = best_sens, specificity = best_spec)
}

# ============================================================
# MULTIVERSE PARAMETERS — FINAL DESIGN
# ============================================================

params <- list(
  # --- Varied parameters ---
  nFactors       = c(4, 6, 8, 10, 12, 15, 20, 25),  # 8 levels (finer around crossover)
  n_respondents  = c(50, 100, 200, 300, 500, 1000),       # 5 levels (added 200, 500)
  pct_careless   = c(0.05, 0.10, 0.25, 0.50),             # 4 levels (same as pilot)
  corProp        = c(0.03, 0.05, 0.07, 0.10, 0.15, 0.20), # 6 levels (finer around optimal)
  z_threshold    = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)  # 11 levels (post-hoc)
)

R_REPLICATIONS  <- 50     # replications per cell
MIN_PAIRS       <- 15     # passed to ReReReRe()

# Fixed parameters
CARELESS_TYPES  <- c(random = 1/3, longstring = 1/3, mixed = 1/3)
CARELESS_LEVELS <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5)  # 50-100% corruption
ITERATIONS      <- 100
ITEMS_CYCLE     <- c(10, 6, 3)  # items per factor, cycling in order
EVAL_THRESHOLD  <- 0.01  # all injected careless = positive (min corruption is 50%)

SEED_BASE       <- 42    # for reproducibility

# Helper: generate nItems vector for a given nFactors
make_nItems <- function(nFactors) {
  rep_len(ITEMS_CYCLE, nFactors)
}

# Benchmark method definitions (name, function, direction)
BENCHMARK_METHODS <- list(
  list(name = "Mahalanobis", fn = compute_mahad,       higher = TRUE),
  list(name = "PersonTotal", fn = compute_ptc,         higher = FALSE),
  list(name = "LongString",  fn = compute_longstring,  higher = TRUE),
  list(name = "IRV",         fn = compute_irv,         higher = FALSE)
)

# ============================================================
# DESIGN SUMMARY
# ============================================================

n_data_conditions <- length(params$nFactors) * length(params$n_respondents) *
                     length(params$pct_careless)
n_corProp <- length(params$corProp)
n_zt <- length(params$z_threshold)
n_rr_calls_per_rep <- n_data_conditions * n_corProp
n_rr_calls_total <- n_rr_calls_per_rep * R_REPLICATIONS
n_result_rows <- n_rr_calls_total * n_zt
n_bench_rows <- n_data_conditions * R_REPLICATIONS * length(BENCHMARK_METHODS)

cat("===============================================================\n")
cat("  FINAL MULTIVERSE DESIGN (z-score + fixed threshold)\n")
cat("===============================================================\n")
cat(sprintf("  nFactors:       %d levels  %s\n", length(params$nFactors),
            paste(params$nFactors, collapse = ", ")))
cat(sprintf("  n_respondents:  %d levels  %s\n", length(params$n_respondents),
            paste(params$n_respondents, collapse = ", ")))
cat(sprintf("  pct_careless:   %d levels  %s\n", length(params$pct_careless),
            paste(params$pct_careless, collapse = ", ")))
cat(sprintf("  corProp:        %d levels  %s\n", n_corProp,
            paste(params$corProp, collapse = ", ")))
cat(sprintf("  z_threshold:    %d levels  (post-hoc, free)\n", n_zt))
cat(sprintf("  Replications:   %d\n", R_REPLICATIONS))
cat("---------------------------------------------------------------\n")
cat(sprintf("  Data conditions per rep:   %d\n", n_data_conditions))
cat(sprintf("  ReReReRe calls per rep:    %d  (x %d corProp)\n", n_rr_calls_per_rep, n_corProp))
cat(sprintf("  Total ReReReRe calls:      %s\n", format(n_rr_calls_total, big.mark = ",")))
cat(sprintf("  Total RR result rows:      %s  (x %d z_threshold)\n",
            format(n_result_rows, big.mark = ","), n_zt))
cat(sprintf("  Benchmark methods:         %d  (Mahalanobis, PTC, LongString, IRV)\n",
            length(BENCHMARK_METHODS)))
cat(sprintf("  Total benchmark rows:      %s  (once per dataset, oracle threshold)\n",
            format(n_bench_rows, big.mark = ",")))
cat(sprintf("  Estimated time:            %.1f-%.1f hours (at 0.5-0.75 s/RR call)\n",
            n_rr_calls_total * 0.5 / 3600, n_rr_calls_total * 0.75 / 3600))
cat("===============================================================\n\n")

# ============================================================
# HELPER: compute classification metrics (same as pilot)
# ============================================================

compute_metrics <- function(z_scores, labels_df, z_threshold_values, eval_threshold) {
  actual <- labels_df$careless_pct >= eval_threshold

  lapply(z_threshold_values, function(zt) {
    flagged <- z_scores <= zt

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

# AUC via Mann-Whitney U (threshold-free)
compute_auc <- function(z_scores, labels_df, eval_threshold) {
  actual <- labels_df$careless_pct >= eval_threshold
  if (any(is.na(z_scores))) return(NA_real_)
  n_pos <- sum(actual)
  n_neg <- sum(!actual)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)

  # Lower z-score = more likely careless, so negate for standard AUC
  scores <- -z_scores
  ranks <- rank(scores, ties.method = "average")
  (sum(ranks[actual]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

# Detection rate by careless pattern type
compute_metrics_by_type <- function(z_scores, labels_df, z_threshold_values, eval_threshold) {
  pattern <- labels_df$pattern
  truly_careless <- labels_df$careless_pct >= eval_threshold

  lapply(z_threshold_values, function(zt) {
    flagged <- z_scores <= zt
    if (sum(truly_careless) == 0) return(NULL)

    types <- unique(pattern[truly_careless])
    types <- types[types != "clean"]
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
# BUILD THE GRID (per-replication)
# ============================================================

data_grid <- expand.grid(
  nFactors = params$nFactors,
  n_respondents = params$n_respondents,
  pct_careless = params$pct_careless,
  stringsAsFactors = FALSE
)
data_grid$total_items <- sapply(data_grid$nFactors, function(nf) sum(make_nItems(nf)))
data_grid$n_lt_p <- data_grid$n_respondents < data_grid$total_items

# ============================================================
# OUTPUT FILES
# ============================================================

results_file    <- "final_multiverse_results.csv"
by_type_file    <- "final_multiverse_results_by_type.csv"
benchmark_file  <- "final_benchmark_results.csv"
checkpoint_file <- "final_multiverse_checkpoint.rds"

# --- Safe write helpers (OneDrive file-lock workaround) ---
safe_write_table <- function(x, file, ..., max_tries = 5) {
  for (i in seq_len(max_tries)) {
    ok <- tryCatch({ write.table(x, file, ...); TRUE },
                   error = function(e) FALSE)
    if (ok) return(invisible(NULL))
    Sys.sleep(0.5 * i)
  }
  stop("Failed to write to ", file, " after ", max_tries, " retries")
}

safe_saveRDS <- function(object, file, max_tries = 5) {
  tmp <- paste0(file, ".tmp")
  for (i in seq_len(max_tries)) {
    ok <- tryCatch({
      saveRDS(object, tmp)
      if (file.exists(file)) file.remove(file)
      file.rename(tmp, file)
      TRUE
    }, error = function(e) FALSE)
    if (ok) return(invisible(NULL))
    Sys.sleep(0.5 * i)
  }
  stop("Failed to save checkpoint to ", file, " after ", max_tries, " retries")
}

# Check for existing checkpoint to resume from
start_rep <- 1
start_row <- 1
if (file.exists(checkpoint_file)) {
  cp <- readRDS(checkpoint_file)
  start_rep <- cp$next_rep
  start_row <- cp$next_row
  cat(sprintf("Resuming from checkpoint: rep %d, row %d of %d\n\n",
              start_rep, start_row, nrow(data_grid)))
}

# Initialize output files if starting completely fresh
if (start_rep == 1 && start_row == 1) {
  header_main <- paste(
    "rep_id", "nFactors", "n_respondents", "total_items", "n_lt_p",
    "pct_careless", "corProp", "z_threshold",
    "TP", "TN", "FP", "FN", "n_truly_careless",
    "accuracy", "sensitivity", "specificity",
    "precision", "f1", "youden", "mcc",
    "false_positive_rate", "false_negative_rate", "n_flagged",
    "mean_z_good", "mean_z_careless", "auc",
    sep = ","
  )
  writeLines(header_main, results_file)

  header_type <- paste(
    "rep_id", "nFactors", "n_respondents", "pct_careless", "corProp", "z_threshold",
    "pattern", "n_cases", "n_detected", "detection_rate", "mean_z",
    sep = ","
  )
  writeLines(header_type, by_type_file)

  header_bench <- paste(
    "rep_id", "nFactors", "n_respondents", "total_items", "n_lt_p",
    "pct_careless", "method",
    "best_mcc", "best_threshold", "sensitivity", "specificity", "auc",
    sep = ","
  )
  writeLines(header_bench, benchmark_file)
}

# ============================================================
# MAIN LOOP: replications x data conditions x corProp
# ============================================================

total_rr_calls <- n_rr_calls_total
completed_calls <- 0

# Count already-completed calls from prior reps
if (start_rep > 1) {
  completed_calls <- (start_rep - 1) * n_rr_calls_per_rep
}
if (start_row > 1) {
  completed_calls <- completed_calls + (start_row - 1) * n_corProp
}

time_start <- Sys.time()

for (rep_id in start_rep:R_REPLICATIONS) {

  cat(sprintf("\n===== REPLICATION %d / %d =====\n", rep_id, R_REPLICATIONS))

  row_start <- if (rep_id == start_rep) start_row else 1

  for (dg_idx in row_start:nrow(data_grid)) {

    dg <- data_grid[dg_idx, ]

    # --- Progress report ---
    elapsed <- as.numeric(difftime(Sys.time(), time_start, units = "mins"))
    if (completed_calls > 0) {
      rate <- elapsed / completed_calls  # min per RR call
      remaining_calls <- total_rr_calls - completed_calls
      eta_min <- rate * remaining_calls
      if (eta_min > 60) {
        eta_str <- sprintf("%.1f hours remaining", eta_min / 60)
      } else {
        eta_str <- sprintf("%.0f min remaining", eta_min)
      }
    } else {
      eta_str <- "estimating..."
    }

    pct_done <- round(100 * completed_calls / total_rr_calls, 1)
    nItems_vec <- make_nItems(dg$nFactors)

    cat(sprintf("[%5.1f%%] Rep %d | Row %d/%d | nF=%d n=%d pct=%.2f | %.1fmin | %s\n",
                pct_done, rep_id, dg_idx, nrow(data_grid),
                dg$nFactors, dg$n_respondents, dg$pct_careless,
                elapsed, eta_str))

    # --- Reproducible seed per cell x replication ---
    # Ensures same dataset for same (rep, nF, n, pct) across runs
    cell_seed <- SEED_BASE + (rep_id - 1) * nrow(data_grid) + dg_idx
    set.seed(cell_seed)

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
      completed_calls <- completed_calls + n_corProp
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
      completed_calls <- completed_calls + n_corProp
      next
    }

    corrupted_data <- injection$data_corrupted
    labels <- injection$labels

    # --- Benchmark methods (once per dataset, before corProp loop) ---
    is_careless <- labels$careless_pct >= EVAL_THRESHOLD

    bench_rows <- vector("list", length(BENCHMARK_METHODS))
    for (bm_idx in seq_along(BENCHMARK_METHODS)) {
      bm <- BENCHMARK_METHODS[[bm_idx]]
      scores <- tryCatch(
        bm$fn(corrupted_data),
        error = function(e) rep(NA_real_, nrow(corrupted_data))
      )
      best <- find_best_threshold(scores, is_careless, higher_is_careless = bm$higher)
      auc_val <- calc_auc_generic(scores, is_careless, higher_is_careless = bm$higher)
      bench_rows[[bm_idx]] <- data.frame(
        rep_id = rep_id,
        nFactors = dg$nFactors,
        n_respondents = dg$n_respondents,
        total_items = dg$total_items,
        n_lt_p = dg$n_lt_p,
        pct_careless = dg$pct_careless,
        method = bm$name,
        best_mcc = best$mcc,
        best_threshold = best$threshold,
        sensitivity = best$sensitivity,
        specificity = best$specificity,
        auc = auc_val,
        stringsAsFactors = FALSE
      )
    }
    bench_df <- do.call(rbind, bench_rows)
    safe_write_table(bench_df, benchmark_file, append = TRUE, sep = ",",
                row.names = FALSE, col.names = FALSE)

    # --- Run ReReReRe for each corProp ---
    for (cp in params$corProp) {

      rr_result <- tryCatch(
        ReReReRe(
          data = corrupted_data,
          corProp = cp,
          cutOff = 0,  # dummy; z_threshold applied post-hoc
          iterations = ITERATIONS,
          min_pairs = MIN_PAIRS,
          align_signs = TRUE,
          progress = FALSE
        ),
        error = function(e) {
          cat("  WARNING: ReReReRe failed (corProp=", cp, "):", conditionMessage(e), "\n")
          return(NULL)
        }
      )

      if (is.null(rr_result)) {
        completed_calls <- completed_calls + 1
        next
      }

      # Compute AUC (threshold-free, once per RR call)
      rr_auc <- compute_auc(rr_result$z_score, labels, EVAL_THRESHOLD)

      # Apply all z_threshold values post-hoc
      metrics <- compute_metrics(rr_result$z_score, labels, params$z_threshold, EVAL_THRESHOLD)
      metrics$rep_id <- rep_id
      metrics$nFactors <- dg$nFactors
      metrics$n_respondents <- dg$n_respondents
      metrics$total_items <- dg$total_items
      metrics$n_lt_p <- dg$n_lt_p
      metrics$pct_careless <- dg$pct_careless
      metrics$corProp <- cp
      metrics$auc <- rr_auc

      # Reorder columns
      metrics <- metrics[, c("rep_id", "nFactors", "n_respondents",
                              "total_items", "n_lt_p",
                              "pct_careless", "corProp", "z_threshold",
                              "TP", "TN", "FP", "FN", "n_truly_careless",
                              "accuracy", "sensitivity", "specificity",
                              "precision", "f1", "youden", "mcc",
                              "false_positive_rate", "false_negative_rate",
                              "n_flagged",
                              "mean_z_good", "mean_z_careless", "auc")]

      # Append to CSV
      safe_write_table(metrics, results_file, append = TRUE, sep = ",",
                  row.names = FALSE, col.names = FALSE)

      # By-type metrics
      type_metrics <- compute_metrics_by_type(rr_result$z_score, labels,
                                               params$z_threshold, EVAL_THRESHOLD)
      if (!is.null(type_metrics) && nrow(type_metrics) > 0) {
        type_metrics$rep_id <- rep_id
        type_metrics$nFactors <- dg$nFactors
        type_metrics$n_respondents <- dg$n_respondents
        type_metrics$pct_careless <- dg$pct_careless
        type_metrics$corProp <- cp

        type_metrics <- type_metrics[, c("rep_id", "nFactors", "n_respondents",
                                          "pct_careless", "corProp", "z_threshold",
                                          "pattern", "n_cases", "n_detected",
                                          "detection_rate", "mean_z")]

        safe_write_table(type_metrics, by_type_file, append = TRUE, sep = ",",
                    row.names = FALSE, col.names = FALSE)
      }

      completed_calls <- completed_calls + 1
    }

    # --- Checkpoint after each data-generation row ---
    safe_saveRDS(list(next_rep = rep_id, next_row = dg_idx + 1,
                      timestamp = Sys.time(), completed_calls = completed_calls),
                 checkpoint_file)
  }

  # After completing a full replication, reset row start
  # and update checkpoint
  safe_saveRDS(list(next_rep = rep_id + 1, next_row = 1,
                    timestamp = Sys.time(), completed_calls = completed_calls),
               checkpoint_file)

  elapsed_rep <- as.numeric(difftime(Sys.time(), time_start, units = "mins"))
  cat(sprintf("\n  Replication %d complete | Total elapsed: %.1f min (%.1f hours)\n",
              rep_id, elapsed_rep, elapsed_rep / 60))
}

# ============================================================
# DONE
# ============================================================

elapsed_total <- difftime(Sys.time(), time_start, units = "hours")

cat("\n===============================================================\n")
cat("  FINAL MULTIVERSE COMPLETE\n")
cat("===============================================================\n")
cat(sprintf("  Total time: %.1f hours\n", as.numeric(elapsed_total)))
cat(sprintf("  ReReReRe calls completed: %s\n", format(completed_calls, big.mark = ",")))
cat(sprintf("  Results: %s\n", results_file))
cat(sprintf("  By-type: %s\n", by_type_file))
cat(sprintf("  Benchmarks: %s\n", benchmark_file))
cat("===============================================================\n")

# Clean up checkpoint
if (file.exists(checkpoint_file)) file.remove(checkpoint_file)
