#### Partial Simulation — weighted_v2 (scale-invariant) on ≤60-item conditions ####
#
# Complements Final_Simulation_Coupled_Partial.R by testing weighted_v2 on the
# same short-questionnaire subset. Same SEED_BASE, same data_grid, same seeding
# formula → byte-identical synthetic data per (rep_id, dg_idx) cell. Results can
# be joined directly for head-to-head comparison.
#
# weighted_v2 uses the scale-invariant within-respondent weighted Pearson
# correlation (fix for the 2026-04-21 Schneider inversion bug). All pairs
# included, weighted by sample |r|. No corProp loop (corProp doesn't apply).
#
# OUTPUTS (separate from the coupled partial run):
#   partial_v2_sim_rr_results.csv
#   partial_v2_sim_rr_by_type.csv
#   partial_v2_sim_checkpoint.rds
#
# EXPECTED RUNTIME: 18,000 calls (360 conditions × 50 reps × 1 scoring).
# Weighted_v2 is slower per call than coupled (uses all pairs), but no corProp
# loop. Estimate: ~6-10 hours.
# ============================================================

library(dplyr)
library(magrittr)
library(lavaan)
library(psych)

tryCatch(
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path)),
  error = function(e) message("Working directory not changed (not in RStudio?)")
)

source("Synthetic_Good_Responses_2.R")
source("Careless_machine_2.R")
# NOTE: we do NOT source ReReReRe.R — weighted_v2 is inlined here to avoid
# any coupling with the live file (which is currently coupled-only).

# ============================================================
# weighted_v2 SCORING (inline, self-contained)
# ============================================================

rowCor_weighted_v2 <- function(A, B, weights, zero_val = 0) {
  sum_w <- sum(weights)
  W <- matrix(weights, nrow(A), ncol(A), byrow = TRUE)
  mean_A <- rowSums(A * W, na.rm = TRUE) / sum_w
  mean_B <- rowSums(B * W, na.rm = TRUE) / sum_w
  A_c <- A - mean_A; B_c <- B - mean_B
  num <- rowSums(W * A_c * B_c, na.rm = TRUE)
  den <- sqrt(rowSums(W * A_c^2, na.rm = TRUE) * rowSums(W * B_c^2, na.rm = TRUE))
  r <- abs(num / den)
  r[is.nan(r) | is.na(r)] <- zero_val
  r
}

score_weighted_v2 <- function(data, iterations = 100, align_signs = TRUE) {
  data <- data[, sapply(data, is.numeric), drop = FALSE]
  N <- nrow(data); J <- ncol(data)
  mat <- as.matrix(data)

  if (align_signs) item_max <- apply(mat, 2, max, na.rm = TRUE)

  rawCorMat <- cor(mat, use = "pairwise.complete.obs")
  corMat <- abs(rawCorMat)
  corMat[upper.tri(corMat, diag = TRUE)] <- NA
  rawCorMat[upper.tri(rawCorMat, diag = TRUE)] <- NA
  if (align_signs) signMat <- sign(rawCorMat)

  pair_idx <- which(!is.na(corMat), arr.ind = TRUE)
  weights <- corMat[pair_idx]
  pair_signs <- sign(rawCorMat[pair_idx])
  k <- nrow(pair_idx)

  A_sel <- mat[, pair_idx[, 1], drop = FALSE]
  B_sel <- mat[, pair_idx[, 2], drop = FALSE]

  if (align_signs) {
    nf <- which(pair_signs < 0)
    if (length(nf) > 0) {
      B_sel[, nf] <- rep(item_max[pair_idx[nf, 2]] + 1, each = N) - B_sel[, nf]
    }
  }

  rowCors <- rowCor_weighted_v2(A_sel, B_sel, weights, zero_val = 0)

  all_RIC <- matrix(NA_real_, nrow = N, ncol = iterations)
  for (i in seq_len(iterations)) {
    idx1 <- sample(J, k, replace = TRUE)
    idx2 <- sample(J, k, replace = TRUE)
    same <- idx1 == idx2
    while (any(same)) {
      idx2[same] <- sample(J, sum(same), replace = TRUE)
      same <- idx1 == idx2
    }
    A_r <- mat[, idx1, drop = FALSE]; B_r <- mat[, idx2, drop = FALSE]
    if (align_signs) {
      ri <- pmax(idx1, idx2); ci <- pmin(idx1, idx2)
      rs <- signMat[cbind(ri, ci)]; rs[is.na(rs)] <- 1
      rn <- which(rs < 0)
      if (length(rn) > 0) {
        B_r[, rn] <- rep(item_max[idx2[rn]] + 1, each = N) - B_r[, rn]
      }
    }
    all_RIC[, i] <- rowCor_weighted_v2(A_r, B_r, weights, zero_val = 0)
  }

  rm <- rowMeans(all_RIC, na.rm = TRUE)
  rs <- apply(all_RIC, 1, sd, na.rm = TRUE)
  z  <- ifelse(rs > 0, (rowCors - rm) / rs, 0)

  data.frame(z_score = z, indCors = rowCors, rand_mean = rm, rand_sd = rs,
             n_pairs = k)
}

# ============================================================
# DESIGN — MUST match Final_Simulation.R for seed congruence
# ============================================================

params <- list(
  nFactors       = c(3, 5, 8, 12, 16, 20),
  items_per_factor = c(3, 5, 7, 10),
  n_respondents  = c(50, 100, 300, 500),
  pct_careless   = c(0.03, 0.05, 0.10, 0.15, 0.25, 0.40),
  z_threshold    = c(0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7,
                     3.0, 3.3, 3.6, 3.9, 4.2, 5.0)
)

R_REPLICATIONS  <- 50
ITERATIONS      <- 100
EVAL_THRESHOLD  <- 0.01
SEED_BASE       <- 2026   # MUST match Final_Simulation.R
MAX_ITEMS       <- 60     # same short-questionnaire threshold

CARELESS_TYPES  <- c(random = 1/3, longstring = 1/3, mixed = 1/3)
CARELESS_LEVELS <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5)

data_grid <- expand.grid(
  nFactors = params$nFactors,
  ipf = params$items_per_factor,
  n_respondents = params$n_respondents,
  pct_careless = params$pct_careless,
  stringsAsFactors = FALSE
)
data_grid$total_items <- data_grid$nFactors * data_grid$ipf
data_grid$n_lt_p <- data_grid$n_respondents < data_grid$total_items

eligible_mask <- data_grid$total_items <= MAX_ITEMS
n_eligible <- sum(eligible_mask)
n_zt <- length(params$z_threshold)
n_rr_calls_per_rep <- n_eligible
n_rr_calls_total <- n_rr_calls_per_rep * R_REPLICATIONS

cat("===============================================================\n")
cat("  PARTIAL SIMULATION — weighted_v2 (<=", MAX_ITEMS, " items)\n", sep = "")
cat("===============================================================\n")
cat(sprintf("  Full grid:                 %d conditions\n", nrow(data_grid)))
cat(sprintf("  Eligible (total_items<=%d): %d conditions\n", MAX_ITEMS, n_eligible))
cat(sprintf("  Replications:              %d\n", R_REPLICATIONS))
cat(sprintf("  Scoring calls per rep:     %s (no corProp loop)\n",
            format(n_rr_calls_per_rep, big.mark = ",")))
cat(sprintf("  Total scoring calls:       %s\n", format(n_rr_calls_total, big.mark = ",")))
cat(sprintf("  Total result rows:         %s  (x %d z_threshold)\n",
            format(n_rr_calls_total * n_zt, big.mark = ","), n_zt))
cat(sprintf("  Estimated time:            ~6-10 hours\n"))
cat("===============================================================\n\n")

# ============================================================
# HELPERS
# ============================================================

compute_metrics <- function(z_scores, labels_df, z_threshold_values, eval_threshold) {
  actual <- labels_df$careless_pct >= eval_threshold
  lapply(z_threshold_values, function(zt) {
    flagged <- z_scores <= zt
    TP <- sum(flagged & actual);  TN <- sum(!flagged & !actual)
    FP <- sum(flagged & !actual); FN <- sum(!flagged & actual)
    sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
    specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
    precision   <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
    f1 <- if (!is.na(precision) && !is.na(sensitivity) && (precision + sensitivity) > 0)
      2 * precision * sensitivity / (precision + sensitivity) else NA_real_
    mcc_denom <- sqrt(as.numeric(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    mcc <- if (mcc_denom > 0) (TP * TN - FP * FN) / mcc_denom else NA_real_
    data.frame(
      z_threshold = zt,
      TP = TP, TN = TN, FP = FP, FN = FN,
      n_truly_careless = sum(actual),
      sensitivity = sensitivity, specificity = specificity,
      precision = precision, f1 = f1, mcc = mcc,
      n_flagged = sum(flagged),
      mean_z_good = mean(z_scores[!actual], na.rm = TRUE),
      mean_z_careless = mean(z_scores[actual], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }) %>% do.call(rbind, .)
}

compute_auc <- function(z_scores, labels_df, eval_threshold) {
  actual <- labels_df$careless_pct >= eval_threshold
  if (any(is.na(z_scores))) return(NA_real_)
  n_pos <- sum(actual); n_neg <- sum(!actual)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  scores <- -z_scores
  ranks <- rank(scores, ties.method = "average")
  (sum(ranks[actual]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

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
        z_threshold = zt, pattern = typ,
        n_cases = sum(mask),
        n_detected = sum(flagged[mask]),
        detection_rate = mean(flagged[mask]),
        mean_z = mean(z_scores[mask], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}

safe_write_table <- function(x, file, ..., max_tries = 5) {
  for (i in seq_len(max_tries)) {
    ok <- tryCatch({ write.table(x, file, ...); TRUE }, error = function(e) FALSE)
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

# ============================================================
# OUTPUT FILES
# ============================================================

results_file    <- "partial_v2_sim_rr_results.csv"
by_type_file    <- "partial_v2_sim_rr_by_type.csv"
checkpoint_file <- "partial_v2_sim_checkpoint.rds"

start_rep <- 1
start_row <- 1
if (file.exists(checkpoint_file)) {
  cp <- readRDS(checkpoint_file)
  start_rep <- cp$next_rep
  start_row <- cp$next_row
  cat(sprintf("Resuming: rep %d, row %d/%d\n\n", start_rep, start_row, nrow(data_grid)))
}

if (start_rep == 1 && start_row == 1) {
  # Header — same schema as Final_Simulation, with "mode_used"="weighted_v2" and
  # corProp=NA (not applicable for weighted_v2)
  header_main <- paste(
    "rep_id", "nFactors", "ipf", "total_items", "n_respondents", "n_lt_p",
    "pct_careless", "corProp", "mode_used", "z_threshold",
    "TP", "TN", "FP", "FN", "n_truly_careless",
    "sensitivity", "specificity", "precision", "f1", "mcc",
    "n_flagged", "mean_z_good", "mean_z_careless", "auc",
    sep = ","
  )
  writeLines(header_main, results_file)

  header_type <- paste(
    "rep_id", "nFactors", "ipf", "total_items", "n_respondents",
    "pct_careless", "corProp", "z_threshold",
    "pattern", "n_cases", "n_detected", "detection_rate", "mean_z",
    sep = ","
  )
  writeLines(header_type, by_type_file)
}

# ============================================================
# MAIN LOOP
# ============================================================

advance_checkpoint <- function(rep_id, dg_idx, n_rows) {
  if (dg_idx < n_rows) list(next_rep = rep_id, next_row = dg_idx + 1)
  else list(next_rep = rep_id + 1, next_row = 1)
}

completed_calls <- 0
if (start_rep > 1) completed_calls <- (start_rep - 1) * n_rr_calls_per_rep

time_start <- Sys.time()

for (rep_id in start_rep:R_REPLICATIONS) {

  cat(sprintf("\n===== REPLICATION %d / %d =====\n", rep_id, R_REPLICATIONS))

  row_start <- if (rep_id == start_rep) start_row else 1

  for (dg_idx in row_start:nrow(data_grid)) {

    dg <- data_grid[dg_idx, ]

    # Skip conditions > MAX_ITEMS
    if (dg$total_items > MAX_ITEMS) {
      safe_saveRDS(advance_checkpoint(rep_id, dg_idx, nrow(data_grid)), checkpoint_file)
      next
    }

    elapsed <- as.numeric(difftime(Sys.time(), time_start, units = "mins"))
    if (completed_calls > 0) {
      rate <- elapsed / completed_calls
      eta_min <- rate * (n_rr_calls_total - completed_calls)
      eta_str <- if (eta_min > 60) sprintf("%.1f h remaining", eta_min / 60)
                 else sprintf("%.0f min remaining", eta_min)
    } else {
      eta_str <- "estimating..."
    }

    pct_done <- round(100 * completed_calls / n_rr_calls_total, 1)
    cat(sprintf("[%5.1f%%] Rep %d | Row %d/%d | nF=%d ipf=%d items=%d n=%d pct=%.2f | %.1fmin | %s\n",
                pct_done, rep_id, dg_idx, nrow(data_grid),
                dg$nFactors, dg$ipf, dg$total_items, dg$n_respondents,
                dg$pct_careless, elapsed, eta_str))

    # --- Seeding MUST match Final_Simulation.R formula ---
    cell_seed <- SEED_BASE + (rep_id - 1) * nrow(data_grid) + dg_idx
    set.seed(cell_seed)

    nItems_vec <- rep(dg$ipf, dg$nFactors)

    clean_data <- tryCatch(
      simulated_good_responses(nConstructs = dg$nFactors, nItems = nItems_vec,
                               n = dg$n_respondents),
      error = function(e) NULL
    )
    if (is.null(clean_data)) {
      completed_calls <- completed_calls + 1
      safe_saveRDS(advance_checkpoint(rep_id, dg_idx, nrow(data_grid)), checkpoint_file)
      next
    }

    injection <- tryCatch(
      inject_careless(data = clean_data, pct_careless = dg$pct_careless,
                      pct_types = CARELESS_TYPES, careless_levels = CARELESS_LEVELS),
      error = function(e) NULL
    )
    if (is.null(injection)) {
      completed_calls <- completed_calls + 1
      safe_saveRDS(advance_checkpoint(rep_id, dg_idx, nrow(data_grid)), checkpoint_file)
      next
    }

    corrupted_data <- injection$data_corrupted
    labels <- injection$labels

    # --- SCORE with weighted_v2 (no corProp loop) ---
    rr_result <- tryCatch(
      score_weighted_v2(corrupted_data, iterations = ITERATIONS, align_signs = TRUE),
      error = function(e) {
        cat("  WARNING: score_weighted_v2 failed:", conditionMessage(e), "\n")
        NULL
      }
    )

    if (is.null(rr_result)) {
      completed_calls <- completed_calls + 1
      safe_saveRDS(advance_checkpoint(rep_id, dg_idx, nrow(data_grid)), checkpoint_file)
      next
    }

    rr_auc <- compute_auc(rr_result$z_score, labels, EVAL_THRESHOLD)
    metrics <- compute_metrics(rr_result$z_score, labels, params$z_threshold, EVAL_THRESHOLD)
    metrics$rep_id <- rep_id
    metrics$nFactors <- dg$nFactors
    metrics$ipf <- dg$ipf
    metrics$total_items <- dg$total_items
    metrics$n_respondents <- dg$n_respondents
    metrics$n_lt_p <- dg$n_lt_p
    metrics$pct_careless <- dg$pct_careless
    metrics$corProp <- NA   # not applicable for weighted_v2
    metrics$mode_used <- "weighted_v2"
    metrics$auc <- rr_auc

    metrics <- metrics[, c("rep_id", "nFactors", "ipf", "total_items",
                           "n_respondents", "n_lt_p",
                           "pct_careless", "corProp", "mode_used", "z_threshold",
                           "TP", "TN", "FP", "FN", "n_truly_careless",
                           "sensitivity", "specificity", "precision", "f1", "mcc",
                           "n_flagged", "mean_z_good", "mean_z_careless", "auc")]

    safe_write_table(metrics, results_file, append = TRUE, sep = ",",
                     row.names = FALSE, col.names = FALSE)

    type_metrics <- compute_metrics_by_type(rr_result$z_score, labels,
                                            params$z_threshold, EVAL_THRESHOLD)
    if (!is.null(type_metrics) && nrow(type_metrics) > 0) {
      type_metrics$rep_id <- rep_id
      type_metrics$nFactors <- dg$nFactors
      type_metrics$ipf <- dg$ipf
      type_metrics$total_items <- dg$total_items
      type_metrics$n_respondents <- dg$n_respondents
      type_metrics$pct_careless <- dg$pct_careless
      type_metrics$corProp <- NA

      type_metrics <- type_metrics[, c("rep_id", "nFactors", "ipf", "total_items",
                                       "n_respondents", "pct_careless", "corProp",
                                       "z_threshold", "pattern", "n_cases",
                                       "n_detected", "detection_rate", "mean_z")]

      safe_write_table(type_metrics, by_type_file, append = TRUE, sep = ",",
                       row.names = FALSE, col.names = FALSE)
    }

    completed_calls <- completed_calls + 1
    safe_saveRDS(advance_checkpoint(rep_id, dg_idx, nrow(data_grid)), checkpoint_file)
  }

  cat(sprintf("\n--- Replication %d complete | Elapsed: %.1f min ---\n",
              rep_id, as.numeric(difftime(Sys.time(), time_start, units = "mins"))))
}

total_time <- as.numeric(difftime(Sys.time(), time_start, units = "mins"))
cat("\n===============================================================\n")
cat("  WEIGHTED_V2 PARTIAL SIMULATION COMPLETE\n")
cat("===============================================================\n")
cat(sprintf("  Total time:     %.1f min (%.1f h)\n", total_time, total_time / 60))
cat(sprintf("  Scoring calls:  %s\n", format(completed_calls, big.mark = ",")))
cat(sprintf("  Results file:   %s\n", results_file))
cat(sprintf("  By-type file:   %s\n", by_type_file))
cat("===============================================================\n")

if (file.exists(checkpoint_file)) file.remove(checkpoint_file)

cat("\nNext step: compare partial_v2_sim_rr_results.csv against the coupled\n")
cat("partial run. Use Compare_Coupled_vs_WeightedV2.R (to be written).\n")
