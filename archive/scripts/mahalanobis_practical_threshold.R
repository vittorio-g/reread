#### Mahalanobis: Practical (chi-square) vs Oracle Threshold ####
# 
# Regenerates the SAME datasets as Final_Multiverse_Runner.R using identical
# seeds, then evaluates Mahalanobis distance with:
#   1. Oracle threshold (post-hoc best MCC from 200 candidates)
#   2. Practical threshold: chi-square(p, alpha=.001) -- standard in literature
#   3. Practical threshold: chi-square(p, alpha=.01)  -- less conservative
#
# This is fast because we skip all ReReReRe permutation steps.
# Only data generation + Mahalanobis computation + evaluation.
#
# Uses IDENTICAL seed logic as Final_Multiverse_Runner.R:
#   cell_seed = SEED_BASE + (rep_id - 1) * nrow(data_grid) + dg_idx

rm(list = ls())

library(dplyr)
library(lavaan)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("Synthetic_Good_Responses_2.R")
source("Careless_machine_2.R")

# ============================================================
# PARAMETERS (must match Final_Multiverse_Runner.R exactly)
# ============================================================

params <- list(
  nFactors      = c(4, 6, 8, 10, 12, 15, 20, 25),
  n_respondents = c(50, 100, 200, 300, 500, 1000),
  pct_careless  = c(0.05, 0.10, 0.25, 0.50)
)

R_REPLICATIONS  <- 50
CARELESS_TYPES  <- c(random = 1/3, longstring = 1/3, mixed = 1/3)
CARELESS_LEVELS <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5)
ITEMS_CYCLE     <- c(10, 6, 3)
EVAL_THRESHOLD  <- 0.01
SEED_BASE       <- 42

make_nItems <- function(nFactors) rep_len(ITEMS_CYCLE, nFactors)

# ============================================================
# MAHALANOBIS FUNCTIONS (same as in runner)
# ============================================================

compute_mahad <- function(data) {
  x <- as.matrix(data)
  center <- colMeans(x)
  cv <- cov(x)
  tryCatch(
    mahalanobis(x, center, cv),
    error = function(e) {
      tryCatch(
        mahalanobis(x, center, cv + diag(1e-3, ncol(x))),
        error = function(e2) rep(NA_real_, nrow(x))
      )
    }
  )
}

calc_mcc <- function(tp, tn, fp, fn) {
  denom <- sqrt(as.numeric(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  if (denom == 0) return(0)
  (tp * tn - fp * fn) / denom
}

calc_auc_generic <- function(scores, is_careless, higher_is_careless = TRUE) {
  if (any(is.na(scores))) return(NA_real_)
  if (!higher_is_careless) scores <- -scores
  n_pos <- sum(is_careless)
  n_neg <- sum(!is_careless)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  ranks <- rank(scores, ties.method = "average")
  (sum(ranks[is_careless]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

find_best_threshold <- function(scores, is_careless, higher_is_careless = TRUE,
                                n_candidates = 200) {
  if (any(is.na(scores))) {
    return(list(mcc = NA_real_, threshold = NA_real_,
                sensitivity = NA_real_, specificity = NA_real_))
  }
  candidates <- quantile(scores, probs = seq(0.005, 0.995, length.out = n_candidates),
                         na.rm = TRUE)
  candidates <- unique(candidates)
  best_mcc <- -Inf; best_thresh <- NA; best_sens <- NA; best_spec <- NA
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

# Evaluate at a FIXED threshold
eval_fixed_threshold <- function(scores, is_careless, threshold, higher_is_careless = TRUE) {
  if (any(is.na(scores))) {
    return(list(mcc = NA_real_, sensitivity = NA_real_, specificity = NA_real_))
  }
  flagged <- if (higher_is_careless) scores >= threshold else scores <= threshold
  tp <- sum(flagged & is_careless); tn <- sum(!flagged & !is_careless)
  fp <- sum(flagged & !is_careless); fn <- sum(!flagged & is_careless)
  mcc <- calc_mcc(tp, tn, fp, fn)
  n_pos <- tp + fn
  sens <- if (n_pos > 0) tp / n_pos else NA_real_
  spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  list(mcc = mcc, sensitivity = sens, specificity = spec)
}

# ============================================================
# BUILD DATA GRID (same as runner)
# ============================================================

data_grid <- expand.grid(
  nFactors = params$nFactors,
  n_respondents = params$n_respondents,
  pct_careless = params$pct_careless,
  stringsAsFactors = FALSE
)
data_grid$total_items <- sapply(data_grid$nFactors, function(nf) sum(make_nItems(nf)))
data_grid$n_lt_p <- data_grid$n_respondents < data_grid$total_items

cat("===============================================================\n")
cat("  MAHALANOBIS: PRACTICAL vs ORACLE THRESHOLD\n")
cat("===============================================================\n")
cat(sprintf("  Data conditions: %d\n", nrow(data_grid)))
cat(sprintf("  Replications:    %d\n", R_REPLICATIONS))
cat(sprintf("  Total datasets:  %d\n", nrow(data_grid) * R_REPLICATIONS))
cat("  Thresholds: oracle, chi-sq(.001), chi-sq(.01)\n")
cat("===============================================================\n\n")

# ============================================================
# OUTPUT FILE
# ============================================================

out_file <- "mahalanobis_practical_results.csv"

header <- paste(
  "rep_id", "nFactors", "n_respondents", "total_items", "pct_careless",
  "auc",
  "oracle_mcc", "oracle_threshold", "oracle_sens", "oracle_spec",
  "chisq001_mcc", "chisq001_threshold", "chisq001_sens", "chisq001_spec",
  "chisq01_mcc", "chisq01_threshold", "chisq01_sens", "chisq01_spec",
  sep = ","
)
writeLines(header, out_file)

# ============================================================
# MAIN LOOP
# ============================================================

total_datasets <- nrow(data_grid) * R_REPLICATIONS
done <- 0
time_start <- Sys.time()

for (rep_id in 1:R_REPLICATIONS) {
  
  cat(sprintf("===== REPLICATION %d / %d =====\n", rep_id, R_REPLICATIONS))
  
  for (dg_idx in 1:nrow(data_grid)) {
    
    dg <- data_grid[dg_idx, ]
    
    # Same seed as Final_Multiverse_Runner.R
    cell_seed <- SEED_BASE + (rep_id - 1) * nrow(data_grid) + dg_idx
    set.seed(cell_seed)
    
    # Generate clean dataset
    nItems_vec <- make_nItems(dg$nFactors)
    clean_data <- tryCatch(
      simulated_good_responses(
        nConstructs = dg$nFactors,
        nItems = nItems_vec,
        n = dg$n_respondents
      ),
      error = function(e) NULL
    )
    if (is.null(clean_data)) { done <- done + 1; next }
    
    # Inject careless
    injection <- tryCatch(
      inject_careless(
        data = clean_data,
        pct_careless = dg$pct_careless,
        pct_types = CARELESS_TYPES,
        careless_levels = CARELESS_LEVELS
      ),
      error = function(e) NULL
    )
    if (is.null(injection)) { done <- done + 1; next }
    
    corrupted_data <- injection$data_corrupted
    labels <- injection$labels
    is_careless <- labels$careless_pct >= EVAL_THRESHOLD
    
    # Compute Mahalanobis D^2
    d2 <- compute_mahad(corrupted_data)
    
    # AUC (threshold-free)
    auc_val <- calc_auc_generic(d2, is_careless, higher_is_careless = TRUE)
    
    # Oracle threshold
    oracle <- find_best_threshold(d2, is_careless, higher_is_careless = TRUE)
    
    # Practical thresholds (chi-square)
    p <- dg$total_items
    chi_001 <- qchisq(0.999, df = p)  # alpha = .001
    chi_01  <- qchisq(0.99,  df = p)  # alpha = .01
    
    prac_001 <- eval_fixed_threshold(d2, is_careless, chi_001, higher_is_careless = TRUE)
    prac_01  <- eval_fixed_threshold(d2, is_careless, chi_01,  higher_is_careless = TRUE)
    
    # Write row
    row <- sprintf("%d,%d,%d,%d,%.2f,%.6f,%.6f,%.4f,%.6f,%.6f,%.6f,%.4f,%.6f,%.6f,%.6f,%.4f,%.6f,%.6f",
      rep_id, dg$nFactors, dg$n_respondents, p, dg$pct_careless,
      auc_val,
      oracle$mcc, oracle$threshold, oracle$sensitivity, oracle$specificity,
      prac_001$mcc, chi_001, prac_001$sensitivity, prac_001$specificity,
      prac_01$mcc, chi_01, prac_01$sensitivity, prac_01$specificity
    )
    cat(row, "\n", file = out_file, append = TRUE)
    
    done <- done + 1
    
    # Progress every 100 datasets
    if (done %% 100 == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), time_start, units = "mins"))
      rate <- elapsed / done
      eta <- rate * (total_datasets - done)
      cat(sprintf("  [%5.1f%%] %d/%d datasets | %.1f min elapsed | ETA %.1f min\n",
                  100 * done / total_datasets, done, total_datasets, elapsed, eta))
    }
  }
}

elapsed_total <- as.numeric(difftime(Sys.time(), time_start, units = "mins"))
cat(sprintf("\nDone! %d datasets in %.1f minutes.\n", done, elapsed_total))
cat(sprintf("Output: %s\n", out_file))

# ============================================================
# QUICK SUMMARY
# ============================================================

cat("\n\n===============================================================\n")
cat("  QUICK SUMMARY\n")
cat("===============================================================\n\n")

results <- read.csv(out_file)

cat("--- Overall ---\n")
cat(sprintf("  Oracle MCC:     %.3f (mean), %.3f (median)\n",
            mean(results$oracle_mcc, na.rm=T), median(results$oracle_mcc, na.rm=T)))
cat(sprintf("  Chi-sq .001 MCC: %.3f (mean), %.3f (median)\n",
            mean(results$chisq001_mcc, na.rm=T), median(results$chisq001_mcc, na.rm=T)))
cat(sprintf("  Chi-sq .01 MCC:  %.3f (mean), %.3f (median)\n",
            mean(results$chisq01_mcc, na.rm=T), median(results$chisq01_mcc, na.rm=T)))
cat(sprintf("  AUC:            %.3f (mean)\n\n", mean(results$auc, na.rm=T)))

cat("--- By nFactors ---\n")
cat(sprintf("%4s  %12s  %12s  %12s  %10s\n", "nF", "Oracle MCC", "ChiSq.001", "ChiSq.01", "AUC"))
cat(strrep("-", 56), "\n")
for (nf in sort(unique(results$nFactors))) {
  sub <- results[results$nFactors == nf, ]
  cat(sprintf("%4d  %12.3f  %12.3f  %12.3f  %10.3f\n",
              nf,
              mean(sub$oracle_mcc, na.rm=T),
              mean(sub$chisq001_mcc, na.rm=T),
              mean(sub$chisq01_mcc, na.rm=T),
              mean(sub$auc, na.rm=T)))
}

cat("\n--- Sensitivity/Specificity at chi-sq .001 ---\n")
cat(sprintf("%4s  %8s  %8s  %12s\n", "nF", "Sens", "Spec", "MCC"))
cat(strrep("-", 38), "\n")
for (nf in sort(unique(results$nFactors))) {
  sub <- results[results$nFactors == nf, ]
  cat(sprintf("%4d  %8.3f  %8.3f  %12.3f\n",
              nf,
              mean(sub$chisq001_sens, na.rm=T),
              mean(sub$chisq001_spec, na.rm=T),
              mean(sub$chisq001_mcc, na.rm=T)))
}

cat("\n--- Oracle MCC loss from using practical threshold ---\n")
cat(sprintf("%4s  %14s  %14s\n", "nF", "Loss (.001)", "Loss (.01)"))
cat(strrep("-", 36), "\n")
for (nf in sort(unique(results$nFactors))) {
  sub <- results[results$nFactors == nf, ]
  oracle_m <- mean(sub$oracle_mcc, na.rm=T)
  chi001_m <- mean(sub$chisq001_mcc, na.rm=T)
  chi01_m  <- mean(sub$chisq01_mcc, na.rm=T)
  cat(sprintf("%4d  %+14.3f  %+14.3f\n", nf, chi001_m - oracle_m, chi01_m - oracle_m))
}

cat("\nDone.\n")
