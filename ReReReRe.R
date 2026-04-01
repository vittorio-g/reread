#### ReReReRe: Permutation-based careless respondent detection ####
#
# Two public functions:
#   ReReReRe()   — standard 2-level method (default, recommended)
#   ReReReRe_F() — per-factor variant using EFA (experimental)
#
# Returns per-respondent:
#   result    - legacy percentile: P(random < coupled)
#   z_score   - (coupled_cor - mean_random) / sd_random
#   indCors   - raw coupled absolute correlation
#   rand_mean - mean of random iteration correlations
#   rand_sd   - sd of random iteration correlations
#   flagged   - binary flag (z_score <= z_threshold)
#   z_threshold_used - the z_threshold used for flagging (useful when auto_z=TRUE)
#   nFactors_detected - number of factors detected by parallel analysis
#   mode_used - "coupled", "weighted", or "efa_d"
#   n_pairs   - number of item pairs used
#
# === ARCHITECTURE (2026-04-01) ===
#
# ReReReRe() — STANDARD (default)
#   mode="auto": weighted (all pairs, |r|-weighted) for ≤60 items,
#                coupled (top-k% pairs by |r|) for >60 items.
#   Validated on 6 real datasets + Johnson IPIP-300 inject-and-detect.
#   Coupled AUC=0.635 on real data (wins 5/6 datasets vs EFA-D's 0.561).
#
# ReReReRe_F() — PER-FACTOR (experimental)
#   Uses EFA to identify within-factor pairs, weights by observed |r|.
#   Wins on simulated data (MCC=0.333 vs Std=0.305 overall).
#   Wins on short questionnaires (<60 items) on real data (Pennycook).
#   Loses on real validation datasets with ground truth labels.
#   Recommended for: secondary analysis, short questionnaires, exploratory use.
#
# === REVISION HISTORY ===
# 2026-03-10b: Added z_score output (percentile has ceiling effect)
# 2026-03-10c: Fixed zero_val for random pairs (1 → 0)
# 2026-03-25:  Added align_signs for reverse-coded items
# 2026-03-26b: Changed align_signs from sign flip to reverse coding
# 2026-03-26c: Added auto_z calibration (total_items based)
# 2026-03-28:  Updated defaults (corProp 0.05→0.03, confirmed z=1.5)
# 2026-03-30:  Added weighted mode for short questionnaires
# 2026-03-31:  Tested EFA-D as default — rejected after real data validation
# 2026-04-01:  Split into ReReReRe() (standard) + ReReReRe_F() (per-factor)
# ===

library(dplyr)
library(magrittr)
library(lavaan)
library(psych)
library(mice)

# --- Calibration lookup: total_items → optimal z_threshold ---
# From 2D calibration simulation (R=30, 10 nF × 6 ipf = 60 cells, n=300, corProp=0.05)
# total_items is the best single predictor: R²=0.476 (vs nF+ipf R²=0.245)
#
# 2D lookup table: nF × ipf → mean optimal z
.cal_2d <- data.frame(
  nF  = c(4,4,4,4,4,4, 6,6,6,6,6,6, 8,8,8,8,8,8, 10,10,10,10,10,10,
          12,12,12,12,12,12, 15,15,15,15,15,15, 18,18,18,18,18,18,
          20,20,20,20,20,20, 25,25,25,25,25,25, 30,30,30,30,30,30),
  ipf = rep(c(3,4,6,8,10,12), 10),
  total = c(12,16,24,32,40,48, 18,24,36,48,60,72, 24,32,48,64,80,96,
            30,40,60,80,100,120, 36,48,72,96,120,144, 45,60,90,120,150,180,
            54,72,108,144,180,216, 60,80,120,160,200,240, 75,100,150,200,250,300,
            90,120,180,240,300,360),
  z   = c(1.05,1.22,0.84,1.01,0.97,0.83, 1.11,0.96,1.05,1.01,0.74,0.80,
          0.98,0.93,0.95,0.70,0.69,0.68, 1.12,0.88,0.66,0.66,0.70,0.67,
          1.04,0.96,0.63,0.76,0.60,0.69, 0.65,0.72,0.53,0.59,0.65,1.09,
          0.56,0.65,0.57,0.78,0.90,1.42, 0.72,0.50,0.50,0.78,1.23,1.77,
          0.73,0.43,0.68,1.31,1.94,2.30, 0.56,0.50,1.03,1.87,2.52,2.78)
)

# Fit LOESS on total_items (best single predictor)
.loess_total <- loess(z ~ total, data = .cal_2d, span = 0.4)

# Simple linear fallback: z = 0.505 + 0.0042 * total_items (R²=0.476)
.linear_coef <- c(intercept = 0.505, slope = 0.0042)

#' Look up calibrated z_threshold from total number of items
#' @param total_items Total items in the questionnaire (integer)
#' @param nf Number of factors (optional, for 2D lookup if available)
#' @return Optimal z_threshold (numeric), clamped to [0.3, 3.5]
get_calibrated_z <- function(total_items, nf = NULL) {
  # Clamp to calibration range
  ti_clamped <- max(min(total_items, 400), 10)

  # Use LOESS prediction
  z_pred <- tryCatch(
    as.numeric(predict(.loess_total, newdata = data.frame(total = ti_clamped))),
    error = function(e) {
      # Fallback to linear model
      .linear_coef["intercept"] + .linear_coef["slope"] * ti_clamped
    }
  )

  # Handle NA from LOESS extrapolation
  if (is.na(z_pred)) {
    z_pred <- .linear_coef["intercept"] + .linear_coef["slope"] * ti_clamped
  }

  # Clamp to reasonable range
  max(0.3, min(3.5, z_pred))
}

# --- Helper: vectorized row-wise absolute correlation ---
# Given two N x k matrices, computes |cor(A[i,], B[i,])| for each row i.
# zero_val: value to return when variance is zero (e.g., longstring respondents)
rowCor_abs <- function(A, B, zero_val = 0) {
  A_mean <- rowMeans(A, na.rm = TRUE)
  B_mean <- rowMeans(B, na.rm = TRUE)
  A_c <- A - A_mean
  B_c <- B - B_mean
  num <- rowSums(A_c * B_c, na.rm = TRUE)
  den <- sqrt(rowSums(A_c^2, na.rm = TRUE) * rowSums(B_c^2, na.rm = TRUE))
  r <- abs(num / den)
  r[is.nan(r) | is.na(r)] <- zero_val
  r
}

# --- Helper: weighted coherence score across all pairs ---
# Instead of computing a single correlation across k pairs (rowCor_abs),
# this computes a weighted average of pair-level standardized agreement.
# Each pair contributes proportionally to its sample-level |r| weight.
#
# A, B: N x k matrices (item values for each pair, after sign alignment)
# weights: vector of length k (|r_sample| for each pair)
# Returns: vector of length N (weighted coherence score per person)
rowCor_weighted <- function(A, B, weights) {
  N <- nrow(A)
  k <- ncol(A)

  # Standardize each column (across respondents)
  A_z <- scale(A, center = TRUE, scale = TRUE)
  B_z <- scale(B, center = TRUE, scale = TRUE)

  # Replace NAs from zero-variance columns
  A_z[is.na(A_z)] <- 0
  B_z[is.na(B_z)] <- 0

  # Cross-products: how well does each person follow each pair's pattern?
  cross_products <- A_z * B_z  # N x k

  # Weight each pair by |r_sample|
  w_matrix <- matrix(weights, nrow = N, ncol = k, byrow = TRUE)
  weighted_sum <- rowSums(cross_products * w_matrix, na.rm = TRUE)
  total_weight <- sum(weights)

  return(weighted_sum / total_weight)
}

# --- Helper: EFA-based within-factor pair selection ---
# Runs parallel analysis + EFA, assigns items to primary factors,
# generates all within-factor pairs, returns pair indices + observed |r| weights.
#
# Returns list(idx_A, idx_B, weights, pair_sign, nF_detected, k) or NULL on failure.
.efa_pairs <- function(mat, align_signs = TRUE, progress = FALSE) {
  N <- nrow(mat)
  p <- ncol(mat)

  # Step 1: Estimate nF via parallel analysis
  pa <- tryCatch({
    suppressMessages(suppressWarnings(
      fa.parallel(mat, fa = "fa", plot = FALSE, n.iter = 20)
    ))
  }, error = function(e) NULL)

  nF_est <- NULL
  if (!is.null(pa) && !is.null(pa$nfact) && pa$nfact >= 1) {
    nF_est <- pa$nfact
  } else {
    # Fallback: Kaiser criterion
    ev <- eigen(cor(mat, use = "pairwise.complete.obs"),
                symmetric = TRUE, only.values = TRUE)$values
    nF_est <- max(1, sum(ev > 1))
  }
  nF_est <- max(1, min(nF_est, floor(p / 2)))

  if (progress) cat(sprintf("  EFA: estimating %d factors...\n", nF_est))

  # Step 2: Run EFA
  efa_result <- tryCatch({
    suppressWarnings(
      fa(mat, nfactors = nF_est, rotate = "oblimin", fm = "minres",
         scores = "none", warnings = FALSE)
    )
  }, error = function(e) {
    # Try with fewer factors
    tryCatch({
      suppressWarnings(
        fa(mat, nfactors = max(1, nF_est - 1), rotate = "oblimin", fm = "minres",
           scores = "none", warnings = FALSE)
      )
    }, error = function(e2) NULL)
  })

  if (is.null(efa_result)) {
    warning("EFA failed — falling back to weighted mode (all pairs).")
    return(NULL)
  }

  # Step 3: Assign items to primary factor
  loadings_mat <- as.matrix(efa_result$loadings[])
  nF_actual <- ncol(loadings_mat)
  primary_factor <- apply(abs(loadings_mat), 1, which.max)

  # Step 4: Generate within-factor pairs
  pair_list <- list()
  for (f in seq_len(nF_actual)) {
    items_f <- which(primary_factor == f)
    if (length(items_f) < 2) next
    combos <- combn(items_f, 2)
    for (ci in seq_len(ncol(combos)))
      pair_list[[length(pair_list) + 1]] <- c(combos[1, ci], combos[2, ci])
  }

  if (length(pair_list) == 0) {
    warning("EFA produced no within-factor pairs — falling back to weighted mode.")
    return(NULL)
  }

  pair_mat <- do.call(rbind, pair_list)
  idx_A <- pair_mat[, 1]
  idx_B <- pair_mat[, 2]

  # Step 5: Weight by observed |r|
  raw_cor_mat <- cor(mat, use = "pairwise.complete.obs")
  pair_r <- numeric(length(idx_A))
  for (j in seq_along(idx_A)) pair_r[j] <- raw_cor_mat[idx_A[j], idx_B[j]]
  pair_abs_r <- abs(pair_r)
  pair_sign <- sign(pair_r)

  # Remove NA pairs
  keep <- !is.na(pair_abs_r)
  idx_A <- idx_A[keep]; idx_B <- idx_B[keep]
  pair_abs_r <- pair_abs_r[keep]; pair_sign <- pair_sign[keep]

  k <- length(idx_A)
  weights <- pair_abs_r
  weights[weights < 1e-6] <- 1e-6

  if (progress) {
    cat(sprintf("  EFA: %d factors, %d within-factor pairs, mean |r|=%.3f\n",
                nF_actual, k, mean(pair_abs_r)))
  }

  list(idx_A = idx_A, idx_B = idx_B, weights = weights,
       pair_sign = pair_sign, nF_detected = nF_est, k = k)
}

ReReReRe <- function(data, #any dataset with questionnaire data
                    corProp=0.03, # proportion of highest correlations (legacy, for coupled mode)
                    cutOff=0.99, # the severity of the evaluation (legacy, for percentile)
                    z_threshold=1.5, # z-score threshold for flagging (or "auto")
                    iterations=100,
                    min_pairs=15, # minimum number of item pairs to use (coupled mode)
                    align_signs=TRUE, # align reverse-coded items using sample correlation signs
                    auto_z=FALSE, # auto-calibrate z_threshold based on total_items
                    mode="auto", # "auto" (default: weighted≤60, coupled>60), "coupled", "weighted", "efa_d"
                    min_r=0.0, # minimum |r| to include a pair in weighted mode (0 = all pairs)
                    progress = F){

  #keep only numeric values
  data <- data[, sapply(data, is.numeric), drop = FALSE]
  N <- nrow(data)
  J <- ncol(data)

  # --- Resolve mode ---
  mode <- match.arg(mode, c("auto", "coupled", "weighted", "efa_d"))
  if (mode == "auto") {
    mode <- if (J <= 60) "weighted" else "coupled"
  }

  if (progress) cat(sprintf("  Mode: %s (%d items)\n", mode, J))

  # --- Auto-calibration based on total_items ---
  nFactors_detected <- NA
  if (auto_z || identical(z_threshold, "auto")) {
    total_items <- J

    # Parallel analysis will be run by EFA anyway, but for non-EFA modes we need it
    if (mode != "efa_d") {
      pa <- tryCatch({
        suppressMessages(suppressWarnings(
          fa.parallel(data, fa = "fa", plot = FALSE, n.iter = 20)
        ))
      }, error = function(e) NULL)
      if (!is.null(pa)) nFactors_detected <- pa$nfact
    }

    z_threshold <- get_calibrated_z(total_items, nf = nFactors_detected)

    if (progress) {
      cat(sprintf("  Auto-calibration: %d items -> z_threshold = %.2f\n",
                  total_items, z_threshold))
    }
  }

  # Ensure z_threshold is numeric
  z_threshold <- as.numeric(z_threshold)

  #convert to matrix for faster operations
  mat <- as.matrix(data)

  # Precompute per-item observed max for reverse coding (align_signs).
  if (align_signs) {
    item_max <- apply(mat, 2, max, na.rm = TRUE)
  }

  # Precompute correlation matrices
  rawCorMat <- cor(mat, use = "pairwise.complete.obs")
  corMat <- abs(rawCorMat)
  corMat[upper.tri(corMat, diag = TRUE)] <- NA
  rawCorMat[upper.tri(rawCorMat, diag = TRUE)] <- NA

  if (align_signs) {
    signMat <- sign(rawCorMat)
  }

  # ====================================================================
  # MODE: EFA-D (default) — EFA within-factor pairs + observed |r| weights
  # ====================================================================
  if (mode == "efa_d") {

    efa_info <- .efa_pairs(mat, align_signs = align_signs, progress = progress)

    if (!is.null(efa_info)) {
      # EFA succeeded — use within-factor pairs
      idx_A <- efa_info$idx_A
      idx_B <- efa_info$idx_B
      weights <- efa_info$weights
      pair_sign <- efa_info$pair_sign
      nFactors_detected <- efa_info$nF_detected
      k <- efa_info$k

      A_coupled <- mat[, idx_A, drop = FALSE]
      B_coupled <- mat[, idx_B, drop = FALSE]

      # Sign alignment using observed correlation sign
      if (align_signs) {
        needs_flip <- which(pair_sign < 0)
        if (length(needs_flip) > 0) {
          flip_max <- item_max[idx_B[needs_flip]]
          B_coupled[, needs_flip] <- rep(flip_max + 1, each = N) - B_coupled[, needs_flip]
        }
      }

      # Weighted coherence score
      rowCors <- rowCor_weighted(A_coupled, B_coupled, weights)

      # Permutation baseline: k random pairs, same weights
      all_RIC <- matrix(NA_real_, nrow = N, ncol = iterations)

      if (progress) pb <- txtProgressBar(min = 0, max = iterations, style = 3)

      for (i in seq_len(iterations)) {
        if (progress) setTxtProgressBar(pb, i)

        rand_idx1 <- sample(J, k, replace = TRUE)
        rand_idx2 <- sample(J, k, replace = TRUE)
        same <- rand_idx1 == rand_idx2
        while (any(same)) {
          rand_idx2[same] <- sample(J, sum(same), replace = TRUE)
          same <- rand_idx1 == rand_idx2
        }

        A_rand <- mat[, rand_idx1, drop = FALSE]
        B_rand <- mat[, rand_idx2, drop = FALSE]

        if (align_signs) {
          ri <- pmax(rand_idx1, rand_idx2)
          ci <- pmin(rand_idx1, rand_idx2)
          rand_signs <- signMat[cbind(ri, ci)]
          rand_signs[is.na(rand_signs)] <- 1
          rand_neg <- which(rand_signs < 0)
          if (length(rand_neg) > 0) {
            rand_flip_max <- item_max[rand_idx2[rand_neg]]
            B_rand[, rand_neg] <- rep(rand_flip_max + 1, each = N) - B_rand[, rand_neg]
          }
        }

        all_RIC[, i] <- rowCor_weighted(A_rand, B_rand, weights)
      }

      if (progress) close(pb)

      mode_used <- "efa_d"

    } else {
      # EFA failed — fallback to weighted (all pairs)
      mode <- "weighted"
      if (progress) cat("  EFA failed, falling back to weighted mode.\n")
    }
  }

  # ====================================================================
  # LEGACY MODE: WEIGHTED (all pairs, |r|-weighted)
  # Also used as fallback when EFA fails
  # ====================================================================
  if (mode == "weighted") {

    pair_idx <- which(!is.na(corMat), arr.ind = TRUE)
    pair_r <- rawCorMat[pair_idx]
    pair_abs_r <- corMat[pair_idx]
    pair_sign <- sign(pair_r)

    keep <- pair_abs_r >= min_r & !is.na(pair_abs_r)
    pair_idx <- pair_idx[keep, , drop = FALSE]
    pair_abs_r <- pair_abs_r[keep]
    pair_sign <- pair_sign[keep]

    k <- nrow(pair_idx)
    weights <- pair_abs_r

    if (progress) {
      cat(sprintf("  Weighted mode: %d pairs, mean |r|=%.3f\n", k, mean(pair_abs_r)))
    }

    A_coupled <- mat[, pair_idx[, 1], drop = FALSE]
    B_coupled <- mat[, pair_idx[, 2], drop = FALSE]

    n_negative <- sum(pair_sign < 0)
    if (align_signs && n_negative > 0) {
      needs_flip <- which(pair_sign < 0)
      flip_max <- item_max[pair_idx[needs_flip, 2]]
      B_coupled[, needs_flip] <- rep(flip_max + 1, each = N) - B_coupled[, needs_flip]
    }

    rowCors <- rowCor_weighted(A_coupled, B_coupled, weights)

    all_RIC <- matrix(NA_real_, nrow = N, ncol = iterations)
    if (progress) pb <- txtProgressBar(min = 0, max = iterations, style = 3)

    for (i in seq_len(iterations)) {
      if (progress) setTxtProgressBar(pb, i)
      rand_idx1 <- sample(J, k, replace = TRUE)
      rand_idx2 <- sample(J, k, replace = TRUE)
      same <- rand_idx1 == rand_idx2
      while (any(same)) { rand_idx2[same] <- sample(J, sum(same), replace = TRUE); same <- rand_idx1 == rand_idx2 }
      A_rand <- mat[, rand_idx1, drop = FALSE]
      B_rand <- mat[, rand_idx2, drop = FALSE]
      if (align_signs) {
        ri <- pmax(rand_idx1, rand_idx2); ci <- pmin(rand_idx1, rand_idx2)
        rand_signs <- signMat[cbind(ri, ci)]; rand_signs[is.na(rand_signs)] <- 1
        rand_neg <- which(rand_signs < 0)
        if (length(rand_neg) > 0) {
          B_rand[, rand_neg] <- rep(item_max[rand_idx2[rand_neg]] + 1, each = N) - B_rand[, rand_neg]
        }
      }
      all_RIC[, i] <- rowCor_weighted(A_rand, B_rand, weights)
    }
    if (progress) close(pb)

    mode_used <- "weighted"
    nFactors_detected <- NA
  }

  # ====================================================================
  # LEGACY MODE: COUPLED (top-k% pairs by |r|)
  # ====================================================================
  if (mode == "coupled") {

    corThreshold <- quantile(corMat, 1 - corProp, na.rm = TRUE)
    couples <- which(corMat >= corThreshold, arr.ind = TRUE)
    k <- nrow(couples)

    if (k < min_pairs) {
      all_cors <- corMat[lower.tri(corMat)]
      all_cors <- all_cors[!is.na(all_cors)]
      n_available <- length(all_cors)
      use_k <- min(min_pairs, n_available)
      if (use_k > 0) {
        corThreshold <- sort(all_cors, decreasing = TRUE)[use_k]
        couples <- which(corMat >= corThreshold, arr.ind = TRUE)
        k <- nrow(couples)
      }
    }

    if (k == 0) {
      warning("No item pairs above threshold. Returning NA.")
      return(data.frame(result = rep(NA, N), indCors = rep(NA, N), flagged = rep(NA, N)))
    }

    coupled_signs <- sign(rawCorMat[couples])
    n_negative <- sum(coupled_signs < 0)

    if (n_negative > 0 && !align_signs) {
      warning(sprintf(
        "ReReReRe: %d of %d coupled pairs (%.0f%%) have negative sample correlations. ",
        n_negative, k, 100 * n_negative / k),
        "Consider setting align_signs=TRUE or reverse-coding items.")
    }

    A_coupled <- mat[, couples[, 1], drop = FALSE]
    B_coupled <- mat[, couples[, 2], drop = FALSE]

    if (align_signs && n_negative > 0) {
      needs_flip <- which(coupled_signs < 0)
      flip_max <- item_max[couples[needs_flip, 2]]
      B_coupled[, needs_flip] <- rep(flip_max + 1, each = N) - B_coupled[, needs_flip]
    }

    rowCors <- rowCor_abs(A_coupled, B_coupled, zero_val = 0)

    all_pairs <- combn(J, 2)
    n_pairs_total <- ncol(all_pairs)
    all_RIC <- matrix(NA_real_, nrow = N, ncol = iterations)
    if (progress) pb <- txtProgressBar(min = 0, max = iterations, style = 3)

    for (i in seq_len(iterations)) {
      if (progress) setTxtProgressBar(pb, i)
      sampled <- sample(n_pairs_total, k)
      idx1 <- all_pairs[1, sampled]; idx2 <- all_pairs[2, sampled]
      A_rand <- mat[, idx1, drop = FALSE]; B_rand <- mat[, idx2, drop = FALSE]
      if (align_signs) {
        ri <- pmax(idx1, idx2); ci <- pmin(idx1, idx2)
        rand_signs <- signMat[cbind(ri, ci)]; rand_signs[is.na(rand_signs)] <- 1
        rand_neg <- which(rand_signs < 0)
        if (length(rand_neg) > 0) {
          B_rand[, rand_neg] <- rep(item_max[idx2[rand_neg]] + 1, each = N) - B_rand[, rand_neg]
        }
      }
      all_RIC[, i] <- rowCor_abs(A_rand, B_rand, zero_val = 0)
    }
    if (progress) close(pb)

    mode_used <- "coupled"
    nFactors_detected <- NA
  }

  corComparedIndex <- rowMeans(all_RIC < rowCors, na.rm = TRUE)

  # --- Z-score: how many SDs is indCors above the respondent's random baseline? ---
  rand_means <- rowMeans(all_RIC, na.rm = TRUE)
  rand_sds   <- apply(all_RIC, 1, sd, na.rm = TRUE)

  z_score <- ifelse(rand_sds > 0,
                    (rowCors - rand_means) / rand_sds,
                    0)

  data.frame(
    result = corComparedIndex,     # legacy percentile score
    z_score = z_score,             # z-score: SDs above random baseline
    indCors = rowCors,             # raw coupled/weighted coherence score
    rand_mean = rand_means,        # mean of random iterations
    rand_sd = rand_sds,            # sd of random iterations
    flagged = z_score <= z_threshold,  # z-score flagging (primary)
    z_threshold_used = z_threshold,    # threshold used
    nFactors_detected = nFactors_detected,  # from parallel analysis
    mode_used = mode_used,         # "efa_d", "coupled", or "weighted"
    n_pairs = k                    # number of item pairs used
  )
}


# ======================================================================
# ReReReRe_F: Per-factor variant (experimental)
# ======================================================================
#
# Uses EFA (parallel analysis + oblimin rotation) to:
#   1. Assign items to factors
#   2. Select only within-factor pairs
#   3. Weight by observed |r|
#
# Advantages over standard ReReReRe:
#   - Wins on simulated data overall (MCC=0.333 vs 0.305)
#   - Wins on short questionnaires (<60 items) even on real data
#   - Theoretically principled: only uses pairs that "should" correlate
#
# Disadvantages:
#   - Loses on real validation datasets with ground truth (AUC 0.561 vs 0.635)
#   - Depends on EFA quality, which degrades with noisy real-world data
#   - Adds computational cost (parallel analysis + EFA)
#
# Use case: secondary/complementary analysis, especially for short instruments.
# ======================================================================

ReReReRe_F <- function(data,
                       z_threshold = 1.5,
                       iterations = 100,
                       align_signs = TRUE,
                       auto_z = FALSE,
                       min_r = 0.0,
                       progress = FALSE) {

  # Delegate to ReReReRe with mode="efa_d"
  ReReReRe(data,
           corProp = 0.03,
           z_threshold = z_threshold,
           iterations = iterations,
           align_signs = align_signs,
           auto_z = auto_z,
           mode = "efa_d",
           min_r = min_r,
           progress = progress)
}
