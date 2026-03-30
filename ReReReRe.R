#### ReReReRe: Permutation-based careless respondent detection ####
#
# Returns per-respondent:
#   result    - legacy percentile: P(random < coupled)
#   z_score   - (coupled_cor - mean_random) / sd_random
#   indCors   - raw coupled absolute correlation
#   rand_mean - mean of random iteration correlations
#   rand_sd   - sd of random iteration correlations
#   flagged   - binary flag (z_score <= z_threshold)
#   z_threshold_used - the z_threshold used for flagging (useful when auto_z=TRUE)
#   nFactors_detected - number of factors detected by parallel analysis (when auto_z=TRUE)
#   mode_used - "coupled" (standard top-k%) or "weighted" (all pairs, |r|-weighted)
#
# === REVISION 2026-03-10b ===
# Added z_score output. The percentile (result) compresses toward 1.0
# with larger questionnaires (ceiling effect). The z-score preserves
# how FAR above random the coupled correlation is, not just WHETHER
# it beats random. With more items, z-scores for good respondents
# INCREASE (better precision → smaller SD → larger z), while careless
# respondents stay near 0. This inverts the ceiling effect.
#
# === REVISION 2026-03-10c ===
# Fixed zero_val for random pairs: changed from 1 to 0.
# Old: zero_val=1 was a percentile-specific hack — it forced
#   randomCor=1 for zero-variance rows so that (1 < indCors) = FALSE,
#   pushing the percentile to 0. But this contaminated rand_mean and
#   rand_sd, poisoning the z-score and delta metrics for partial
#   longstring respondents (bimodal random iteration distributions).
# New: zero_val=0 for both coupled and random. Undefined correlation
#   = 0 (no evidence of linear relationship). Percentile still works:
#   pure longstring → coupled=0, random=0, mean(0<0)=0 → flagged.
#
# === REVISION 2026-03-25 ===
# Added align_signs parameter (default TRUE) to fix reverse-coded item
# cancellation. Problem: the sample-level correlation matrix uses abs()
# to select coupled pairs, so pairs with negative correlations (e.g.,
# positively-keyed x negatively-keyed items within the same factor)
# are correctly identified as "strongly associated." But at the
# individual level, rowCor_abs computes a SINGLE correlation across
# all k pairs. Pairs with positive and negative sample-level correlations
# contribute opposite signs, cancelling each other out INSIDE cor()
# before abs() is applied. Result: indCors ≈ 0 for everyone.
#
# Fix: use the sign of the raw (signed) sample correlation to align
# all pairs before computing individual-level correlations. For each
# pair where the sample-level correlation is negative, we reverse-code
# the B column using (max+1-x) rather than multiplying by -1. Simple
# sign flip (-1*x) distorts the centering in rowCor_abs: flipped values
# fall outside the original Likert range, creating disproportionate
# deviations from the mean and giving negative-correlation pairs more
# weight in the individual-level correlation. Proper reverse coding
# keeps all values in the original range, ensuring balanced contribution
# from all pairs. Assumes max(observed) = max(scale), which holds for
# Likert data with reasonable sample sizes. Applied to both coupled
# and random pairs.
#
# === REVISION 2026-03-26b ===
# Changed align_signs from sign flip (*-1) to proper reverse coding
# (max+1-x). Sign flip distorted the individual-level correlation by
# giving negative-correlation pairs disproportionate weight after
# centering. Reverse coding preserves the value range and centering
# balance. Uses per-item observed max as scale max.
#
# === REVISION 2026-03-26c ===
# Added auto_z parameter for automatic z_threshold calibration.
# v1 (nF-only): LOESS on nF from 1D calibration. R²=0.245.
# v2 (total_items, CURRENT): uses total_items = ncol(data) as predictor.
#   Simple linear: z = 0.505 + 0.0042 * total_items (R²=0.476).
#   2D LOESS surface over nF × ipf for when parallel analysis is available.
#   total_items alone explains 2x more variance than nF+ipf separately.
#   Auto-z closes 38% of the gap between fixed z=1.5 and oracle.
#   Calibrated on: 10 nF (4-30) × 6 ipf (3-12) = 60 cells, 30 reps each.
#
# === REVISION 2026-03-28 ===
# Updated defaults based on full multiverse (5,670 RR calls, 85,050 rows):
#   corProp: 0.05 → 0.03 (MCC 0.352 vs 0.327, fewer but more selective pairs)
#   z_threshold: 1.5 (confirmed as robust universal default)
#   auto_z: now uses 2D calibration (total_items based)
# Variable importance (eta²): ipf=31.1%, nF=26.1%, pct=3.6%, corProp=1.7%, n=1.0%
#
# === REVISION 2026-03-30 ===
# Added weighted mode for short questionnaires (≤60 items).
# Standard ReReReRe selects top-k% pairs by |r|; weighted mode uses ALL pairs
# but weights each by its sample-level |r|. This maximizes information when
# few high-correlation pairs are available.
#
# mode parameter: "auto" (default), "coupled", "weighted"
#   "auto"    → weighted if total_items ≤ 60, coupled otherwise
#   "coupled" → standard top-k% pair selection (original algorithm)
#   "weighted"→ all pairs weighted by |r_sample|
#
# Simulation results (18 conditions, nF=4-20, ipf=3-10, 3 reps):
#   ≤30 items:  weighted +37% MCC over standard
#   30-60 items: weighted +24% MCC
#   60-100 items: standard −8% better
#   >100 items:  standard −16% better
# The crossover is at ~60 items, justifying the auto-switch threshold.
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

ReReReRe <- function(data, #any dataset with questionnaire data
                    corProp=0.03, # proportion of highest correlations (0.03 beats 0.05 in multiverse)
                    cutOff=0.99, # the severity of the evaluation (legacy, for percentile)
                    z_threshold=1.5, # z-score threshold for flagging (or "auto")
                    iterations=100,
                    min_pairs=15, # minimum number of item pairs to use
                    align_signs=TRUE, # align reverse-coded items using sample correlation signs
                    auto_z=FALSE, # auto-calibrate z_threshold based on total_items
                    mode="auto", # "auto", "coupled" (standard top-k%), "weighted" (all pairs, |r|-weighted)
                    min_r=0.0, # minimum |r| to include a pair in weighted mode (0 = all pairs)
                    progress = F){

  #keep only numeric values
  data <- data[, sapply(data, is.numeric), drop = FALSE]
  N <- nrow(data)
  J <- ncol(data)

  # --- Resolve mode: auto selects weighted (≤60 items) or coupled (>60 items) ---
  mode <- match.arg(mode, c("auto", "coupled", "weighted"))
  use_weighted <- FALSE
  if (mode == "auto") {
    use_weighted <- (J <= 60)
    if (progress) {
      cat(sprintf("  Mode: auto → %s (%d items %s 60)\n",
                  ifelse(use_weighted, "weighted", "coupled"),
                  J, ifelse(use_weighted, "<=", ">")))
    }
  } else if (mode == "weighted") {
    use_weighted <- TRUE
  }

  # --- Auto-calibration based on total_items ---
  # Uses total_items (ncol) as the primary predictor for optimal z_threshold.
  # Optionally also runs parallel analysis to estimate nF for the output.
  # total_items explains 2x more variance than nF+ipf separately (R²=0.476).
  nFactors_detected <- NA
  if (auto_z || identical(z_threshold, "auto")) {
    total_items <- J  # J = ncol(data), already computed above

    # Try parallel analysis to estimate nF (for output, and as secondary info)
    pa <- tryCatch({
      suppressMessages(suppressWarnings(
        fa.parallel(data, fa = "fa", plot = FALSE, n.iter = 20)
      ))
    }, error = function(e) NULL)

    if (!is.null(pa)) {
      nFactors_detected <- pa$nfact
    }

    # Calibrate z from total_items (the best predictor)
    z_threshold <- get_calibrated_z(total_items, nf = nFactors_detected)

    if (progress) {
      cat(sprintf("  Auto-calibration: %d items, nF=%s -> z_threshold = %.2f\n",
                  total_items,
                  ifelse(is.na(nFactors_detected), "?", as.character(nFactors_detected)),
                  z_threshold))
    }
  } else if (is.character(z_threshold) && z_threshold == "auto") {
    auto_z <- TRUE
  }

  # Ensure z_threshold is numeric at this point
  z_threshold <- as.numeric(z_threshold)

  #convert to matrix for faster operations
  mat <- as.matrix(data)

  # Precompute per-item observed max for reverse coding (align_signs).
  # Assumes max(observed) = max(scale), valid for Likert data with n >= ~50.
  if (align_signs) {
    item_max <- apply(mat, 2, max, na.rm = TRUE)
  }

  #### COMPUTING COUPLED CORRELATION ####

  #computing correlations: raw (signed) and absolute
  rawCorMat <- cor(mat, use = "pairwise.complete.obs")
  corMat <- abs(rawCorMat)

  #substituting upper triangle (diagonal included) with NAs
  corMat[upper.tri(corMat, diag = TRUE)] <- NA
  rawCorMat[upper.tri(rawCorMat, diag = TRUE)] <- NA

  # Precompute the full sign matrix for random pair reverse-coding
  if (align_signs) {
    signMat <- sign(rawCorMat)
  }

  if (use_weighted) {
    # ====================================================================
    # WEIGHTED MODE: use ALL pairs, weighted by |r_sample|
    # Better for short questionnaires (≤60 items) where few high-|r| pairs
    # exist. Every pair contributes, but strong correlations count more.
    # ====================================================================

    # Get all pairs from lower triangle
    pair_idx <- which(!is.na(corMat), arr.ind = TRUE)
    pair_r <- rawCorMat[pair_idx]
    pair_abs_r <- corMat[pair_idx]
    pair_sign <- sign(pair_r)

    # Filter by min_r (default 0 = all pairs)
    keep <- pair_abs_r >= min_r & !is.na(pair_abs_r)
    pair_idx <- pair_idx[keep, , drop = FALSE]
    pair_r <- pair_r[keep]
    pair_abs_r <- pair_abs_r[keep]
    pair_sign <- pair_sign[keep]

    k <- nrow(pair_idx)
    weights <- pair_abs_r

    if (progress) {
      cat(sprintf("  Weighted mode: using %d pairs (min_r=%.2f), mean |r|=%.3f\n",
                  k, min_r, mean(pair_abs_r)))
    }

    # Build A and B matrices
    A_coupled <- mat[, pair_idx[, 1], drop = FALSE]
    B_coupled <- mat[, pair_idx[, 2], drop = FALSE]

    # Sign alignment: reverse-code B for negative-correlation pairs
    n_negative <- sum(pair_sign < 0)
    if (align_signs && n_negative > 0) {
      needs_flip <- which(pair_sign < 0)
      flip_max <- item_max[pair_idx[needs_flip, 2]]
      B_coupled[, needs_flip] <- rep(flip_max + 1, each = N) - B_coupled[, needs_flip]
    }

    # Weighted coherence score (coupled)
    rowCors <- rowCor_weighted(A_coupled, B_coupled, weights)

    # --- Random permutation baseline (weighted) ---
    all_RIC <- matrix(NA_real_, nrow = N, ncol = iterations)

    if (progress) pb <- txtProgressBar(min = 0, max = iterations, style = 3)

    for (i in seq_len(iterations)) {
      if (progress) setTxtProgressBar(pb, i)

      # Sample k random pairs (ensuring different items)
      rand_idx1 <- sample(J, k, replace = TRUE)
      rand_idx2 <- sample(J, k, replace = TRUE)
      same <- rand_idx1 == rand_idx2
      while (any(same)) {
        rand_idx2[same] <- sample(J, sum(same), replace = TRUE)
        same <- rand_idx1 == rand_idx2
      }

      A_rand <- mat[, rand_idx1, drop = FALSE]
      B_rand <- mat[, rand_idx2, drop = FALSE]

      # Sign alignment for random pairs
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

      # Use SAME weights as coupled (keeps comparison fair)
      all_RIC[, i] <- rowCor_weighted(A_rand, B_rand, weights)
    }

    if (progress) close(pb)

  } else {
    # ====================================================================
    # COUPLED MODE (standard): select top-k% pairs by |r|
    # Better for longer questionnaires (>60 items) where high-|r| pairs
    # provide a cleaner signal than including weak pairs as noise.
    # ====================================================================

    #getting the threshold correlation based on the proportion of highest correlation we decided to include.
    corThreshold <- quantile(corMat, 1 - corProp, na.rm = TRUE)

    #getting the row and col indices with values higher than threshold
    couples <- which(corMat >= corThreshold, arr.ind = TRUE)
    k <- nrow(couples)

    # Enforce minimum number of pairs for stable individual-level correlations.
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

    # --- Sign alignment for reverse-coded items ---
    coupled_signs <- sign(rawCorMat[couples])
    n_negative <- sum(coupled_signs < 0)

    if (n_negative > 0 && !align_signs) {
      warning(sprintf(
        "ReReReRe: %d of %d coupled pairs (%.0f%%) have negative sample correlations. ",
        n_negative, k, 100 * n_negative / k),
        "This typically indicates reverse-coded items that will cancel each other ",
        "in the individual-level correlation. Consider setting align_signs=TRUE ",
        "or reverse-coding items before running ReReReRe.")
    }

    # Extract N x k matrices for the coupled pairs
    A_coupled <- mat[, couples[, 1], drop = FALSE]
    B_coupled <- mat[, couples[, 2], drop = FALSE]

    # Apply sign alignment: reverse-code B columns where sample correlation is negative
    if (align_signs && n_negative > 0) {
      needs_flip <- which(coupled_signs < 0)
      flip_max <- item_max[couples[needs_flip, 2]]
      B_coupled[, needs_flip] <- rep(flip_max + 1, each = N) - B_coupled[, needs_flip]
    }

    # zero_val=0: for zero-variance rows (longstring), return 0 so they are easily spotted
    rowCors <- rowCor_abs(A_coupled, B_coupled, zero_val = 0)

    #### RANDOM PERMUTATION ITERATIONS ####

    all_pairs <- combn(J, 2) # 2 x C(J,2) matrix
    n_pairs <- ncol(all_pairs)

    all_RIC <- matrix(NA_real_, nrow = N, ncol = iterations)

    if (progress) pb <- txtProgressBar(min = 0, max = iterations, style = 3)

    for (i in seq_len(iterations)){
      if (progress) setTxtProgressBar(pb, i)

      sampled <- sample(n_pairs, k)
      idx1 <- all_pairs[1, sampled]
      idx2 <- all_pairs[2, sampled]

      A_rand <- mat[, idx1, drop = FALSE]
      B_rand <- mat[, idx2, drop = FALSE]

      if (align_signs) {
        ri <- pmax(idx1, idx2)
        ci <- pmin(idx1, idx2)
        rand_signs <- signMat[cbind(ri, ci)]
        rand_signs[is.na(rand_signs)] <- 1
        rand_neg <- which(rand_signs < 0)
        if (length(rand_neg) > 0) {
          rand_flip_max <- item_max[idx2[rand_neg]]
          B_rand[, rand_neg] <- rep(rand_flip_max + 1, each = N) - B_rand[, rand_neg]
        }
      }

      all_RIC[, i] <- rowCor_abs(A_rand, B_rand, zero_val = 0)
    }

    if (progress) close(pb)
  }

  corComparedIndex <- rowMeans(all_RIC < rowCors, na.rm = TRUE)

  # --- Z-score: how many SDs is indCors above the respondent's random baseline? ---
  # mean and sd of random iterations per respondent (row-wise)
  rand_means <- rowMeans(all_RIC, na.rm = TRUE)
  rand_sds   <- apply(all_RIC, 1, sd, na.rm = TRUE)

  # z = (coupled_cor - mean_random) / sd_random
  # When sd is 0 (all random iterations identical, very rare), return 0
  z_score <- ifelse(rand_sds > 0,
                    (rowCors - rand_means) / rand_sds,
                    0)

  mode_used <- ifelse(use_weighted, "weighted", "coupled")

  data.frame(
    result = corComparedIndex,     # legacy percentile score
    z_score = z_score,             # z-score: SDs above random baseline
    indCors = rowCors,             # raw coupled/weighted coherence score
    rand_mean = rand_means,        # mean of random iterations
    rand_sd = rand_sds,            # sd of random iterations
    flagged = z_score <= z_threshold,  # z-score flagging (primary)
    z_threshold_used = z_threshold,    # threshold used (useful when auto_z=TRUE)
    nFactors_detected = nFactors_detected,  # from parallel analysis (NA if not auto)
    mode_used = mode_used          # "coupled" or "weighted"
  )
}
