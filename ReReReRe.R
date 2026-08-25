#### ReReReRe: Permutation-based careless respondent detection ####
#
# Two public functions:
#   ReReReRe()   — standard 2-level method (default, recommended)
#   ReReReRe_F() — per-factor variant using EFA (experimental)
#
# Returns per-respondent:
#   result       - legacy percentile: P(random < coupled)
#   z_score      - (coupled_cor - mean_random) / sd_random, MINUS variance penalty if enabled
#   z_score_raw  - z_score WITHOUT variance penalty (always reported for reference)
#   indCors      - raw coupled absolute correlation
#   rand_mean    - mean of random iteration correlations
#   rand_sd      - sd of random iteration correlations
#   flagged      - binary flag (z_score <= z_threshold)
#   z_threshold_used - the z_threshold used for flagging (useful when auto_z=TRUE)
#   nFactors_detected - number of factors detected by parallel analysis
#   mode_used    - "coupled", "weighted", or "efa_d"
#   n_pairs      - number of item pairs used
#   variance_penalty_used - TRUE/FALSE
#
# === VARIANCE PENALTY (2026-04-21) ===
#
# Optional flag to detect straight-lining and acquiescent response patterns,
# which the base z-score alone identifies only weakly. When variance_penalty=TRUE:
#
#   z_score = z_score_raw - vp_alpha * exp(-sd_respondent / vp_beta)
#
# A respondent with zero within-person SD (pure straight-liner) receives the full
# penalty (vp_alpha), pushing their z-score strongly negative. Respondents with
# normal response variance are essentially unaffected.
#
# Tuned on simulated data (α=3, β=0.5): +0.053 MCC on pure_straight,
# +0.030 on acquiescent, 0.000 on random/longstring (no damage).
# Defaults: vp_alpha=3.0, vp_beta=0.5.
#
# === SPLIT-HALF (2026-04-21) ===
#
# Available in ReReReRe_SplitHalf.R as score_split_half(). Splits the top-k%
# coupled pairs into halves across B random splits and scores each half
# separately. Two aggregation strategies:
#   aggregation="min":  returns min z across all 2*B halves — catches partial
#                        carelessness aggressively (sensitivity 0.97-0.99 on
#                        random/longstring/mixed at corruption>50%) but flags
#                        many low-corruption respondents (specificity 0.46).
#                        Better used as a continuous "degree of carelessness"
#                        score than a binary flag.
#   aggregation="mean": returns mean z across halves — best binary classifier
#                        overall (MCC=0.589 vs std=0.567, std+VP=0.582).
#
# Use when partial carelessness (60-80% corruption) is the main concern.
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
# 2026-04-27:  Added rescale parameter (default "proportion") for mixed-scale
#              robustness. Generalized align_signs reverse-coding from
#              (item_max + 1) - x to (item_max + item_min) - x so it works on
#              any rescaled scale. Identical to legacy on raw 1..K Likert.
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
       pair_sign = pair_sign, nF_detected = nF_est, k = k,
       primary_factor = primary_factor, nF_actual = nF_actual)
}

# --- Helper: sample cross-factor pairs (idea 2 from buone_idee.md) ---
# Samples k pairs where the two items belong to DIFFERENT factors.
# This ensures the permutation baseline is purely cross-factor noise.
#
# primary_factor: integer vector of length J (factor assignment per item)
# k: number of pairs to sample
# J: total number of items
# Returns list(idx1, idx2)
.sample_cross_factor_pairs <- function(primary_factor, k, J) {
  # Precompute cross-factor items for each factor
  factors <- unique(primary_factor)
  if (length(factors) < 2) {
    # Only 1 factor: can't do cross-factor, fall back to unrestricted
    idx1 <- sample(J, k, replace = TRUE)
    idx2 <- sample(J, k, replace = TRUE)
    same <- idx1 == idx2
    while (any(same)) { idx2[same] <- sample(J, sum(same), replace = TRUE); same <- idx1 == idx2 }
    return(list(idx1 = idx1, idx2 = idx2))
  }

  cross_items <- lapply(factors, function(f) which(primary_factor != f))
  names(cross_items) <- as.character(factors)

  idx1 <- sample(J, k, replace = TRUE)
  idx2 <- integer(k)

  for (i in seq_len(k)) {
    pool <- cross_items[[as.character(primary_factor[idx1[i]])]]
    idx2[i] <- pool[sample.int(length(pool), 1)]
  }

  list(idx1 = idx1, idx2 = idx2)
}

# --- Helper: compute per-factor z-scores ---
# Runs EFA, computes z-score per factor, returns matrix + metadata.
# Used by per-factor variant functions.
#
# Returns list(factor_z_matrix, primary_factor, nF_detected, n_factors_used, efa_info)
.compute_per_factor_z <- function(mat, iterations = 50, align_signs = TRUE,
                                   min_items_per_factor = 3,
                                   cross_factor_baseline = FALSE,
                                   external_efa = NULL) {
  N <- nrow(mat); p <- ncol(mat)
  item_max <- if (align_signs) apply(mat, 2, max, na.rm = TRUE) else NULL

  # Step 1: EFA (or use external)
  if (!is.null(external_efa)) {
    primary_factor <- external_efa$primary_factor
    nF_est <- external_efa$nF_detected
    nF_actual <- length(unique(primary_factor))
  } else {
    pa <- tryCatch({
      suppressMessages(suppressWarnings(
        fa.parallel(mat, fa = "fa", plot = FALSE, n.iter = 20)
      ))
    }, error = function(e) NULL)

    nF_est <- if (!is.null(pa) && !is.null(pa$nfact) && pa$nfact >= 1) pa$nfact
              else max(1, sum(eigen(cor(mat, use="pairwise.complete.obs"),
                                     symmetric=TRUE, only.values=TRUE)$values > 1))
    nF_est <- max(1, min(nF_est, floor(p / 2)))

    efa_result <- tryCatch({
      suppressWarnings(fa(mat, nfactors = nF_est, rotate = "oblimin", fm = "minres",
                          scores = "none", warnings = FALSE))
    }, error = function(e) {
      tryCatch({
        suppressWarnings(fa(mat, nfactors = max(1, nF_est-1), rotate = "oblimin",
                            fm = "minres", scores = "none", warnings = FALSE))
      }, error = function(e2) NULL)
    })

    if (is.null(efa_result)) return(NULL)

    loadings_mat <- as.matrix(efa_result$loadings[])
    nF_actual <- ncol(loadings_mat)
    primary_factor <- apply(abs(loadings_mat), 1, which.max)
  }

  # Precompute
  raw_cor_mat <- cor(mat, use = "pairwise.complete.obs")
  signMat <- sign(raw_cor_mat)
  signMat[upper.tri(signMat, diag = TRUE)] <- NA

  # Step 2: Per-factor computation
  factor_z_matrix <- matrix(NA_real_, nrow = N, ncol = nF_actual)
  factors_used <- 0

  for (f in seq_len(nF_actual)) {
    items_f <- which(primary_factor == f)
    if (length(items_f) < min_items_per_factor) next

    combos <- combn(items_f, 2)
    idx_A <- combos[1, ]; idx_B <- combos[2, ]
    k_f <- length(idx_A)

    pair_r <- numeric(k_f)
    for (j in seq_len(k_f)) pair_r[j] <- raw_cor_mat[idx_A[j], idx_B[j]]
    pair_abs_r <- abs(pair_r); pair_sign <- sign(pair_r)

    keep <- !is.na(pair_abs_r)
    if (sum(keep) < 2) next
    idx_A <- idx_A[keep]; idx_B <- idx_B[keep]
    pair_abs_r <- pair_abs_r[keep]; pair_sign <- pair_sign[keep]
    k_f <- length(idx_A)
    weights <- pair_abs_r; weights[weights < 1e-6] <- 1e-6

    A_f <- mat[, idx_A, drop = FALSE]; B_f <- mat[, idx_B, drop = FALSE]
    if (align_signs) {
      nf <- which(pair_sign < 0)
      if (length(nf) > 0) B_f[, nf] <- rep(item_max[idx_B[nf]] + 1, each = N) - B_f[, nf]
    }

    coupled_f <- rowCor_weighted(A_f, B_f, weights)

    # Permutation baseline
    rand_f <- matrix(NA_real_, N, iterations)
    for (iter in seq_len(iterations)) {
      if (cross_factor_baseline && length(unique(primary_factor)) >= 2) {
        cf <- .sample_cross_factor_pairs(primary_factor, k_f, p)
        ri1 <- cf$idx1; ri2 <- cf$idx2
      } else {
        ri1 <- sample(p, k_f, replace = TRUE); ri2 <- sample(p, k_f, replace = TRUE)
        same <- ri1 == ri2
        while (any(same)) { ri2[same] <- sample(p, sum(same), replace = TRUE); same <- ri1 == ri2 }
      }
      Ar <- mat[, ri1, drop = FALSE]; Br <- mat[, ri2, drop = FALSE]
      if (align_signs) {
        ri <- pmax(ri1, ri2); ci <- pmin(ri1, ri2)
        rs <- signMat[cbind(ri, ci)]; rs[is.na(rs)] <- 1
        rn <- which(rs < 0)
        if (length(rn) > 0) Br[, rn] <- rep(item_max[ri2[rn]] + 1, each = N) - Br[, rn]
      }
      rand_f[, iter] <- rowCor_weighted(Ar, Br, weights)
    }

    rm_f <- rowMeans(rand_f, na.rm = TRUE)
    rsd_f <- apply(rand_f, 1, sd, na.rm = TRUE); rsd_f[rsd_f == 0] <- 1e-10
    factor_z_matrix[, f] <- (coupled_f - rm_f) / rsd_f
    factors_used <- factors_used + 1
  }

  if (factors_used == 0) return(NULL)

  valid_cols <- which(colSums(!is.na(factor_z_matrix)) > 0)
  fz <- factor_z_matrix[, valid_cols, drop = FALSE]

  list(factor_z_matrix = fz,
       primary_factor = primary_factor,
       nF_detected = nF_est,
       n_factors_used = factors_used)
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
                    cross_factor_baseline=FALSE, # use cross-factor only random pairs (requires EFA)
                    variance_penalty=FALSE, # penalize z-score for low within-respondent SD (detects straight-lining/acquiescence)
                    vp_alpha=3.0, # variance-penalty strength (only used if variance_penalty=TRUE)
                    vp_beta=0.5,  # variance-penalty decay (only used if variance_penalty=TRUE)
                    rescale = c("proportion", "none", "minmax", "zscore"), # per-item rescaling for mixed-scale robustness
                    progress = F){

  rescale <- match.arg(rescale)

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

  # --- Per-item rescaling (mixed Likert scales) ---
  # Pearson correlation is invariant to identical affine transforms applied to
  # all items, so per-item rescaling has zero cost on homogeneous-scale data
  # (every column gets divided by the same constant). On heterogeneous scales
  # (e.g. some items 1-5, some 1-7) it removes the bias toward longer scales
  # in the row-wise correlation used by the coupled mode.
  #
  # "proportion" (default): x / max(x). Cheapest and provably non-destructive
  #   on uniform-scale Likert data. Recommended for general use.
  # "minmax": (x - min) / (max - min). Maps every item to [0, 1].
  # "zscore": (x - mean) / sd. Most aggressive; removes mean/spread differences.
  # "none": legacy behaviour (raw Likert values).
  if (rescale != "none" && J >= 1) {
    item_min_pre <- apply(mat, 2, min, na.rm = TRUE)
    item_max_pre <- apply(mat, 2, max, na.rm = TRUE)

    if (rescale == "proportion") {
      denom <- item_max_pre
      denom[!is.finite(denom) | denom == 0] <- 1
      mat <- sweep(mat, 2, denom, "/")
    } else if (rescale == "minmax") {
      rng <- item_max_pre - item_min_pre
      rng[!is.finite(rng) | rng == 0] <- 1
      mat <- sweep(mat, 2, item_min_pre, "-")
      mat <- sweep(mat, 2, rng, "/")
    } else if (rescale == "zscore") {
      mat <- scale(mat, center = TRUE, scale = TRUE)
      attr(mat, "scaled:center") <- NULL
      attr(mat, "scaled:scale")  <- NULL
      mat[is.nan(mat)] <- 0
    }
  }

  # Mixed-scale heuristic warning when user opts out of rescaling.
  if (rescale == "none" && J >= 2) {
    rng_check <- apply(mat, 2, function(v) diff(range(v, na.rm = TRUE)))
    rng_check <- rng_check[is.finite(rng_check) & rng_check > 0]
    if (length(rng_check) >= 2 && (max(rng_check) / min(rng_check)) >= 1.5) {
      warning("ReReReRe: items appear to span heterogeneous response scales ",
              "(max/min item range >= 1.5). Consider rescale='proportion'.")
    }
  }

  # Precompute per-item observed min/max for reverse coding (align_signs).
  # Computed AFTER rescaling so the reverse formula (item_max + item_min - x)
  # is correct on the rescaled scale. On raw Likert with min=1, this reduces
  # to the legacy (item_max + 1) - x.
  if (align_signs) {
    item_max <- apply(mat, 2, max, na.rm = TRUE)
    item_min <- apply(mat, 2, min, na.rm = TRUE)
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
          flip_min <- item_min[idx_B[needs_flip]]
          B_coupled[, needs_flip] <- rep(flip_max + flip_min, each = N) - B_coupled[, needs_flip]
        }
      }

      # Weighted coherence score
      rowCors <- rowCor_weighted(A_coupled, B_coupled, weights)

      # Permutation baseline: k random pairs, same weights
      all_RIC <- matrix(NA_real_, nrow = N, ncol = iterations)

      if (progress) pb <- txtProgressBar(min = 0, max = iterations, style = 3)

      # Get primary_factor for cross-factor baseline
      efa_primary_factor <- efa_info$primary_factor

      for (i in seq_len(iterations)) {
        if (progress) setTxtProgressBar(pb, i)

        if (cross_factor_baseline && length(unique(efa_primary_factor)) >= 2) {
          cf <- .sample_cross_factor_pairs(efa_primary_factor, k, J)
          rand_idx1 <- cf$idx1; rand_idx2 <- cf$idx2
        } else {
          rand_idx1 <- sample(J, k, replace = TRUE)
          rand_idx2 <- sample(J, k, replace = TRUE)
          same <- rand_idx1 == rand_idx2
          while (any(same)) {
            rand_idx2[same] <- sample(J, sum(same), replace = TRUE)
            same <- rand_idx1 == rand_idx2
          }
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
            rand_flip_min <- item_min[rand_idx2[rand_neg]]
            B_rand[, rand_neg] <- rep(rand_flip_max + rand_flip_min, each = N) - B_rand[, rand_neg]
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

  # --- Cross-factor baseline: get factor assignments if needed ---
  cf_factor <- NULL
  if (cross_factor_baseline && mode %in% c("weighted", "coupled")) {
    efa_for_cf <- .efa_pairs(mat, align_signs = align_signs, progress = FALSE)
    if (!is.null(efa_for_cf)) {
      cf_factor <- efa_for_cf$primary_factor
      nFactors_detected <- efa_for_cf$nF_detected
    } else {
      cross_factor_baseline <- FALSE  # EFA failed, disable
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
      flip_min <- item_min[pair_idx[needs_flip, 2]]
      B_coupled[, needs_flip] <- rep(flip_max + flip_min, each = N) - B_coupled[, needs_flip]
    }

    rowCors <- rowCor_weighted(A_coupled, B_coupled, weights)

    all_RIC <- matrix(NA_real_, nrow = N, ncol = iterations)
    if (progress) pb <- txtProgressBar(min = 0, max = iterations, style = 3)

    for (i in seq_len(iterations)) {
      if (progress) setTxtProgressBar(pb, i)
      if (cross_factor_baseline && !is.null(cf_factor) && length(unique(cf_factor)) >= 2) {
        cf <- .sample_cross_factor_pairs(cf_factor, k, J)
        rand_idx1 <- cf$idx1; rand_idx2 <- cf$idx2
      } else {
        rand_idx1 <- sample(J, k, replace = TRUE)
        rand_idx2 <- sample(J, k, replace = TRUE)
        same <- rand_idx1 == rand_idx2
        while (any(same)) { rand_idx2[same] <- sample(J, sum(same), replace = TRUE); same <- rand_idx1 == rand_idx2 }
      }
      A_rand <- mat[, rand_idx1, drop = FALSE]
      B_rand <- mat[, rand_idx2, drop = FALSE]
      if (align_signs) {
        ri <- pmax(rand_idx1, rand_idx2); ci <- pmin(rand_idx1, rand_idx2)
        rand_signs <- signMat[cbind(ri, ci)]; rand_signs[is.na(rand_signs)] <- 1
        rand_neg <- which(rand_signs < 0)
        if (length(rand_neg) > 0) {
          rb_idx <- rand_idx2[rand_neg]
          B_rand[, rand_neg] <- rep(item_max[rb_idx] + item_min[rb_idx], each = N) - B_rand[, rand_neg]
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
      flip_min <- item_min[couples[needs_flip, 2]]
      B_coupled[, needs_flip] <- rep(flip_max + flip_min, each = N) - B_coupled[, needs_flip]
    }

    rowCors <- rowCor_abs(A_coupled, B_coupled, zero_val = 0)

    all_pairs <- combn(J, 2)
    n_pairs_total <- ncol(all_pairs)
    all_RIC <- matrix(NA_real_, nrow = N, ncol = iterations)
    if (progress) pb <- txtProgressBar(min = 0, max = iterations, style = 3)

    # Precompute cross-factor pair pool for coupled mode
    cf_cross_pairs <- NULL
    if (cross_factor_baseline && !is.null(cf_factor) && length(unique(cf_factor)) >= 2) {
      # Filter all_pairs to only cross-factor
      pf1 <- cf_factor[all_pairs[1, ]]; pf2 <- cf_factor[all_pairs[2, ]]
      cf_mask <- pf1 != pf2
      cf_cross_pairs <- all_pairs[, cf_mask, drop = FALSE]
    }

    for (i in seq_len(iterations)) {
      if (progress) setTxtProgressBar(pb, i)
      if (!is.null(cf_cross_pairs) && ncol(cf_cross_pairs) >= k) {
        sampled <- sample(ncol(cf_cross_pairs), k)
        idx1 <- cf_cross_pairs[1, sampled]; idx2 <- cf_cross_pairs[2, sampled]
      } else {
        sampled <- sample(n_pairs_total, k)
        idx1 <- all_pairs[1, sampled]; idx2 <- all_pairs[2, sampled]
      }
      A_rand <- mat[, idx1, drop = FALSE]; B_rand <- mat[, idx2, drop = FALSE]
      if (align_signs) {
        ri <- pmax(idx1, idx2); ci <- pmin(idx1, idx2)
        rand_signs <- signMat[cbind(ri, ci)]; rand_signs[is.na(rand_signs)] <- 1
        rand_neg <- which(rand_signs < 0)
        if (length(rand_neg) > 0) {
          rb_idx <- idx2[rand_neg]
          B_rand[, rand_neg] <- rep(item_max[rb_idx] + item_min[rb_idx], each = N) - B_rand[, rand_neg]
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

  z_score_raw <- ifelse(rand_sds > 0,
                        (rowCors - rand_means) / rand_sds,
                        0)

  # --- Optional variance penalty (detects straight-lining / acquiescence) ---
  # Rationale: a respondent whose within-person SD is near zero (e.g. all 7s)
  # produces rowCor = 0 for any pair, so z_score collapses to 0 regardless of
  # actual carelessness. Subtracting alpha * exp(-sd_resp / beta) pushes such
  # respondents to strongly negative z, increasing detection confidence.
  # Defaults (alpha=3, beta=0.5) tuned on simulated data (2026-04-21):
  # +0.053 MCC on pure_straight, +0.030 on acquiescent, 0.000 on random/longstring.
  if (isTRUE(variance_penalty)) {
    sd_resp <- apply(as.matrix(data), 1, sd, na.rm = TRUE)
    sd_resp[is.na(sd_resp)] <- 0
    penalty <- vp_alpha * exp(-sd_resp / max(vp_beta, 1e-6))
    z_score <- z_score_raw - penalty
  } else {
    z_score <- z_score_raw
  }

  data.frame(
    result = corComparedIndex,     # legacy percentile score
    z_score = z_score,             # z-score (with variance_penalty applied if requested)
    z_score_raw = z_score_raw,     # original z-score without penalty (for reference)
    indCors = rowCors,             # raw coupled/weighted coherence score
    rand_mean = rand_means,        # mean of random iterations
    rand_sd = rand_sds,            # sd of random iterations
    flagged = z_score <= z_threshold,  # z-score flagging (primary)
    z_threshold_used = z_threshold,    # threshold used
    nFactors_detected = nFactors_detected,  # from parallel analysis
    mode_used = mode_used,         # "efa_d", "coupled", or "weighted"
    n_pairs = k,                   # number of item pairs used
    variance_penalty_used = isTRUE(variance_penalty),
    rescale_used = rescale         # per-item rescaling applied
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
                       cross_factor_baseline = FALSE,
                       rescale = c("proportion", "none", "minmax", "zscore"),
                       progress = FALSE) {

  rescale <- match.arg(rescale)

  # Delegate to ReReReRe with mode="efa_d"
  ReReReRe(data,
           corProp = 0.03,
           z_threshold = z_threshold,
           iterations = iterations,
           align_signs = align_signs,
           auto_z = auto_z,
           mode = "efa_d",
           min_r = min_r,
           cross_factor_baseline = cross_factor_baseline,
           rescale = rescale,
           progress = progress)
}


# ======================================================================
# ReReReRe_F variants (experimental)
# ======================================================================

#' Variant A: Per-factor z, flagging by proportion of factors below threshold
#'
#' @param factor_z_threshold z-score threshold per factor (below = "failed")
#' @param flag_thresholds proportion of failed factors to flag (vector for post-hoc sweep)
#' @return data.frame with prop_low (proportion of factors with z < factor_z_threshold),
#'         mean_z, var_z, and flag columns for each flag_threshold
ReReReRe_F_proplow <- function(data, iterations = 50, align_signs = TRUE,
                                factor_z_threshold = 1.0,
                                cross_factor_baseline = FALSE) {
  data <- data[, sapply(data, is.numeric), drop = FALSE]
  mat <- as.matrix(data); N <- nrow(mat)

  pf <- .compute_per_factor_z(mat, iterations = iterations, align_signs = align_signs,
                                cross_factor_baseline = cross_factor_baseline)
  if (is.null(pf)) return(NULL)

  fz <- pf$factor_z_matrix
  prop_low <- rowMeans(fz < factor_z_threshold, na.rm = TRUE)
  mean_z <- rowMeans(fz, na.rm = TRUE)
  var_z <- apply(fz, 1, var, na.rm = TRUE)

  data.frame(prop_low = prop_low, mean_z = mean_z, var_z = var_z,
             nF_detected = pf$nF_detected, n_factors_used = pf$n_factors_used)
}


#' Variant B: Per-factor z, using mean and variance as features
#'
#' Careless respondents have LOW mean_z and potentially HIGH or LOW var_z.
#' Combined score: mean_z - lambda * sqrt(var_z)
ReReReRe_F_meanvar <- function(data, iterations = 50, align_signs = TRUE,
                                lambda = 1.0,
                                cross_factor_baseline = FALSE) {
  data <- data[, sapply(data, is.numeric), drop = FALSE]
  mat <- as.matrix(data); N <- nrow(mat)

  pf <- .compute_per_factor_z(mat, iterations = iterations, align_signs = align_signs,
                                cross_factor_baseline = cross_factor_baseline)
  if (is.null(pf)) return(NULL)

  fz <- pf$factor_z_matrix
  mean_z <- rowMeans(fz, na.rm = TRUE)
  var_z <- apply(fz, 1, var, na.rm = TRUE)
  combined <- mean_z - lambda * sqrt(pmax(var_z, 0))

  data.frame(mean_z = mean_z, var_z = var_z, combined = combined,
             nF_detected = pf$nF_detected, n_factors_used = pf$n_factors_used)
}


#' Variant C: Iterative EFA — run once, clean, re-EFA, re-score
#'
#' Round 1: EFA on full sample → z-scores → flag at lenient threshold
#' Round 2: EFA on clean sample → re-score ALL respondents with clean structure
ReReReRe_F_iterative <- function(data, iterations = 50, align_signs = TRUE,
                                  initial_z_threshold = 2.0,
                                  cross_factor_baseline = FALSE) {
  data <- data[, sapply(data, is.numeric), drop = FALSE]
  mat <- as.matrix(data); N <- nrow(mat)

  # Round 1: full sample EFA
  rr1 <- ReReReRe(data, corProp = 0.03, z_threshold = initial_z_threshold,
                   iterations = iterations, align_signs = align_signs,
                   mode = "efa_d", cross_factor_baseline = cross_factor_baseline)

  flagged_r1 <- rr1$z_score <= initial_z_threshold
  n_flagged <- sum(flagged_r1)

  # Guard: don't remove too many (keep at least 50% of sample)
  if (n_flagged >= N * 0.5 || n_flagged == 0) {
    # No improvement possible, return round 1 results
    return(rr1)
  }

  # Round 2: EFA on clean subsample
  mat_clean <- mat[!flagged_r1, , drop = FALSE]

  # Run EFA on clean data
  pa2 <- tryCatch({
    suppressMessages(suppressWarnings(
      fa.parallel(mat_clean, fa = "fa", plot = FALSE, n.iter = 20)
    ))
  }, error = function(e) NULL)

  nF2 <- if (!is.null(pa2) && !is.null(pa2$nfact) && pa2$nfact >= 1) pa2$nfact
         else max(1, sum(eigen(cor(mat_clean, use="pairwise.complete.obs"),
                                symmetric=TRUE, only.values=TRUE)$values > 1))
  nF2 <- max(1, min(nF2, floor(ncol(mat_clean) / 2)))

  efa2 <- tryCatch({
    suppressWarnings(fa(mat_clean, nfactors = nF2, rotate = "oblimin", fm = "minres",
                        scores = "none", warnings = FALSE))
  }, error = function(e) NULL)

  if (is.null(efa2)) return(rr1)  # fallback to round 1

  # Use clean EFA to score ALL respondents
  loadings2 <- as.matrix(efa2$loadings[])
  primary_factor2 <- apply(abs(loadings2), 1, which.max)

  # Re-run ReReReRe on full data using the clean factor structure
  # We can't directly pass external_efa to ReReReRe, so we use .compute_per_factor_z
  # with external_efa to get per-factor z-scores, then aggregate as mean_z

  # But actually we want the global EFA-D score with clean factor pairs
  # Rebuild pairs from clean EFA structure applied to full correlation matrix
  raw_cor_full <- cor(mat, use = "pairwise.complete.obs")
  item_max <- apply(mat, 2, max, na.rm = TRUE)
  signMat <- sign(raw_cor_full)
  signMat[upper.tri(signMat, diag = TRUE)] <- NA

  # Generate within-factor pairs from clean EFA
  nF2_actual <- ncol(loadings2)
  pair_list <- list()
  for (f in seq_len(nF2_actual)) {
    items_f <- which(primary_factor2 == f)
    if (length(items_f) < 2) next
    combos <- combn(items_f, 2)
    for (ci in seq_len(ncol(combos)))
      pair_list[[length(pair_list) + 1]] <- c(combos[1, ci], combos[2, ci])
  }

  if (length(pair_list) == 0) return(rr1)

  pair_mat <- do.call(rbind, pair_list)
  idx_A <- pair_mat[, 1]; idx_B <- pair_mat[, 2]
  pair_r <- numeric(length(idx_A))
  for (j in seq_along(idx_A)) pair_r[j] <- raw_cor_full[idx_A[j], idx_B[j]]
  pair_abs_r <- abs(pair_r); pair_sign <- sign(pair_r)
  keep <- !is.na(pair_abs_r)
  idx_A <- idx_A[keep]; idx_B <- idx_B[keep]
  pair_abs_r <- pair_abs_r[keep]; pair_sign <- pair_sign[keep]
  k <- length(idx_A)
  weights <- pair_abs_r; weights[weights < 1e-6] <- 1e-6

  A_coupled <- mat[, idx_A, drop = FALSE]; B_coupled <- mat[, idx_B, drop = FALSE]
  if (align_signs) {
    nf <- which(pair_sign < 0)
    if (length(nf) > 0) B_coupled[, nf] <- rep(item_max[idx_B[nf]] + 1, each = N) - B_coupled[, nf]
  }

  rowCors <- rowCor_weighted(A_coupled, B_coupled, weights)

  # Permutation baseline
  J <- ncol(mat)
  all_RIC <- matrix(NA_real_, nrow = N, ncol = iterations)
  for (i in seq_len(iterations)) {
    if (cross_factor_baseline && length(unique(primary_factor2)) >= 2) {
      cf <- .sample_cross_factor_pairs(primary_factor2, k, J)
      ri1 <- cf$idx1; ri2 <- cf$idx2
    } else {
      ri1 <- sample(J, k, replace = TRUE); ri2 <- sample(J, k, replace = TRUE)
      same <- ri1 == ri2
      while (any(same)) { ri2[same] <- sample(J, sum(same), replace = TRUE); same <- ri1 == ri2 }
    }
    Ar <- mat[, ri1, drop = FALSE]; Br <- mat[, ri2, drop = FALSE]
    if (align_signs) {
      ri <- pmax(ri1, ri2); ci <- pmin(ri1, ri2)
      rs <- signMat[cbind(ri, ci)]; rs[is.na(rs)] <- 1
      rn <- which(rs < 0)
      if (length(rn) > 0) Br[, rn] <- rep(item_max[ri2[rn]] + 1, each = N) - Br[, rn]
    }
    all_RIC[, i] <- rowCor_weighted(Ar, Br, weights)
  }

  rand_means <- rowMeans(all_RIC, na.rm = TRUE)
  rand_sds <- apply(all_RIC, 1, sd, na.rm = TRUE)
  z_score <- ifelse(rand_sds > 0, (rowCors - rand_means) / rand_sds, 0)

  data.frame(
    result = rowMeans(all_RIC < rowCors, na.rm = TRUE),
    z_score = z_score,
    indCors = rowCors,
    rand_mean = rand_means,
    rand_sd = rand_sds,
    flagged = z_score <= rr1$z_threshold_used[1],
    z_threshold_used = rr1$z_threshold_used[1],
    nFactors_detected = nF2,
    mode_used = "efa_d_iterative",
    n_pairs = k
  )
}

## Official name (2026-07-05): the method is now "ReReRe" (3 Re, matching
## re-re.re). ReReRe() is the canonical name; ReReReRe() is kept as a
## backward-compatible alias so existing scripts keep working.
ReReRe   <- ReReReRe
ReReRe_F <- ReReReRe_F
