# ReReReRe_SplitHalf.R
# Split-half z-score to detect partial carelessness.
#
# Rationale:
#   - A fully-attentive respondent is coherent across ALL pair subsets.
#   - A fully-careless respondent is incoherent across ALL pair subsets.
#   - A PARTIALLY careless respondent (e.g. 60-70% corrupted) has some
#     preserved pairs and some corrupted pairs. The mean z across all pairs
#     is moderate, but SOME random split concentrates corrupted pairs in one
#     half, producing a very low z for that half.
#
# score_split_half() returns the MIN z-score across B random splits x 2 halves.
# Partial carelessness shows up as low min-z even when mean-z is moderate.

score_split_half <- function(data,
                             corProp = 0.03,
                             iterations = 100,
                             align_signs = TRUE,
                             min_pairs = 15,
                             n_splits = 5,
                             aggregation = "min",   # "min", "mean", "median"
                             seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  data <- data[, sapply(data, is.numeric), drop = FALSE]
  N <- nrow(data); J <- ncol(data)
  mat <- as.matrix(data)
  item_max <- apply(mat, 2, max, na.rm = TRUE)

  # --- Correlation matrix ---
  rawCorMat <- cor(mat, use = "pairwise.complete.obs")
  corMat <- abs(rawCorMat)
  corMat[upper.tri(corMat, diag = TRUE)] <- NA
  rawCorMat[upper.tri(rawCorMat, diag = TRUE)] <- NA
  signMat <- sign(rawCorMat)

  # --- Find top-k% coupled pairs ---
  corThreshold <- quantile(corMat, 1 - corProp, na.rm = TRUE)
  couples <- which(corMat >= corThreshold, arr.ind = TRUE)
  k <- nrow(couples)
  if (k < min_pairs) {
    all_cors <- corMat[lower.tri(corMat)]
    all_cors <- all_cors[!is.na(all_cors)]
    use_k <- min(min_pairs, length(all_cors))
    corThreshold <- sort(all_cors, decreasing = TRUE)[use_k]
    couples <- which(corMat >= corThreshold, arr.ind = TRUE)
    k <- nrow(couples)
  }
  if (k < 4) stop("Not enough coupled pairs for split-half (need >=4)")

  # --- Align signs on coupled pairs ---
  coupled_signs <- sign(rawCorMat[couples])
  A_coupled <- mat[, couples[, 1], drop = FALSE]
  B_coupled <- mat[, couples[, 2], drop = FALSE]
  n_negative <- sum(coupled_signs < 0)
  if (align_signs && n_negative > 0) {
    needs_flip <- which(coupled_signs < 0)
    flip_max <- item_max[couples[needs_flip, 2]]
    B_coupled[, needs_flip] <- rep(flip_max + 1, each = N) - B_coupled[, needs_flip]
  }

  # --- rowCor helper (same as main) ---
  rowCor_abs_local <- function(A, B) {
    A_mean <- rowMeans(A, na.rm = TRUE); B_mean <- rowMeans(B, na.rm = TRUE)
    A_c <- A - A_mean; B_c <- B - B_mean
    num <- rowSums(A_c * B_c, na.rm = TRUE)
    den <- sqrt(rowSums(A_c^2, na.rm = TRUE) * rowSums(B_c^2, na.rm = TRUE))
    r <- abs(num / den); r[is.nan(r) | is.na(r)] <- 0; r
  }

  all_pairs <- combn(J, 2)
  n_pairs_total <- ncol(all_pairs)
  k_half <- k %/% 2

  # Iterations per half — reduce so total random-pair calls ~= standard RR cost
  iter_per_half <- max(30, iterations %/% max(1, n_splits))

  # Each split generates 2 halves -> z_scores_all is N x (n_splits * 2)
  z_scores_all <- matrix(NA_real_, nrow = N, ncol = n_splits * 2)
  half_col <- 0

  for (s in seq_len(n_splits)) {
    perm <- sample(k)
    h_idx_list <- list(perm[1:k_half], perm[(k_half + 1):(2 * k_half)])

    for (h in 1:2) {
      half_col <- half_col + 1
      hidx <- h_idx_list[[h]]
      kh <- length(hidx)

      # Observed on half
      rowCors_h <- rowCor_abs_local(
        A_coupled[, hidx, drop = FALSE],
        B_coupled[, hidx, drop = FALSE]
      )

      # Random baseline for half
      all_RIC_h <- matrix(NA_real_, nrow = N, ncol = iter_per_half)
      for (i in seq_len(iter_per_half)) {
        sampled <- sample(n_pairs_total, kh)
        idx1 <- all_pairs[1, sampled]; idx2 <- all_pairs[2, sampled]
        A_r <- mat[, idx1, drop = FALSE]; B_r <- mat[, idx2, drop = FALSE]
        if (align_signs) {
          ri <- pmax(idx1, idx2); ci <- pmin(idx1, idx2)
          rs <- signMat[cbind(ri, ci)]; rs[is.na(rs)] <- 1
          rneg <- which(rs < 0)
          if (length(rneg) > 0)
            B_r[, rneg] <- rep(item_max[idx2[rneg]] + 1, each = N) - B_r[, rneg]
        }
        all_RIC_h[, i] <- rowCor_abs_local(A_r, B_r)
      }

      rm_h <- rowMeans(all_RIC_h, na.rm = TRUE)
      rs_h <- apply(all_RIC_h, 1, sd, na.rm = TRUE)
      z_scores_all[, half_col] <- ifelse(rs_h > 0,
                                         (rowCors_h - rm_h) / rs_h, 0)
    }
  }

  # Aggregation
  z_final <- switch(aggregation,
    "min"    = apply(z_scores_all, 1, min,    na.rm = TRUE),
    "mean"   = apply(z_scores_all, 1, mean,   na.rm = TRUE),
    "median" = apply(z_scores_all, 1, median, na.rm = TRUE),
    stop("aggregation must be 'min', 'mean', or 'median'")
  )

  list(
    z_score       = z_final,
    z_min         = apply(z_scores_all, 1, min,  na.rm = TRUE),
    z_mean        = apply(z_scores_all, 1, mean, na.rm = TRUE),
    z_max         = apply(z_scores_all, 1, max,  na.rm = TRUE),
    z_per_half    = z_scores_all,
    k = k, n_splits = n_splits, iter_per_half = iter_per_half
  )
}
