# Iterative coupled ReReReRe: re-selects pairs using the correlation matrix
# computed on the "cleaner" subset (after removing first-pass suspects).
#
# Idea:
#   Pass 1: standard RR on full data -> z1.
#   Trim: remove bottom-X% (most suspect).
#   Pass 2: compute cor matrix on clean subset -> top-k% pairs from clean cor.
#           Score ALL respondents using these pairs.
#
# Useful when contamination by careless corrupts the sample correlation matrix
# (high careless rate, especially with random respondents).

score_iter_coupled <- function(data,
                                corProp = 0.03,
                                iterations = 100,
                                align_signs = TRUE,
                                trim_pct = 0.20,
                                min_pairs = 15,
                                first_z_threshold = "quantile") {
  data <- data[, sapply(data, is.numeric), drop = FALSE]
  N <- nrow(data); J <- ncol(data); mat <- as.matrix(data)

  # ---- Pass 1: standard RR ----
  rr1 <- ReReReRe(data, corProp = corProp, iterations = iterations,
                  align_signs = align_signs, mode = "coupled",
                  min_pairs = min_pairs, variance_penalty = FALSE)
  z1 <- rr1$z_score

  # ---- Trim: respondents at the bottom trim_pct of z1 flagged as suspects ----
  thr <- quantile(z1, trim_pct, na.rm = TRUE)
  clean_mask <- z1 > thr
  if (sum(clean_mask) < 50 || sum(clean_mask) == N) {
    # Degenerate trim, return pass 1
    return(list(z_score = z1, pass = 1))
  }

  mat_clean <- mat[clean_mask, , drop = FALSE]

  # ---- Pass 2: select pairs from CLEAN correlation matrix, score ALL ----
  rawCorMat <- cor(mat_clean, use = "pairwise.complete.obs")
  corMat <- abs(rawCorMat)
  corMat[upper.tri(corMat, diag = TRUE)] <- NA
  rawCorMat[upper.tri(rawCorMat, diag = TRUE)] <- NA
  signMat <- sign(rawCorMat)
  item_max <- apply(mat, 2, max, na.rm = TRUE)

  corThreshold <- quantile(corMat, 1 - corProp, na.rm = TRUE)
  couples <- which(corMat >= corThreshold, arr.ind = TRUE)
  k <- nrow(couples)
  if (k < min_pairs) {
    all_cors <- corMat[lower.tri(corMat)]; all_cors <- all_cors[!is.na(all_cors)]
    use_k <- min(min_pairs, length(all_cors))
    corThreshold <- sort(all_cors, decreasing = TRUE)[use_k]
    couples <- which(corMat >= corThreshold, arr.ind = TRUE)
    k <- nrow(couples)
  }
  if (k == 0) return(list(z_score = z1, pass = 1))

  coupled_signs <- sign(rawCorMat[couples])
  A_c <- mat[, couples[, 1], drop = FALSE]
  B_c <- mat[, couples[, 2], drop = FALSE]
  if (align_signs) {
    nn <- which(coupled_signs < 0)
    if (length(nn) > 0)
      B_c[, nn] <- rep(item_max[couples[nn, 2]] + 1, each = N) - B_c[, nn]
  }

  # rowCor_abs inline
  rowCor_abs_local <- function(A, B) {
    Am <- rowMeans(A, na.rm = TRUE); Bm <- rowMeans(B, na.rm = TRUE)
    Ac <- A - Am; Bc <- B - Bm
    num <- rowSums(Ac * Bc, na.rm = TRUE)
    den <- sqrt(rowSums(Ac^2, na.rm = TRUE) * rowSums(Bc^2, na.rm = TRUE))
    r <- abs(num / den); r[is.nan(r) | is.na(r)] <- 0; r
  }

  rowCors <- rowCor_abs_local(A_c, B_c)

  # Permutation baseline using ALL items (as per standard RR)
  all_pairs <- combn(J, 2)
  n_pairs_total <- ncol(all_pairs)
  all_RIC <- matrix(NA_real_, nrow = N, ncol = iterations)
  for (i in seq_len(iterations)) {
    sampled <- sample(n_pairs_total, k)
    idx1 <- all_pairs[1, sampled]; idx2 <- all_pairs[2, sampled]
    A_r <- mat[, idx1, drop = FALSE]; B_r <- mat[, idx2, drop = FALSE]
    if (align_signs) {
      ri <- pmax(idx1, idx2); ci <- pmin(idx1, idx2)
      rs <- signMat[cbind(ri, ci)]; rs[is.na(rs)] <- 1
      rn <- which(rs < 0)
      if (length(rn) > 0)
        B_r[, rn] <- rep(item_max[idx2[rn]] + 1, each = N) - B_r[, rn]
    }
    all_RIC[, i] <- rowCor_abs_local(A_r, B_r)
  }
  rm_v <- rowMeans(all_RIC, na.rm = TRUE)
  rs_v <- apply(all_RIC, 1, sd, na.rm = TRUE)
  z2 <- ifelse(rs_v > 0, (rowCors - rm_v) / rs_v, 0)

  list(z_score = z2, pass = 2, z_pass1 = z1, trim_threshold = thr, k = k)
}
