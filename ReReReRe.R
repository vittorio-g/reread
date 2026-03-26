#### ReReReRe: Permutation-based careless respondent detection ####
#
# Returns per-respondent:
#   result    - legacy percentile: P(random < coupled)
#   z_score   - (coupled_cor - mean_random) / sd_random
#   indCors   - raw coupled absolute correlation
#   rand_mean - mean of random iteration correlations
#   rand_sd   - sd of random iteration correlations
#   flagged   - legacy binary flag (result <= cutOff)
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
# ===

library(dplyr)
library(magrittr)
library(lavaan)
library(psych)
library(mice)

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

ReReReRe <- function(data, #any dataset with questionnaire data
                    corProp=0.05, # the proportion of highest correlations to use
                    cutOff=0.99, # the severity of the evaluation
                    iterations=100,
                    min_pairs=15, # minimum number of item pairs to use
                    align_signs=TRUE, # align reverse-coded items using sample correlation signs
                    progress = F){

  #keep only numeric values
  data <- data[, sapply(data, is.numeric), drop = FALSE]
  N <- nrow(data)
  J <- ncol(data)

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

  #getting the threshold correlation based on the proportion of highest correlation we decided to include.
  corThreshold <- quantile(corMat, 1 - corProp, na.rm = TRUE)

  #getting the row and col indices with values higher than threshold
  couples <- which(corMat >= corThreshold, arr.ind = TRUE)
  k <- nrow(couples)

  # Enforce minimum number of pairs for stable individual-level correlations.
  # With fewer than ~10-15 pairs, cor() across k points is extremely noisy
  # (k=2 always gives r=±1, k=3-5 is barely meaningful).
  # If corProp yields too few pairs, we lower the threshold to include more.
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
  # Get the sign of each coupled pair's sample-level correlation.
  # Used to flip B columns so all pairs contribute the same direction
  # at the individual level, preventing cancellation in rowCor_abs.
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

  # Precompute the full sign matrix for random pair reverse-coding
  # (only the lower triangle is populated, matching rawCorMat)
  if (align_signs) {
    signMat <- sign(rawCorMat)
  }

  #getting the correlation for each individual (VECTORIZED)
  # Extract N x k matrices for the coupled pairs
  A_coupled <- mat[, couples[, 1], drop = FALSE]
  B_coupled <- mat[, couples[, 2], drop = FALSE]

  # Apply sign alignment: reverse-code B columns where sample correlation is negative
  if (align_signs && n_negative > 0) {
    # For each negatively-correlated coupled pair, reverse-code the B item
    # using (max+1-x). This keeps values in the original Likert range,
    # preventing centering distortion that sign flip (*-1) would cause.
    needs_flip <- which(coupled_signs < 0)
    # Get the max of each B item that needs flipping
    flip_max <- item_max[couples[needs_flip, 2]]
    # Reverse code: (max+1) - x, applied column-wise
    B_coupled[, needs_flip] <- rep(flip_max + 1, each = N) - B_coupled[, needs_flip]
  }

  # zero_val=0: for zero-variance rows (longstring), return 0 so they are easily spotted
  rowCors <- rowCor_abs(A_coupled, B_coupled, zero_val = 0)

  #### RANDOM PERMUTATION ITERATIONS ####

  #all possible column pair indices
  all_pairs <- combn(J, 2) # 2 x C(J,2) matrix
  n_pairs <- ncol(all_pairs)

  all_RIC <- matrix(NA_real_, nrow = N, ncol = iterations)

  if(progress == TRUE){
    pb <- txtProgressBar(min = 0, max = iterations, style = 3)
  }

  #computing random correlations
  for (i in seq_len(iterations)){

    if(progress == TRUE) setTxtProgressBar(pb, i)

    #Sampling k random column pairs (same number as coupled pairs)
    sampled <- sample(n_pairs, k)
    idx1 <- all_pairs[1, sampled]
    idx2 <- all_pairs[2, sampled]

    #computing the random individual correlation (VECTORIZED)
    A_rand <- mat[, idx1, drop = FALSE]
    B_rand <- mat[, idx2, drop = FALSE]

    # Apply sign alignment to random pairs too (reverse coding, same logic)
    if (align_signs) {
      # Look up sign from signMat. signMat has lower triangle populated
      # (row > col), so ensure correct indexing.
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

    # zero_val=0: undefined correlation → 0 (no evidence of relationship).
    all_RIC[, i] <- rowCor_abs(A_rand, B_rand, zero_val = 0)
  }

  if(progress == TRUE) close(pb)

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

  data.frame(
    result = corComparedIndex,     # legacy percentile score
    z_score = z_score,             # z-score: SDs above random baseline
    indCors = rowCors,             # raw coupled correlation
    rand_mean = rand_means,        # mean of random iterations
    rand_sd = rand_sds,            # sd of random iterations
    flagged = corComparedIndex <= cutOff  # legacy flagging
  )
}
