###############################################################################
# Test: EFA Option D — EFA for pair SELECTION + observed |r| for WEIGHTING
#
# Difference from Option A:
#   A: within-factor pairs, weight = |lambda_i * lambda_j| (loading product)
#   D: within-factor pairs, weight = |r_observed|           (empirical cor)
#
# The idea: EFA tells us WHICH pairs matter (within-factor = theoretically
# justified), but the observed |r| tells us HOW MUCH they matter (empirical,
# no model assumptions on weight magnitude). If the EFA assigns items wrong,
# the observed |r| compensates (low |r| for cross-factor pairs = low weight).
#
# Also tests: sign alignment from observed r (not loading signs).
###############################################################################

library(lavaan)
library(psych)
library(dplyr)

source("Synthetic_Good_Responses_2.R")
source("Careless_machine_2.R")
source("ReReReRe.R")

# ============================================================================
# HELPERS
# ============================================================================

evaluate_mcc <- function(predicted, actual) {
  tp <- sum(predicted == 1 & actual == 1)
  tn <- sum(predicted == 0 & actual == 0)
  fp <- sum(predicted == 1 & actual == 0)
  fn <- sum(predicted == 0 & actual == 1)
  denom <- as.double(tp + fp) * as.double(tp + fn) *
           as.double(tn + fp) * as.double(tn + fn)
  if (denom == 0) return(0)
  return((as.double(tp) * tn - as.double(fp) * fn) / sqrt(denom))
}

find_best_z <- function(z_scores, labels) {
  best_mcc <- -1; best_z <- NA
  for (z in seq(0.1, 3.0, 0.1)) {
    pred <- as.integer(z_scores <= z)
    mcc <- evaluate_mcc(pred, labels)
    if (mcc > best_mcc) { best_mcc <- mcc; best_z <- z }
  }
  return(list(z = best_z, mcc = best_mcc))
}

# ============================================================================
# OPTION D: EFA for pair selection + |r_observed| for weighting
# ============================================================================

ReReReRe_efa_D <- function(data, iterations = 50, align_signs = TRUE,
                            nf_override = NULL, min_r = 0.0) {
  data <- data[, sapply(data, is.numeric), drop = FALSE]
  mat <- as.matrix(data)
  N <- nrow(mat)
  p <- ncol(mat)

  item_max <- apply(mat, 2, max, na.rm = TRUE)

  # --- Step 1: Estimate nF ---
  nF_est <- nf_override
  if (is.null(nF_est)) {
    pa <- tryCatch({
      suppressMessages(suppressWarnings(
        fa.parallel(mat, fa = "fa", plot = FALSE, n.iter = 20)
      ))
    }, error = function(e) NULL)

    if (!is.null(pa) && !is.null(pa$nfact) && pa$nfact >= 1) {
      nF_est <- pa$nfact
    } else {
      ev <- eigen(cor(mat, use = "pairwise.complete.obs"),
                  symmetric = TRUE, only.values = TRUE)$values
      nF_est <- max(1, sum(ev > 1))
    }
  }
  nF_est <- max(1, min(nF_est, floor(p / 2)))

  # --- Step 2: Run EFA ---
  efa_result <- tryCatch({
    suppressWarnings(
      fa(mat, nfactors = nF_est, rotate = "oblimin", fm = "minres",
         scores = "none", warnings = FALSE)
    )
  }, error = function(e) {
    tryCatch({
      suppressWarnings(
        fa(mat, nfactors = max(1, nF_est - 1), rotate = "oblimin", fm = "minres",
           scores = "none", warnings = FALSE)
      )
    }, error = function(e2) NULL)
  })

  if (is.null(efa_result)) {
    warning("EFA failed.")
    return(data.frame(z_score = rep(NA_real_, N), indCors = rep(NA_real_, N),
                      rand_mean = rep(NA_real_, N), rand_sd = rep(NA_real_, N),
                      nF_detected = nF_est, n_pairs = NA_integer_))
  }

  # --- Step 3: Assign items to primary factor ---
  loadings_mat <- as.matrix(efa_result$loadings[])
  nF_actual <- ncol(loadings_mat)
  primary_factor <- apply(abs(loadings_mat), 1, which.max)

  # --- Step 4: Generate within-factor pairs ---
  pair_list <- list()
  for (f in seq_len(nF_actual)) {
    items_f <- which(primary_factor == f)
    if (length(items_f) < 2) next
    combos <- combn(items_f, 2)
    for (ci in seq_len(ncol(combos))) {
      pair_list[[length(pair_list) + 1]] <- c(combos[1, ci], combos[2, ci])
    }
  }

  if (length(pair_list) == 0) {
    warning("No within-factor pairs.")
    return(data.frame(z_score = rep(NA_real_, N), indCors = rep(NA_real_, N),
                      rand_mean = rep(NA_real_, N), rand_sd = rep(NA_real_, N),
                      nF_detected = nF_est, n_pairs = 0L))
  }

  pair_mat <- do.call(rbind, pair_list)
  idx_A <- pair_mat[, 1]
  idx_B <- pair_mat[, 2]

  # --- Step 5: OPTION D DIFFERENCE — weight by OBSERVED |r|, not loading product ---
  raw_cor_mat <- cor(mat, use = "pairwise.complete.obs")

  # Get observed correlation for each within-factor pair
  pair_r <- numeric(length(idx_A))
  for (j in seq_along(idx_A)) {
    pair_r[j] <- raw_cor_mat[idx_A[j], idx_B[j]]
  }
  pair_abs_r <- abs(pair_r)
  pair_sign <- sign(pair_r)

  # Filter by min_r
  keep <- pair_abs_r >= min_r & !is.na(pair_abs_r)
  if (sum(keep) == 0) {
    warning("No pairs above min_r threshold.")
    return(data.frame(z_score = rep(NA_real_, N), indCors = rep(NA_real_, N),
                      rand_mean = rep(NA_real_, N), rand_sd = rep(NA_real_, N),
                      nF_detected = nF_est, n_pairs = 0L))
  }

  idx_A <- idx_A[keep]
  idx_B <- idx_B[keep]
  pair_r <- pair_r[keep]
  pair_abs_r <- pair_abs_r[keep]
  pair_sign <- pair_sign[keep]

  k <- length(idx_A)
  weights <- pair_abs_r
  weights[weights < 1e-6] <- 1e-6

  # Build A and B matrices
  A_coupled <- mat[, idx_A, drop = FALSE]
  B_coupled <- mat[, idx_B, drop = FALSE]

  # Sign alignment using OBSERVED correlation sign (not loading sign)
  if (align_signs) {
    needs_flip <- which(pair_sign < 0)
    if (length(needs_flip) > 0) {
      flip_max <- item_max[idx_B[needs_flip]]
      B_coupled[, needs_flip] <- rep(flip_max + 1, each = N) - B_coupled[, needs_flip]
    }
  }

  # Weighted coherence score
  indCors <- rowCor_weighted(A_coupled, B_coupled, weights)

  # --- Permutation baseline ---
  signMat <- sign(raw_cor_mat)
  signMat[upper.tri(signMat, diag = TRUE)] <- NA

  rand_scores <- matrix(NA_real_, N, iterations)
  for (iter in seq_len(iterations)) {
    rand_idx1 <- sample(p, k, replace = TRUE)
    rand_idx2 <- sample(p, k, replace = TRUE)
    same <- rand_idx1 == rand_idx2
    while (any(same)) {
      rand_idx2[same] <- sample(p, sum(same), replace = TRUE)
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

    rand_scores[, iter] <- rowCor_weighted(A_rand, B_rand, weights)
  }

  rand_mean <- rowMeans(rand_scores, na.rm = TRUE)
  rand_sd <- apply(rand_scores, 1, sd, na.rm = TRUE)
  rand_sd[rand_sd == 0] <- 1e-10
  z_score <- (indCors - rand_mean) / rand_sd

  data.frame(z_score = z_score, indCors = indCors,
             rand_mean = rand_mean, rand_sd = rand_sd,
             nF_detected = nF_est, n_pairs = k)
}

# ============================================================================
# OPTION A (for comparison): EFA selection + loading weights
# ============================================================================

ReReReRe_efa_A <- function(data, iterations = 50, align_signs = TRUE,
                            nf_override = NULL) {
  data <- data[, sapply(data, is.numeric), drop = FALSE]
  mat <- as.matrix(data)
  N <- nrow(mat)
  p <- ncol(mat)
  item_max <- apply(mat, 2, max, na.rm = TRUE)

  nF_est <- nf_override
  if (is.null(nF_est)) {
    pa <- tryCatch({
      suppressMessages(suppressWarnings(
        fa.parallel(mat, fa = "fa", plot = FALSE, n.iter = 20)
      ))
    }, error = function(e) NULL)
    if (!is.null(pa) && !is.null(pa$nfact) && pa$nfact >= 1) {
      nF_est <- pa$nfact
    } else {
      ev <- eigen(cor(mat, use = "pairwise.complete.obs"),
                  symmetric = TRUE, only.values = TRUE)$values
      nF_est <- max(1, sum(ev > 1))
    }
  }
  nF_est <- max(1, min(nF_est, floor(p / 2)))

  efa_result <- tryCatch({
    suppressWarnings(fa(mat, nfactors = nF_est, rotate = "oblimin", fm = "minres",
                        scores = "none", warnings = FALSE))
  }, error = function(e) {
    tryCatch({
      suppressWarnings(fa(mat, nfactors = max(1, nF_est - 1), rotate = "oblimin",
                          fm = "minres", scores = "none", warnings = FALSE))
    }, error = function(e2) NULL)
  })

  if (is.null(efa_result)) {
    return(data.frame(z_score = rep(NA_real_, N), indCors = rep(NA_real_, N),
                      rand_mean = rep(NA_real_, N), rand_sd = rep(NA_real_, N),
                      nF_detected = nF_est, n_pairs = NA_integer_))
  }

  loadings_mat <- as.matrix(efa_result$loadings[])
  nF_actual <- ncol(loadings_mat)
  primary_factor <- apply(abs(loadings_mat), 1, which.max)
  primary_loading <- numeric(p)
  for (i in seq_len(p)) primary_loading[i] <- loadings_mat[i, primary_factor[i]]

  pair_list <- list()
  for (f in seq_len(nF_actual)) {
    items_f <- which(primary_factor == f)
    if (length(items_f) < 2) next
    combos <- combn(items_f, 2)
    for (ci in seq_len(ncol(combos)))
      pair_list[[length(pair_list) + 1]] <- c(combos[1, ci], combos[2, ci])
  }
  if (length(pair_list) == 0) {
    return(data.frame(z_score = rep(NA_real_, N), indCors = rep(NA_real_, N),
                      rand_mean = rep(NA_real_, N), rand_sd = rep(NA_real_, N),
                      nF_detected = nF_est, n_pairs = 0L))
  }

  pair_mat <- do.call(rbind, pair_list)
  idx_A <- pair_mat[, 1]; idx_B <- pair_mat[, 2]
  k <- length(idx_A)

  # Option A: weight = |loading product|
  loading_product <- primary_loading[idx_A] * primary_loading[idx_B]
  weights <- abs(loading_product)
  weights[weights < 1e-6] <- 1e-6
  pair_agree <- sign(loading_product)

  A_coupled <- mat[, idx_A, drop = FALSE]
  B_coupled <- mat[, idx_B, drop = FALSE]

  if (align_signs) {
    needs_flip <- which(pair_agree < 0)
    if (length(needs_flip) > 0) {
      flip_max <- item_max[idx_B[needs_flip]]
      B_coupled[, needs_flip] <- rep(flip_max + 1, each = N) - B_coupled[, needs_flip]
    }
  }

  indCors <- rowCor_weighted(A_coupled, B_coupled, weights)

  raw_cor_mat <- cor(mat, use = "pairwise.complete.obs")
  signMat <- sign(raw_cor_mat)
  signMat[upper.tri(signMat, diag = TRUE)] <- NA

  rand_scores <- matrix(NA_real_, N, iterations)
  for (iter in seq_len(iterations)) {
    rand_idx1 <- sample(p, k, replace = TRUE)
    rand_idx2 <- sample(p, k, replace = TRUE)
    same <- rand_idx1 == rand_idx2
    while (any(same)) { rand_idx2[same] <- sample(p, sum(same), replace = TRUE); same <- rand_idx1 == rand_idx2 }
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
    rand_scores[, iter] <- rowCor_weighted(A_rand, B_rand, weights)
  }

  rand_mean <- rowMeans(rand_scores, na.rm = TRUE)
  rand_sd <- apply(rand_scores, 1, sd, na.rm = TRUE)
  rand_sd[rand_sd == 0] <- 1e-10
  z_score <- (indCors - rand_mean) / rand_sd

  data.frame(z_score = z_score, indCors = indCors,
             rand_mean = rand_mean, rand_sd = rand_sd,
             nF_detected = nF_est, n_pairs = k)
}

# ============================================================================
# WEIGHTED (inline, for comparison)
# ============================================================================

ReReReRe_weighted_inline <- function(data, iterations = 50, align_signs = TRUE) {
  data <- data[, sapply(data, is.numeric), drop = FALSE]
  mat <- as.matrix(data)
  N <- nrow(mat); p <- ncol(mat)
  item_max <- apply(mat, 2, max, na.rm = TRUE)

  cor_matrix <- cor(mat, use = "pairwise.complete.obs")
  pair_idx <- which(lower.tri(cor_matrix), arr.ind = TRUE)
  pair_r <- cor_matrix[pair_idx]
  pair_abs_r <- abs(pair_r)
  keep <- !is.na(pair_abs_r)
  pair_idx <- pair_idx[keep, , drop = FALSE]
  pair_r <- pair_r[keep]; pair_abs_r <- pair_abs_r[keep]
  pair_sign <- sign(pair_r)

  k <- nrow(pair_idx); weights <- pair_abs_r
  A_coupled <- mat[, pair_idx[, 1], drop = FALSE]
  B_coupled <- mat[, pair_idx[, 2], drop = FALSE]

  if (align_signs) {
    nf <- which(pair_sign < 0)
    if (length(nf) > 0) {
      B_coupled[, nf] <- rep(item_max[pair_idx[nf, 2]] + 1, each = N) - B_coupled[, nf]
    }
  }

  indCors <- rowCor_weighted(A_coupled, B_coupled, weights)

  signMat <- sign(cor_matrix); signMat[upper.tri(signMat, diag = TRUE)] <- NA
  rand_scores <- matrix(NA_real_, N, iterations)
  for (iter in seq_len(iterations)) {
    ri1 <- sample(p, k, replace = TRUE); ri2 <- sample(p, k, replace = TRUE)
    same <- ri1 == ri2
    while (any(same)) { ri2[same] <- sample(p, sum(same), replace = TRUE); same <- ri1 == ri2 }
    Ar <- mat[, ri1, drop = FALSE]; Br <- mat[, ri2, drop = FALSE]
    if (align_signs) {
      ri <- pmax(ri1, ri2); ci <- pmin(ri1, ri2)
      rs <- signMat[cbind(ri, ci)]; rs[is.na(rs)] <- 1
      rn <- which(rs < 0)
      if (length(rn) > 0) Br[, rn] <- rep(item_max[ri2[rn]] + 1, each = N) - Br[, rn]
    }
    rand_scores[, iter] <- rowCor_weighted(Ar, Br, weights)
  }

  rm <- rowMeans(rand_scores, na.rm = TRUE)
  rsd <- apply(rand_scores, 1, sd, na.rm = TRUE); rsd[rsd == 0] <- 1e-10
  data.frame(z_score = (indCors - rm) / rsd, indCors = indCors, rand_mean = rm, rand_sd = rsd)
}

# ============================================================================
# SIMULATION: Compare FOUR methods
# ============================================================================

cat("\n============================================================\n")
cat("OPTION D: EFA selection + |r| weights vs A, Weighted, Standard\n")
cat("============================================================\n\n")

nF_levels <- c(4, 6, 8, 10, 15, 20)
ipf_levels <- c(3, 6, 10)
n_reps <- 3
N <- 300
pct_careless <- 0.10
iterations <- 50

results <- data.frame()
total_cond <- length(nF_levels) * length(ipf_levels) * n_reps
counter <- 0

for (nF in nF_levels) {
  for (ipf in ipf_levels) {
    total_items <- nF * ipf
    for (rep_i in seq_len(n_reps)) {
      counter <- counter + 1
      seed_val <- 42000 + (counter - 1)
      set.seed(seed_val)

      cat(sprintf("[%d/%d] nF=%d, ipf=%d (%d items), rep=%d ... ",
                  counter, total_cond, nF, ipf, total_items, rep_i))

      good_data <- tryCatch(simulated_good_responses(nF, rep(ipf, nF), N),
                            error = function(e) { cat("SKIP\n"); NULL })
      if (is.null(good_data)) next

      inj <- inject_careless(good_data, pct_careless, seed = seed_val + 10000)
      corrupted <- inj$data_corrupted
      labels_vec <- as.integer(inj$labels$careless)

      # Standard
      r_std <- tryCatch(ReReReRe(corrupted, corProp=0.03, iterations=iterations,
                                  align_signs=TRUE, mode="coupled"), error=function(e) NULL)
      b_std <- if (!is.null(r_std)) find_best_z(r_std$z_score, labels_vec) else list(z=NA, mcc=NA)

      # Weighted
      r_wt <- tryCatch(ReReReRe_weighted_inline(corrupted, iterations=iterations),
                        error=function(e) NULL)
      b_wt <- if (!is.null(r_wt)) find_best_z(r_wt$z_score, labels_vec) else list(z=NA, mcc=NA)

      # EFA Option A (loading weights)
      r_efaA <- tryCatch(ReReReRe_efa_A(corrupted, iterations=iterations),
                          error=function(e) NULL)
      b_efaA <- if (!is.null(r_efaA)) find_best_z(r_efaA$z_score, labels_vec) else list(z=NA, mcc=NA)

      # EFA Option D (observed |r| weights)
      r_efaD <- tryCatch(ReReReRe_efa_D(corrupted, iterations=iterations),
                          error=function(e) NULL)
      b_efaD <- if (!is.null(r_efaD)) find_best_z(r_efaD$z_score, labels_vec) else list(z=NA, mcc=NA)

      nF_det <- if (!is.null(r_efaD)) r_efaD$nF_detected[1] else NA
      k_D <- if (!is.null(r_efaD)) r_efaD$n_pairs[1] else NA

      cat(sprintf("std=%.3f wt=%.3f efaA=%.3f efaD=%.3f (nF=%s, k=%s)\n",
                  ifelse(is.na(b_std$mcc), 0, b_std$mcc),
                  ifelse(is.na(b_wt$mcc), 0, b_wt$mcc),
                  ifelse(is.na(b_efaA$mcc), 0, b_efaA$mcc),
                  ifelse(is.na(b_efaD$mcc), 0, b_efaD$mcc),
                  ifelse(is.na(nF_det), "?", nF_det),
                  ifelse(is.na(k_D), "?", k_D)))

      results <- rbind(results, data.frame(
        nF=nF, ipf=ipf, total_items=total_items, rep=rep_i,
        mcc_std=b_std$mcc, z_std=b_std$z,
        mcc_wt=b_wt$mcc, z_wt=b_wt$z,
        mcc_efaA=b_efaA$mcc, z_efaA=b_efaA$z,
        mcc_efaD=b_efaD$mcc, z_efaD=b_efaD$z,
        nF_detected=nF_det, n_pairs_D=k_D
      ))
    }
  }
}

# ============================================================================
# SUMMARY
# ============================================================================

results$item_bin <- cut(results$total_items, breaks=c(0,30,60,100,200,Inf),
                        labels=c("<30","30-60","60-100","100-200",">200"))

cat("\n\n============================================================\n")
cat("SUMMARY BY ITEM BIN\n")
cat("============================================================\n\n")

bin_tab <- results %>%
  group_by(item_bin) %>%
  summarise(n=n(),
            std=mean(mcc_std, na.rm=T), wt=mean(mcc_wt, na.rm=T),
            efaA=mean(mcc_efaA, na.rm=T), efaD=mean(mcc_efaD, na.rm=T),
            D_vs_A=mean(mcc_efaD - mcc_efaA, na.rm=T),
            D_vs_std=mean(mcc_efaD - mcc_std, na.rm=T),
            D_vs_wt=mean(mcc_efaD - mcc_wt, na.rm=T),
            .groups="drop")

cat(sprintf("%-8s %3s %7s %7s %7s %7s %8s %8s %8s\n",
    "Items", "n", "Std", "Wt", "EFA-A", "EFA-D", "D-A", "D-Std", "D-Wt"))
cat(paste(rep("-", 80), collapse=""), "\n")
for (i in seq_len(nrow(bin_tab))) {
  r <- bin_tab[i,]
  cat(sprintf("%-8s %3d %7.3f %7.3f %7.3f %7.3f %+8.3f %+8.3f %+8.3f\n",
      r$item_bin, r$n, r$std, r$wt, r$efaA, r$efaD, r$D_vs_A, r$D_vs_std, r$D_vs_wt))
}

cat("\n\n============================================================\n")
cat("SUMMARY BY nF\n")
cat("============================================================\n\n")

nf_tab <- results %>%
  group_by(nF) %>%
  summarise(std=mean(mcc_std,na.rm=T), wt=mean(mcc_wt,na.rm=T),
            efaA=mean(mcc_efaA,na.rm=T), efaD=mean(mcc_efaD,na.rm=T), .groups="drop")
print(as.data.frame(nf_tab), digits=3)

cat("\n\n============================================================\n")
cat("SUMMARY BY nF x ipf\n")
cat("============================================================\n\n")

cell_tab <- results %>%
  group_by(nF, ipf, total_items) %>%
  summarise(std=mean(mcc_std,na.rm=T), wt=mean(mcc_wt,na.rm=T),
            efaA=mean(mcc_efaA,na.rm=T), efaD=mean(mcc_efaD,na.rm=T),
            nF_det=mean(nF_detected,na.rm=T), .groups="drop") %>%
  arrange(total_items)
print(as.data.frame(cell_tab), digits=3)

cat("\n\n============================================================\n")
cat("OVERALL MEANS\n")
cat("============================================================\n\n")

cat(sprintf("Standard:       %.3f\n", mean(results$mcc_std, na.rm=T)))
cat(sprintf("Weighted:       %.3f\n", mean(results$mcc_wt, na.rm=T)))
cat(sprintf("EFA Option A:   %.3f (loading weights)\n", mean(results$mcc_efaA, na.rm=T)))
cat(sprintf("EFA Option D:   %.3f (observed |r| weights)\n", mean(results$mcc_efaD, na.rm=T)))

cat("\n\nBest method per bin:\n")
for (i in seq_len(nrow(bin_tab))) {
  r <- bin_tab[i,]
  vals <- c(Std=r$std, Wt=r$wt, EFA_A=r$efaA, EFA_D=r$efaD)
  best <- names(which.max(vals))
  cat(sprintf("  %s: %s (%.3f)\n", r$item_bin, best, max(vals)))
}

write.csv(results, "test_efa_optionD_results.csv", row.names=FALSE)
cat("\n\nResults saved to test_efa_optionD_results.csv\n")
cat("\n=== DONE ===\n")
