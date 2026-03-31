###############################################################################
# Test: EFA-based ReReReRe (Option A) vs Standard vs Weighted
#
# EFA-based pair selection:
#   1. Parallel analysis -> estimate nF
#   2. EFA with oblimin rotation -> assign items to primary factor
#   3. Generate ALL within-factor pairs
#   4. Weight each pair by |lambda_i * lambda_j|
#   5. Sign alignment via loading product sign
#   6. Permutation baseline: random pairs from ALL items (same k, same weights)
#
# Compares THREE methods:
#   - Standard ReReReRe (top-k% pairs)
#   - Weighted ReReReRe (all pairs weighted by |r|)
#   - EFA-based ReReReRe (within-factor pairs weighted by loading products)
###############################################################################

library(lavaan)
library(psych)
library(dplyr)

source("Synthetic_Good_Responses_2.R")
source("Careless_machine_2.R")
source("ReReReRe.R")

# ============================================================================
# HELPER: rowCor_weighted (same as in ReReReRe.R)
# ============================================================================

# Already defined in ReReReRe.R, but kept here for clarity:
# rowCor_weighted <- function(A, B, weights) { ... }

# ============================================================================
# HELPER: MCC and oracle best z
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
  best_mcc <- -1
  best_z <- NA
  for (z in seq(0.1, 3.0, 0.1)) {
    pred <- as.integer(z_scores <= z)
    mcc <- evaluate_mcc(pred, labels)
    if (mcc > best_mcc) {
      best_mcc <- mcc
      best_z <- z
    }
  }
  return(list(z = best_z, mcc = best_mcc))
}

# ============================================================================
# WEIGHTED ReReReRe (all pairs, |r|-weighted) -- inline implementation
# ============================================================================

ReReReRe_weighted_inline <- function(data, iterations = 50, align_signs = TRUE) {
  data <- data[, sapply(data, is.numeric), drop = FALSE]
  mat <- as.matrix(data)
  N <- nrow(mat)
  p <- ncol(mat)

  # Precompute item max for reverse coding
  item_max <- apply(mat, 2, max, na.rm = TRUE)

  # Sample-level correlations
  cor_matrix <- cor(mat, use = "pairwise.complete.obs")
  raw_cor <- cor_matrix  # signed

  # Get ALL unique pairs from lower triangle
  pair_idx <- which(lower.tri(cor_matrix), arr.ind = TRUE)
  pair_r <- raw_cor[pair_idx]
  pair_abs_r <- abs(pair_r)
  pair_sign <- sign(pair_r)

  # Remove NA pairs
  keep <- !is.na(pair_abs_r)
  pair_idx <- pair_idx[keep, , drop = FALSE]
  pair_r <- pair_r[keep]
  pair_abs_r <- pair_abs_r[keep]
  pair_sign <- pair_sign[keep]

  k <- nrow(pair_idx)
  weights <- pair_abs_r

  # Build A and B matrices
  A_coupled <- mat[, pair_idx[, 1], drop = FALSE]
  B_coupled <- mat[, pair_idx[, 2], drop = FALSE]

  # Sign alignment: reverse-code B for negative-correlation pairs
  if (align_signs) {
    needs_flip <- which(pair_sign < 0)
    if (length(needs_flip) > 0) {
      flip_max <- item_max[pair_idx[needs_flip, 2]]
      B_coupled[, needs_flip] <- rep(flip_max + 1, each = N) - B_coupled[, needs_flip]
    }
  }

  # Weighted coherence score (coupled)
  indCors <- rowCor_weighted(A_coupled, B_coupled, weights)

  # Precompute sign matrix for random pair alignment
  signMat <- sign(raw_cor)
  signMat[upper.tri(signMat, diag = TRUE)] <- NA

  # Permutation baseline
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

  # Z-scores
  rand_mean <- rowMeans(rand_scores, na.rm = TRUE)
  rand_sd <- apply(rand_scores, 1, sd, na.rm = TRUE)
  rand_sd[rand_sd == 0] <- 1e-10
  z_score <- (indCors - rand_mean) / rand_sd

  data.frame(
    z_score = z_score,
    indCors = indCors,
    rand_mean = rand_mean,
    rand_sd = rand_sd
  )
}

# ============================================================================
# EFA-BASED ReReReRe (Option A)
# ============================================================================

ReReReRe_efa <- function(data, iterations = 50, align_signs = TRUE,
                         nf_override = NULL) {
  # data: numeric dataframe/matrix (rows=respondents, cols=items)
  # iterations: number of random permutations
  # align_signs: use loading signs for reverse-coded items
  # nf_override: if given, skip parallel analysis and use this nF

  data <- data[, sapply(data, is.numeric), drop = FALSE]
  mat <- as.matrix(data)
  N <- nrow(mat)
  p <- ncol(mat)

  # Precompute item max for reverse coding
  item_max <- apply(mat, 2, max, na.rm = TRUE)

  # --- Step 1: Estimate nF via parallel analysis ---
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
      # Fallback: use eigenvalue > 1 rule
      ev <- eigen(cor(mat, use = "pairwise.complete.obs"),
                  symmetric = TRUE, only.values = TRUE)$values
      nF_est <- max(1, sum(ev > 1))
      warning(sprintf("Parallel analysis failed, using eigenvalue>1 rule: nF=%d", nF_est))
    }
  }

  # Clamp nF to sensible range
  nF_est <- max(1, min(nF_est, floor(p / 2)))

  # --- Step 2: Run EFA ---
  efa_result <- tryCatch({
    suppressWarnings(
      fa(mat, nfactors = nF_est, rotate = "oblimin", fm = "minres",
         scores = "none", warnings = FALSE)
    )
  }, error = function(e) {
    # Try with fewer factors
    tryCatch({
      nF_reduced <- max(1, nF_est - 1)
      suppressWarnings(
        fa(mat, nfactors = nF_reduced, rotate = "oblimin", fm = "minres",
           scores = "none", warnings = FALSE)
      )
    }, error = function(e2) NULL)
  })

  if (is.null(efa_result)) {
    warning("EFA failed completely. Returning NA z-scores.")
    return(data.frame(
      z_score = rep(NA_real_, N),
      indCors = rep(NA_real_, N),
      rand_mean = rep(NA_real_, N),
      rand_sd = rep(NA_real_, N),
      nF_detected = nF_est,
      n_pairs = NA_integer_
    ))
  }

  # --- Step 3: Assign items to primary factor ---
  loadings_mat <- as.matrix(efa_result$loadings[])  # p x nF matrix
  nF_actual <- ncol(loadings_mat)

  # For each item: which factor has the highest |loading|?
  primary_factor <- apply(abs(loadings_mat), 1, which.max)  # vector of length p
  primary_loading <- numeric(p)
  for (i in seq_len(p)) {
    primary_loading[i] <- loadings_mat[i, primary_factor[i]]  # SIGNED loading
  }

  # --- Step 4: Generate ALL within-factor pairs ---
  pair_list <- list()
  for (f in seq_len(nF_actual)) {
    items_in_factor <- which(primary_factor == f)
    n_items_f <- length(items_in_factor)
    if (n_items_f < 2) next

    # All unique pairs within this factor
    combos <- combn(items_in_factor, 2)
    for (c_idx in seq_len(ncol(combos))) {
      i_item <- combos[1, c_idx]
      j_item <- combos[2, c_idx]
      pair_list[[length(pair_list) + 1]] <- c(i_item, j_item)
    }
  }

  if (length(pair_list) == 0) {
    warning("No within-factor pairs found. All factors have <2 items.")
    return(data.frame(
      z_score = rep(NA_real_, N),
      indCors = rep(NA_real_, N),
      rand_mean = rep(NA_real_, N),
      rand_sd = rep(NA_real_, N),
      nF_detected = nF_est,
      n_pairs = 0L
    ))
  }

  # Convert to index matrix (k x 2)
  pair_mat <- do.call(rbind, pair_list)
  k <- nrow(pair_mat)
  idx_A <- pair_mat[, 1]
  idx_B <- pair_mat[, 2]

  # --- Step 5: Compute weights and sign alignment ---
  # Weight = |lambda_i * lambda_j| (product of primary loadings)
  # Sign alignment: if lambda_i * lambda_j < 0, items disagree -> reverse-code B
  loading_product <- primary_loading[idx_A] * primary_loading[idx_B]
  weights <- abs(loading_product)

  # Ensure no zero weights (very weak items)
  weights[weights < 1e-6] <- 1e-6

  # Sign of loading product determines alignment
  pair_agree <- sign(loading_product)  # +1 = same direction, -1 = reverse-coded

  # Build A and B matrices
  A_coupled <- mat[, idx_A, drop = FALSE]
  B_coupled <- mat[, idx_B, drop = FALSE]

  # Apply sign alignment: reverse-code B where loading product is negative
  if (align_signs) {
    needs_flip <- which(pair_agree < 0)
    if (length(needs_flip) > 0) {
      flip_max <- item_max[idx_B[needs_flip]]
      B_coupled[, needs_flip] <- rep(flip_max + 1, each = N) - B_coupled[, needs_flip]
    }
  }

  # --- Step 6: Compute weighted coherence score ---
  indCors <- rowCor_weighted(A_coupled, B_coupled, weights)

  # --- Step 7: Permutation baseline ---
  # Sample random pairs from ALL items (not restricted to within-factor)
  # Same number k, same weights vector (to keep comparison fair)

  # Precompute sign matrix for random pair alignment
  raw_cor <- cor(mat, use = "pairwise.complete.obs")
  signMat <- sign(raw_cor)
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

    # Sign alignment for random pairs using sample correlation signs
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
    rand_scores[, iter] <- rowCor_weighted(A_rand, B_rand, weights)
  }

  # --- Step 8: Z-scores ---
  rand_mean <- rowMeans(rand_scores, na.rm = TRUE)
  rand_sd <- apply(rand_scores, 1, sd, na.rm = TRUE)
  rand_sd[rand_sd == 0] <- 1e-10
  z_score <- (indCors - rand_mean) / rand_sd

  data.frame(
    z_score = z_score,
    indCors = indCors,
    rand_mean = rand_mean,
    rand_sd = rand_sd,
    nF_detected = nF_est,
    n_pairs = k
  )
}


# ============================================================================
# SIMULATION: Compare three methods
# ============================================================================

cat("\n============================================================\n")
cat("EFA-BASED vs WEIGHTED vs STANDARD ReReReRe\n")
cat("============================================================\n\n")

# Test conditions
nF_levels <- c(4, 6, 8, 10, 15, 20)
ipf_levels <- c(3, 6, 10)
n_reps <- 3
N <- 300
pct_careless <- 0.10
iterations <- 50

# Results storage
results <- data.frame()

total_conditions <- length(nF_levels) * length(ipf_levels) * n_reps
cond_counter <- 0

for (nF in nF_levels) {
  for (ipf in ipf_levels) {
    total_items <- nF * ipf
    for (rep_i in seq_len(n_reps)) {
      cond_counter <- cond_counter + 1
      seed_val <- 42000 + (cond_counter - 1)
      set.seed(seed_val)

      cat(sprintf("[%d/%d] nF=%d, ipf=%d (%d items), rep=%d ... ",
                  cond_counter, total_conditions, nF, ipf, total_items, rep_i))

      # --- Generate data ---
      good_data <- tryCatch({
        simulated_good_responses(nF, rep(ipf, nF), N)
      }, error = function(e) {
        cat(sprintf("SKIP (data gen failed: %s)\n", e$message))
        NULL
      })
      if (is.null(good_data)) next

      inj <- inject_careless(good_data, pct_careless, seed = seed_val + 10000)
      corrupted <- inj$data_corrupted
      labels_vec <- as.integer(inj$labels$careless)

      # --- Method 1: Standard ReReReRe ---
      res_std <- tryCatch({
        ReReReRe(corrupted, corProp = 0.03, iterations = iterations,
                 align_signs = TRUE, mode = "coupled")
      }, error = function(e) NULL)

      if (!is.null(res_std)) {
        best_std <- find_best_z(res_std$z_score, labels_vec)
      } else {
        best_std <- list(z = NA, mcc = NA)
      }

      # --- Method 2: Weighted ReReReRe (all pairs, |r|-weighted) ---
      res_wt <- tryCatch({
        ReReReRe_weighted_inline(corrupted, iterations = iterations, align_signs = TRUE)
      }, error = function(e) NULL)

      if (!is.null(res_wt)) {
        best_wt <- find_best_z(res_wt$z_score, labels_vec)
      } else {
        best_wt <- list(z = NA, mcc = NA)
      }

      # --- Method 3: EFA-based ReReReRe ---
      res_efa <- tryCatch({
        ReReReRe_efa(corrupted, iterations = iterations, align_signs = TRUE)
      }, error = function(e) NULL)

      if (!is.null(res_efa)) {
        best_efa <- find_best_z(res_efa$z_score, labels_vec)
        nF_detected <- res_efa$nF_detected[1]
        n_efa_pairs <- res_efa$n_pairs[1]
      } else {
        best_efa <- list(z = NA, mcc = NA)
        nF_detected <- NA
        n_efa_pairs <- NA
      }

      cat(sprintf("std=%.3f wt=%.3f efa=%.3f (nF_det=%s, k=%s)\n",
                  ifelse(is.na(best_std$mcc), 0, best_std$mcc),
                  ifelse(is.na(best_wt$mcc), 0, best_wt$mcc),
                  ifelse(is.na(best_efa$mcc), 0, best_efa$mcc),
                  ifelse(is.na(nF_detected), "?", as.character(nF_detected)),
                  ifelse(is.na(n_efa_pairs), "?", as.character(n_efa_pairs))))

      results <- rbind(results, data.frame(
        nF = nF,
        ipf = ipf,
        total_items = total_items,
        rep = rep_i,
        mcc_standard = best_std$mcc,
        z_standard = best_std$z,
        mcc_weighted = best_wt$mcc,
        z_weighted = best_wt$z,
        mcc_efa = best_efa$mcc,
        z_efa = best_efa$z,
        nF_detected = nF_detected,
        n_efa_pairs = n_efa_pairs,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# ============================================================================
# SUMMARY TABLES
# ============================================================================

cat("\n\n============================================================\n")
cat("RESULTS BY nF x ipf\n")
cat("============================================================\n\n")

summary_by_cell <- results %>%
  group_by(nF, ipf, total_items) %>%
  summarise(
    mean_std = mean(mcc_standard, na.rm = TRUE),
    mean_wt  = mean(mcc_weighted, na.rm = TRUE),
    mean_efa = mean(mcc_efa, na.rm = TRUE),
    mean_nF_det = mean(nF_detected, na.rm = TRUE),
    mean_k_efa = mean(n_efa_pairs, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(total_items)

print(as.data.frame(summary_by_cell), digits = 3, row.names = FALSE)

# --- Summary by item range bins ---
cat("\n\n============================================================\n")
cat("SUMMARY BY TOTAL ITEMS BINS\n")
cat("============================================================\n\n")

results$item_bin <- cut(results$total_items,
                        breaks = c(0, 30, 60, 100, 200, Inf),
                        labels = c("<30", "30-60", "60-100", "100-200", ">200"),
                        right = TRUE)

summary_by_bin <- results %>%
  group_by(item_bin) %>%
  summarise(
    n_cells = n(),
    mean_std = mean(mcc_standard, na.rm = TRUE),
    mean_wt  = mean(mcc_weighted, na.rm = TRUE),
    mean_efa = mean(mcc_efa, na.rm = TRUE),
    efa_vs_std = mean(mcc_efa - mcc_standard, na.rm = TRUE),
    efa_vs_wt  = mean(mcc_efa - mcc_weighted, na.rm = TRUE),
    .groups = "drop"
  )

cat("Mean oracle MCC by total items bin:\n\n")
print(as.data.frame(summary_by_bin), digits = 3, row.names = FALSE)

# --- Summary by nF ---
cat("\n\n============================================================\n")
cat("SUMMARY BY nF (averaged over ipf)\n")
cat("============================================================\n\n")

summary_by_nF <- results %>%
  group_by(nF) %>%
  summarise(
    mean_std = mean(mcc_standard, na.rm = TRUE),
    mean_wt  = mean(mcc_weighted, na.rm = TRUE),
    mean_efa = mean(mcc_efa, na.rm = TRUE),
    .groups = "drop"
  )

print(as.data.frame(summary_by_nF), digits = 3, row.names = FALSE)

# --- Overall ---
cat("\n\n============================================================\n")
cat("OVERALL MEANS\n")
cat("============================================================\n\n")

cat(sprintf("Standard ReReReRe:  %.3f\n", mean(results$mcc_standard, na.rm = TRUE)))
cat(sprintf("Weighted ReReReRe:  %.3f\n", mean(results$mcc_weighted, na.rm = TRUE)))
cat(sprintf("EFA-based ReReReRe: %.3f\n", mean(results$mcc_efa, na.rm = TRUE)))

# --- Best method per bin ---
cat("\n\nBest method per item bin:\n")
for (i in seq_len(nrow(summary_by_bin))) {
  row <- summary_by_bin[i, ]
  vals <- c(std = row$mean_std, wt = row$mean_wt, efa = row$mean_efa)
  best <- names(which.max(vals))
  cat(sprintf("  %s items: %s (%.3f)\n", row$item_bin, best, max(vals)))
}

# Save results
write.csv(results, "test_efa_results.csv", row.names = FALSE)
cat("\n\nResults saved to test_efa_results.csv\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
