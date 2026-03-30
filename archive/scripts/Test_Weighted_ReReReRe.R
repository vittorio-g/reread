###############################################################################
# Test: Weighted ReReReRe vs Standard ReReReRe
# Instead of selecting top-k% pairs, use ALL pairs weighted by |r_sample|
###############################################################################

library(lavaan)
library(psych)

source("Synthetic_Good_Responses_2.R")
source("Careless_machine_2.R")

# ============================================================================
# WEIGHTED rowCor — each pair weighted by sample-level |r|
# ============================================================================

rowCor_weighted <- function(A, B, weights) {
  # A, B: N x k matrices (item values for each pair)
  # weights: vector of length k (|r_sample| for each pair)
  # Returns: vector of length N (weighted coherence score per person)

  N <- nrow(A)
  k <- ncol(A)

  # Center each column (item) by its mean
  A_centered <- scale(A, center=TRUE, scale=TRUE)  # standardize
  B_centered <- scale(B, center=TRUE, scale=TRUE)

  # Replace NAs from zero-variance columns
  A_centered[is.na(A_centered)] <- 0
  B_centered[is.na(B_centered)] <- 0

  # For each person: weighted sum of (z_A * z_B) across pairs
  # This is a weighted average of pair-level "agreement"
  cross_products <- A_centered * B_centered  # N x k

  # Weight each pair
  w_matrix <- matrix(weights, nrow=N, ncol=k, byrow=TRUE)

  weighted_sum <- rowSums(cross_products * w_matrix, na.rm=TRUE)
  total_weight <- sum(weights)

  return(weighted_sum / total_weight)
}

# ============================================================================
# WEIGHTED ReReReRe
# ============================================================================

ReReReRe_weighted <- function(data, iterations=100, align_signs=TRUE,
                               min_r=0.0, z_threshold=NULL) {
  # min_r: minimum |r| to include a pair (0 = use all pairs)
  # Everything else same logic as standard ReReReRe

  data <- as.data.frame(data)
  N <- nrow(data)
  p <- ncol(data)

  # 1. Sample-level correlation matrix
  cor_matrix <- cor(data, use="pairwise.complete.obs")

  # 2. Get ALL pairs and their |r|
  pair_idx <- which(upper.tri(cor_matrix), arr.ind=TRUE)
  pair_r <- cor_matrix[pair_idx]
  pair_abs_r <- abs(pair_r)
  pair_sign <- sign(pair_r)

  # Filter by min_r
  keep <- pair_abs_r >= min_r & !is.na(pair_abs_r)
  pair_idx <- pair_idx[keep, ]
  pair_r <- pair_r[keep]
  pair_abs_r <- pair_abs_r[keep]
  pair_sign <- pair_sign[keep]

  k <- nrow(pair_idx)
  cat(sprintf("  Weighted mode: using %d pairs (min_r=%.2f), mean |r|=%.3f\n",
      k, min_r, mean(pair_abs_r)))

  # 3. Build A and B matrices
  A <- as.matrix(data[, pair_idx[,1]])
  B <- as.matrix(data[, pair_idx[,2]])

  # 4. Align signs if needed
  if (align_signs) {
    neg_pairs <- pair_sign < 0
    if (any(neg_pairs)) {
      for (j in which(neg_pairs)) {
        col_max <- max(B[, j], na.rm=TRUE)
        B[, j] <- (col_max + 1) - B[, j]
      }
      # After alignment, update weights (use original |r|)
    }
  }

  # 5. Compute weighted coherence score (coupled)
  weights <- pair_abs_r  # weight = sample-level |r|
  indCors <- rowCor_weighted(A, B, weights)

  # 6. Permutation baseline
  rand_scores <- matrix(NA, N, iterations)

  for (iter in 1:iterations) {
    # Random pairs: shuffle column assignments
    rand_idx1 <- sample(p, k, replace=TRUE)
    rand_idx2 <- sample(p, k, replace=TRUE)

    # Ensure different items in each pair
    same <- rand_idx1 == rand_idx2
    while (any(same)) {
      rand_idx2[same] <- sample(p, sum(same), replace=TRUE)
      same <- rand_idx1 == rand_idx2
    }

    A_rand <- as.matrix(data[, rand_idx1])
    B_rand <- as.matrix(data[, rand_idx2])

    # Get random pair correlations for weighting
    rand_r <- numeric(k)
    for (j in 1:k) {
      rand_r[j] <- abs(cor(A_rand[,j], B_rand[,j], use="pairwise.complete.obs"))
    }
    rand_r[is.na(rand_r)] <- 0

    # Align signs for random pairs too
    if (align_signs) {
      for (j in 1:k) {
        raw_r <- cor(A_rand[,j], B_rand[,j], use="pairwise.complete.obs")
        if (!is.na(raw_r) && raw_r < 0) {
          col_max <- max(B_rand[,j], na.rm=TRUE)
          B_rand[,j] <- (col_max + 1) - B_rand[,j]
        }
      }
    }

    # Use SAME weights as coupled (not random pair weights)
    # This keeps the comparison fair — coupled pairs are weighted by their
    # known high correlation, random pairs by the same weights
    rand_scores[, iter] <- rowCor_weighted(A_rand, B_rand, weights)
  }

  # 7. Z-scores
  rand_mean <- rowMeans(rand_scores, na.rm=TRUE)
  rand_sd <- apply(rand_scores, 1, sd, na.rm=TRUE)
  rand_sd[rand_sd == 0] <- 1e-10

  z_score <- (indCors - rand_mean) / rand_sd

  # 8. Flag
  if (is.null(z_threshold)) {
    total_items <- p
    z_threshold <- max(0.3, min(3.5, 0.505 + 0.0042 * total_items))
    cat(sprintf("  Auto z_threshold: %.2f\n", z_threshold))
  }

  flagged <- as.integer(z_score <= z_threshold)

  return(data.frame(
    z_score = z_score,
    indCors = indCors,
    rand_mean = rand_mean,
    flagged = flagged,
    z_threshold_used = z_threshold
  ))
}


# ============================================================================
# STANDARD ReReReRe (for comparison)
# ============================================================================
source("ReReReRe.R")


# ============================================================================
# TEST: Compare on simulated data
# ============================================================================

cat("\n============================================================\n")
cat("WEIGHTED vs STANDARD ReReReRe — SIMULATION COMPARISON\n")
cat("============================================================\n")

evaluate_mcc <- function(predicted, actual) {
  tp <- sum(predicted == 1 & actual == 1)
  tn <- sum(predicted == 0 & actual == 0)
  fp <- sum(predicted == 1 & actual == 0)
  fn <- sum(predicted == 0 & actual == 1)
  denom <- as.double(tp+fp) * as.double(tp+fn) * as.double(tn+fp) * as.double(tn+fn)
  if (denom == 0) return(0)
  return((as.double(tp)*tn - as.double(fp)*fn) / sqrt(denom))
}

find_best_z <- function(z_scores, labels) {
  best_mcc <- -1; best_z <- NA
  for (z in seq(0.1, 3.0, 0.1)) {
    pred <- as.integer(z_scores <= z)
    mcc <- evaluate_mcc(pred, labels)
    if (mcc > best_mcc) { best_mcc <- mcc; best_z <- z }
  }
  return(list(z=best_z, mcc=best_mcc))
}

# Test conditions
conditions <- expand.grid(
  nF = c(4, 6, 8, 10, 15, 20),
  ipf = c(3, 6, 10),
  stringsAsFactors = FALSE
)

n <- 300
pct_careless <- 0.10
reps <- 3

results <- data.frame()

for (i in 1:nrow(conditions)) {
  nF <- conditions$nF[i]
  ipf <- conditions$ipf[i]
  total_items <- nF * ipf

  cat(sprintf("\n--- nF=%d, ipf=%d (total=%d items) ---\n", nF, ipf, total_items))

  for (rep in 1:reps) {
    set.seed(1000*i + rep)

    # Generate data
    tryCatch({
      good_data <- simulated_good_responses(nF, rep(ipf, nF), n)
      inj <- inject_careless(good_data, pct_careless)
      data <- inj$data
      labels <- as.integer(inj$labels$careless)

      # --- Standard ReReReRe ---
      rr_std <- ReReReRe(data, corProp=0.03, iterations=50, align_signs=TRUE)
      best_std <- find_best_z(rr_std$z_score, labels)

      # Fixed z=1.5
      mcc_std_15 <- evaluate_mcc(as.integer(rr_std$z_score <= 1.5), labels)

      # --- Weighted ReReReRe ---
      rr_w0 <- ReReReRe_weighted(data, iterations=50, align_signs=TRUE,
                                   min_r=0.0, z_threshold=1.5)
      mcc_w0_15 <- evaluate_mcc(rr_w0$flagged, labels)
      best_w0 <- find_best_z(rr_w0$z_score, labels)

      # Weighted with min_r=0.10 (skip very weak pairs)
      rr_w10 <- ReReReRe_weighted(data, iterations=50, align_signs=TRUE,
                                    min_r=0.10, z_threshold=1.5)
      mcc_w10_15 <- evaluate_mcc(rr_w10$flagged, labels)
      best_w10 <- find_best_z(rr_w10$z_score, labels)

      results <- rbind(results, data.frame(
        nF=nF, ipf=ipf, total_items=total_items, rep=rep,
        std_mcc_oracle=best_std$mcc, std_mcc_z15=mcc_std_15,
        w0_mcc_oracle=best_w0$mcc, w0_mcc_z15=mcc_w0_15,
        w10_mcc_oracle=best_w10$mcc, w10_mcc_z15=mcc_w10_15
      ))

      cat(sprintf("  rep %d: std=%.3f(%.3f) w_all=%.3f(%.3f) w_10=%.3f(%.3f) [oracle(z15)]\n",
          rep, best_std$mcc, mcc_std_15, best_w0$mcc, mcc_w0_15, best_w10$mcc, mcc_w10_15))

    }, error = function(e) {
      cat(sprintf("  rep %d: ERROR: %s\n", rep, e$message))
    })
  }
}

# ============================================================================
# SUMMARY TABLE
# ============================================================================

cat("\n\n============================================================\n")
cat("SUMMARY: Mean MCC by condition\n")
cat("============================================================\n")

library(dplyr)

summary_tab <- results %>%
  group_by(nF, ipf, total_items) %>%
  summarise(
    std_oracle = mean(std_mcc_oracle),
    std_z15 = mean(std_mcc_z15),
    w_all_oracle = mean(w0_mcc_oracle),
    w_all_z15 = mean(w0_mcc_z15),
    w_10_oracle = mean(w10_mcc_oracle),
    w_10_z15 = mean(w10_mcc_z15),
    .groups="drop"
  ) %>%
  arrange(total_items)

cat("\nnF  ipf  items | Standard      | Weighted(all) | Weighted(r>.10)\n")
cat("               | oracle  z=1.5 | oracle  z=1.5 | oracle  z=1.5\n")
cat("---------------------------------------------------------------\n")
for (j in 1:nrow(summary_tab)) {
  r <- summary_tab[j, ]
  cat(sprintf("%2d  %2d  %4d  | %5.3f  %5.3f | %5.3f  %5.3f | %5.3f  %5.3f\n",
      r$nF, r$ipf, r$total_items,
      r$std_oracle, r$std_z15,
      r$w_all_oracle, r$w_all_z15,
      r$w_10_oracle, r$w_10_z15))
}

# Overall means
cat("\n--- Overall means ---\n")
cat(sprintf("Standard:        oracle=%.3f, z=1.5=%.3f\n",
    mean(results$std_mcc_oracle), mean(results$std_mcc_z15)))
cat(sprintf("Weighted(all):   oracle=%.3f, z=1.5=%.3f\n",
    mean(results$w0_mcc_oracle), mean(results$w0_mcc_z15)))
cat(sprintf("Weighted(r>.10): oracle=%.3f, z=1.5=%.3f\n",
    mean(results$w10_mcc_oracle), mean(results$w10_mcc_z15)))

# By total_items bins
cat("\n--- By total items ---\n")
results$items_bin <- cut(results$total_items, breaks=c(0, 30, 60, 100, 200, 400),
                          labels=c("<30", "30-60", "60-100", "100-200", ">200"))
bin_summary <- results %>%
  group_by(items_bin) %>%
  summarise(
    n = n(),
    std = mean(std_mcc_oracle),
    w_all = mean(w0_mcc_oracle),
    w_10 = mean(w10_mcc_oracle),
    .groups="drop"
  )
print(bin_summary)

# Save results
write.csv(results, "archive/weighted_rerere_comparison.csv", row.names=FALSE)

cat("\n=== DONE ===\n")
