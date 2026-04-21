# Simulation_WeightedV2_vs_Top.R
# Confronto veloce: weighted_v2 vs i top performer dalla sim All Variants
# (std, efa_d, iterative). 45 condizioni × 3 reps.
#
# Setup identico a Simulation_AllVariants.R: N=300, 15% careless,
# corruzione 10-100%, GT = >50%. Runtime atteso: ~45-70 minuti.

suppressPackageStartupMessages({
  library(dplyr)
})
source("Synthetic_Good_Responses_2.R")
source("Careless_machine_2.R")
source("ReReReRe.R")

# ============================================================
# weighted_v2: correlazione pesata propriamente
# ============================================================
rowCor_weighted_v2 <- function(A, B, weights, zero_val = 0) {
  sum_w <- sum(weights)
  W <- matrix(weights, nrow(A), ncol(A), byrow = TRUE)
  mean_A <- rowSums(A * W, na.rm = TRUE) / sum_w
  mean_B <- rowSums(B * W, na.rm = TRUE) / sum_w
  A_c <- A - mean_A; B_c <- B - mean_B
  num <- rowSums(W * A_c * B_c, na.rm = TRUE)
  den <- sqrt(rowSums(W * A_c^2, na.rm = TRUE) * rowSums(W * B_c^2, na.rm = TRUE))
  r <- abs(num / den); r[is.nan(r) | is.na(r)] <- zero_val; r
}

score_weighted_v2 <- function(data, iterations = 100, align_signs = TRUE) {
  data <- data[, sapply(data, is.numeric), drop = FALSE]
  N <- nrow(data); J <- ncol(data); mat <- as.matrix(data)
  if (align_signs) item_max <- apply(mat, 2, max, na.rm = TRUE)
  rawCorMat <- cor(mat, use = "pairwise.complete.obs")
  corMat <- abs(rawCorMat)
  corMat[upper.tri(corMat, diag = TRUE)] <- NA
  rawCorMat[upper.tri(rawCorMat, diag = TRUE)] <- NA
  if (align_signs) signMat <- sign(rawCorMat)
  pair_idx <- which(!is.na(corMat), arr.ind = TRUE)
  weights <- corMat[pair_idx]
  pair_signs <- sign(rawCorMat[pair_idx])
  k <- nrow(pair_idx)
  A_sel <- mat[, pair_idx[, 1], drop = FALSE]
  B_sel <- mat[, pair_idx[, 2], drop = FALSE]
  if (align_signs) {
    nf <- which(pair_signs < 0)
    if (length(nf) > 0)
      B_sel[, nf] <- rep(item_max[pair_idx[nf, 2]] + 1, each = N) - B_sel[, nf]
  }
  rowCors <- rowCor_weighted_v2(A_sel, B_sel, weights, 0)
  all_RIC <- matrix(NA_real_, nrow = N, ncol = iterations)
  for (i in seq_len(iterations)) {
    idx1 <- sample(J, k, replace = TRUE); idx2 <- sample(J, k, replace = TRUE)
    same <- idx1 == idx2
    while (any(same)) { idx2[same] <- sample(J, sum(same), replace = TRUE); same <- idx1 == idx2 }
    A_r <- mat[, idx1, drop = FALSE]; B_r <- mat[, idx2, drop = FALSE]
    if (align_signs) {
      ri <- pmax(idx1, idx2); ci <- pmin(idx1, idx2)
      rs <- signMat[cbind(ri, ci)]; rs[is.na(rs)] <- 1
      rn <- which(rs < 0)
      if (length(rn) > 0)
        B_r[, rn] <- rep(item_max[idx2[rn]] + 1, each = N) - B_r[, rn]
    }
    all_RIC[, i] <- rowCor_weighted_v2(A_r, B_r, weights, 0)
  }
  rm <- rowMeans(all_RIC, na.rm = TRUE)
  rs <- apply(all_RIC, 1, sd, na.rm = TRUE)
  ifelse(rs > 0, (rowCors - rm) / rs, 0)
}

# ============================================================
# Griglia: stessa di Simulation_AllVariants.R ma 3 reps
# ============================================================
NF_LEVELS  <- c(4, 8, 12, 20, 30)
IPF_LEVELS <- c(3, 6, 10)
REPS <- 3
N_FIXED <- 300
PCT_CARELESS <- 0.15
ITER <- 100
CORPROP <- 0.03
SEED_BASE <- 2026

# Corruzione realistica 10-100%, GT solo >50%
CARELESS_LEVELS <- seq(0.1, 1.0, 0.1)
CARELESS_TYPES <- c(random = 1/3, longstring = 1/3, mixed = 1/3)
GT_CUTOFF <- 0.50

Z_GRID <- seq(0.1, 3.0, 0.2)

conditions <- expand.grid(nF = NF_LEVELS, ipf = IPF_LEVELS,
                          stringsAsFactors = FALSE)
conditions$total_items <- conditions$nF * conditions$ipf
TOTAL_CONDITIONS <- nrow(conditions) * REPS

cat("============================================================\n")
cat(" WEIGHTED_V2 vs TOP METHODS\n")
cat("============================================================\n")
cat(sprintf("%d conditions × %d reps = %d runs\n",
            nrow(conditions), REPS, TOTAL_CONDITIONS))
cat("Methods: std, efa_d, iterative, weighted_v2\n\n")

# ============================================================
# Helpers
# ============================================================
compute_mcc_at <- function(z, lab, zt) {
  flag <- z <= zt
  TP <- sum(flag & lab); TN <- sum(!flag & !lab)
  FP <- sum(flag & !lab); FN <- sum(!flag & lab)
  den <- sqrt(as.numeric(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  if (den == 0) 0 else (TP*TN - FP*FN) / den
}
oracle_mcc <- function(z, lab) {
  max(sapply(Z_GRID, function(zt) compute_mcc_at(z, lab, zt)), na.rm = TRUE)
}

# ============================================================
# Loop
# ============================================================
results <- list(); idx <- 0
t0 <- Sys.time()

for (ci in seq_len(nrow(conditions))) {
  nF <- conditions$nF[ci]
  ipf <- conditions$ipf[ci]
  total <- conditions$total_items[ci]

  for (rep_id in seq_len(REPS)) {
    idx <- idx + 1
    seed_i <- SEED_BASE + (rep_id - 1) * 13 + ci
    set.seed(seed_i)

    clean_data <- simulated_good_responses(nConstructs = nF,
                                           nItems = rep(ipf, nF),
                                           n = N_FIXED)
    inj <- inject_careless(clean_data,
                           pct_careless = PCT_CARELESS,
                           pct_types = CARELESS_TYPES,
                           careless_levels = CARELESS_LEVELS)
    labels <- inj$labels$careless_pct > GT_CUTOFF
    n_careless_gt <- sum(labels)
    n_noise <- sum(inj$labels$careless_pct > 0 & inj$labels$careless_pct <= GT_CUTOFF)

    cat(sprintf("[%d/%d] nF=%d ipf=%d (%d items) rep=%d | GT:%d/noise:%d | ",
                idx, TOTAL_CONDITIONS, nF, ipf, total, rep_id,
                n_careless_gt, n_noise))

    data_mat <- inj$data_corrupted

    # --- 4 methods ---
    z_std <- ReReReRe(data_mat, corProp = CORPROP, iterations = ITER,
                      align_signs = TRUE, mode = "auto")$z_score
    z_efa <- ReReReRe(data_mat, iterations = ITER, align_signs = TRUE,
                      mode = "efa_d")$z_score
    z_iter <- ReReReRe_F_iterative(data_mat, iterations = ITER,
                                    align_signs = TRUE,
                                    initial_z_threshold = 1.0)$z_score
    z_v2 <- score_weighted_v2(data_mat, iterations = ITER, align_signs = TRUE)

    # Oracle MCC for each
    mcc_std  <- oracle_mcc(z_std, labels)
    mcc_efa  <- oracle_mcc(z_efa, labels)
    mcc_iter <- oracle_mcc(z_iter, labels)
    mcc_v2   <- oracle_mcc(z_v2, labels)

    results[[idx]] <- data.frame(
      nF = nF, ipf = ipf, total_items = total, rep = rep_id,
      std = round(mcc_std, 3),
      efa_d = round(mcc_efa, 3),
      iterative = round(mcc_iter, 3),
      weighted_v2 = round(mcc_v2, 3)
    )

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    eta <- if (idx > 1) (elapsed / idx) * (TOTAL_CONDITIONS - idx) else NA
    cat(sprintf("std=%.3f efa=%.3f iter=%.3f v2=%.3f | %.1f min ETA %.0f\n",
                mcc_std, mcc_efa, mcc_iter, mcc_v2, elapsed, eta))
  }
}

# ============================================================
# Aggregate & save
# ============================================================
df <- do.call(rbind, results)
write.csv(df, "sim_weighted_v2_results.csv", row.names = FALSE)

cat("\n============================================================\n")
cat(" RESULTS\n")
cat("============================================================\n\n")

cat("=== OVERALL MEAN MCC ===\n")
overall <- df %>%
  summarise(across(c(std, efa_d, iterative, weighted_v2),
                   ~ round(mean(.x, na.rm = TRUE), 3)))
print(as.data.frame(overall), row.names = FALSE)

cat("\n=== BY ITEMS BIN ===\n")
df_bin <- df %>%
  mutate(bin = cut(total_items,
                   breaks = c(0, 30, 60, 100, 200, 1000),
                   labels = c("<30", "30-60", "60-100", "100-200", ">200"))) %>%
  group_by(bin) %>%
  summarise(across(c(std, efa_d, iterative, weighted_v2),
                   ~ round(mean(.x, na.rm = TRUE), 3)),
            .groups = "drop")
print(as.data.frame(df_bin), row.names = FALSE)

cat("\n=== BY nF ===\n")
df_nF <- df %>%
  group_by(nF) %>%
  summarise(across(c(std, efa_d, iterative, weighted_v2),
                   ~ round(mean(.x, na.rm = TRUE), 3)),
            .groups = "drop")
print(as.data.frame(df_nF), row.names = FALSE)

cat("\n=== WINNERS per cell ===\n")
win_df <- df %>%
  group_by(nF, ipf, total_items) %>%
  summarise(std = round(mean(std, na.rm = TRUE), 3),
            efa_d = round(mean(efa_d, na.rm = TRUE), 3),
            iterative = round(mean(iterative, na.rm = TRUE), 3),
            weighted_v2 = round(mean(weighted_v2, na.rm = TRUE), 3),
            .groups = "drop") %>%
  rowwise() %>%
  mutate(winner = c("std","efa_d","iterative","weighted_v2")[
    which.max(c(std, efa_d, iterative, weighted_v2))]) %>%
  ungroup()
print(as.data.frame(win_df), row.names = FALSE)

cat(sprintf("\nTotal runtime: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("Results saved: sim_weighted_v2_results.csv\n")
cat("=== DONE ===\n")
