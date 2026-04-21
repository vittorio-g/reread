#### Pilot: weighted_v2 vs coupled on simulated short questionnaires ####
#
# Goal: quickly verify weighted_v2 does NOT crash on simulated data before
# committing to a full re-sim.
#
# Decision rule (rough):
#   - PASS if weighted_v2 MCC within 0.05 of coupled on most cells AND no AUC inversion
#   - FAIL if weighted_v2 drops ≥0.10 MCC on multiple cells, or AUC < 0.5 anywhere
#
# Small grid: 3 nF × 3 ipf × 1 N × 1 pct × 10 reps = 90 cells per method.
# Runtime: ~5-15 min.
# ============================================================

setwd("C:/Users/User/OneDrive - CNR/Claude/ReReReRe")
library(dplyr); library(pROC)
source("Synthetic_Good_Responses_2.R")
source("Careless_machine_2.R")
source("ReReReRe.R")   # coupled-only; used for the baseline

# ============================================================
# weighted_v2 inline (same as Test_Weighted_V2 and the full partial script)
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
# Pilot grid
# ============================================================

grid <- expand.grid(
  nF  = c(3, 8, 16),     # low, mid, high within ≤60 items zone
  ipf = c(3, 5, 10),     # low, mid, high items/factor
  stringsAsFactors = FALSE
)
# Filter to ≤60 items only
grid$total <- grid$nF * grid$ipf
grid <- grid[grid$total <= 60, ]
N_FIXED   <- 300
PCT       <- 0.10
REPS      <- 10
ITER      <- 100
SEED_BASE <- 2026
MAX_ITEMS <- 60

CARELESS_TYPES  <- c(random = 1/3, longstring = 1/3, mixed = 1/3)
CARELESS_LEVELS <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5)
Z_GRID <- seq(0.1, 3.0, 0.1)

cat("============================================================\n")
cat(" PILOT weighted_v2 vs coupled — simulated data sanity check\n")
cat("============================================================\n")
cat(sprintf("Grid: %d cells (≤%d items) × %d reps = %d scoring calls per method\n",
            nrow(grid), MAX_ITEMS, REPS, nrow(grid) * REPS))
print(grid[, c("nF", "ipf", "total")], row.names = FALSE)

# ============================================================
# Helpers
# ============================================================
compute_mcc <- function(z, lab, zt) {
  flagged <- z <= zt
  TP <- sum(flagged & lab); TN <- sum(!flagged & !lab)
  FP <- sum(flagged & !lab); FN <- sum(!flagged & lab)
  den <- sqrt(as.numeric(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  if (den == 0) 0 else (TP*TN - FP*FN) / den
}

best_mcc <- function(z, lab) {
  mccs <- sapply(Z_GRID, function(zt) compute_mcc(z, lab, zt))
  list(mcc = max(mccs, na.rm = TRUE), z = Z_GRID[which.max(mccs)])
}

auc_z <- function(z, lab) {
  tryCatch(as.numeric(auc(roc(lab, z, direction = ">", quiet = TRUE))),
           error = function(e) NA)
}

# ============================================================
# Main loop
# ============================================================
results <- list(); idx <- 0
t0 <- Sys.time()

for (g in seq_len(nrow(grid))) {
  nF <- grid$nF[g]; ipf <- grid$ipf[g]; total <- nF * ipf
  cat(sprintf("\n--- nF=%d ipf=%d (items=%d) ---\n", nF, ipf, total))

  for (rep_id in 1:REPS) {
    cell_seed <- SEED_BASE + (rep_id - 1) * 576 + g   # arbitrary, just consistent
    set.seed(cell_seed)

    clean_data <- simulated_good_responses(nConstructs = nF,
                                            nItems = rep(ipf, nF),
                                            n = N_FIXED)
    inj <- inject_careless(clean_data, pct_careless = PCT,
                           pct_types = CARELESS_TYPES,
                           careless_levels = CARELESS_LEVELS)
    labels <- inj$labels$careless_pct >= 0.01

    # coupled (current default)
    rr_c <- ReReReRe(inj$data_corrupted, corProp = 0.03, iterations = ITER,
                     align_signs = TRUE)
    z_c <- rr_c$z_score

    # weighted_v2
    z_v2 <- score_weighted_v2(inj$data_corrupted, iterations = ITER, align_signs = TRUE)

    # Metrics
    b_c  <- best_mcc(z_c,  labels); a_c  <- auc_z(z_c,  labels)
    b_v2 <- best_mcc(z_v2, labels); a_v2 <- auc_z(z_v2, labels)
    m15_c  <- compute_mcc(z_c,  labels, 1.5)
    m15_v2 <- compute_mcc(z_v2, labels, 1.5)

    idx <- idx + 1
    results[[idx]] <- data.frame(
      nF = nF, ipf = ipf, total = total, rep = rep_id,
      AUC_coupled = round(a_c, 3),   AUC_v2     = round(a_v2, 3),
      MCC15_coupled = round(m15_c, 3), MCC15_v2 = round(m15_v2, 3),
      MCCoracle_coupled = round(b_c$mcc, 3), oracle_z_c = b_c$z,
      MCCoracle_v2      = round(b_v2$mcc, 3), oracle_z_v2 = b_v2$z
    )
  }

  sub <- do.call(rbind, results[(idx - REPS + 1):idx])
  cat(sprintf("  AUC:         coupled %.3f ± %.3f   v2 %.3f ± %.3f\n",
      mean(sub$AUC_coupled), sd(sub$AUC_coupled),
      mean(sub$AUC_v2),      sd(sub$AUC_v2)))
  cat(sprintf("  MCC (z=1.5): coupled %.3f ± %.3f   v2 %.3f ± %.3f\n",
      mean(sub$MCC15_coupled), sd(sub$MCC15_coupled),
      mean(sub$MCC15_v2),      sd(sub$MCC15_v2)))
  cat(sprintf("  MCC oracle:  coupled %.3f            v2 %.3f\n",
      mean(sub$MCCoracle_coupled), mean(sub$MCCoracle_v2)))
}

elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("\nElapsed: %.1f min\n", elapsed))

# ============================================================
# Aggregate and decide
# ============================================================
df <- do.call(rbind, results)
write.csv(df, "pilot_weighted_v2_results.csv", row.names = FALSE)

agg <- df %>%
  group_by(nF, ipf, total) %>%
  summarise(across(c(AUC_coupled, AUC_v2, MCC15_coupled, MCC15_v2,
                     MCCoracle_coupled, MCCoracle_v2),
                   ~ round(mean(.), 3)),
            .groups = "drop")
agg$AUC_delta_v2_minus_c   <- round(agg$AUC_v2 - agg$AUC_coupled, 3)
agg$MCCoracle_delta         <- round(agg$MCCoracle_v2 - agg$MCCoracle_coupled, 3)
agg$MCC15_delta             <- round(agg$MCC15_v2 - agg$MCC15_coupled, 3)

cat("\n===========================================================\n")
cat(" AGGREGATE BY CELL (mean over", REPS, "reps)\n")
cat("===========================================================\n")
print(as.data.frame(agg), row.names = FALSE)

# ============================================================
# Verdict
# ============================================================
cat("\n===========================================================\n")
cat(" VERDICT\n")
cat("===========================================================\n")

n_cells         <- nrow(agg)
n_auc_inverted  <- sum(agg$AUC_v2 < 0.5)
n_bad_mcc       <- sum(agg$MCCoracle_delta < -0.10)
n_close_mcc     <- sum(abs(agg$MCCoracle_delta) <= 0.05)
n_v2_wins       <- sum(agg$MCCoracle_delta > 0.05)

cat(sprintf("Cells: %d\n", n_cells))
cat(sprintf("  AUC_v2 < 0.5 (inversion):           %d / %d\n", n_auc_inverted, n_cells))
cat(sprintf("  MCC oracle drop > 0.10 vs coupled:  %d / %d\n", n_bad_mcc, n_cells))
cat(sprintf("  MCC oracle within ±0.05 of coupled: %d / %d\n", n_close_mcc, n_cells))
cat(sprintf("  MCC oracle beats coupled by >0.05:  %d / %d\n", n_v2_wins, n_cells))

verdict <- "UNCLEAR"
if (n_auc_inverted == 0 && n_bad_mcc == 0) {
  verdict <- "PASS — safe to proceed to full partial weighted_v2 simulation"
} else if (n_auc_inverted > 0 || n_bad_mcc >= 3) {
  verdict <- "FAIL — weighted_v2 degrades on simulated data; do NOT run full sim"
}

cat(sprintf("\nVERDICT: %s\n", verdict))
cat("\nResults saved: pilot_weighted_v2_results.csv\n")
cat("=== DONE ===\n")
