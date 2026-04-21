#### Test weighted_v2 (scale-invariant) on real validation datasets ####
#
# Proposed fix: replace rowCor_weighted (column-standardized cross-product sum,
# NOT scale-invariant) with rowCor_weighted_v2 (weighted Pearson correlation
# WITHIN-respondent, scale-invariant).
#
# Compares on Schroeders / Schneider / Niessen:
#   coupled        — current validated method (top-|r|, rowCor_abs)
#   weighted_v2    — NEW fix (all pairs weighted, within-respondent correlation)
#   weighted_OLD   — broken version (all pairs, column-z cross products) — sanity check
#
# Expected: weighted_v2 should
#   - Match weighted_OLD AUC on datasets where OLD didn't invert
#   - RECOVER the signal on Schneider (OLD: 0.304, target: >0.7)
#   - Handle straight-liners correctly (via scale invariance)

setwd("C:/Users/User/OneDrive - CNR/Claude/ReReReRe")
library(dplyr)
library(pROC)
library(haven)

# ============================================================
# SCORING HELPERS (all self-contained)
# ============================================================

rowCor_abs <- function(A, B, zero_val = 0) {
  A_mean <- rowMeans(A, na.rm = TRUE); B_mean <- rowMeans(B, na.rm = TRUE)
  A_c <- A - A_mean; B_c <- B - B_mean
  num <- rowSums(A_c * B_c, na.rm = TRUE)
  den <- sqrt(rowSums(A_c^2, na.rm = TRUE) * rowSums(B_c^2, na.rm = TRUE))
  r <- abs(num / den)
  r[is.nan(r) | is.na(r)] <- zero_val
  r
}

# NEW: weighted Pearson correlation within-respondent (scale-invariant)
rowCor_weighted_v2 <- function(A, B, weights, zero_val = 0) {
  sum_w <- sum(weights)
  W <- matrix(weights, nrow(A), ncol(A), byrow = TRUE)
  mean_A <- rowSums(A * W, na.rm = TRUE) / sum_w
  mean_B <- rowSums(B * W, na.rm = TRUE) / sum_w
  A_c <- A - mean_A; B_c <- B - mean_B
  num <- rowSums(W * A_c * B_c, na.rm = TRUE)
  den <- sqrt(rowSums(W * A_c^2, na.rm = TRUE) * rowSums(W * B_c^2, na.rm = TRUE))
  r <- abs(num / den)
  r[is.nan(r) | is.na(r)] <- zero_val
  r
}

# OLD broken: column-z cross-product sum (NOT scale-invariant within respondent)
rowCor_weighted_OLD <- function(A, B, weights) {
  N <- nrow(A); k <- ncol(A)
  A_z <- scale(A, center = TRUE, scale = TRUE)
  B_z <- scale(B, center = TRUE, scale = TRUE)
  A_z[is.na(A_z)] <- 0; B_z[is.na(B_z)] <- 0
  cp <- A_z * B_z
  W <- matrix(weights, N, k, byrow = TRUE)
  rowSums(cp * W, na.rm = TRUE) / sum(weights)
}

# ============================================================
# UNIFIED SCORING: all three methods in one function
# ============================================================

score <- function(data, method, corProp = 0.03, iterations = 100,
                  min_pairs = 15, align_signs = TRUE, seed = 42) {
  set.seed(seed)
  data <- data[, sapply(data, is.numeric), drop = FALSE]
  N <- nrow(data); J <- ncol(data)
  mat <- as.matrix(data)

  if (align_signs) item_max <- apply(mat, 2, max, na.rm = TRUE)

  rawCorMat <- cor(mat, use = "pairwise.complete.obs")
  corMat <- abs(rawCorMat)
  corMat[upper.tri(corMat, diag = TRUE)] <- NA
  rawCorMat[upper.tri(rawCorMat, diag = TRUE)] <- NA
  if (align_signs) signMat <- sign(rawCorMat)

  # Select pairs
  if (method == "coupled") {
    thr <- quantile(corMat, 1 - corProp, na.rm = TRUE)
    pair_idx <- which(corMat >= thr, arr.ind = TRUE)
    if (nrow(pair_idx) < min_pairs) {
      all_cors <- corMat[lower.tri(corMat)]
      all_cors <- all_cors[!is.na(all_cors)]
      thr <- sort(all_cors, decreasing = TRUE)[min(min_pairs, length(all_cors))]
      pair_idx <- which(corMat >= thr, arr.ind = TRUE)
    }
    weights <- NULL
  } else {
    pair_idx <- which(!is.na(corMat), arr.ind = TRUE)
    weights <- corMat[pair_idx]
  }

  k <- nrow(pair_idx)
  pair_signs <- sign(rawCorMat[pair_idx])

  A_sel <- mat[, pair_idx[, 1], drop = FALSE]
  B_sel <- mat[, pair_idx[, 2], drop = FALSE]

  if (align_signs) {
    nf <- which(pair_signs < 0)
    if (length(nf) > 0) {
      B_sel[, nf] <- rep(item_max[pair_idx[nf, 2]] + 1, each = N) - B_sel[, nf]
    }
  }

  score_fn <- switch(method,
                     coupled     = function(A, B) rowCor_abs(A, B, zero_val = 0),
                     weighted_v2 = function(A, B) rowCor_weighted_v2(A, B, weights, zero_val = 0),
                     weighted_OLD = function(A, B) rowCor_weighted_OLD(A, B, weights))

  rowCors <- score_fn(A_sel, B_sel)

  # Permutation baseline
  all_RIC <- matrix(NA_real_, nrow = N, ncol = iterations)
  all_pairs_comb <- combn(J, 2)
  n_pairs_total <- ncol(all_pairs_comb)

  for (i in seq_len(iterations)) {
    if (method == "coupled") {
      samp <- sample(n_pairs_total, k)
      idx1 <- all_pairs_comb[1, samp]; idx2 <- all_pairs_comb[2, samp]
    } else {
      idx1 <- sample(J, k, replace = TRUE)
      idx2 <- sample(J, k, replace = TRUE)
      same <- idx1 == idx2
      while (any(same)) {
        idx2[same] <- sample(J, sum(same), replace = TRUE)
        same <- idx1 == idx2
      }
    }
    A_r <- mat[, idx1, drop = FALSE]; B_r <- mat[, idx2, drop = FALSE]
    if (align_signs) {
      ri <- pmax(idx1, idx2); ci <- pmin(idx1, idx2)
      rs <- signMat[cbind(ri, ci)]; rs[is.na(rs)] <- 1
      rn <- which(rs < 0)
      if (length(rn) > 0) {
        B_r[, rn] <- rep(item_max[idx2[rn]] + 1, each = N) - B_r[, rn]
      }
    }
    all_RIC[, i] <- score_fn(A_r, B_r)
  }

  rm <- rowMeans(all_RIC, na.rm = TRUE)
  rs <- apply(all_RIC, 1, sd, na.rm = TRUE)
  z  <- ifelse(rs > 0, (rowCors - rm) / rs, 0)

  list(indCors = rowCors, rand_mean = rm, rand_sd = rs, z_score = z, k = k)
}

# ============================================================
# EVAL HELPERS
# ============================================================

compute_mcc <- function(tp, tn, fp, fn) {
  tp <- as.double(tp); tn <- as.double(tn); fp <- as.double(fp); fn <- as.double(fn)
  den <- sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  if (is.na(den) || den == 0) return(0)
  (tp*tn - fp*fn) / den
}

eval_z <- function(z, labels, z_thr) {
  flagged <- as.integer(z <= z_thr)
  tp <- sum(flagged==1 & labels==1); tn <- sum(flagged==0 & labels==0)
  fp <- sum(flagged==1 & labels==0); fn <- sum(flagged==0 & labels==1)
  list(mcc = compute_mcc(tp, tn, fp, fn),
       sens = if ((tp+fn)>0) tp/(tp+fn) else NA,
       spec = if ((tn+fp)>0) tn/(tn+fp) else NA)
}

find_best_z <- function(z, labels, zr = seq(0.1, 3.0, 0.1)) {
  mccs <- sapply(zr, function(t) eval_z(z, labels, t)$mcc)
  list(mcc = max(mccs, na.rm = TRUE), z = zr[which.max(mccs)])
}

compute_auc <- function(score, labels) {
  tryCatch(as.numeric(auc(roc(labels, score, direction = ">", quiet = TRUE))),
           error = function(e) NA)
}

run_all_methods <- function(items, labels, dataset_name, seed = 42) {
  cat(sprintf("\n=== %s (N=%d, items=%d, careless=%d = %.1f%%) ===\n",
              dataset_name, nrow(items), ncol(items),
              sum(labels), 100*mean(labels)))

  methods <- c("coupled", "weighted_v2", "weighted_OLD")
  out <- lapply(methods, function(m) {
    r <- score(items, method = m, seed = seed)
    auc_v <- compute_auc(r$z_score, labels)
    e15 <- eval_z(r$z_score, labels, 1.5)
    best <- find_best_z(r$z_score, labels)

    # Also inspect: indCors / rand_mean good vs careless (diagnostic)
    g <- labels == 0; c_ <- labels == 1
    d_indCors   <- mean(r$indCors[g])    - mean(r$indCors[c_])
    d_rand_mean <- mean(r$rand_mean[g])  - mean(r$rand_mean[c_])
    d_lift <- (mean(r$indCors[g])  - mean(r$rand_mean[g])) -
              (mean(r$indCors[c_]) - mean(r$rand_mean[c_]))

    data.frame(
      Method = m, k = r$k,
      AUC = round(auc_v, 3),
      MCC_z15 = round(e15$mcc, 3),
      Sens_z15 = round(e15$sens, 3),
      Spec_z15 = round(e15$spec, 3),
      MCC_oracle = round(best$mcc, 3),
      Oracle_z = best$z,
      Delta_indCors_goodVSbad = round(d_indCors, 3),
      Delta_lift_goodVSbad = round(d_lift, 3)
    )
  })
  do.call(rbind, out)
}

# ============================================================
# DATASETS
# ============================================================

cat("\n=========================================================\n")
cat("  TEST weighted_v2 (scale-invariant) on real datasets\n")
cat("=========================================================\n")

all_results <- list()

# ---- Schroeders ----
d1 <- read.csv("external_datasets/01_Schroeders_2022/data_mod_resp.csv", sep=";")
lab1 <- d1$Careless
it1 <- d1 %>% select(starts_with("HE"))
ok1 <- complete.cases(it1)
it1 <- it1[ok1,]; lab1 <- lab1[ok1]
all_results[[1]] <- cbind(Dataset = "Schroeders2022",
                          run_all_methods(it1, lab1, "Schroeders 2022"))

# ---- Schneider ----
d2 <- read.csv("external_datasets/06_Schneider_QoL/carersp.csv")
lab2 <- d2$c01
cols2 <- grep("^(dep|pain|cog|fat)\\d", names(d2), value = TRUE)
it2 <- d2[, cols2]
ok2 <- complete.cases(it2) & !is.na(lab2)
it2 <- it2[ok2,]; lab2 <- lab2[ok2]
all_results[[2]] <- cbind(Dataset = "Schneider_QoL",
                          run_all_methods(it2, lab2, "Schneider QoL"))

# ---- Niessen ----
d3 <- read_sav("external_datasets/08_Niessen_2016/Raw_data.sav")
d3 <- d3 %>% filter(Use_me == 1)
lab3 <- as.integer(d3$Conditon)
cols3 <- grep("^[ECANO]\\d+$", names(d3), value = TRUE)
it3 <- as.data.frame(d3[, cols3])
it3 <- data.frame(lapply(it3, as.numeric))
ok3 <- complete.cases(it3)
it3 <- it3[ok3,]; lab3 <- lab3[ok3]
all_results[[3]] <- cbind(Dataset = "Niessen2016",
                          run_all_methods(it3, lab3, "Niessen 2016"))

# ============================================================
# SUMMARY
# ============================================================

cat("\n\n=========================================================\n")
cat("  SUMMARY TABLE\n")
cat("=========================================================\n\n")

summary_df <- do.call(rbind, all_results)
print(summary_df, row.names = FALSE)
write.csv(summary_df, "test_weighted_v2_results.csv", row.names = FALSE)

cat("\nSaved: test_weighted_v2_results.csv\n")

# Headline comparison per dataset
cat("\n--- Headline: AUC by method ---\n")
for (ds in unique(summary_df$Dataset)) {
  sub <- summary_df[summary_df$Dataset == ds, ]
  cat(sprintf("%-20s  coupled=%.3f  weighted_v2=%.3f  weighted_OLD=%.3f\n",
              ds, sub$AUC[sub$Method == "coupled"],
              sub$AUC[sub$Method == "weighted_v2"],
              sub$AUC[sub$Method == "weighted_OLD"]))
}

cat("\n=== DONE ===\n")
