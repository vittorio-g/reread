#### Diagnose: Why does weighted mode invert on Schneider? ####
# Schneider: 31 items × 4 domains, 6.6% careless (latent class)
# Observed: RR_auto (weighted) AUC=0.304, RR_coupled AUC=0.730
# Expected: both > 0.5 in the same direction

setwd("C:/Users/User/OneDrive - CNR/Claude/ReReReRe")
library(dplyr)
library(pROC)
source("ReReReRe.R")

# --- Load Schneider ---
d2 <- read.csv("external_datasets/06_Schneider_QoL/carersp.csv")
labels <- d2$c01
item_cols <- grep("^(dep|pain|cog|fat)\\d", names(d2), value=TRUE)
items <- d2[, item_cols]
ok <- complete.cases(items) & !is.na(labels)
items <- items[ok,]; labels <- labels[ok]

N <- nrow(items); J <- ncol(items)
cat(sprintf("N=%d items=%d careless=%d\n\n", N, J, sum(labels)))

# --- Domain map ---
dom <- sub("\\d+$", "", item_cols)
cat("Domain counts:\n"); print(table(dom))

# --- Item correlation matrix structure ---
R <- cor(items, use="pairwise.complete.obs")
cat("\nCorrelation summary (lower triangle):\n")
r_low <- R[lower.tri(R)]
cat(sprintf("  min=%.3f  max=%.3f  median=%.3f  mean=%.3f\n",
  min(r_low), max(r_low), median(r_low), mean(r_low)))
cat(sprintf("  |r| < 0.1: %d (%.0f%%)\n", sum(abs(r_low)<0.1), 100*mean(abs(r_low)<0.1)))
cat(sprintf("  r < 0 (negative): %d (%.0f%%)\n", sum(r_low<0), 100*mean(r_low<0)))

# Within-domain vs between-domain
pair_mat <- which(lower.tri(R), arr.ind=TRUE)
dom_A <- dom[pair_mat[,1]]; dom_B <- dom[pair_mat[,2]]
within <- dom_A == dom_B
cat(sprintf("\nWithin-domain pairs: %d, mean |r|=%.3f\n", sum(within), mean(abs(r_low[within]))))
cat(sprintf("Between-domain pairs: %d, mean |r|=%.3f\n", sum(!within), mean(abs(r_low[!within]))))

# --- Run weighted and coupled ---
set.seed(42)
cat("\n--- Running WEIGHTED mode ---\n")
rr_w <- ReReReRe(items, iterations=100, align_signs=TRUE, mode="weighted")

set.seed(42)
cat("--- Running COUPLED mode ---\n")
rr_c <- ReReReRe(items, iterations=100, align_signs=TRUE, mode="coupled", corProp=0.03)

# --- Compare indCors / rand_mean / z_score distributions by label ---
compare_dists <- function(rr, name, lab) {
  good <- lab == 0; car <- lab == 1
  cat(sprintf("\n[%s]\n", name))
  cat(sprintf("  indCors:   good mean=%.4f sd=%.4f | careless mean=%.4f sd=%.4f | Δ=%.4f\n",
      mean(rr$indCors[good]), sd(rr$indCors[good]),
      mean(rr$indCors[car]),  sd(rr$indCors[car]),
      mean(rr$indCors[good]) - mean(rr$indCors[car])))
  cat(sprintf("  rand_mean: good mean=%.4f sd=%.4f | careless mean=%.4f sd=%.4f | Δ=%.4f\n",
      mean(rr$rand_mean[good]), sd(rr$rand_mean[good]),
      mean(rr$rand_mean[car]),  sd(rr$rand_mean[car]),
      mean(rr$rand_mean[good]) - mean(rr$rand_mean[car])))
  cat(sprintf("  rand_sd:   good mean=%.4f | careless mean=%.4f\n",
      mean(rr$rand_sd[good]), mean(rr$rand_sd[car])))
  cat(sprintf("  z_score:   good mean=%.3f sd=%.3f | careless mean=%.3f sd=%.3f | Δ=%.3f\n",
      mean(rr$z_score[good]), sd(rr$z_score[good]),
      mean(rr$z_score[car]),  sd(rr$z_score[car]),
      mean(rr$z_score[good]) - mean(rr$z_score[car])))
  auc_gt <- tryCatch(as.numeric(auc(roc(lab, rr$z_score, direction=">", quiet=TRUE))), error=function(e) NA)
  auc_lt <- tryCatch(as.numeric(auc(roc(lab, rr$z_score, direction="<", quiet=TRUE))), error=function(e) NA)
  cat(sprintf("  AUC z_score: direction '>' = %.3f   direction '<' = %.3f\n", auc_gt, auc_lt))
  auc_i_gt <- tryCatch(as.numeric(auc(roc(lab, rr$indCors, direction=">", quiet=TRUE))), error=function(e) NA)
  cat(sprintf("  AUC indCors: direction '>' = %.3f\n", auc_i_gt))
}
compare_dists(rr_w, "WEIGHTED", labels)
compare_dists(rr_c, "COUPLED",  labels)

# --- Examine the weight distribution in weighted mode ---
# Replicate the pair-selection logic from weighted mode
mat <- as.matrix(items)
rawCor <- cor(mat, use="pairwise.complete.obs")
absCor <- abs(rawCor); absCor[upper.tri(absCor, diag=TRUE)] <- NA
rawCor_lt <- rawCor; rawCor_lt[upper.tri(rawCor_lt, diag=TRUE)] <- NA

pair_idx <- which(!is.na(absCor), arr.ind=TRUE)
pair_abs_r <- absCor[pair_idx]
pair_sign  <- sign(rawCor_lt[pair_idx])

cat(sprintf("\nWeighted-mode pair pool: %d pairs\n", length(pair_abs_r)))
cat(sprintf("  |r| quartiles: %s\n",
    paste(sprintf("%.3f", quantile(pair_abs_r, c(.25,.5,.75,.9,.99))), collapse=" / ")))
cat(sprintf("  Sum of all weights: %.2f\n", sum(pair_abs_r)))
cat(sprintf("  Top 3%% weight share: %.1f%%\n",
    100 * sum(sort(pair_abs_r, decreasing=TRUE)[1:round(length(pair_abs_r)*.03)]) / sum(pair_abs_r)))
cat(sprintf("  Negative-sign pairs: %d (%.1f%%)\n",
    sum(pair_sign < 0), 100*mean(pair_sign < 0)))

# --- Key diagnostic: what happens in weighted mode permutation? ---
# In weighted mode, random pairs are drawn from J items, but weights are
# assigned in order (weights[i] applies to random pair i, not to that pair's actual |r|).
# Let's see if this mismatch creates a biased baseline.

set.seed(99)
k <- nrow(pair_idx)
weights <- pair_abs_r
cat(sprintf("\nRe-running permutation diagnostic (k=%d)...\n", k))

# Pick one iteration manually and inspect
rand_idx1 <- sample(J, k, replace=TRUE)
rand_idx2 <- sample(J, k, replace=TRUE)
same <- rand_idx1 == rand_idx2
while (any(same)) { rand_idx2[same] <- sample(J, sum(same), replace=TRUE); same <- rand_idx1 == rand_idx2 }

# For each random pair, what's its ACTUAL |r| vs its ASSIGNED weight?
actual_r <- sapply(seq_len(k), function(j) abs(rawCor[rand_idx1[j], rand_idx2[j]]))
assigned_w <- weights
cat(sprintf("  Random pairs: actual mean|r|=%.3f  assigned mean weight=%.3f\n",
    mean(actual_r, na.rm=TRUE), mean(assigned_w)))
cat(sprintf("  Correlation between actual and assigned: %.3f\n",
    cor(actual_r, assigned_w, use="complete.obs")))

# --- What if weights were aligned with actual random-pair |r|? ---
# Simulate: compute weighted coherence with weights = actual |r| of each random pair
# (This is what a CORRECT null would do)

A_rand <- mat[, rand_idx1, drop=FALSE]
B_rand <- mat[, rand_idx2, drop=FALSE]

A_z <- scale(A_rand, center=TRUE, scale=TRUE)
B_z <- scale(B_rand, center=TRUE, scale=TRUE)
A_z[is.na(A_z)] <- 0; B_z[is.na(B_z)] <- 0
cp <- A_z * B_z

# Version 1: current algorithm (weights = coupled weights, misassigned)
w_current <- matrix(assigned_w, N, k, byrow=TRUE)
score_current <- rowSums(cp * w_current) / sum(assigned_w)

# Version 2: corrected (weights = actual |r| of random pair)
w_corrected <- matrix(actual_r, N, k, byrow=TRUE)
score_corrected <- rowSums(cp * w_corrected, na.rm=TRUE) / sum(actual_r, na.rm=TRUE)

cat("\nOne permutation iteration baseline comparison:\n")
cat(sprintf("  Current (misassigned weights): good=%.4f  careless=%.4f  diff=%.4f\n",
    mean(score_current[labels==0]), mean(score_current[labels==1]),
    mean(score_current[labels==0]) - mean(score_current[labels==1])))
cat(sprintf("  Corrected (true |r| weights):  good=%.4f  careless=%.4f  diff=%.4f\n",
    mean(score_corrected[labels==0]), mean(score_corrected[labels==1]),
    mean(score_corrected[labels==0]) - mean(score_corrected[labels==1])))

# --- Coupled indCors vs weighted indCors: what's different? ---
cat("\nLow-|r| pair contribution (weighted includes all 465 pairs, coupled only top ~14):\n")

# Take good respondent #1 and careless respondent #1
good_id <- which(labels==0)[1]; car_id <- which(labels==1)[1]
for (rid in c(good_id, car_id)) {
  # Row-level cross products weighted
  lab <- if (labels[rid]==0) "GOOD" else "CARELESS"
  z_A <- scale(mat[, pair_idx[,1], drop=FALSE], center=TRUE, scale=TRUE)[rid, ]
  z_B <- scale(mat[, pair_idx[,2], drop=FALSE], center=TRUE, scale=TRUE)[rid, ]
  cp_row <- z_A * z_B * sign(rawCor_lt[pair_idx])  # align signs
  # Sort pairs by |r|
  ord <- order(pair_abs_r, decreasing=TRUE)
  top14 <- ord[1:14]    # coupled set
  bottom_rest <- ord[15:length(ord)]  # excluded by coupled

  cat(sprintf("\n  Respondent %d (%s):\n", rid, lab))
  cat(sprintf("    TOP-14 pairs  (|r| %.2f-%.2f): mean cp=%.3f, sum(cp*w)=%.2f\n",
      pair_abs_r[top14[14]], pair_abs_r[top14[1]],
      mean(cp_row[top14], na.rm=TRUE),
      sum(cp_row[top14] * pair_abs_r[top14], na.rm=TRUE)))
  cat(sprintf("    REST pairs    (|r| %.2f-%.2f): mean cp=%.3f, sum(cp*w)=%.2f\n",
      min(pair_abs_r[bottom_rest]), max(pair_abs_r[bottom_rest]),
      mean(cp_row[bottom_rest], na.rm=TRUE),
      sum(cp_row[bottom_rest] * pair_abs_r[bottom_rest], na.rm=TRUE)))
}

cat("\n=== DONE ===\n")
