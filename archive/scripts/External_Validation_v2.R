#### External Validation v2 — Updated defaults (2026-03-28) ####
# Uses new ReReReRe defaults: corProp=0.03, tests both auto_z and fixed z=1.5
# Compares against practical Mahalanobis (chi-sq .001)

setwd("C:/Users/vitto/Desktop/ReReReRe")
library(dplyr)
library(pROC)
source("ReReReRe.R")

# ============================================================
# Helpers
# ============================================================
compute_mcc <- function(tp, tn, fp, fn) {
  tp <- as.double(tp); tn <- as.double(tn)
  fp <- as.double(fp); fn <- as.double(fn)
  num <- (tp * tn) - (fp * fn)
  den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  if (is.na(den) || den == 0) return(0)
  num / den
}

evaluate <- function(predicted, true_label) {
  tp <- sum(predicted == 1 & true_label == 1, na.rm=TRUE)
  tn <- sum(predicted == 0 & true_label == 0, na.rm=TRUE)
  fp <- sum(predicted == 1 & true_label == 0, na.rm=TRUE)
  fn <- sum(predicted == 0 & true_label == 1, na.rm=TRUE)
  mcc <- compute_mcc(tp, tn, fp, fn)
  sens <- if ((tp+fn)>0) tp/(tp+fn) else NA
  spec <- if ((tn+fp)>0) tn/(tn+fp) else NA
  list(MCC=mcc, Sens=sens, Spec=spec, TP=tp, TN=tn, FP=fp, FN=fn)
}

compute_auc <- function(score, true_label, direction) {
  tryCatch({
    roc_obj <- roc(true_label, score, direction=direction, quiet=TRUE)
    as.numeric(auc(roc_obj))
  }, error = function(e) NA)
}

run_mah <- function(data_items, alpha=0.001) {
  mat <- as.matrix(data_items)
  p <- ncol(mat)
  complete_idx <- complete.cases(mat)
  center <- colMeans(mat[complete_idx,], na.rm=TRUE)
  cov_mat <- cov(mat[complete_idx,], use="pairwise.complete.obs")
  if (any(is.na(cov_mat)) || det(cov_mat) < 1e-10) cov_mat <- cov_mat + diag(0.01, p)
  d2 <- mahalanobis(mat, center, cov_mat)
  threshold <- qchisq(1-alpha, df=p)
  list(d2=d2, flagged=d2>threshold)
}

# Find best z for RR (oracle)
find_best_z <- function(z_scores, labels, z_range=seq(0.1, 3.0, 0.1)) {
  best_mcc <- -Inf; best_z <- 1.5
  for (z in z_range) {
    flagged <- as.integer(z_scores <= z)
    ev <- evaluate(flagged, labels)
    if (ev$MCC > best_mcc) { best_mcc <- ev$MCC; best_z <- z }
  }
  list(z=best_z, mcc=best_mcc)
}

cat("\n====================================================\n")
cat("   EXTERNAL VALIDATION v2 (2026-03-28)\n")
cat("   New defaults: corProp=0.03\n")
cat("   Tests: auto_z, fixed z=1.5, oracle best z\n")
cat("   Benchmark: Mahalanobis chi-sq(.001)\n")
cat("====================================================\n\n")

results <- list()

# ============================================================
# DATASET 1: Schroeders 2022 (HEXACO-60, experimental)
# ============================================================
cat("=== 1. Schroeders et al. 2022 ===\n")
cat("  60 items, ~6 factors (HEXACO), experimental induction\n")

d1 <- read.csv("external_datasets/01_Schroeders_2022/data_mod_resp.csv", sep=";")
d1_label <- d1$Careless
d1_items <- d1 %>% select(starts_with("HE"))
complete1 <- complete.cases(d1_items)
d1_items <- d1_items[complete1,]; d1_label <- d1_label[complete1]
cat(sprintf("  N=%d, items=%d, careless=%d (%.1f%%)\n",
    nrow(d1_items), ncol(d1_items), sum(d1_label), mean(d1_label)*100))

# Run ReReReRe with new defaults
rr1 <- ReReReRe(d1_items, corProp=0.03, iterations=100, align_signs=TRUE)
rr1_auto <- ReReReRe(d1_items, corProp=0.03, iterations=100, align_signs=TRUE, auto_z=TRUE)
mah1 <- run_mah(d1_items)

# Evaluate at multiple thresholds
best1 <- find_best_z(rr1$z_score, d1_label)
ev1_15 <- evaluate(as.integer(rr1$z_score <= 1.5), d1_label)
ev1_auto <- evaluate(rr1_auto$flagged, d1_label)
ev1_mah <- evaluate(mah1$flagged, d1_label)

auc1_rr <- compute_auc(rr1$z_score, d1_label, ">")
auc1_mah <- compute_auc(mah1$d2, d1_label, "<")

cat(sprintf("  RR z=1.5:    MCC=%.3f  Sens=%.3f  Spec=%.3f\n", ev1_15$MCC, ev1_15$Sens, ev1_15$Spec))
cat(sprintf("  RR auto-z:   MCC=%.3f  Sens=%.3f  Spec=%.3f  (z=%.2f)\n",
    ev1_auto$MCC, ev1_auto$Sens, ev1_auto$Spec, rr1_auto$z_threshold_used[1]))
cat(sprintf("  RR oracle:   MCC=%.3f  (z=%.1f)\n", best1$mcc, best1$z))
cat(sprintf("  Mahalanobis: MCC=%.3f  Sens=%.3f  Spec=%.3f\n", ev1_mah$MCC, ev1_mah$Sens, ev1_mah$Spec))
cat(sprintf("  AUC: RR=%.3f, Mah=%.3f\n", auc1_rr, auc1_mah))

results[[1]] <- data.frame(
  Dataset="Schroeders2022", Items=ncol(d1_items), nF_approx=6, N=nrow(d1_items),
  Pct_careless=round(mean(d1_label)*100,1), GT_type="Experimental",
  RR_AUC=round(auc1_rr,3), Mah_AUC=round(auc1_mah,3),
  RR_MCC_z15=round(ev1_15$MCC,3), RR_MCC_auto=round(ev1_auto$MCC,3),
  RR_MCC_oracle=round(best1$mcc,3), RR_oracle_z=best1$z,
  Mah_MCC=round(ev1_mah$MCC,3), Auto_z=round(rr1_auto$z_threshold_used[1],2))

# ============================================================
# DATASET 2: Schneider QoL (31 items, 4 factors, latent class)
# ============================================================
cat("\n=== 2. Schneider et al. (QoL) ===\n")
cat("  31 items, 4 factors, latent class ground truth\n")

d2 <- read.csv("external_datasets/06_Schneider_QoL/carersp.csv")
d2_label <- d2$c01
item_cols2 <- grep("^(dep|pain|cog|fat)\\d", names(d2), value=TRUE)
d2_items <- d2[, item_cols2]
complete2 <- complete.cases(d2_items) & !is.na(d2_label)
d2_items <- d2_items[complete2,]; d2_label <- d2_label[complete2]
cat(sprintf("  N=%d, items=%d, careless=%d (%.1f%%)\n",
    nrow(d2_items), ncol(d2_items), sum(d2_label), mean(d2_label)*100))

rr2 <- ReReReRe(d2_items, corProp=0.03, iterations=100, align_signs=TRUE)
rr2_auto <- ReReReRe(d2_items, corProp=0.03, iterations=100, align_signs=TRUE, auto_z=TRUE)
mah2 <- run_mah(d2_items)

best2 <- find_best_z(rr2$z_score, d2_label)
ev2_15 <- evaluate(as.integer(rr2$z_score <= 1.5), d2_label)
ev2_auto <- evaluate(rr2_auto$flagged, d2_label)
ev2_mah <- evaluate(mah2$flagged, d2_label)

auc2_rr <- compute_auc(rr2$z_score, d2_label, ">")
auc2_mah <- compute_auc(mah2$d2, d2_label, "<")

cat(sprintf("  RR z=1.5:    MCC=%.3f  Sens=%.3f  Spec=%.3f\n", ev2_15$MCC, ev2_15$Sens, ev2_15$Spec))
cat(sprintf("  RR auto-z:   MCC=%.3f  Sens=%.3f  Spec=%.3f  (z=%.2f)\n",
    ev2_auto$MCC, ev2_auto$Sens, ev2_auto$Spec, rr2_auto$z_threshold_used[1]))
cat(sprintf("  RR oracle:   MCC=%.3f  (z=%.1f)\n", best2$mcc, best2$z))
cat(sprintf("  Mahalanobis: MCC=%.3f  Sens=%.3f  Spec=%.3f\n", ev2_mah$MCC, ev2_mah$Sens, ev2_mah$Spec))
cat(sprintf("  AUC: RR=%.3f, Mah=%.3f\n", auc2_rr, auc2_mah))

results[[2]] <- data.frame(
  Dataset="Schneider_QoL", Items=ncol(d2_items), nF_approx=4, N=nrow(d2_items),
  Pct_careless=round(mean(d2_label)*100,1), GT_type="Latent class",
  RR_AUC=round(auc2_rr,3), Mah_AUC=round(auc2_mah,3),
  RR_MCC_z15=round(ev2_15$MCC,3), RR_MCC_auto=round(ev2_auto$MCC,3),
  RR_MCC_oracle=round(best2$mcc,3), RR_oracle_z=best2$z,
  Mah_MCC=round(ev2_mah$MCC,3), Auto_z=round(rr2_auto$z_threshold_used[1],2))

# ============================================================
# DATASET 3: Niessen 2016 (IPIP-100, 5 factors, speed)
# ============================================================
cat("\n=== 3. Niessen et al. 2016 ===\n")
cat("  100 items, 5 factors (Big Five), speed manipulation\n")

library(haven)
d3 <- read_sav("external_datasets/08_Niessen_2016/Raw_data.sav")
d3 <- d3 %>% filter(Use_me == 1)
d3_label <- as.integer(d3$Conditon)  # 1 = speed condition (careless)
# Items are prefixed E/C/A/N/O with numbers 1-100 (interleaved Big Five)
item_cols3 <- grep("^[ECANO]\\d+$", names(d3), value=TRUE)
d3_items <- as.data.frame(d3[, item_cols3])
d3_items <- data.frame(lapply(d3_items, as.numeric))
complete3 <- complete.cases(d3_items)
d3_items <- d3_items[complete3,]; d3_label <- d3_label[complete3]
cat(sprintf("  N=%d, items=%d, careless=%d (%.1f%%)\n",
    nrow(d3_items), ncol(d3_items), sum(d3_label), mean(d3_label)*100))

rr3 <- ReReReRe(d3_items, corProp=0.03, iterations=100, align_signs=TRUE)
rr3_auto <- ReReReRe(d3_items, corProp=0.03, iterations=100, align_signs=TRUE, auto_z=TRUE)
mah3 <- run_mah(d3_items)

best3 <- find_best_z(rr3$z_score, d3_label)
ev3_15 <- evaluate(as.integer(rr3$z_score <= 1.5), d3_label)
ev3_auto <- evaluate(rr3_auto$flagged, d3_label)
ev3_mah <- evaluate(mah3$flagged, d3_label)

auc3_rr <- compute_auc(rr3$z_score, d3_label, ">")
auc3_mah <- compute_auc(mah3$d2, d3_label, "<")

cat(sprintf("  RR z=1.5:    MCC=%.3f  Sens=%.3f  Spec=%.3f\n", ev3_15$MCC, ev3_15$Sens, ev3_15$Spec))
cat(sprintf("  RR auto-z:   MCC=%.3f  Sens=%.3f  Spec=%.3f  (z=%.2f)\n",
    ev3_auto$MCC, ev3_auto$Sens, ev3_auto$Spec, rr3_auto$z_threshold_used[1]))
cat(sprintf("  RR oracle:   MCC=%.3f  (z=%.1f)\n", best3$mcc, best3$z))
cat(sprintf("  Mahalanobis: MCC=%.3f  Sens=%.3f  Spec=%.3f\n", ev3_mah$MCC, ev3_mah$Sens, ev3_mah$Spec))
cat(sprintf("  AUC: RR=%.3f, Mah=%.3f\n", auc3_rr, auc3_mah))

results[[3]] <- data.frame(
  Dataset="Niessen2016", Items=ncol(d3_items), nF_approx=5, N=nrow(d3_items),
  Pct_careless=round(mean(d3_label)*100,1), GT_type="Experimental (speed)",
  RR_AUC=round(auc3_rr,3), Mah_AUC=round(auc3_mah,3),
  RR_MCC_z15=round(ev3_15$MCC,3), RR_MCC_auto=round(ev3_auto$MCC,3),
  RR_MCC_oracle=round(best3$mcc,3), RR_oracle_z=best3$z,
  Mah_MCC=round(ev3_mah$MCC,3), Auto_z=round(rr3_auto$z_threshold_used[1],2))

# ============================================================
# SUMMARY
# ============================================================
cat("\n\n====================================================\n")
cat("   SUMMARY TABLE\n")
cat("====================================================\n\n")

summary_df <- do.call(rbind, results)
write.csv(summary_df, "archive/external_validation_v2.csv", row.names=FALSE)

cat(sprintf("%-15s %4s %3s %4s %5s | %5s %5s | %5s %5s %5s %4s | %5s\n",
    "Dataset", "Items", "nF", "N", "Pct%",
    "RR_AUC", "Mah_AUC",
    "z=1.5", "auto", "oracl", "oz",
    "Mah"))
cat(paste(rep("-", 90), collapse=""), "\n")

for (i in 1:nrow(summary_df)) {
  r <- summary_df[i,]
  cat(sprintf("%-15s %4d %3d %4d %5.1f | %5.3f %5.3f | %5.3f %5.3f %5.3f %4.1f | %5.3f\n",
      r$Dataset, r$Items, r$nF_approx, r$N, r$Pct_careless,
      r$RR_AUC, r$Mah_AUC,
      r$RR_MCC_z15, r$RR_MCC_auto, r$RR_MCC_oracle, r$RR_oracle_z,
      r$Mah_MCC))
}

cat("\nSaved results to archive/external_validation_v2.csv\n")
cat("\n=== DONE ===\n")
