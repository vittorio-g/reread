#### Goldammer 2024: External Validation on Experimental Data ####
# 5 studies with experimental manipulation (honest vs instructed careless)
# BFI-2 facets, 60 items, 6-point scale, ~15 facets

setwd("C:/Users/vitto/Desktop/ReReReRe")
library(pROC)
source("ReReReRe.R")

compute_mcc <- function(tp, tn, fp, fn) {
  tp <- as.double(tp); tn <- as.double(tn); fp <- as.double(fp); fn <- as.double(fn)
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
  list(MCC=compute_mcc(tp, tn, fp, fn), Sens=tp/(tp+fn), Spec=tn/(tn+fp),
       TP=tp, TN=tn, FP=fp, FN=fn)
}

compute_auc <- function(score, true_label, direction) {
  tryCatch({
    as.numeric(auc(roc(true_label, score, direction=direction, quiet=TRUE)))
  }, error = function(e) NA)
}

run_mah <- function(data_items, alpha=0.001) {
  mat <- as.matrix(data_items)
  p <- ncol(mat)
  center <- colMeans(mat, na.rm=TRUE)
  cov_mat <- cov(mat, use="pairwise.complete.obs")
  if (any(is.na(cov_mat)) || det(cov_mat) < 1e-10) cov_mat <- cov_mat + diag(0.01, p)
  d2 <- mahalanobis(mat, center, cov_mat)
  threshold <- qchisq(1-alpha, df=p)
  list(d2=d2, flagged=d2>threshold)
}

find_best_z <- function(z_scores, labels, z_range=seq(0.1, 3.0, 0.1)) {
  best_mcc <- -Inf; best_z <- 1.5
  for (z in z_range) {
    ev <- evaluate(as.integer(z_scores <= z), labels)
    if (ev$MCC > best_mcc) { best_mcc <- ev$MCC; best_z <- z }
  }
  list(z=best_z, mcc=best_mcc)
}

cat("\n====================================================\n")
cat("   GOLDAMMER 2024: EXPERIMENTAL CARELESS DETECTION\n")
cat("   BFI-2, 60 items, 6-point scale, ~15 facets\n")
cat("   corProp=0.03 (new default)\n")
cat("====================================================\n\n")

results <- list()

# ============================================================
# Helper to process one study
# ============================================================
process_study <- function(file, study_name, gt_col="careless_all") {
  cat(sprintf("\n=== %s ===\n", study_name))

  d <- read.csv(file)
  cat(sprintf("  N=%d, cols=%d\n", nrow(d), ncol(d)))

  # Get ground truth
  labels <- as.integer(d[[gt_col]])

  # Get item columns (original items, not R-coded)
  item_cols <- grep("^[ecano]_[a-z]+[0-9]+$", names(d), value=TRUE)
  items <- d[, item_cols]

  # Remove incomplete
  complete <- complete.cases(items) & !is.na(labels)
  items <- items[complete, ]
  labels <- labels[complete]

  cat(sprintf("  Items=%d, N=%d (complete), Careless=%d (%.1f%%)\n",
      ncol(items), nrow(items), sum(labels), mean(labels)*100))

  # Run ReReReRe
  rr <- ReReReRe(items, corProp=0.03, iterations=100, align_signs=TRUE)
  rr_auto <- ReReReRe(items, corProp=0.03, iterations=100, align_signs=TRUE, auto_z=TRUE)
  mah <- run_mah(items)

  # Evaluate
  best <- find_best_z(rr$z_score, labels)
  ev_15 <- evaluate(as.integer(rr$z_score <= 1.5), labels)
  ev_auto <- evaluate(rr_auto$flagged, labels)
  ev_mah <- evaluate(mah$flagged, labels)

  auc_rr <- compute_auc(rr$z_score, labels, ">")
  auc_mah <- compute_auc(mah$d2, labels, "<")

  cat(sprintf("  RR z=1.5:    MCC=%.3f  Sens=%.3f  Spec=%.3f\n", ev_15$MCC, ev_15$Sens, ev_15$Spec))
  cat(sprintf("  RR auto-z:   MCC=%.3f  Sens=%.3f  Spec=%.3f  (z=%.2f)\n",
      ev_auto$MCC, ev_auto$Sens, ev_auto$Spec, rr_auto$z_threshold_used[1]))
  cat(sprintf("  RR oracle:   MCC=%.3f  (z=%.1f)\n", best$mcc, best$z))
  cat(sprintf("  Mahalanobis: MCC=%.3f  Sens=%.3f  Spec=%.3f\n", ev_mah$MCC, ev_mah$Sens, ev_mah$Spec))
  cat(sprintf("  AUC: RR=%.3f, Mah=%.3f\n", auc_rr, auc_mah))

  data.frame(
    Dataset=study_name, Items=ncol(items), N=nrow(items),
    Pct_careless=round(mean(labels)*100,1), GT_type="Experimental",
    RR_AUC=round(auc_rr,3), Mah_AUC=round(auc_mah,3),
    RR_MCC_z15=round(ev_15$MCC,3), RR_MCC_auto=round(ev_auto$MCC,3),
    RR_MCC_oracle=round(best$mcc,3), RR_oracle_z=best$z,
    Mah_MCC=round(ev_mah$MCC,3), Auto_z=round(rr_auto$z_threshold_used[1],2))
}

# ============================================================
# Run all studies
# ============================================================
basedir <- "external_datasets/13_Goldammer_2024"

# Study 1: BFI-2 with bidirectional keying (60 items)
results[[1]] <- process_study(file.path(basedir, "Study_1.csv"), "Gold_S1_BFI2_bidir")

# Study 2: IPIP-5F with unidirectional keying (60 items)
results[[2]] <- process_study(file.path(basedir, "Study_2.csv"), "Gold_S2_IPIP_unidir")

# Study 3: BFI-2 longitudinal (use t1 only for simplicity)
cat("\n=== Goldammer Study 3 (t1 only) ===\n")
d3 <- read.csv(file.path(basedir, "Study_3.csv"))
# condition: 0=honest, 1=33% careless, 2=100% careless → label any careless
labels3 <- as.integer(d3$condition > 0)
# Get t1 items only
item_cols_t1 <- grep("^[ecano]_[a-z]+[0-9]+_t1$", names(d3), value=TRUE)
if (length(item_cols_t1) == 0) {
  # Try without _t1 suffix
  item_cols_t1 <- grep("^[ecano]_[a-z]+[0-9]+$", names(d3), value=TRUE)
}
cat(sprintf("  t1 items found: %d\n", length(item_cols_t1)))

if (length(item_cols_t1) >= 20) {
  items3 <- d3[, item_cols_t1]
  complete3 <- complete.cases(items3) & !is.na(labels3)
  items3 <- items3[complete3, ]; labels3 <- labels3[complete3]
  cat(sprintf("  N=%d, items=%d, careless=%d\n", nrow(items3), ncol(items3), sum(labels3)))

  rr3 <- ReReReRe(items3, corProp=0.03, iterations=100, align_signs=TRUE)
  rr3_auto <- ReReReRe(items3, corProp=0.03, iterations=100, align_signs=TRUE, auto_z=TRUE)
  mah3 <- run_mah(items3)

  best3 <- find_best_z(rr3$z_score, labels3)
  ev3_15 <- evaluate(as.integer(rr3$z_score <= 1.5), labels3)
  ev3_auto <- evaluate(rr3_auto$flagged, labels3)
  ev3_mah <- evaluate(mah3$flagged, labels3)
  auc3_rr <- compute_auc(rr3$z_score, labels3, ">")
  auc3_mah <- compute_auc(mah3$d2, labels3, "<")

  cat(sprintf("  RR z=1.5: MCC=%.3f | auto-z: MCC=%.3f (z=%.2f) | oracle: MCC=%.3f (z=%.1f) | Mah: MCC=%.3f\n",
      ev3_15$MCC, ev3_auto$MCC, rr3_auto$z_threshold_used[1], best3$mcc, best3$z, ev3_mah$MCC))
  cat(sprintf("  AUC: RR=%.3f, Mah=%.3f\n", auc3_rr, auc3_mah))

  results[[3]] <- data.frame(
    Dataset="Gold_S3_longitudinal", Items=ncol(items3), N=nrow(items3),
    Pct_careless=round(mean(labels3)*100,1), GT_type="Experimental",
    RR_AUC=round(auc3_rr,3), Mah_AUC=round(auc3_mah,3),
    RR_MCC_z15=round(ev3_15$MCC,3), RR_MCC_auto=round(ev3_auto$MCC,3),
    RR_MCC_oracle=round(best3$mcc,3), RR_oracle_z=best3$z,
    Mah_MCC=round(ev3_mah$MCC,3), Auto_z=round(rr3_auto$z_threshold_used[1],2))
}

# ============================================================
# Also test: careless_33 (partial) vs careless_100 (full) separately
# ============================================================
cat("\n\n=== BREAKDOWN: 33% vs 100% careless (Study 1) ===\n")
d1 <- read.csv(file.path(basedir, "Study_1.csv"))
item_cols <- grep("^[ecano]_[a-z]+[0-9]+$", names(d1), value=TRUE)
items1 <- d1[, item_cols]

# 100% careless vs honest
labels_100 <- ifelse(d1$condition == 2, 1, ifelse(d1$condition == 0, 0, NA))
ok <- complete.cases(items1) & !is.na(labels_100)
rr_100 <- ReReReRe(items1[ok,], corProp=0.03, iterations=100, align_signs=TRUE)
best_100 <- find_best_z(rr_100$z_score, labels_100[ok])
auc_100 <- compute_auc(rr_100$z_score, labels_100[ok], ">")
cat(sprintf("  100%% careless vs honest: AUC=%.3f, oracle MCC=%.3f (z=%.1f)\n",
    auc_100, best_100$mcc, best_100$z))

# 33% careless vs honest
labels_33 <- ifelse(d1$condition == 1, 1, ifelse(d1$condition == 0, 0, NA))
ok33 <- complete.cases(items1) & !is.na(labels_33)
rr_33 <- ReReReRe(items1[ok33,], corProp=0.03, iterations=100, align_signs=TRUE)
best_33 <- find_best_z(rr_33$z_score, labels_33[ok33])
auc_33 <- compute_auc(rr_33$z_score, labels_33[ok33], ">")
cat(sprintf("  33%% careless vs honest:  AUC=%.3f, oracle MCC=%.3f (z=%.1f)\n",
    auc_33, best_33$mcc, best_33$z))

# ============================================================
# SUMMARY
# ============================================================
cat("\n\n====================================================\n")
cat("   SUMMARY: GOLDAMMER 2024 VALIDATION\n")
cat("====================================================\n\n")

summary_df <- do.call(rbind, results)
write.csv(summary_df, "archive/goldammer_validation.csv", row.names=FALSE)

for (i in 1:nrow(summary_df)) {
  r <- summary_df[i,]
  cat(sprintf("%-25s Items=%d N=%3d Pct=%.0f%% | AUC: RR=%.3f Mah=%.3f | MCC: z1.5=%.3f auto=%.3f oracle=%.3f(z=%.1f) Mah=%.3f\n",
      r$Dataset, r$Items, r$N, r$Pct_careless,
      r$RR_AUC, r$Mah_AUC,
      r$RR_MCC_z15, r$RR_MCC_auto, r$RR_MCC_oracle, r$RR_oracle_z, r$Mah_MCC))
}

cat("\n=== DONE ===\n")
