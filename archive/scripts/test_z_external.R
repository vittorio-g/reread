#### Test z thresholds on external datasets ####
source("ReReReRe.R")
library(pROC)
library(haven)

# Helper: compute MCC from z_scores, ground_truth, threshold
compute_metrics <- function(z_scores, gt, z_thresh) {
  flagged <- as.integer(z_scores <= z_thresh)
  TP <- sum(flagged == 1 & gt == 1)
  TN <- sum(flagged == 0 & gt == 0)
  FP <- sum(flagged == 1 & gt == 0)
  FN <- sum(flagged == 0 & gt == 1)
  denom <- as.double(TP+FP) * (TP+FN) * (TN+FP) * (TN+FN)
  mcc <- if (denom == 0) 0 else (TP*TN - FP*FN) / sqrt(denom)
  sens <- if (TP+FN > 0) TP/(TP+FN) else 0
  spec <- if (TN+FP > 0) TN/(TN+FP) else 0
  list(mcc=mcc, sens=sens, spec=spec, TP=TP, TN=TN, FP=FP, FN=FN)
}

z_thresholds <- seq(0, 3.5, by=0.5)

# ==============================================================
# 1. SCHNEIDER QoL
# ==============================================================
cat("=== SCHNEIDER QoL ===\n")
d3 <- read.csv("external_datasets/06_Schneider_QoL/carersp.csv")
item_cols3 <- grep("^(dep|pain|cog|fat)\\d", names(d3), value=TRUE)
d3_complete <- d3[complete.cases(d3[,item_cols3]),]
d3_items <- d3_complete[, item_cols3]
gt3 <- d3_complete$c01  # 1 = careless

cat(sprintf("N=%d, items=%d, careless=%d (%.1f%%)\n",
    nrow(d3_items), ncol(d3_items), sum(gt3==1), 100*mean(gt3==1)))

# Item ranges
mat3 <- as.matrix(d3_items)
cat("Item ranges:\n")
for (col in item_cols3) {
  cat(sprintf("  %s: %d-%d (unique: %d)\n", col, min(mat3[,col]), max(mat3[,col]),
      length(unique(mat3[,col]))))
}

rr3 <- ReReReRe(d3_items, corProp=0.05, iterations=100, min_pairs=15, align_signs=TRUE)
auc3 <- as.numeric(roc(gt3, rr3$z_score, direction=">", quiet=TRUE)$auc)
cat(sprintf("AUC = %.3f\n", auc3))

cat(sprintf("%-6s %-8s %-8s %-8s %-6s %-6s %-6s %-6s\n",
    "z", "MCC", "Sens", "Spec", "TP", "TN", "FP", "FN"))
for (z in z_thresholds) {
  m <- compute_metrics(rr3$z_score, gt3, z)
  cat(sprintf("%-6.1f %-8.3f %-8.3f %-8.3f %-6d %-6d %-6d %-6d\n",
      z, m$mcc, m$sens, m$spec, m$TP, m$TN, m$FP, m$FN))
}

# Diagnostics
cat(sprintf("\nz_score good:    mean=%.2f  sd=%.2f\n",
    mean(rr3$z_score[gt3==0]), sd(rr3$z_score[gt3==0])))
cat(sprintf("z_score careless: mean=%.2f  sd=%.2f\n",
    mean(rr3$z_score[gt3==1]), sd(rr3$z_score[gt3==1])))
cat(sprintf("indCors good:    mean=%.4f\n", mean(rr3$indCors[gt3==0])))
cat(sprintf("indCors careless: mean=%.4f\n", mean(rr3$indCors[gt3==1])))
cat(sprintf("rand_mean good:    mean=%.4f\n", mean(rr3$rand_mean[gt3==0])))
cat(sprintf("rand_mean careless: mean=%.4f\n", mean(rr3$rand_mean[gt3==1])))

# Negative coupled pairs?
rawCorMat <- cor(mat3, use="pairwise.complete.obs")
absCorMat <- abs(rawCorMat)
absCorMat[upper.tri(absCorMat, diag=TRUE)] <- NA
rawCorMat[upper.tri(rawCorMat, diag=TRUE)] <- NA
thresh3 <- quantile(absCorMat, 1-0.05, na.rm=TRUE)
couples3 <- which(absCorMat >= thresh3, arr.ind=TRUE)
signs3 <- sign(rawCorMat[couples3])
cat(sprintf("Coupled pairs: %d total, %d negative (%.1f%%)\n",
    length(signs3), sum(signs3<0), 100*mean(signs3<0)))

# ==============================================================
# 2. SCHROEDERS HEXACO
# ==============================================================
cat("\n=== SCHROEDERS HEXACO ===\n")
d1 <- read.csv2("external_datasets/01_Schroeders_2022/data_mod_resp.csv")
item_cols1 <- grep("^HE01_", names(d1), value=TRUE)
d1_complete <- d1[complete.cases(d1[,item_cols1]),]
d1_items <- d1_complete[, item_cols1]
gt1 <- d1_complete$Careless  # 0/1

cat(sprintf("N=%d, items=%d, careless=%d (%.1f%%)\n",
    nrow(d1_items), ncol(d1_items), sum(gt1==1), 100*mean(gt1==1)))

# Item ranges
mat1 <- as.matrix(d1_items)
ranges1 <- apply(mat1, 2, range)
cat(sprintf("Item range: %d-%d across all items\n", min(ranges1[1,]), max(ranges1[2,])))
n_neg_pairs1 <- NA

rr1 <- ReReReRe(d1_items, corProp=0.05, iterations=100, min_pairs=15, align_signs=TRUE)
auc1 <- as.numeric(roc(gt1, rr1$z_score, direction=">", quiet=TRUE)$auc)
cat(sprintf("AUC = %.3f\n", auc1))

cat(sprintf("%-6s %-8s %-8s %-8s\n", "z", "MCC", "Sens", "Spec"))
for (z in z_thresholds) {
  m <- compute_metrics(rr1$z_score, gt1, z)
  cat(sprintf("%-6.1f %-8.3f %-8.3f %-8.3f\n", z, m$mcc, m$sens, m$spec))
}

cat(sprintf("\nz_score good:    mean=%.2f  sd=%.2f\n",
    mean(rr1$z_score[gt1==0]), sd(rr1$z_score[gt1==0])))
cat(sprintf("z_score careless: mean=%.2f  sd=%.2f\n",
    mean(rr1$z_score[gt1==1]), sd(rr1$z_score[gt1==1])))

rawCorMat1 <- cor(mat1, use="pairwise.complete.obs")
absCorMat1 <- abs(rawCorMat1)
absCorMat1[upper.tri(absCorMat1, diag=TRUE)] <- NA
rawCorMat1[upper.tri(rawCorMat1, diag=TRUE)] <- NA
thresh1 <- quantile(absCorMat1, 1-0.05, na.rm=TRUE)
couples1 <- which(absCorMat1 >= thresh1, arr.ind=TRUE)
signs1 <- sign(rawCorMat1[couples1])
cat(sprintf("Coupled pairs: %d total, %d negative (%.1f%%)\n",
    length(signs1), sum(signs1<0), 100*mean(signs1<0)))

# ==============================================================
# 3. NIESSEN IPIP
# ==============================================================
cat("\n=== NIESSEN IPIP ===\n")
d2_raw <- read_sav("external_datasets/08_Niessen_2016/Raw_data.sav")
d2 <- d2_raw[d2_raw$Use_me == 1, ]
item_cols2 <- grep("^[ECANO][0-9]+$", names(d2), value=TRUE)
d2_complete <- d2[complete.cases(d2[, item_cols2]), ]
d2_items <- as.data.frame(d2_complete[, item_cols2])
gt2 <- as.integer(d2_complete$Conditon)  # 0/1

cat(sprintf("N=%d, items=%d, careless=%d (%.1f%%)\n",
    nrow(d2_items), ncol(d2_items), sum(gt2==1), 100*mean(gt2==1)))

mat2 <- as.matrix(d2_items)
ranges2 <- apply(mat2, 2, range)
cat(sprintf("Item range: %d-%d across all items\n", min(ranges2[1,]), max(ranges2[2,])))

rr2 <- ReReReRe(d2_items, corProp=0.05, iterations=100, min_pairs=15, align_signs=TRUE)
auc2 <- as.numeric(roc(gt2, rr2$z_score, direction=">", quiet=TRUE)$auc)
cat(sprintf("AUC = %.3f\n", auc2))

cat(sprintf("%-6s %-8s %-8s %-8s\n", "z", "MCC", "Sens", "Spec"))
for (z in z_thresholds) {
  m <- compute_metrics(rr2$z_score, gt2, z)
  cat(sprintf("%-6.1f %-8.3f %-8.3f %-8.3f\n", z, m$mcc, m$sens, m$spec))
}

cat(sprintf("\nz_score good:    mean=%.2f  sd=%.2f\n",
    mean(rr2$z_score[gt2==0]), sd(rr2$z_score[gt2==0])))
cat(sprintf("z_score careless: mean=%.2f  sd=%.2f\n",
    mean(rr2$z_score[gt2==1]), sd(rr2$z_score[gt2==1])))

rawCorMat2 <- cor(mat2, use="pairwise.complete.obs")
absCorMat2 <- abs(rawCorMat2)
absCorMat2[upper.tri(absCorMat2, diag=TRUE)] <- NA
rawCorMat2[upper.tri(rawCorMat2, diag=TRUE)] <- NA
thresh2 <- quantile(absCorMat2, 1-0.05, na.rm=TRUE)
couples2 <- which(absCorMat2 >= thresh2, arr.ind=TRUE)
signs2 <- sign(rawCorMat2[couples2])
cat(sprintf("Coupled pairs: %d total, %d negative (%.1f%%)\n",
    length(signs2), sum(signs2<0), 100*mean(signs2<0)))

# ==============================================================
# SUMMARY
# ==============================================================
cat("\n=== SUMMARY (reverse coding align_signs) ===\n")
cat(sprintf("%-20s %-8s %-10s %-10s %-10s %-10s\n",
    "Dataset", "AUC", "MCC@1.0", "MCC@1.5", "MCC@2.0", "Best_z"))
for (info in list(
  list(name="Schneider", auc=auc3, rr=rr3, gt=gt3),
  list(name="Schroeders", auc=auc1, rr=rr1, gt=gt1),
  list(name="Niessen", auc=auc2, rr=rr2, gt=gt2)
)) {
  m10 <- compute_metrics(info$rr$z_score, info$gt, 1.0)
  m15 <- compute_metrics(info$rr$z_score, info$gt, 1.5)
  m20 <- compute_metrics(info$rr$z_score, info$gt, 2.0)
  best_mcc <- -1; best_z <- NA
  for (z in seq(0, 3.5, by=0.5)) {
    m <- compute_metrics(info$rr$z_score, info$gt, z)
    if (m$mcc > best_mcc) { best_mcc <- m$mcc; best_z <- z }
  }
  cat(sprintf("%-20s %-8.3f %-10.3f %-10.3f %-10.3f z=%.1f (%.3f)\n",
      info$name, info$auc, m10$mcc, m15$mcc, m20$mcc, best_z, best_mcc))
}

cat("\nDone.\n")
