#### Test auto_z on real datasets ####
# Compares ReReReRe with auto_z vs fixed z vs Mahalanobis
# on the 3 validated external datasets.
#
# 2026-03-26

library(dplyr)
library(pROC)
library(haven)

setwd("C:/Users/vitto/Desktop/ReReReRe")
source("ReReReRe.R")

# --- Helpers ---
compute_mcc <- function(tp, tn, fp, fn) {
  tp <- as.double(tp); tn <- as.double(tn)
  fp <- as.double(fp); fn <- as.double(fn)
  num <- (tp * tn) - (fp * fn)
  den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  if (is.na(den) || den == 0) return(0)
  num / den
}

evaluate <- function(predicted_flag, true_label) {
  tp <- sum(predicted_flag == 1 & true_label == 1, na.rm = TRUE)
  tn <- sum(predicted_flag == 0 & true_label == 0, na.rm = TRUE)
  fp <- sum(predicted_flag == 1 & true_label == 0, na.rm = TRUE)
  fn <- sum(predicted_flag == 0 & true_label == 1, na.rm = TRUE)
  mcc <- compute_mcc(tp, tn, fp, fn)
  sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA
  spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA
  list(MCC = mcc, Sens = sens, Spec = spec, TP = tp, FP = fp, FN = fn, TN = tn)
}

compute_auc <- function(score, true_label, direction) {
  tryCatch({
    roc_obj <- roc(true_label, score, direction = direction, quiet = TRUE)
    as.numeric(auc(roc_obj))
  }, error = function(e) NA)
}

run_mahalanobis <- function(data_items, alpha = 0.001) {
  mat <- as.matrix(data_items)
  p <- ncol(mat)
  complete_idx <- complete.cases(mat)
  center <- colMeans(mat[complete_idx, ], na.rm = TRUE)
  cov_mat <- cov(mat[complete_idx, ], use = "pairwise.complete.obs")
  if (any(is.na(cov_mat)) || det(cov_mat) < 1e-10) cov_mat <- cov_mat + diag(0.01, p)
  d2 <- mahalanobis(mat, center, cov_mat)
  threshold <- qchisq(1 - alpha, df = p)
  list(d2 = d2, flagged = d2 > threshold)
}

cat("\n============================================================\n")
cat("   AUTO-Z VALIDATION ON REAL DATASETS\n")
cat("   Comparing: auto_z | fixed z=1.5 | fixed z=2.0 | Mahalanobis\n")
cat("============================================================\n\n")

results <- list()

# ============================================================
# DATASET 1: Schneider QoL (4 factors, 31 items)
# ============================================================
cat("=== Schneider QoL (31 items, ~4 factors) ===\n")
d <- read.csv("Dataset/dataset_6/carersp.csv")
label <- d$c01
items <- d[, grep("^(dep|pain|cog|fat)\\d", names(d))]
cc <- complete.cases(items) & !is.na(label)
items <- items[cc, ]; label <- label[cc]
cat(sprintf("  N=%d, %d careless (%.1f%%)\n", length(label), sum(label), mean(label)*100))

# Run auto_z
cat("  Running auto_z...\n")
rr_auto <- ReReReRe(items, auto_z = TRUE, progress = FALSE)
cat(sprintf("  -> nF detected: %d, z_threshold: %.3f\n",
            rr_auto$nFactors_detected[1], rr_auto$z_threshold_used[1]))

# Run fixed z=1.5
rr_15 <- ReReReRe(items, z_threshold = 1.5, progress = FALSE)
# Run fixed z=2.0
rr_20 <- ReReReRe(items, z_threshold = 2.0, progress = FALSE)
# Mahalanobis
mah <- run_mahalanobis(items)

# Evaluate all
e_auto <- evaluate(rr_auto$flagged, label)
e_15   <- evaluate(rr_15$flagged, label)
e_20   <- evaluate(rr_20$flagged, label)
e_mah  <- evaluate(mah$flagged, label)

auc_rr  <- compute_auc(rr_auto$z_score, label, ">")
auc_mah <- compute_auc(mah$d2, label, "<")

cat(sprintf("  auto_z (z=%.2f): MCC=%.3f  Sens=%.3f  Spec=%.3f\n",
            rr_auto$z_threshold_used[1], e_auto$MCC, e_auto$Sens, e_auto$Spec))
cat(sprintf("  fixed z=1.5:     MCC=%.3f  Sens=%.3f  Spec=%.3f\n", e_15$MCC, e_15$Sens, e_15$Spec))
cat(sprintf("  fixed z=2.0:     MCC=%.3f  Sens=%.3f  Spec=%.3f\n", e_20$MCC, e_20$Sens, e_20$Spec))
cat(sprintf("  Mahalanobis:     MCC=%.3f  Sens=%.3f  Spec=%.3f\n", e_mah$MCC, e_mah$Sens, e_mah$Spec))
cat(sprintf("  AUC: RR=%.3f, Mah=%.3f\n\n", auc_rr, auc_mah))

results[["Schneider"]] <- list(nF = rr_auto$nFactors_detected[1],
                                z_auto = rr_auto$z_threshold_used[1],
                                mcc_auto = e_auto$MCC, mcc_15 = e_15$MCC,
                                mcc_20 = e_20$MCC, mcc_mah = e_mah$MCC,
                                auc_rr = auc_rr, auc_mah = auc_mah)


# ============================================================
# DATASET 2: Schroeders HEXACO (60 items, ~6 factors)
# ============================================================
cat("=== Schroeders HEXACO (60 items, ~6 factors) ===\n")
d <- read.csv("Dataset/dataset_1/data_mod_resp.csv", sep = ";")
label <- d$Careless
items <- d %>% select(starts_with("HE"))
cc <- complete.cases(items) & !is.na(label)
items <- items[cc, ]; label <- label[cc]
cat(sprintf("  N=%d, %d careless (%.1f%%)\n", length(label), sum(label), mean(label)*100))

cat("  Running auto_z...\n")
rr_auto <- ReReReRe(items, auto_z = TRUE, progress = FALSE)
cat(sprintf("  -> nF detected: %d, z_threshold: %.3f\n",
            rr_auto$nFactors_detected[1], rr_auto$z_threshold_used[1]))

rr_15 <- ReReReRe(items, z_threshold = 1.5, progress = FALSE)
rr_20 <- ReReReRe(items, z_threshold = 2.0, progress = FALSE)
mah <- run_mahalanobis(items)

e_auto <- evaluate(rr_auto$flagged, label)
e_15   <- evaluate(rr_15$flagged, label)
e_20   <- evaluate(rr_20$flagged, label)
e_mah  <- evaluate(mah$flagged, label)

auc_rr  <- compute_auc(rr_auto$z_score, label, ">")
auc_mah <- compute_auc(mah$d2, label, "<")

cat(sprintf("  auto_z (z=%.2f): MCC=%.3f  Sens=%.3f  Spec=%.3f\n",
            rr_auto$z_threshold_used[1], e_auto$MCC, e_auto$Sens, e_auto$Spec))
cat(sprintf("  fixed z=1.5:     MCC=%.3f  Sens=%.3f  Spec=%.3f\n", e_15$MCC, e_15$Sens, e_15$Spec))
cat(sprintf("  fixed z=2.0:     MCC=%.3f  Sens=%.3f  Spec=%.3f\n", e_20$MCC, e_20$Sens, e_20$Spec))
cat(sprintf("  Mahalanobis:     MCC=%.3f  Sens=%.3f  Spec=%.3f\n", e_mah$MCC, e_mah$Sens, e_mah$Spec))
cat(sprintf("  AUC: RR=%.3f, Mah=%.3f\n\n", auc_rr, auc_mah))

results[["Schroeders"]] <- list(nF = rr_auto$nFactors_detected[1],
                                 z_auto = rr_auto$z_threshold_used[1],
                                 mcc_auto = e_auto$MCC, mcc_15 = e_15$MCC,
                                 mcc_20 = e_20$MCC, mcc_mah = e_mah$MCC,
                                 auc_rr = auc_rr, auc_mah = auc_mah)


# ============================================================
# DATASET 3: Niessen IPIP-100 (100 items, 5 factors)
# ============================================================
cat("=== Niessen IPIP-100 (100 items, 5 factors) ===\n")
d <- read_sav("Dataset/dataset_8/Raw_data.sav")
cat(sprintf("  Raw: %d rows\n", nrow(d)))

# Filter Use_me == 1 if column exists
if ("Use_me" %in% names(d)) {
  d <- d[d$Use_me == 1, ]
  cat(sprintf("  After Use_me filter: %d rows\n", nrow(d)))
}

# Items: E1-E100
item_cols <- paste0("E", 1:100)
item_cols <- item_cols[item_cols %in% names(d)]
items <- as.data.frame(d[, item_cols])
items <- data.frame(lapply(items, as.numeric))

# Ground truth: Conditon (note the typo in original data)
cond_col <- if ("Conditon" %in% names(d)) "Conditon" else if ("Condition" %in% names(d)) "Condition" else NULL
if (!is.null(cond_col)) {
  label <- as.numeric(d[[cond_col]])
} else {
  stop("Cannot find condition column")
}

cc <- complete.cases(items) & !is.na(label)
items <- items[cc, ]; label <- label[cc]
cat(sprintf("  N=%d, %d careless (%.1f%%)\n", length(label), sum(label), mean(label)*100))

cat("  Running auto_z...\n")
rr_auto <- ReReReRe(items, auto_z = TRUE, progress = FALSE)
cat(sprintf("  -> nF detected: %d, z_threshold: %.3f\n",
            rr_auto$nFactors_detected[1], rr_auto$z_threshold_used[1]))

rr_15 <- ReReReRe(items, z_threshold = 1.5, progress = FALSE)
rr_20 <- ReReReRe(items, z_threshold = 2.0, progress = FALSE)
mah <- run_mahalanobis(items)

e_auto <- evaluate(rr_auto$flagged, label)
e_15   <- evaluate(rr_15$flagged, label)
e_20   <- evaluate(rr_20$flagged, label)
e_mah  <- evaluate(mah$flagged, label)

auc_rr  <- compute_auc(rr_auto$z_score, label, ">")
auc_mah <- compute_auc(mah$d2, label, "<")

cat(sprintf("  auto_z (z=%.2f): MCC=%.3f  Sens=%.3f  Spec=%.3f\n",
            rr_auto$z_threshold_used[1], e_auto$MCC, e_auto$Sens, e_auto$Spec))
cat(sprintf("  fixed z=1.5:     MCC=%.3f  Sens=%.3f  Spec=%.3f\n", e_15$MCC, e_15$Sens, e_15$Spec))
cat(sprintf("  fixed z=2.0:     MCC=%.3f  Sens=%.3f  Spec=%.3f\n", e_20$MCC, e_20$Sens, e_20$Spec))
cat(sprintf("  Mahalanobis:     MCC=%.3f  Sens=%.3f  Spec=%.3f\n", e_mah$MCC, e_mah$Sens, e_mah$Spec))
cat(sprintf("  AUC: RR=%.3f, Mah=%.3f\n\n", auc_rr, auc_mah))

results[["Niessen"]] <- list(nF = rr_auto$nFactors_detected[1],
                              z_auto = rr_auto$z_threshold_used[1],
                              mcc_auto = e_auto$MCC, mcc_15 = e_15$MCC,
                              mcc_20 = e_20$MCC, mcc_mah = e_mah$MCC,
                              auc_rr = auc_rr, auc_mah = auc_mah)


# ============================================================
# SUMMARY
# ============================================================
cat("\n============================================================\n")
cat("   SUMMARY TABLE\n")
cat("============================================================\n\n")
cat(sprintf("%-12s %4s %6s %10s %10s %10s %10s %8s %8s\n",
            "Dataset", "nF", "z_auto", "MCC_auto", "MCC_z1.5", "MCC_z2.0", "MCC_Mah", "AUC_RR", "AUC_Mah"))
cat(paste(rep("-", 90), collapse=""), "\n")
for (nm in names(results)) {
  r <- results[[nm]]
  cat(sprintf("%-12s %4d %6.2f %10.3f %10.3f %10.3f %10.3f %8.3f %8.3f\n",
              nm, r$nF, r$z_auto, r$mcc_auto, r$mcc_15, r$mcc_20, r$mcc_mah,
              r$auc_rr, r$auc_mah))
}

cat("\n=== DONE ===\n")
