#### External Validation of ReReReRe on Real Datasets ####
# Compares ReReReRe (fixed defaults) vs practical Mahalanobis (chi-sq .001)
# on datasets with known/labeled careless respondents.
#
# 2026-03-25

library(dplyr)
library(pROC)

# Source the algorithm
source("ReReReRe.R")

# ============================================================
# Helper functions
# ============================================================

compute_mcc <- function(tp, tn, fp, fn) {
  tp <- as.double(tp); tn <- as.double(tn)
  fp <- as.double(fp); fn <- as.double(fn)
  num <- (tp * tn) - (fp * fn)
  den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  if (is.na(den) || den == 0) return(0)
  num / den
}

evaluate_method <- function(predicted_flag, true_label) {
  # true_label: 1 = careless, 0 = diligent
  # predicted_flag: TRUE/1 = flagged as careless
  tp <- sum(predicted_flag == 1 & true_label == 1, na.rm = TRUE)
  tn <- sum(predicted_flag == 0 & true_label == 0, na.rm = TRUE)
  fp <- sum(predicted_flag == 1 & true_label == 0, na.rm = TRUE)
  fn <- sum(predicted_flag == 0 & true_label == 1, na.rm = TRUE)

  mcc <- compute_mcc(tp, tn, fp, fn)
  sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA
  spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA
  ppv  <- if ((tp + fp) > 0) tp / (tp + fp) else NA
  f1   <- if (!is.na(sens) && !is.na(ppv) && (sens + ppv) > 0) 2 * sens * ppv / (sens + ppv) else NA

  data.frame(
    N = length(true_label),
    N_careless = sum(true_label == 1, na.rm = TRUE),
    Pct_careless = round(mean(true_label == 1, na.rm = TRUE) * 100, 1),
    TP = tp, TN = tn, FP = fp, FN = fn,
    MCC = round(mcc, 3),
    Sensitivity = round(sens, 3),
    Specificity = round(spec, 3),
    PPV = round(ppv, 3),
    F1 = round(f1, 3)
  )
}

run_mahalanobis_practical <- function(data_items, alpha = 0.001) {
  # Practical Mahalanobis: flag if D^2 > chi-square(p, alpha)
  # This is the standard approach from Meade & Craig (2012)
  mat <- as.matrix(data_items)
  p <- ncol(mat)

  # Handle missing values: use complete cases for covariance
  complete_idx <- complete.cases(mat)
  if (sum(complete_idx) < p + 1) {
    warning("Too few complete cases for Mahalanobis")
    return(rep(NA, nrow(mat)))
  }

  center <- colMeans(mat[complete_idx, ], na.rm = TRUE)
  cov_mat <- cov(mat[complete_idx, ], use = "pairwise.complete.obs")

  # Regularize if singular
  if (any(is.na(cov_mat)) || det(cov_mat) < 1e-10) {
    cov_mat <- cov_mat + diag(0.01, p)
  }

  d2 <- mahalanobis(mat, center, cov_mat)
  threshold <- qchisq(1 - alpha, df = p)

  cat(sprintf("  Mahalanobis: p=%d, chi-sq threshold=%.1f, flagged=%d/%d\n",
              p, threshold, sum(d2 > threshold, na.rm = TRUE), length(d2)))

  list(d2 = d2, flagged = d2 > threshold, threshold = threshold)
}

run_rererere <- function(data_items, corProp = 0.05, z_threshold = 2.0) {
  cat(sprintf("  ReReReRe: %d respondents x %d items, corProp=%.2f\n",
              nrow(data_items), ncol(data_items), corProp))

  rr <- ReReReRe(data_items, corProp = corProp, iterations = 100,
                 min_pairs = 15, align_signs = TRUE)
  rr$flagged_z <- rr$z_score <= z_threshold

  cat(sprintf("  ReReReRe: flagged=%d/%d (z <= %.1f)\n",
              sum(rr$flagged_z, na.rm = TRUE), nrow(rr), z_threshold))
  rr
}

# ============================================================
# Also compute AUC for each method
# ============================================================
compute_auc <- function(score, true_label, direction) {
  # pROC direction convention:
  #   ">" means controls have HIGHER values than cases
  #   "<" means controls have LOWER values than cases
  # For z_score: good respondents (controls) have HIGH z → direction=">"
  # For D^2:    good respondents (controls) have LOW D^2 → direction="<"
  tryCatch({
    roc_obj <- roc(true_label, score, direction = direction, quiet = TRUE)
    as.numeric(auc(roc_obj))
  }, error = function(e) NA)
}


cat("\n====================================================\n")
cat("   EXTERNAL VALIDATION: ReReReRe vs Mahalanobis\n")
cat("   Fixed defaults: corProp=0.05, z_threshold=2.0\n")
cat("   Mahalanobis: chi-square(p, alpha=.001)\n")
cat("====================================================\n\n")

results_all <- list()

# ============================================================
# DATASET 1: Schroeders et al. 2022
# ============================================================
cat("=== Dataset 1: Schroeders et al. 2022 ===\n")
cat("  Experimental induction: instructed careless vs diligent\n")

d1 <- read.csv("external_datasets/01_Schroeders_2022/data_mod_resp.csv", sep = ";")
cat(sprintf("  Raw: %d rows, %d cols\n", nrow(d1), ncol(d1)))

# Ground truth: Careless column (1 = careless, 0 = diligent)
d1_label <- d1$Careless

# Items: all columns starting with HE
d1_items <- d1 %>% select(starts_with("HE"))
cat(sprintf("  Items: %d\n", ncol(d1_items)))

# Remove rows with all NA
complete <- complete.cases(d1_items)
d1_items <- d1_items[complete, ]
d1_label <- d1_label[complete]
cat(sprintf("  After removing incomplete: %d rows\n", nrow(d1_items)))

# Run methods
rr1 <- run_rererere(d1_items)
mah1 <- run_mahalanobis_practical(d1_items)

# Evaluate
cat("\n  --- Results ---\n")
rr1_eval <- evaluate_method(rr1$flagged_z, d1_label)
mah1_eval <- evaluate_method(mah1$flagged, d1_label)

rr1_auc <- compute_auc(rr1$z_score, d1_label, ">")
mah1_auc <- compute_auc(mah1$d2, d1_label, "<")

rr1_eval$AUC <- round(rr1_auc, 3)
mah1_eval$AUC <- round(mah1_auc, 3)

rr1_eval$Method <- "ReReReRe"
mah1_eval$Method <- "Mahalanobis"
rr1_eval$Dataset <- "Schroeders2022"
mah1_eval$Dataset <- "Schroeders2022"

cat("  ReReReRe:    "); cat(sprintf("MCC=%.3f  Sens=%.3f  Spec=%.3f  AUC=%.3f\n", rr1_eval$MCC, rr1_eval$Sensitivity, rr1_eval$Specificity, rr1_auc))
cat("  Mahalanobis: "); cat(sprintf("MCC=%.3f  Sens=%.3f  Spec=%.3f  AUC=%.3f\n", mah1_eval$MCC, mah1_eval$Sensitivity, mah1_eval$Specificity, mah1_auc))

results_all <- c(results_all, list(rr1_eval, mah1_eval))


# ============================================================
# DATASET 2: Brühlmann et al. 2020
# ============================================================
cat("\n=== Dataset 2: Brühlmann et al. 2020 ===\n")
cat("  CrowdFlower crowdsourced data with multiple CR indicators\n")

# Read the semicolon-delimited, CR-terminated, Latin-1 file
d2_raw <- read.csv("external_datasets/05_Bruhlmann_2020/data_anon.csv",
                    sep = ";", fileEncoding = "latin1",
                    stringsAsFactors = FALSE)
cat(sprintf("  Raw: %d rows, %d cols\n", nrow(d2_raw), ncol(d2_raw)))

# Ground truth options:
# - v_IRI: instructed response item (wrong answer = careless)
# - v_Bogus_Item: bogus item
# - quality: CrowdFlower quality rating (low/high?)
# - v_Serious: self-reported seriousness
# Let's use multiple indicators. First check what they look like.
cat(sprintf("  v_IRI values: %s\n", paste(sort(unique(d2_raw$v_IRI)), collapse=", ")))
cat(sprintf("  v_Bogus_Item values: %s\n", paste(sort(unique(d2_raw$v_Bogus_Item)), collapse=", ")))
cat(sprintf("  quality values: %s\n", paste(sort(unique(d2_raw$quality)), collapse=", ")))
cat(sprintf("  v_Serious values: %s\n", paste(sort(unique(d2_raw$v_Serious)), collapse=", ")))

# Items: PANAS (v_PANAS_1-10), AttrakDiff (v_AT_*), NeedFulfillment (v_NF_*),
# TechCommitment (v_TC_*), VisualAesthetics (v_VAWI_*), BFI-44 (v_BFI_*)
item_cols <- grep("^v_(PANAS|AT_|NF_|TC_F|VAWI_|BFI_)", names(d2_raw), value = TRUE)
# Exclude v_IRI and v_Bogus_Item (they're embedded attention checks, not real items)
item_cols <- setdiff(item_cols, c("v_IRI", "v_Bogus_Item"))
d2_items <- d2_raw[, item_cols]
cat(sprintf("  Items: %d\n", ncol(d2_items)))

# Convert to numeric (some might be character)
d2_items <- data.frame(lapply(d2_items, as.numeric))

# Ground truth: Use bogus item failure as primary indicator
# v_Bogus_Item: What's the correct answer? Typically these are items like
# "Please select 'Agree' for this item" — the wrong answer = careless
# Let's use a composite: fail IRI OR fail bogus item OR low self-reported seriousness
# But first, we need to know the correct answers. Let's check the R script.
# For now, use the quality column from CrowdFlower as primary ground truth
d2_label <- ifelse(d2_raw$quality == "low", 1, 0)
cat(sprintf("  Ground truth (quality=low): %d careless / %d total\n",
            sum(d2_label == 1, na.rm = TRUE), length(d2_label)))

# If quality doesn't differentiate well, try IRI
# IRI is typically: correct answer = specific value, wrong = careless
# Let's check: for a 7-point scale, IRI correct is often 1 or 7
# We'll use quality as primary, and also report IRI-based results

# Remove incomplete
complete2 <- complete.cases(d2_items) & !is.na(d2_label)
d2_items_c <- d2_items[complete2, ]
d2_label_c <- d2_label[complete2]
cat(sprintf("  After removing incomplete: %d rows\n", nrow(d2_items_c)))

# Run methods
rr2 <- run_rererere(d2_items_c)
mah2 <- run_mahalanobis_practical(d2_items_c)

# Evaluate
cat("\n  --- Results (quality-based label) ---\n")
rr2_eval <- evaluate_method(rr2$flagged_z, d2_label_c)
mah2_eval <- evaluate_method(mah2$flagged, d2_label_c)

rr2_auc <- compute_auc(rr2$z_score, d2_label_c, ">")
mah2_auc <- compute_auc(mah2$d2, d2_label_c, "<")

rr2_eval$AUC <- round(rr2_auc, 3)
mah2_eval$AUC <- round(mah2_auc, 3)

rr2_eval$Method <- "ReReReRe"
mah2_eval$Method <- "Mahalanobis"
rr2_eval$Dataset <- "Bruhlmann2020"
mah2_eval$Dataset <- "Bruhlmann2020"

cat("  ReReReRe:    "); cat(sprintf("MCC=%.3f  Sens=%.3f  Spec=%.3f  AUC=%.3f\n", rr2_eval$MCC, rr2_eval$Sensitivity, rr2_eval$Specificity, rr2_auc))
cat("  Mahalanobis: "); cat(sprintf("MCC=%.3f  Sens=%.3f  Spec=%.3f  AUC=%.3f\n", mah2_eval$MCC, mah2_eval$Sensitivity, mah2_eval$Specificity, mah2_auc))

results_all <- c(results_all, list(rr2_eval, mah2_eval))


# ============================================================
# DATASET 3: Schneider et al. — Quality of Life
# ============================================================
cat("\n=== Dataset 3: Schneider et al. (QoL) ===\n")
cat("  Internet-based QoL assessments, latent class for CR\n")

d3 <- read.csv("external_datasets/06_Schneider_QoL/carersp.csv")
cat(sprintf("  Raw: %d rows, %d cols\n", nrow(d3), ncol(d3)))

# Ground truth: c01 (1 = careless class, 0 = attentive class)
d3_label <- d3$c01

# Items: dep01-dep08, pain01-pain08(?), cog01-cog08, fat01-fat07
item_cols3 <- grep("^(dep|pain|cog|fat)\\d", names(d3), value = TRUE)
d3_items <- d3[, item_cols3]
cat(sprintf("  Items: %d\n", ncol(d3_items)))

# Remove incomplete
complete3 <- complete.cases(d3_items) & !is.na(d3_label)
d3_items_c <- d3_items[complete3, ]
d3_label_c <- d3_label[complete3]
cat(sprintf("  After removing incomplete: %d rows, %d careless\n",
            nrow(d3_items_c), sum(d3_label_c == 1)))

# Run methods
rr3 <- run_rererere(d3_items_c)
mah3 <- run_mahalanobis_practical(d3_items_c)

# Evaluate
cat("\n  --- Results ---\n")
rr3_eval <- evaluate_method(rr3$flagged_z, d3_label_c)
mah3_eval <- evaluate_method(mah3$flagged, d3_label_c)

rr3_auc <- compute_auc(rr3$z_score, d3_label_c, ">")
mah3_auc <- compute_auc(mah3$d2, d3_label_c, "<")

rr3_eval$AUC <- round(rr3_auc, 3)
mah3_eval$AUC <- round(mah3_auc, 3)

rr3_eval$Method <- "ReReReRe"
mah3_eval$Method <- "Mahalanobis"
rr3_eval$Dataset <- "Schneider_QoL"
mah3_eval$Dataset <- "Schneider_QoL"

cat("  ReReReRe:    "); cat(sprintf("MCC=%.3f  Sens=%.3f  Spec=%.3f  AUC=%.3f\n", rr3_eval$MCC, rr3_eval$Sensitivity, rr3_eval$Specificity, rr3_auc))
cat("  Mahalanobis: "); cat(sprintf("MCC=%.3f  Sens=%.3f  Spec=%.3f  AUC=%.3f\n", mah3_eval$MCC, mah3_eval$Sensitivity, mah3_eval$Specificity, mah3_auc))

results_all <- c(results_all, list(rr3_eval, mah3_eval))


# ============================================================
# DATASET 4: Niessen et al. 2016 (SPSS)
# ============================================================
cat("\n=== Dataset 4: Niessen et al. 2016 ===\n")
cat("  Experimental speed manipulation in web questionnaire\n")

if (!requireNamespace("haven", quietly = TRUE)) {
  cat("  Installing haven for SPSS files...\n")
  install.packages("haven", repos = "https://cloud.r-project.org", quiet = TRUE)
}
library(haven)

d4_raw <- read_sav("external_datasets/08_Niessen_2016/Raw_data.sav")
cat(sprintf("  Raw: %d rows, %d cols\n", nrow(d4_raw), ncol(d4_raw)))
cat(sprintf("  Column names: %s\n", paste(names(d4_raw), collapse=", ")))

# We need to identify: (1) item columns, (2) ground truth column
# Look for a condition/group column
# Print unique values of likely candidate columns
for (col in names(d4_raw)) {
  u <- length(unique(d4_raw[[col]]))
  if (u <= 10 && u >= 2) {
    cat(sprintf("  %s: %d unique values: %s\n", col, u,
                paste(sort(unique(as.character(d4_raw[[col]]))), collapse=", ")))
  }
}

# Try to identify item columns (those with many unique values on a Likert scale)
potential_items <- sapply(d4_raw, function(x) {
  u <- sort(unique(na.omit(as.numeric(x))))
  length(u) >= 3 && length(u) <= 10 && min(u) >= 0 && max(u) <= 7
})
item_names4 <- names(d4_raw)[potential_items]
cat(sprintf("  Potential item columns (3-10 unique vals, range 0-7): %d\n", length(item_names4)))
if (length(item_names4) > 0) {
  cat(sprintf("  First 10: %s\n", paste(head(item_names4, 10), collapse=", ")))
}

# We'll attempt to run if we can identify items and ground truth
# For now, report what we find and attempt a run if possible

if (length(item_names4) >= 10) {
  # Look for a condition column — typically named something like "condition", "group", "speed"
  cond_candidates <- grep("cond|group|speed|instruct|care|type", names(d4_raw), value = TRUE, ignore.case = TRUE)
  cat(sprintf("  Condition column candidates: %s\n", paste(cond_candidates, collapse=", ")))

  # If we found a clear condition column, use it
  if (length(cond_candidates) > 0) {
    cond_col <- cond_candidates[1]
    d4_label <- as.numeric(d4_raw[[cond_col]])
    cat(sprintf("  Using '%s' as ground truth. Values: %s\n", cond_col,
                paste(sort(unique(d4_label)), collapse=", ")))

    d4_items <- as.data.frame(d4_raw[, item_names4])
    d4_items <- data.frame(lapply(d4_items, as.numeric))

    complete4 <- complete.cases(d4_items) & !is.na(d4_label)
    d4_items_c <- d4_items[complete4, ]
    d4_label_c <- d4_label[complete4]

    rr4 <- run_rererere(d4_items_c)
    mah4 <- run_mahalanobis_practical(d4_items_c)

    rr4_eval <- evaluate_method(rr4$flagged_z, d4_label_c)
    mah4_eval <- evaluate_method(mah4$flagged, d4_label_c)

    rr4_auc <- compute_auc(rr4$z_score, d4_label_c, ">")
    mah4_auc <- compute_auc(mah4$d2, d4_label_c, "<")

    rr4_eval$AUC <- round(rr4_auc, 3)
    mah4_eval$AUC <- round(mah4_auc, 3)

    rr4_eval$Method <- "ReReReRe"
    mah4_eval$Method <- "Mahalanobis"
    rr4_eval$Dataset <- "Niessen2016"
    mah4_eval$Dataset <- "Niessen2016"

    cat("\n  --- Results ---\n")
    cat("  ReReReRe:    "); cat(sprintf("MCC=%.3f  Sens=%.3f  Spec=%.3f  AUC=%.3f\n", rr4_eval$MCC, rr4_eval$Sensitivity, rr4_eval$Specificity, rr4_auc))
    cat("  Mahalanobis: "); cat(sprintf("MCC=%.3f  Sens=%.3f  Spec=%.3f  AUC=%.3f\n", mah4_eval$MCC, mah4_eval$Sensitivity, mah4_eval$Specificity, mah4_auc))

    results_all <- c(results_all, list(rr4_eval, mah4_eval))
  } else {
    cat("  Could not identify condition column. Skipping.\n")
    cat("  All column names: ", paste(names(d4_raw), collapse=", "), "\n")
  }
} else {
  cat("  Could not identify sufficient item columns. Skipping.\n")
}


# ============================================================
# SUMMARY TABLE
# ============================================================
cat("\n\n====================================================\n")
cat("   SUMMARY: All Datasets\n")
cat("====================================================\n\n")

results_df <- bind_rows(results_all)
results_df <- results_df %>% select(Dataset, Method, N, N_careless, Pct_careless,
                                     MCC, AUC, Sensitivity, Specificity, PPV, F1)

# Print nicely
print(results_df, row.names = FALSE)

# Save
write.csv(results_df, "external_validation_results.csv", row.names = FALSE)
cat("\nResults saved to external_validation_results.csv\n")

# Quick comparison
cat("\n--- Head-to-head ---\n")
for (ds in unique(results_df$Dataset)) {
  rr_row <- results_df[results_df$Dataset == ds & results_df$Method == "ReReReRe", ]
  mah_row <- results_df[results_df$Dataset == ds & results_df$Method == "Mahalanobis", ]

  mcc_winner <- if (rr_row$MCC > mah_row$MCC) "ReReReRe" else "Mahalanobis"
  auc_winner <- if (!is.na(rr_row$AUC) && !is.na(mah_row$AUC) && rr_row$AUC > mah_row$AUC) "ReReReRe" else "Mahalanobis"

  cat(sprintf("  %s: MCC winner=%s (%.3f vs %.3f), AUC winner=%s (%.3f vs %.3f)\n",
              ds, mcc_winner, rr_row$MCC, mah_row$MCC,
              auc_winner, rr_row$AUC, mah_row$AUC))
}
