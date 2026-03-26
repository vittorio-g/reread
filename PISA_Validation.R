#### PISA 2018 External Validation ####
# ReReReRe vs practical Mahalanobis on PISA 2018 background questionnaire
# Ground truth: screen-time-based C/IER weights (Ulitzsch et al. 2023)
#
# Strategy: work per-country to avoid country-level mean differences
# inflating correlations. Start with Italy (ITA) — good complete-case
# rate and Marcello's home country.

library(haven)
library(mclust)
library(pROC)
source("ReReReRe.R")

# ============================================================
# 1. Build item column names from scale mapping
# ============================================================
mapping <- read.csv("external_datasets/11_PISA_2018/PISA_2018_scale_mapping.csv",
                    stringsAsFactors = FALSE)

# Parse item patterns like "ST097Q01TA-ST097Q05TA" into actual column names
parse_items <- function(pattern, n_items) {
  # Extract prefix (e.g., "ST097Q") and suffix (e.g., "TA")
  parts <- strsplit(pattern, "-")[[1]]
  first <- parts[1]
  # Extract the prefix before the number, and the suffix after
  prefix <- sub("(ST\\d+Q)\\d+(\\w+)", "\\1", first)
  suffix <- sub("(ST\\d+Q)\\d+(\\w+)", "\\2", first)
  sprintf("%s%02d%s", prefix, 1:n_items, suffix)
}

# Build item list per scale
scale_items <- list()
all_item_cols <- c()
for (i in 1:nrow(mapping)) {
  items <- parse_items(mapping$item_pattern[i], mapping$n_items[i])
  scale_items[[mapping$scale_code[i]]] <- items
  all_item_cols <- c(all_item_cols, items)
}

cat(sprintf("Total item columns: %d across %d scales\n",
            length(all_item_cols), length(scale_items)))

# Exclude binary (2-option) scales for RR analysis
binary_scales <- mapping$scale_code[mapping$n_options == 2]
cat("Binary scales (excluded):", paste(binary_scales, collapse=", "), "\n")

likert_scales <- mapping$scale_code[mapping$n_options > 2]
likert_items <- unlist(scale_items[likert_scales])
cat(sprintf("Likert items: %d across %d scales\n",
            length(likert_items), length(likert_scales)))

# Timing column names: STnnn_TT
timing_cols <- paste0(mapping$scale_code, "_TT")

# ============================================================
# 2. Load data (selective columns to manage memory)
# ============================================================
cat("\nLoading PISA student questionnaire (selective columns)...\n")
# Read only the columns we need: ID, country, + item columns
needed_stu <- c("CNTSTUID", "CNT", all_item_cols)

# haven can read selected columns
stu <- read_sav("external_datasets/11_PISA_2018/STU/CY07_MSU_STU_QQQ.sav",
                col_select = any_of(needed_stu))
cat(sprintf("STU loaded: %d rows x %d cols\n", nrow(stu), ncol(stu)))

# Check which item columns we actually got
found_items <- intersect(all_item_cols, names(stu))
missing_items <- setdiff(all_item_cols, names(stu))
cat(sprintf("Found %d/%d item columns\n", length(found_items), length(all_item_cols)))
if (length(missing_items) > 0) {
  cat("Missing:", paste(head(missing_items, 20), collapse=", "), "\n")
}

found_likert <- intersect(likert_items, names(stu))
cat(sprintf("Found Likert items: %d/%d\n", length(found_likert), length(likert_items)))

cat("\nLoading PISA timing data...\n")
needed_tim <- c("CNTSTUID", "CNT", timing_cols)
tim <- read_sav("external_datasets/11_PISA_2018/TIM/CY07_MSU_STU_TIM.sav",
                col_select = any_of(needed_tim))
cat(sprintf("TIM loaded: %d rows x %d cols\n", nrow(tim), ncol(tim)))

found_timing <- intersect(timing_cols, names(tim))
cat(sprintf("Found %d/%d timing columns\n", length(found_timing), length(timing_cols)))

# ============================================================
# 3. Filter to target country and merge
# ============================================================
target_country <- "ITA"
cat(sprintf("\n=== Filtering to %s ===\n", target_country))

stu_c <- stu[stu$CNT == target_country, ]
tim_c <- tim[tim$CNT == target_country, ]
cat(sprintf("STU %s: %d rows\n", target_country, nrow(stu_c)))
cat(sprintf("TIM %s: %d rows\n", target_country, nrow(tim_c)))

# Merge on CNTSTUID
merged <- merge(stu_c, tim_c, by="CNTSTUID", suffixes=c("", ".tim"))
cat(sprintf("Merged: %d rows\n", nrow(merged)))

# Complete cases on Likert items only
item_mat <- as.data.frame(lapply(merged[, found_likert], as.numeric))
cc <- complete.cases(item_mat)
cat(sprintf("Complete Likert cases: %d (%.1f%%)\n", sum(cc), 100*mean(cc)))

merged_cc <- merged[cc, ]
item_mat_cc <- item_mat[cc, ]

# Free memory
rm(stu, tim, stu_c, tim_c)
gc()

# ============================================================
# 3b. Longstring pre-screening (shared for RR and Mah)
# ============================================================
# RR is broken by within-factor contiguous longstrings: straight-lining
# inflates coupled correlations, making careless respondents look GOOD.
# Standard longstring detection (max consecutive identical values) is
# applied as a shared first screen before both methods.

cat("\n=== Longstring pre-screening ===\n")
longstring_max <- apply(as.matrix(item_mat_cc), 1, function(x) {
  r <- rle(x)
  max(r$lengths)
})

cat(sprintf("Longstring distribution: mean=%.1f, median=%d, max=%d\n",
    mean(longstring_max), median(longstring_max), max(longstring_max)))
cat(sprintf("Quantiles: 90%%=%d  95%%=%d  99%%=%d\n",
    quantile(longstring_max, 0.90),
    quantile(longstring_max, 0.95),
    quantile(longstring_max, 0.99)))

# Threshold: flag if max consecutive identical > items_in_largest_scale
# For PISA, largest Likert scale has 9 items (ST186). Use 10 as threshold.
# This is conservative: a run of 10+ identical values across scale boundaries
# is very unlikely from attentive responding.
ls_threshold <- 10
ls_flagged <- longstring_max > ls_threshold
cat(sprintf("Longstring flagged (>%d): %d / %d (%.1f%%)\n",
    ls_threshold, sum(ls_flagged), length(ls_flagged), 100*mean(ls_flagged)))

# Save labels before filtering (for combined evaluation later)
pre_ls_n <- nrow(item_mat_cc)

# Remove longstring respondents from both data and merged (for C/IER weights)
merged_cc <- merged_cc[!ls_flagged, ]
item_mat_cc <- item_mat_cc[!ls_flagged, ]
cat(sprintf("After longstring removal: %d respondents\n", nrow(item_mat_cc)))

# ============================================================
# 4. Compute C/IER attentiveness weights per scale
# ============================================================
cat("\n=== Computing C/IER weights ===\n")

# For each scale with timing data, run Ulitzsch's method
# (simplified: no sampling weights for our purposes)
att_weights_per_scale <- list()
scales_processed <- 0

for (sc in mapping$scale_code) {
  tt_col <- paste0(sc, "_TT")
  if (!(tt_col %in% names(merged_cc))) next

  n_items_sc <- mapping$n_items[mapping$scale_code == sc]
  time_raw <- as.numeric(merged_cc[[tt_col]])

  # Skip if too many NAs
  valid <- !is.na(time_raw) & time_raw > 0
  if (sum(valid) < 50) {
    cat(sprintf("  %s: skipped (%d valid times)\n", sc, sum(valid)))
    next
  }

  # Log geometric mean per item (Ulitzsch's formula)
  time_log <- log((time_raw[valid] / 1000) ^ (1/n_items_sc))

  # Gaussian mixture decomposition
  tryCatch({
    BIC_res <- mclustBIC(time_log, G=1:9, modelNames="V")
    best_G <- which.max(BIC_res[,1])

    if (best_G == 1) {
      # Single component — no C/IER separation possible
      att <- rep(1.0, length(time_raw))
      att[!valid] <- NA
    } else {
      mod <- Mclust(time_log, G=best_G, modelNames="V")
      # C/IER component = smallest mean (fastest)
      indC <- which.min(mod$parameters$mean)
      # Attentiveness = 1 - P(C/IER)
      att_valid <- 1 - mod$z[, indC]
      att <- rep(NA_real_, length(time_raw))
      att[valid] <- att_valid
    }

    att_weights_per_scale[[sc]] <- att
    scales_processed <- scales_processed + 1
    if (scales_processed %% 10 == 0) cat(sprintf("  Processed %d scales...\n", scales_processed))
  }, error = function(e) {
    cat(sprintf("  %s: ERROR — %s\n", sc, e$message))
  })
}

cat(sprintf("Scales with attentiveness weights: %d\n", length(att_weights_per_scale)))

# Aggregate: mean attentiveness weight per respondent across scales
att_matrix <- do.call(cbind, att_weights_per_scale)
mean_att <- rowMeans(att_matrix, na.rm=TRUE)

cat(sprintf("\nMean attentiveness: %.3f (sd=%.3f)\n", mean(mean_att, na.rm=TRUE), sd(mean_att, na.rm=TRUE)))
qq <- quantile(mean_att, c(0.01, 0.05, 0.10, 0.25, 0.50), na.rm=TRUE)
cat(sprintf("Quantiles: 1%%=%.3f  5%%=%.3f  10%%=%.3f  25%%=%.3f  50%%=%.3f\n",
    qq[1], qq[2], qq[3], qq[4], qq[5]))

# Binary label: mean_att < 0.5 → careless
careless_label <- as.integer(mean_att < 0.5)
cat(sprintf("Careless (att < 0.5): %d / %d (%.1f%%)\n",
    sum(careless_label, na.rm=TRUE), sum(!is.na(careless_label)),
    100*mean(careless_label, na.rm=TRUE)))

# Also try stricter threshold
careless_strict <- as.integer(mean_att < 0.3)
cat(sprintf("Careless (att < 0.3): %d / %d (%.1f%%)\n",
    sum(careless_strict, na.rm=TRUE), sum(!is.na(careless_strict)),
    100*mean(careless_strict, na.rm=TRUE)))

# ============================================================
# 5. Run ReReReRe
# ============================================================
cat("\n=== Running ReReReRe ===\n")
cat(sprintf("Data: %d respondents x %d Likert items\n", nrow(item_mat_cc), ncol(item_mat_cc)))

rr <- ReReReRe(item_mat_cc, corProp=0.05, iterations=100, min_pairs=15, align_signs=TRUE)
rr$flagged_z <- rr$z_score <= 2.0

cat(sprintf("RR flagged (z<=2.0): %d / %d (%.1f%%)\n",
    sum(rr$flagged_z), nrow(item_mat_cc), 100*mean(rr$flagged_z)))

cat(sprintf("\nz_score summary: mean=%.2f, median=%.2f, sd=%.2f\n",
    mean(rr$z_score), median(rr$z_score), sd(rr$z_score)))

# ============================================================
# 6. Run Mahalanobis (practical threshold)
# ============================================================
cat("\n=== Running Mahalanobis ===\n")
mat_m <- as.matrix(item_mat_cc)
center <- colMeans(mat_m)
cov_mat <- cov(mat_m)
p <- ncol(mat_m)

# Regularize if needed
if (any(is.na(cov_mat)) || det(cov_mat) < 1e-10) {
  cat("Regularizing covariance matrix\n")
  cov_mat <- cov_mat + diag(0.01, p)
}

d2 <- mahalanobis(mat_m, center, cov_mat)
mah_thresh <- qchisq(0.999, df=p)
mah_flag <- d2 > mah_thresh

cat(sprintf("Mah threshold (chi-sq, p=%d): %.1f\n", p, mah_thresh))
cat(sprintf("Mah flagged: %d / %d (%.1f%%)\n",
    sum(mah_flag), length(mah_flag), 100*mean(mah_flag)))

# ============================================================
# 7. Evaluate against C/IER weights
# ============================================================
cat("\n=== Evaluation ===\n")

# --- AUC against CONTINUOUS attentiveness (most fair) ---
# Higher z_score = more attentive, higher att_weight = more attentive
# So we compute correlation between z_score and att_weight
cor_rr_att <- cor(rr$z_score, mean_att, use="complete.obs")
cor_mah_att <- cor(d2, mean_att, use="complete.obs")
cat(sprintf("Correlation with attentiveness:\n"))
cat(sprintf("  RR z_score:  r = %.3f (should be positive)\n", cor_rr_att))
cat(sprintf("  Mah D^2:     r = %.3f (should be negative)\n", cor_mah_att))

# --- AUC against binary label (att < 0.5) ---
compute_mcc <- function(pred, true) {
  tp <- sum(pred & true); tn <- sum(!pred & !true)
  fp <- sum(pred & !true); fn <- sum(!pred & true)
  num <- as.double(tp)*tn - as.double(fp)*fn
  den <- sqrt(as.double(tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  if (den == 0) return(0)
  num / den
}

valid_idx <- !is.na(careless_label)

cat(sprintf("\n--- Binary label: att < 0.5 (%d careless) ---\n",
    sum(careless_label[valid_idx])))

rr_auc <- tryCatch(
  as.numeric(auc(roc(careless_label[valid_idx], rr$z_score[valid_idx],
                     direction=">", quiet=TRUE))),
  error = function(e) NA)

mah_auc <- tryCatch(
  as.numeric(auc(roc(careless_label[valid_idx], d2[valid_idx],
                     direction="<", quiet=TRUE))),
  error = function(e) NA)

rr_mcc <- compute_mcc(rr$flagged_z[valid_idx], careless_label[valid_idx]==1)
mah_mcc <- compute_mcc(mah_flag[valid_idx], careless_label[valid_idx]==1)

rr_sens <- sum(rr$flagged_z[valid_idx] & careless_label[valid_idx]==1) /
           sum(careless_label[valid_idx]==1)
rr_spec <- sum(!rr$flagged_z[valid_idx] & careless_label[valid_idx]==0) /
           sum(careless_label[valid_idx]==0)
mah_sens <- sum(mah_flag[valid_idx] & careless_label[valid_idx]==1) /
            sum(careless_label[valid_idx]==1)
mah_spec <- sum(!mah_flag[valid_idx] & careless_label[valid_idx]==0) /
            sum(careless_label[valid_idx]==0)

cat(sprintf("ReReReRe:    MCC=%.3f  AUC=%.3f  Sens=%.3f  Spec=%.3f\n",
    rr_mcc, rr_auc, rr_sens, rr_spec))
cat(sprintf("Mahalanobis: MCC=%.3f  AUC=%.3f  Sens=%.3f  Spec=%.3f\n",
    mah_mcc, mah_auc, mah_sens, mah_spec))

# --- Also try stricter label (att < 0.3) ---
valid_s <- !is.na(careless_strict)
cat(sprintf("\n--- Binary label: att < 0.3 (%d careless) ---\n",
    sum(careless_strict[valid_s])))

rr_auc_s <- tryCatch(
  as.numeric(auc(roc(careless_strict[valid_s], rr$z_score[valid_s],
                     direction=">", quiet=TRUE))),
  error = function(e) NA)
mah_auc_s <- tryCatch(
  as.numeric(auc(roc(careless_strict[valid_s], d2[valid_s],
                     direction="<", quiet=TRUE))),
  error = function(e) NA)

cat(sprintf("ReReReRe:    AUC=%.3f\n", rr_auc_s))
cat(sprintf("Mahalanobis: AUC=%.3f\n", mah_auc_s))

# --- z_score by attentiveness groups ---
cat("\n=== z_score by attentiveness quartiles ===\n")
att_q <- cut(mean_att, breaks=quantile(mean_att, c(0, 0.10, 0.25, 0.50, 0.75, 1.0), na.rm=TRUE),
             include.lowest=TRUE, labels=c("0-10%", "10-25%", "25-50%", "50-75%", "75-100%"))
for (lev in levels(att_q)) {
  idx <- att_q == lev & !is.na(att_q)
  cat(sprintf("  %s (n=%d): z mean=%.2f, median=%.2f\n",
      lev, sum(idx), mean(rr$z_score[idx]), median(rr$z_score[idx])))
}

cat("\nDone.\n")
