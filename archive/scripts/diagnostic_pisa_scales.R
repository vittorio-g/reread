#### PISA Scale Diagnostic ####
# Check per-item observed min/max and test 0-1 rescaling

library(haven)
library(mclust)
library(pROC)
source("ReReReRe.R")

# === 1. Load and filter (same as PISA_Validation.R) ===
mapping <- read.csv("external_datasets/11_PISA_2018/PISA_2018_scale_mapping.csv",
                    stringsAsFactors = FALSE)

parse_items <- function(pattern, n_items) {
  parts <- strsplit(pattern, "-")[[1]]
  first <- parts[1]
  prefix <- sub("(ST\\d+Q)\\d+(\\w+)", "\\1", first)
  suffix <- sub("(ST\\d+Q)\\d+(\\w+)", "\\2", first)
  sprintf("%s%02d%s", prefix, 1:n_items, suffix)
}

scale_items <- list()
all_item_cols <- c()
for (i in 1:nrow(mapping)) {
  items <- parse_items(mapping$item_pattern[i], mapping$n_items[i])
  scale_items[[mapping$scale_code[i]]] <- items
  all_item_cols <- c(all_item_cols, items)
}

likert_scales <- mapping$scale_code[mapping$n_options > 2]
likert_items <- unlist(scale_items[likert_scales])

needed_stu <- c("CNTSTUID", "CNT", all_item_cols)
stu <- read_sav("external_datasets/11_PISA_2018/STU/CY07_MSU_STU_QQQ.sav",
                col_select = any_of(needed_stu))
needed_tim <- c("CNTSTUID", "CNT", paste0(mapping$scale_code, "_TT"))
tim <- read_sav("external_datasets/11_PISA_2018/TIM/CY07_MSU_STU_TIM.sav",
                col_select = any_of(needed_tim))

found_likert <- intersect(likert_items, names(stu))
stu_c <- stu[stu$CNT == "ITA", ]
tim_c <- tim[tim$CNT == "ITA", ]
merged <- merge(stu_c, tim_c, by="CNTSTUID", suffixes=c("", ".tim"))
item_mat <- as.data.frame(lapply(merged[, found_likert], as.numeric))
cc <- complete.cases(item_mat)
merged_cc <- merged[cc, ]
item_mat_cc <- item_mat[cc, ]
rm(stu, tim, stu_c, tim_c); gc()

# Longstring filter
longstring_max <- apply(as.matrix(item_mat_cc), 1, function(x) max(rle(x)$lengths))
ls_flagged <- longstring_max > 10
merged_cc <- merged_cc[!ls_flagged, ]
item_mat_cc <- item_mat_cc[!ls_flagged, ]
cat(sprintf("After filtering: %d respondents x %d items\n\n", nrow(item_mat_cc), ncol(item_mat_cc)))

# === 2. Per-item scale diagnostics ===
cat("=== PER-ITEM OBSERVED RANGES ===\n\n")
mat <- as.matrix(item_mat_cc)
item_min <- apply(mat, 2, min, na.rm=TRUE)
item_max <- apply(mat, 2, max, na.rm=TRUE)
item_nuniq <- apply(mat, 2, function(x) length(unique(x)))

# Map items to scales and expected n_options
item_to_scale <- character(length(found_likert))
item_expected_max <- integer(length(found_likert))
names(item_to_scale) <- found_likert
names(item_expected_max) <- found_likert
for (sc in likert_scales) {
  sc_items <- intersect(scale_items[[sc]], found_likert)
  item_to_scale[sc_items] <- sc
  item_expected_max[sc_items] <- mapping$n_options[mapping$scale_code == sc]
}

# Summary by scale
cat(sprintf("%-12s %-8s %-6s %-6s %-6s %-8s %-6s %s\n",
    "Scale", "Expected", "ObsMin", "ObsMax", "Unique", "Match?", "Items", "Problem items"))
cat(paste(rep("-", 90), collapse=""), "\n")

problem_items <- c()
for (sc in likert_scales) {
  sc_items <- intersect(scale_items[[sc]], found_likert)
  if (length(sc_items) == 0) next
  exp_max <- mapping$n_options[mapping$scale_code == sc]
  obs_mins <- item_min[sc_items]
  obs_maxs <- item_max[sc_items]
  obs_uniqs <- item_nuniq[sc_items]

  all_match <- all(obs_mins == 1) && all(obs_maxs == exp_max)
  probs <- sc_items[obs_maxs != exp_max | obs_mins != 1]
  if (length(probs) > 0) problem_items <- c(problem_items, probs)

  cat(sprintf("%-12s %-8d %-6s %-6s %-6s %-8s %-6d %s\n",
      sc, exp_max,
      paste(unique(obs_mins), collapse="/"),
      paste(unique(obs_maxs), collapse="/"),
      paste(range(obs_uniqs), collapse="-"),
      ifelse(all_match, "OK", "MISMATCH"),
      length(sc_items),
      ifelse(length(probs) > 0, paste(probs, collapse=", "), "")))
}

if (length(problem_items) > 0) {
  cat(sprintf("\n%d problem items found:\n", length(problem_items)))
  for (it in problem_items) {
    vals <- sort(unique(mat[, it]))
    cat(sprintf("  %s (scale %s, expected 1-%d): observed values = %s\n",
        it, item_to_scale[it], item_expected_max[it],
        paste(vals, collapse=", ")))
  }
} else {
  cat("\nAll items match expected ranges.\n")
}

# === 3. Compute C/IER weights (same as PISA_Validation.R) ===
cat("\n=== Computing C/IER weights ===\n")
att_weights_per_scale <- list()
for (sc in mapping$scale_code) {
  tt_col <- paste0(sc, "_TT")
  if (!(tt_col %in% names(merged_cc))) next
  n_items_sc <- mapping$n_items[mapping$scale_code == sc]
  time_raw <- as.numeric(merged_cc[[tt_col]])
  valid <- !is.na(time_raw) & time_raw > 0
  if (sum(valid) < 50) next
  time_log <- log((time_raw[valid] / 1000) ^ (1/n_items_sc))
  tryCatch({
    BIC_res <- mclustBIC(time_log, G=1:9, modelNames="V")
    best_G <- which.max(BIC_res[,1])
    if (best_G == 1) {
      att <- rep(1.0, length(time_raw)); att[!valid] <- NA
    } else {
      mod <- Mclust(time_log, G=best_G, modelNames="V")
      indC <- which.min(mod$parameters$mean)
      att_valid <- 1 - mod$z[, indC]
      att <- rep(NA_real_, length(time_raw)); att[valid] <- att_valid
    }
    att_weights_per_scale[[sc]] <- att
  }, error = function(e) NULL)
}
att_matrix <- do.call(cbind, att_weights_per_scale)
mean_att <- rowMeans(att_matrix, na.rm=TRUE)
careless_label <- as.integer(mean_att < 0.5)
valid_idx <- !is.na(careless_label)

compute_mcc <- function(pred, true) {
  tp <- sum(pred & true); tn <- sum(!pred & !true)
  fp <- sum(pred & !true); fn <- sum(!pred & true)
  num <- as.double(tp)*tn - as.double(fp)*fn
  den <- sqrt(as.double(tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  if (den == 0) return(0)
  num / den
}

# === 4. Run RR on RAW data (baseline — already done but repeat for comparison) ===
cat("\n=== RR on RAW data (no rescaling) ===\n")
rr_raw <- ReReReRe(item_mat_cc, corProp=0.05, iterations=100, min_pairs=15, align_signs=TRUE)

rr_raw_auc <- as.numeric(auc(roc(careless_label[valid_idx], rr_raw$z_score[valid_idx],
                                  direction=">", quiet=TRUE)))
cat(sprintf("RAW:  AUC=%.3f  z_mean=%.2f  z_sd=%.2f\n",
    rr_raw_auc, mean(rr_raw$z_score), sd(rr_raw$z_score)))

# === 5. Rescale to 0-1 per item using THEORETICAL range ===
cat("\n=== RR on 0-1 rescaled data (theoretical range) ===\n")
item_mat_01 <- item_mat_cc
for (it in found_likert) {
  exp_max <- item_expected_max[it]
  if (is.na(exp_max) || exp_max <= 1) next
  # Theoretical range: 1 to exp_max → rescale to 0-1
  item_mat_01[[it]] <- (item_mat_01[[it]] - 1) / (exp_max - 1)
}

# Verify rescaling
cat("After rescaling — sample of ranges:\n")
mat01 <- as.matrix(item_mat_01)
for (sc in likert_scales[1:5]) {
  sc_items <- intersect(scale_items[[sc]], found_likert)
  if (length(sc_items) == 0) next
  mins01 <- apply(mat01[, sc_items, drop=FALSE], 2, min)
  maxs01 <- apply(mat01[, sc_items, drop=FALSE], 2, max)
  cat(sprintf("  %s (was 1-%d): now %.2f-%.2f\n",
      sc, mapping$n_options[mapping$scale_code == sc],
      min(mins01), max(maxs01)))
}

rr_01 <- ReReReRe(item_mat_01, corProp=0.05, iterations=100, min_pairs=15, align_signs=TRUE)

rr_01_auc <- as.numeric(auc(roc(careless_label[valid_idx], rr_01$z_score[valid_idx],
                                 direction=">", quiet=TRUE)))
cat(sprintf("0-1:  AUC=%.3f  z_mean=%.2f  z_sd=%.2f\n",
    rr_01_auc, mean(rr_01$z_score), sd(rr_01$z_score)))

# === 6. Compare across z thresholds ===
cat("\n=== MCC across z thresholds ===\n")
cat(sprintf("%-8s %-10s %-10s %-10s %-10s\n", "z_thr", "RAW_MCC", "01_MCC", "RAW_Sens", "01_Sens"))
for (z in seq(0, 5, 0.5)) {
  flag_raw <- rr_raw$z_score[valid_idx] <= z
  flag_01  <- rr_01$z_score[valid_idx] <= z
  mcc_raw <- compute_mcc(flag_raw, careless_label[valid_idx]==1)
  mcc_01  <- compute_mcc(flag_01, careless_label[valid_idx]==1)
  sens_raw <- sum(flag_raw & careless_label[valid_idx]==1) / sum(careless_label[valid_idx]==1)
  sens_01  <- sum(flag_01 & careless_label[valid_idx]==1) / sum(careless_label[valid_idx]==1)
  cat(sprintf("%-8.1f %-10.3f %-10.3f %-10.3f %-10.3f\n", z, mcc_raw, mcc_01, sens_raw, sens_01))
}

# === 7. Mahalanobis for reference ===
cat("\n=== Mahalanobis reference ===\n")
mat_m <- as.matrix(item_mat_cc)
cov_mat <- cov(mat_m)
if (det(cov_mat) < 1e-10) cov_mat <- cov_mat + diag(0.01, ncol(mat_m))
d2 <- mahalanobis(mat_m, colMeans(mat_m), cov_mat)
mah_auc <- as.numeric(auc(roc(careless_label[valid_idx], d2[valid_idx], direction="<", quiet=TRUE)))
cat(sprintf("Mah AUC=%.3f\n", mah_auc))

cat("\nDone.\n")
