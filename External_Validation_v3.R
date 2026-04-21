#### External Validation v3 (2026-04-21) ####
# Re-validation after algorithm changes (2026-03-30 to 2026-04-01):
#   - mode="auto" default: weighted (≤60 items), coupled (>60 items)
#   - Split into ReReReRe() + ReReReRe_F()
#
# Tests:
#   - ReReReRe()        (auto mode, default)
#   - ReReReRe(mode="coupled")  (force coupled, pre-fix comparison)
#   - ReReReRe_F()      (per-factor EFA variant)
#   - Mahalanobis chi-sq(.001) practical threshold
#
# Datasets: Schroeders, Schneider, Niessen (3 with ground truth)

setwd("C:/Users/User/OneDrive - CNR/Claude/ReReReRe")
library(dplyr)
library(pROC)
library(haven)
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

evaluate <- function(predicted, true_label) {
  tp <- sum(predicted == 1 & true_label == 1, na.rm=TRUE)
  tn <- sum(predicted == 0 & true_label == 0, na.rm=TRUE)
  fp <- sum(predicted == 1 & true_label == 0, na.rm=TRUE)
  fn <- sum(predicted == 0 & true_label == 1, na.rm=TRUE)
  list(MCC=compute_mcc(tp, tn, fp, fn),
       Sens = if ((tp+fn)>0) tp/(tp+fn) else NA,
       Spec = if ((tn+fp)>0) tn/(tn+fp) else NA)
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
    if (!is.na(ev$MCC) && ev$MCC > best_mcc) { best_mcc <- ev$MCC; best_z <- z }
  }
  list(z=best_z, mcc=best_mcc)
}

# Run a single method and collect metrics
run_method <- function(items, labels, method_name, rr_call) {
  rr <- rr_call(items)
  auc_val <- compute_auc(rr$z_score, labels, ">")
  best <- find_best_z(rr$z_score, labels)
  ev15 <- evaluate(as.integer(rr$z_score <= 1.5), labels)
  ev_auto <- if ("flagged" %in% names(rr)) evaluate(as.integer(rr$flagged), labels)
             else list(MCC=NA, Sens=NA, Spec=NA)
  mode_used <- if ("mode_used" %in% names(rr)) rr$mode_used[1] else NA
  data.frame(Method=method_name, Mode=mode_used,
             AUC=round(auc_val,3),
             MCC_z15=round(ev15$MCC,3), Sens_z15=round(ev15$Sens,3), Spec_z15=round(ev15$Spec,3),
             MCC_oracle=round(best$mcc,3), Oracle_z=best$z)
}

cat("\n====================================================\n")
cat("   EXTERNAL VALIDATION v3 (2026-04-21)\n")
cat("   NEW ARCHITECTURE: ReReReRe() auto-switch\n")
cat("   + ReReReRe_F() per-factor + Mahalanobis\n")
cat("====================================================\n\n")

all_results <- list()

# ============================================================
# DATASET 1: Schroeders 2022 (60 items, ~6 factors, experimental)
# ============================================================
cat("=== 1. Schroeders 2022 (60 items, ~6 factors) ===\n")
d1 <- read.csv("external_datasets/01_Schroeders_2022/data_mod_resp.csv", sep=";")
d1_label <- d1$Careless
d1_items <- d1 %>% select(starts_with("HE"))
ok1 <- complete.cases(d1_items)
d1_items <- d1_items[ok1,]; d1_label <- d1_label[ok1]
cat(sprintf("N=%d items=%d careless=%d (%.1f%%)\n",
    nrow(d1_items), ncol(d1_items), sum(d1_label), mean(d1_label)*100))

set.seed(42)
r1 <- rbind(
  run_method(d1_items, d1_label, "RR_auto",
             function(x) ReReReRe(x, iterations=100, align_signs=TRUE)),
  run_method(d1_items, d1_label, "RR_coupled",
             function(x) ReReReRe(x, iterations=100, align_signs=TRUE, mode="coupled")),
  run_method(d1_items, d1_label, "RR_F",
             function(x) ReReReRe_F(x, iterations=100, align_signs=TRUE))
)
mah1 <- run_mah(d1_items)
ev_mah1 <- evaluate(as.integer(mah1$flagged), d1_label)
auc_mah1 <- compute_auc(mah1$d2, d1_label, "<")
r1 <- rbind(r1, data.frame(Method="Mahalanobis", Mode="chi-sq.001",
             AUC=round(auc_mah1,3),
             MCC_z15=round(ev_mah1$MCC,3), Sens_z15=round(ev_mah1$Sens,3), Spec_z15=round(ev_mah1$Spec,3),
             MCC_oracle=NA, Oracle_z=NA))
r1 <- cbind(Dataset="Schroeders2022", r1)
print(r1); cat("\n")
all_results[[1]] <- r1

# ============================================================
# DATASET 2: Schneider QoL (31 items, 4 factors, latent class)
# ============================================================
cat("=== 2. Schneider QoL (31 items, 4 factors) ===\n")
d2 <- read.csv("external_datasets/06_Schneider_QoL/carersp.csv")
d2_label <- d2$c01
item_cols2 <- grep("^(dep|pain|cog|fat)\\d", names(d2), value=TRUE)
d2_items <- d2[, item_cols2]
ok2 <- complete.cases(d2_items) & !is.na(d2_label)
d2_items <- d2_items[ok2,]; d2_label <- d2_label[ok2]
cat(sprintf("N=%d items=%d careless=%d (%.1f%%)\n",
    nrow(d2_items), ncol(d2_items), sum(d2_label), mean(d2_label)*100))

set.seed(42)
r2 <- rbind(
  run_method(d2_items, d2_label, "RR_auto",
             function(x) ReReReRe(x, iterations=100, align_signs=TRUE)),
  run_method(d2_items, d2_label, "RR_coupled",
             function(x) ReReReRe(x, iterations=100, align_signs=TRUE, mode="coupled")),
  run_method(d2_items, d2_label, "RR_F",
             function(x) ReReReRe_F(x, iterations=100, align_signs=TRUE))
)
mah2 <- run_mah(d2_items)
ev_mah2 <- evaluate(as.integer(mah2$flagged), d2_label)
auc_mah2 <- compute_auc(mah2$d2, d2_label, "<")
r2 <- rbind(r2, data.frame(Method="Mahalanobis", Mode="chi-sq.001",
             AUC=round(auc_mah2,3),
             MCC_z15=round(ev_mah2$MCC,3), Sens_z15=round(ev_mah2$Sens,3), Spec_z15=round(ev_mah2$Spec,3),
             MCC_oracle=NA, Oracle_z=NA))
r2 <- cbind(Dataset="Schneider_QoL", r2)
print(r2); cat("\n")
all_results[[2]] <- r2

# ============================================================
# DATASET 3: Niessen 2016 (100 items, 5 factors, speed)
# ============================================================
cat("=== 3. Niessen 2016 (100 items, 5 factors) ===\n")
d3 <- read_sav("external_datasets/08_Niessen_2016/Raw_data.sav")
d3 <- d3 %>% filter(Use_me == 1)
d3_label <- as.integer(d3$Conditon)
item_cols3 <- grep("^[ECANO]\\d+$", names(d3), value=TRUE)
d3_items <- as.data.frame(d3[, item_cols3])
d3_items <- data.frame(lapply(d3_items, as.numeric))
ok3 <- complete.cases(d3_items)
d3_items <- d3_items[ok3,]; d3_label <- d3_label[ok3]
cat(sprintf("N=%d items=%d careless=%d (%.1f%%)\n",
    nrow(d3_items), ncol(d3_items), sum(d3_label), mean(d3_label)*100))

set.seed(42)
r3 <- rbind(
  run_method(d3_items, d3_label, "RR_auto",
             function(x) ReReReRe(x, iterations=100, align_signs=TRUE)),
  run_method(d3_items, d3_label, "RR_coupled",
             function(x) ReReReRe(x, iterations=100, align_signs=TRUE, mode="coupled")),
  run_method(d3_items, d3_label, "RR_weighted",
             function(x) ReReReRe(x, iterations=100, align_signs=TRUE, mode="weighted")),
  run_method(d3_items, d3_label, "RR_F",
             function(x) ReReReRe_F(x, iterations=100, align_signs=TRUE))
)
mah3 <- run_mah(d3_items)
ev_mah3 <- evaluate(as.integer(mah3$flagged), d3_label)
auc_mah3 <- compute_auc(mah3$d2, d3_label, "<")
r3 <- rbind(r3, data.frame(Method="Mahalanobis", Mode="chi-sq.001",
             AUC=round(auc_mah3,3),
             MCC_z15=round(ev_mah3$MCC,3), Sens_z15=round(ev_mah3$Sens,3), Spec_z15=round(ev_mah3$Spec,3),
             MCC_oracle=NA, Oracle_z=NA))
r3 <- cbind(Dataset="Niessen2016", r3)
print(r3); cat("\n")
all_results[[3]] <- r3

# ============================================================
# SUMMARY TABLE
# ============================================================
cat("\n====================================================\n")
cat("   SUMMARY\n")
cat("====================================================\n\n")
summary_df <- do.call(rbind, all_results)
write.csv(summary_df, "external_validation_v3.csv", row.names=FALSE)
print(summary_df)

cat("\nSaved: external_validation_v3.csv\n")
cat("\n=== DONE ===\n")
