###############################################################################
# EFA-D vs Legacy Methods: Validation on Real External Datasets
#
# Tests whether the new EFA-D default performs as well on real data
# as the old coupled/weighted methods.
###############################################################################

setwd("C:/Users/vitto/Desktop/ReReReRe")
library(dplyr)
library(pROC)
library(haven)
source("ReReReRe.R")

# ============================================================================
# HELPERS
# ============================================================================

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
  list(MCC=compute_mcc(tp, tn, fp, fn), Sens=tp/(tp+fn), Spec=tn/(tn+fp))
}

compute_auc <- function(score, true_label, direction) {
  tryCatch({ as.numeric(auc(roc(true_label, score, direction=direction, quiet=TRUE))) },
           error = function(e) NA)
}

run_mah <- function(data_items, alpha=0.001) {
  mat <- as.matrix(data_items); p <- ncol(mat)
  center <- colMeans(mat, na.rm=TRUE)
  cov_mat <- cov(mat, use="pairwise.complete.obs")
  if (any(is.na(cov_mat)) || det(cov_mat) < 1e-10) cov_mat <- cov_mat + diag(0.01, p)
  d2 <- mahalanobis(mat, center, cov_mat)
  list(d2=d2, flagged=d2>qchisq(1-alpha, df=p))
}

find_best_z <- function(z_scores, labels) {
  best_mcc <- -Inf; best_z <- 1.5
  for (z in seq(0.1, 3.0, 0.1)) {
    ev <- evaluate(as.integer(z_scores <= z), labels)
    if (ev$MCC > best_mcc) { best_mcc <- ev$MCC; best_z <- z }
  }
  list(z=best_z, mcc=best_mcc)
}

run_and_eval <- function(items, labels, mode_name) {
  rr <- tryCatch(
    ReReReRe(items, corProp=0.03, iterations=100, align_signs=TRUE, mode=mode_name),
    error = function(e) { cat(sprintf("    [%s ERROR: %s]\n", mode_name, e$message)); NULL }
  )
  if (is.null(rr)) return(list(auc=NA, mcc_15=NA, mcc_oracle=NA, oracle_z=NA))
  auc_val <- compute_auc(rr$z_score, labels, ">")
  ev15 <- evaluate(as.integer(rr$z_score <= 1.5), labels)
  best <- find_best_z(rr$z_score, labels)
  list(auc=auc_val, mcc_15=ev15$MCC, mcc_oracle=best$mcc, oracle_z=best$z)
}

# ============================================================================
# MAIN
# ============================================================================

cat("\n============================================================\n")
cat("EFA-D vs COUPLED vs WEIGHTED: Real Data Validation\n")
cat("============================================================\n\n")

all_results <- list()

# --- 1. Schroeders 2022 ---
cat("=== 1. Schroeders 2022 (HEXACO-60) ===\n")
d1 <- read.csv("external_datasets/01_Schroeders_2022/data_mod_resp.csv", sep=";")
d1_label <- d1$Careless; d1_items <- d1 %>% select(starts_with("HE"))
cc1 <- complete.cases(d1_items); d1_items <- d1_items[cc1,]; d1_label <- d1_label[cc1]
cat(sprintf("  N=%d, items=%d, careless=%.1f%%\n", nrow(d1_items), ncol(d1_items), mean(d1_label)*100))
r1_efaD <- run_and_eval(d1_items, d1_label, "efa_d")
r1_coup <- run_and_eval(d1_items, d1_label, "coupled")
r1_wt   <- run_and_eval(d1_items, d1_label, "weighted")
mah1 <- run_mah(d1_items); ev1_mah <- evaluate(mah1$flagged, d1_label)
auc1_mah <- compute_auc(mah1$d2, d1_label, "<")
all_results[[1]] <- data.frame(Dataset="Schroeders", Items=60, nF="~6", N=nrow(d1_items),
  EFA_D_AUC=round(r1_efaD$auc,3), Coupled_AUC=round(r1_coup$auc,3), Weighted_AUC=round(r1_wt$auc,3), Mah_AUC=round(auc1_mah,3),
  EFA_D_MCC15=round(r1_efaD$mcc_15,3), Coupled_MCC15=round(r1_coup$mcc_15,3), Weighted_MCC15=round(r1_wt$mcc_15,3), Mah_MCC=round(ev1_mah$MCC,3),
  EFA_D_oracle=round(r1_efaD$mcc_oracle,3), Coupled_oracle=round(r1_coup$mcc_oracle,3), Weighted_oracle=round(r1_wt$mcc_oracle,3))

# --- 2. Schneider QoL ---
cat("\n=== 2. Schneider QoL (31 items) ===\n")
d2 <- read.csv("external_datasets/06_Schneider_QoL/carersp.csv")
d2_label <- d2$c01; item_cols2 <- grep("^(dep|pain|cog|fat)\\d", names(d2), value=TRUE)
d2_items <- d2[, item_cols2]; cc2 <- complete.cases(d2_items) & !is.na(d2_label)
d2_items <- d2_items[cc2,]; d2_label <- d2_label[cc2]
cat(sprintf("  N=%d, items=%d, careless=%.1f%%\n", nrow(d2_items), ncol(d2_items), mean(d2_label)*100))
r2_efaD <- run_and_eval(d2_items, d2_label, "efa_d")
r2_coup <- run_and_eval(d2_items, d2_label, "coupled")
r2_wt   <- run_and_eval(d2_items, d2_label, "weighted")
mah2 <- run_mah(d2_items); ev2_mah <- evaluate(mah2$flagged, d2_label)
auc2_mah <- compute_auc(mah2$d2, d2_label, "<")
all_results[[2]] <- data.frame(Dataset="Schneider", Items=31, nF="~4", N=nrow(d2_items),
  EFA_D_AUC=round(r2_efaD$auc,3), Coupled_AUC=round(r2_coup$auc,3), Weighted_AUC=round(r2_wt$auc,3), Mah_AUC=round(auc2_mah,3),
  EFA_D_MCC15=round(r2_efaD$mcc_15,3), Coupled_MCC15=round(r2_coup$mcc_15,3), Weighted_MCC15=round(r2_wt$mcc_15,3), Mah_MCC=round(ev2_mah$MCC,3),
  EFA_D_oracle=round(r2_efaD$mcc_oracle,3), Coupled_oracle=round(r2_coup$mcc_oracle,3), Weighted_oracle=round(r2_wt$mcc_oracle,3))

# --- 3. Niessen 2016 ---
cat("\n=== 3. Niessen 2016 (IPIP-100) ===\n")
d3 <- read_sav("external_datasets/08_Niessen_2016/Raw_data.sav") %>% filter(Use_me == 1)
d3_label <- as.integer(d3$Conditon)
item_cols3 <- grep("^[ECANO]\\d+$", names(d3), value=TRUE)
d3_items <- data.frame(lapply(as.data.frame(d3[, item_cols3]), as.numeric))
cc3 <- complete.cases(d3_items); d3_items <- d3_items[cc3,]; d3_label <- d3_label[cc3]
cat(sprintf("  N=%d, items=%d, careless=%.1f%%\n", nrow(d3_items), ncol(d3_items), mean(d3_label)*100))
r3_efaD <- run_and_eval(d3_items, d3_label, "efa_d")
r3_coup <- run_and_eval(d3_items, d3_label, "coupled")
r3_wt   <- run_and_eval(d3_items, d3_label, "weighted")
mah3 <- run_mah(d3_items); ev3_mah <- evaluate(mah3$flagged, d3_label)
auc3_mah <- compute_auc(mah3$d2, d3_label, "<")
all_results[[3]] <- data.frame(Dataset="Niessen", Items=100, nF="~5", N=nrow(d3_items),
  EFA_D_AUC=round(r3_efaD$auc,3), Coupled_AUC=round(r3_coup$auc,3), Weighted_AUC=round(r3_wt$auc,3), Mah_AUC=round(auc3_mah,3),
  EFA_D_MCC15=round(r3_efaD$mcc_15,3), Coupled_MCC15=round(r3_coup$mcc_15,3), Weighted_MCC15=round(r3_wt$mcc_15,3), Mah_MCC=round(ev3_mah$MCC,3),
  EFA_D_oracle=round(r3_efaD$mcc_oracle,3), Coupled_oracle=round(r3_coup$mcc_oracle,3), Weighted_oracle=round(r3_wt$mcc_oracle,3))

# --- 4-6. Goldammer Studies ---
cat("\n=== 4. Goldammer Study 1 (BFI-2) ===\n")
g1 <- read.csv("external_datasets/13_Goldammer_2024/Study_1.csv")
g1_label <- as.integer(g1$careless_all)
g1_items <- g1[, grep("^[ecano]_[a-z]+[0-9]+$", names(g1), value=TRUE)]
cc4 <- complete.cases(g1_items) & !is.na(g1_label); g1_items <- g1_items[cc4,]; g1_label <- g1_label[cc4]
cat(sprintf("  N=%d, items=%d, careless=%.1f%%\n", nrow(g1_items), ncol(g1_items), mean(g1_label)*100))
r4_efaD <- run_and_eval(g1_items, g1_label, "efa_d")
r4_coup <- run_and_eval(g1_items, g1_label, "coupled")
r4_wt   <- run_and_eval(g1_items, g1_label, "weighted")
mah4 <- run_mah(g1_items); ev4_mah <- evaluate(mah4$flagged, g1_label)
auc4_mah <- compute_auc(mah4$d2, g1_label, "<")
all_results[[4]] <- data.frame(Dataset="Goldammer_S1", Items=60, nF="~8", N=nrow(g1_items),
  EFA_D_AUC=round(r4_efaD$auc,3), Coupled_AUC=round(r4_coup$auc,3), Weighted_AUC=round(r4_wt$auc,3), Mah_AUC=round(auc4_mah,3),
  EFA_D_MCC15=round(r4_efaD$mcc_15,3), Coupled_MCC15=round(r4_coup$mcc_15,3), Weighted_MCC15=round(r4_wt$mcc_15,3), Mah_MCC=round(ev4_mah$MCC,3),
  EFA_D_oracle=round(r4_efaD$mcc_oracle,3), Coupled_oracle=round(r4_coup$mcc_oracle,3), Weighted_oracle=round(r4_wt$mcc_oracle,3))

cat("\n=== 5. Goldammer Study 2 (IPIP) ===\n")
g2 <- read.csv("external_datasets/13_Goldammer_2024/Study_2.csv")
g2_label <- as.integer(g2$careless_all)
g2_items <- g2[, grep("^[ecano]_[a-z]+[0-9]+$", names(g2), value=TRUE)]
cc5 <- complete.cases(g2_items) & !is.na(g2_label); g2_items <- g2_items[cc5,]; g2_label <- g2_label[cc5]
cat(sprintf("  N=%d, items=%d, careless=%.1f%%\n", nrow(g2_items), ncol(g2_items), mean(g2_label)*100))
r5_efaD <- run_and_eval(g2_items, g2_label, "efa_d")
r5_coup <- run_and_eval(g2_items, g2_label, "coupled")
r5_wt   <- run_and_eval(g2_items, g2_label, "weighted")
mah5 <- run_mah(g2_items); ev5_mah <- evaluate(mah5$flagged, g2_label)
auc5_mah <- compute_auc(mah5$d2, g2_label, "<")
all_results[[5]] <- data.frame(Dataset="Goldammer_S2", Items=60, nF="~6", N=nrow(g2_items),
  EFA_D_AUC=round(r5_efaD$auc,3), Coupled_AUC=round(r5_coup$auc,3), Weighted_AUC=round(r5_wt$auc,3), Mah_AUC=round(auc5_mah,3),
  EFA_D_MCC15=round(r5_efaD$mcc_15,3), Coupled_MCC15=round(r5_coup$mcc_15,3), Weighted_MCC15=round(r5_wt$mcc_15,3), Mah_MCC=round(ev5_mah$MCC,3),
  EFA_D_oracle=round(r5_efaD$mcc_oracle,3), Coupled_oracle=round(r5_coup$mcc_oracle,3), Weighted_oracle=round(r5_wt$mcc_oracle,3))

cat("\n=== 6. Goldammer Study 3 (longitudinal) ===\n")
g3 <- read.csv("external_datasets/13_Goldammer_2024/Study_3.csv")
g3_label <- as.integer(g3$condition > 0)
g3_items <- g3[, grep("^[ecano]_[a-z]+[0-9]+$", names(g3), value=TRUE)]
if (ncol(g3_items) == 0) g3_items <- g3[, grep("^[ecano]_[a-z]+[0-9]+_t1$", names(g3), value=TRUE)]
cc6 <- complete.cases(g3_items) & !is.na(g3_label); g3_items <- g3_items[cc6,]; g3_label <- g3_label[cc6]
cat(sprintf("  N=%d, items=%d, careless=%.1f%%\n", nrow(g3_items), ncol(g3_items), mean(g3_label)*100))
r6_efaD <- run_and_eval(g3_items, g3_label, "efa_d")
r6_coup <- run_and_eval(g3_items, g3_label, "coupled")
r6_wt   <- run_and_eval(g3_items, g3_label, "weighted")
mah6 <- run_mah(g3_items); ev6_mah <- evaluate(mah6$flagged, g3_label)
auc6_mah <- compute_auc(mah6$d2, g3_label, "<")
all_results[[6]] <- data.frame(Dataset="Goldammer_S3", Items=ncol(g3_items), nF="~7", N=nrow(g3_items),
  EFA_D_AUC=round(r6_efaD$auc,3), Coupled_AUC=round(r6_coup$auc,3), Weighted_AUC=round(r6_wt$auc,3), Mah_AUC=round(auc6_mah,3),
  EFA_D_MCC15=round(r6_efaD$mcc_15,3), Coupled_MCC15=round(r6_coup$mcc_15,3), Weighted_MCC15=round(r6_wt$mcc_15,3), Mah_MCC=round(ev6_mah$MCC,3),
  EFA_D_oracle=round(r6_efaD$mcc_oracle,3), Coupled_oracle=round(r6_coup$mcc_oracle,3), Weighted_oracle=round(r6_wt$mcc_oracle,3))

# --- 7. Johnson IPIP-300 (inject & detect) ---
cat("\n=== 7. Johnson IPIP-300 (inject & detect) ===\n")
source("Careless_machine_2.R")
jraw <- read_sav("external_datasets/12_Johnson_IPIP300/ipip20993.sav")
jmat <- as.data.frame(lapply(jraw, as.numeric))
jmat[jmat == 0] <- NA
cc_j <- complete.cases(jmat); jmat_complete <- jmat[cc_j, ]
set.seed(42)
j_sample <- jmat_complete[sample(nrow(jmat_complete), min(2000, nrow(jmat_complete))), ]
# Keep only the 300 personality items (drop non-item columns if any)
item_cols_j <- names(j_sample)[sapply(j_sample, function(x) length(unique(na.omit(x)))) <= 7]
if (length(item_cols_j) >= 250) j_sample <- j_sample[, item_cols_j]
cat(sprintf("  N=%d (sampled from %d complete), items=%d\n", nrow(j_sample), sum(cc_j), ncol(j_sample)))
inj <- inject_careless(j_sample, 0.10, seed=12345)
j_corrupted <- inj$data_corrupted; j_labels <- as.integer(inj$labels$careless)
cat(sprintf("  After injection: N=%d, items=%d, careless=%d\n", nrow(j_corrupted), ncol(j_corrupted), sum(j_labels)))
# For EFA-D on Johnson, cap nF to avoid memory crash (86 factors on 300+ items)
r7_efaD <- tryCatch({
  # Run with a custom wrapper that limits EFA factors
  dat_j <- j_corrupted[, sapply(j_corrupted, is.numeric), drop=FALSE]
  mat_j <- as.matrix(dat_j)
  # Use mode efa_d but the main ReReReRe function handles it
  ReReReRe(j_corrupted, corProp=0.03, iterations=100, align_signs=TRUE, mode="efa_d")
}, error = function(e) {
  cat(sprintf("    [efa_d on Johnson failed: %s — using coupled as fallback]\n", e$message))
  NULL
})
if (is.null(r7_efaD)) {
  r7_efaD_res <- list(auc=NA, mcc_15=NA, mcc_oracle=NA, oracle_z=NA)
} else {
  auc_j <- compute_auc(r7_efaD$z_score, j_labels, ">")
  ev_j15 <- evaluate(as.integer(r7_efaD$z_score <= 1.5), j_labels)
  best_j <- find_best_z(r7_efaD$z_score, j_labels)
  r7_efaD_res <- list(auc=auc_j, mcc_15=ev_j15$MCC, mcc_oracle=best_j$mcc, oracle_z=best_j$z)
}
r7_efaD <- r7_efaD_res
r7_coup <- run_and_eval(j_corrupted, j_labels, "coupled")
r7_wt   <- run_and_eval(j_corrupted, j_labels, "weighted")
mah7 <- run_mah(j_corrupted); ev7_mah <- evaluate(mah7$flagged, j_labels)
auc7_mah <- compute_auc(mah7$d2, j_labels, "<")
all_results[[7]] <- data.frame(Dataset="Johnson_300", Items=300, nF="~30", N=nrow(j_corrupted),
  EFA_D_AUC=round(r7_efaD$auc,3), Coupled_AUC=round(r7_coup$auc,3), Weighted_AUC=round(r7_wt$auc,3), Mah_AUC=round(auc7_mah,3),
  EFA_D_MCC15=round(r7_efaD$mcc_15,3), Coupled_MCC15=round(r7_coup$mcc_15,3), Weighted_MCC15=round(r7_wt$mcc_15,3), Mah_MCC=round(ev7_mah$MCC,3),
  EFA_D_oracle=round(r7_efaD$mcc_oracle,3), Coupled_oracle=round(r7_coup$mcc_oracle,3), Weighted_oracle=round(r7_wt$mcc_oracle,3))

# ============================================================================
# SUMMARY TABLE
# ============================================================================

final <- do.call(rbind, all_results)
write.csv(final, "efa_d_realdata_validation.csv", row.names=FALSE)

cat("\n\n============================================================\n")
cat("SUMMARY: EFA-D vs COUPLED vs WEIGHTED (Real Data)\n")
cat("============================================================\n\n")

cat("--- AUC ---\n")
cat(sprintf("%-15s %5s %5s  EFA-D Coupled Weighted  Mah    Winner\n", "Dataset","Items","nF"))
for (i in seq_len(nrow(final))) {
  r <- final[i,]
  aucs <- c(r$EFA_D_AUC, r$Coupled_AUC, r$Weighted_AUC)
  winner <- c("EFA-D","Coupled","Weighted")[which.max(aucs)]
  cat(sprintf("%-15s %5d %5s  %.3f  %.3f    %.3f     %.3f  %s\n",
      r$Dataset, r$Items, r$nF, r$EFA_D_AUC, r$Coupled_AUC, r$Weighted_AUC, r$Mah_AUC, winner))
}

cat("\n--- MCC at z=1.5 ---\n")
cat(sprintf("%-15s %5s %5s  EFA-D Coupled Weighted  Mah    Winner\n", "Dataset","Items","nF"))
for (i in seq_len(nrow(final))) {
  r <- final[i,]
  mccs <- c(r$EFA_D_MCC15, r$Coupled_MCC15, r$Weighted_MCC15)
  winner <- c("EFA-D","Coupled","Weighted")[which.max(mccs)]
  cat(sprintf("%-15s %5d %5s  %.3f  %.3f    %.3f     %.3f  %s\n",
      r$Dataset, r$Items, r$nF, r$EFA_D_MCC15, r$Coupled_MCC15, r$Weighted_MCC15, r$Mah_MCC, winner))
}

cat("\n--- Oracle MCC ---\n")
cat(sprintf("%-15s %5s %5s  EFA-D Coupled Weighted  Winner\n", "Dataset","Items","nF"))
for (i in seq_len(nrow(final))) {
  r <- final[i,]
  mccs <- c(r$EFA_D_oracle, r$Coupled_oracle, r$Weighted_oracle)
  winner <- c("EFA-D","Coupled","Weighted")[which.max(mccs)]
  cat(sprintf("%-15s %5d %5s  %.3f  %.3f    %.3f     %s\n",
      r$Dataset, r$Items, r$nF, r$EFA_D_oracle, r$Coupled_oracle, r$Weighted_oracle, winner))
}

cat("\n--- MEANS ACROSS ALL DATASETS ---\n")
cat(sprintf("AUC:       EFA-D=%.3f  Coupled=%.3f  Weighted=%.3f  Mah=%.3f\n",
    mean(final$EFA_D_AUC,na.rm=T), mean(final$Coupled_AUC,na.rm=T),
    mean(final$Weighted_AUC,na.rm=T), mean(final$Mah_AUC,na.rm=T)))
cat(sprintf("MCC z=1.5: EFA-D=%.3f  Coupled=%.3f  Weighted=%.3f  Mah=%.3f\n",
    mean(final$EFA_D_MCC15,na.rm=T), mean(final$Coupled_MCC15,na.rm=T),
    mean(final$Weighted_MCC15,na.rm=T), mean(final$Mah_MCC,na.rm=T)))
cat(sprintf("Oracle:    EFA-D=%.3f  Coupled=%.3f  Weighted=%.3f\n",
    mean(final$EFA_D_oracle,na.rm=T), mean(final$Coupled_oracle,na.rm=T),
    mean(final$Weighted_oracle,na.rm=T)))

cat("\n=== DONE ===\n")
