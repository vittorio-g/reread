###############################################################################
# Pennycook & Rand — Method comparison with max 20% flagging
###############################################################################
setwd("C:/Users/vitto/Desktop/ReReReRe")
library(dplyr); library(tidyr); library(psych)
source("ReReReRe.R")

cat("\n============================================================\n")
cat("PENNYCOOK & RAND — MAX 20% FLAGGING\n")
cat("============================================================\n\n")

d1 <- read.csv("Dataset/dataset_1/Pennycook & Rand (Study 1).csv")
d2 <- read.csv("Dataset/dataset_1/Pennycook & Rand (Study 2).csv")

acc_items1 <- c(paste0("Fake", 1:15, "_Accurate"), paste0("Real", 1:15, "_Accurate"))
items1 <- d1[, acc_items1]; cc1 <- complete.cases(items1)
items1_c <- items1[cc1, ]; d1_c <- d1[cc1, ]

acc_items2 <- c(paste0("Fake", 1:12, "_Accurate"), paste0("Real", 1:12, "_Accurate"))
items2 <- d2[, acc_items2]; d2$CRT <- d2$CRT_ACC; cc2 <- complete.cases(items2)
items2_c <- items2[cc2, ]; d2_c <- d2[cc2, ]

cat(sprintf("Study 1: N=%d, %d items\n", nrow(items1_c), ncol(items1_c)))
cat(sprintf("Study 2: N=%d, %d items\n", nrow(items2_c), ncol(items2_c)))

# Run all methods
cat("\n--- Running methods ---\n")
run_m <- function(items, mode_name) {
  cat(sprintf("  %s... ", mode_name))
  rr <- ReReReRe(items, corProp=0.03, iterations=100, align_signs=TRUE, mode=mode_name)
  cat("done\n"); rr
}

rr1_c <- run_m(items1_c, "coupled");  rr1_w <- run_m(items1_c, "weighted"); rr1_e <- run_m(items1_c, "efa_d")
rr2_c <- run_m(items2_c, "coupled");  rr2_w <- run_m(items2_c, "weighted"); rr2_e <- run_m(items2_c, "efa_d")

# Find z threshold that flags exactly ~20%
find_z_for_pct <- function(z_scores, target_pct=0.20) {
  quantile(z_scores, probs=target_pct)
}

cat("\n--- Z thresholds for 20% flagging ---\n")
methods_s1 <- list(coupled=rr1_c, weighted=rr1_w, efa_d=rr1_e)
methods_s2 <- list(coupled=rr2_c, weighted=rr2_w, efa_d=rr2_e)

for (pct in c(0.05, 0.10, 0.15, 0.20)) {
  cat(sprintf("\n  Target: %.0f%% flagged\n", pct*100))
  for (mn in names(methods_s1)) {
    z1 <- find_z_for_pct(methods_s1[[mn]]$z_score, pct)
    z2 <- find_z_for_pct(methods_s2[[mn]]$z_score, pct)
    cat(sprintf("    %-10s S1: z<=%.3f  S2: z<=%.3f\n", mn, z1, z2))
  }
}

# ============================================================================
# Evaluate at 5%, 10%, 15%, 20% flagging
# ============================================================================

eval_at_pct <- function(d, z_scores, target_pct, items_mat) {
  z_cut <- quantile(z_scores, probs=target_pct)
  flagged <- z_scores <= z_cut
  n_flag <- sum(flagged)
  clean <- d[!flagged, ]

  if (nrow(clean) < 30) return(NULL)

  r_all   <- cor(d$CRT, d$Discernment, use="complete.obs")
  r_clean <- cor(clean$CRT, clean$Discernment, use="complete.obs")

  reg_all   <- lm(Discernment ~ Age + Sex + Education + Conserv + CRT, data=d)
  reg_clean <- lm(Discernment ~ Age + Sex + Education + Conserv + CRT, data=clean)

  # Per-group discernment correlations
  disc_cols <- intersect(c("L_Discernment", "C_Discernment", "N_Discernment"), names(d))

  cl_all <- d[d$ClintonTrump==1,]; tr_all <- d[d$ClintonTrump==2,]
  cl_cln <- clean[clean$ClintonTrump==1,]; tr_cln <- clean[clean$ClintonTrump==2,]

  # Mean |r| across all discernment correlations (both groups)
  all_rs <- clean_rs <- c()
  for (dc in disc_cols) {
    if (nrow(cl_all)>10) all_rs <- c(all_rs, cor(cl_all$CRT, cl_all[[dc]], use="complete.obs"))
    if (nrow(tr_all)>10) all_rs <- c(all_rs, cor(tr_all$CRT, tr_all[[dc]], use="complete.obs"))
    if (nrow(cl_cln)>10) clean_rs <- c(clean_rs, cor(cl_cln$CRT, cl_cln[[dc]], use="complete.obs"))
    if (nrow(tr_cln)>10) clean_rs <- c(clean_rs, cor(tr_cln$CRT, tr_cln[[dc]], use="complete.obs"))
  }

  # Response variability of flagged
  flag_resp_sd <- mean(apply(items_mat[flagged, ], 1, sd, na.rm=TRUE))
  clean_resp_sd <- mean(apply(items_mat[!flagged, ], 1, sd, na.rm=TRUE))

  data.frame(
    pct_target = target_pct,
    n_flagged = n_flag,
    pct_actual = mean(flagged),
    z_cutoff = z_cut,
    r_all = r_all, r_clean = r_clean, dr = r_clean - r_all,
    R2_all = summary(reg_all)$r.squared,
    R2_clean = summary(reg_clean)$r.squared,
    dR2 = summary(reg_clean)$r.squared - summary(reg_all)$r.squared,
    beta_crt_all = coef(reg_all)["CRT"],
    beta_crt_clean = coef(reg_clean)["CRT"],
    mean_disc_r_all = mean(all_rs), mean_disc_r_clean = mean(clean_rs),
    d_mean_disc_r = mean(clean_rs) - mean(all_rs),
    flag_resp_sd = flag_resp_sd, clean_resp_sd = clean_resp_sd
  )
}

cat("\n\n============================================================\n")
cat("RESULTS — CAPPED AT 5%, 10%, 15%, 20% FLAGGING\n")
cat("============================================================\n")

target_pcts <- c(0.05, 0.10, 0.15, 0.20)

for (si in 1:2) {
  d_full <- if (si==1) d1_c else d2_c
  items_mat <- if (si==1) items1_c else items2_c
  ms <- if (si==1) methods_s1 else methods_s2

  cat(sprintf("\n\n========== STUDY %d (N=%d) ==========\n", si, nrow(d_full)))

  # Paper baseline
  r_base <- cor(d_full$CRT, d_full$Discernment, use="complete.obs")
  reg_base <- lm(Discernment ~ Age + Sex + Education + Conserv + CRT, data=d_full)
  cat(sprintf("Baseline: r(CRT,Disc)=%.3f, R²=%.3f, beta(CRT)=%.3f\n",
      r_base, summary(reg_base)$r.squared, coef(reg_base)["CRT"]))

  cat(sprintf("\n%-10s %5s  %5s  z_cut  r_all  r_cln  Δr      R²_all R²_cln  ΔR²     β_all  β_cln   SD_flg SD_cln\n",
      "Method", "%tgt", "N_flg"))

  for (mn in names(ms)) {
    for (pct in target_pcts) {
      res <- eval_at_pct(d_full, ms[[mn]]$z_score, pct, items_mat)
      if (is.null(res)) next
      cat(sprintf("%-10s %4.0f%%  %5d  %5.2f  %.3f  %.3f  %+.3f   %.3f  %.3f  %+.3f   %.3f  %.3f   %.2f   %.2f\n",
          mn, pct*100, res$n_flagged, res$z_cutoff,
          res$r_all, res$r_clean, res$dr,
          res$R2_all, res$R2_clean, res$dR2,
          res$beta_crt_all, res$beta_crt_clean,
          res$flag_resp_sd, res$clean_resp_sd))
    }
    cat("\n")
  }
}

# ============================================================================
# Detailed: Table 1/3 at 20% flagging
# ============================================================================

cat("\n\n============================================================\n")
cat("TABLE 1/3 CORRELATIONS — 20% FLAGGING\n")
cat("============================================================\n")

for (si in 1:2) {
  d_full <- if (si==1) d1_c else d2_c
  ms <- if (si==1) methods_s1 else methods_s2
  disc_cols <- intersect(c("L_Discernment","C_Discernment","N_Discernment"), names(d_full))

  cat(sprintf("\n===== STUDY %d =====\n", si))

  # Baseline
  cat("\nBASELINE (all respondents):\n")
  for (grp in c("Clinton","Trump")) {
    g <- d_full[d_full$ClintonTrump == ifelse(grp=="Clinton",1,2), ]
    cat(sprintf("  %-10s (N=%d):", grp, nrow(g)))
    for (dc in disc_cols) cat(sprintf("  %s=%+.3f", gsub("_Discernment","",dc), cor(g$CRT, g[[dc]], use="complete.obs")))
    cat("\n")
  }

  for (mn in names(ms)) {
    z_cut <- quantile(ms[[mn]]$z_score, 0.20)
    clean <- d_full[ms[[mn]]$z_score > z_cut, ]
    n_flag <- sum(ms[[mn]]$z_score <= z_cut)

    cat(sprintf("\n%s (flagged %d, z<=%.2f):\n", toupper(mn), n_flag, z_cut))
    for (grp in c("Clinton","Trump")) {
      g <- clean[clean$ClintonTrump == ifelse(grp=="Clinton",1,2), ]
      cat(sprintf("  %-10s (N=%d):", grp, nrow(g)))
      for (dc in disc_cols) cat(sprintf("  %s=%+.3f", gsub("_Discernment","",dc), cor(g$CRT, g[[dc]], use="complete.obs")))
      cat("\n")
    }
  }
}

# ============================================================================
# Profile of the 20% flagged
# ============================================================================

cat("\n\n============================================================\n")
cat("PROFILE OF 20%% FLAGGED\n")
cat("============================================================\n")

for (si in 1:2) {
  d_full <- if (si==1) d1_c else d2_c
  items_mat <- if (si==1) items1_c else items2_c
  ms <- if (si==1) methods_s1 else methods_s2

  cat(sprintf("\n===== STUDY %d =====\n", si))
  cat(sprintf("%-10s  %6s %6s %6s %6s %6s %6s %6s\n",
      "Method", "CRT_f", "CRT_c", "Disc_f", "Disc_c", "SD_f", "SD_c", "Agree"))

  for (mn in names(ms)) {
    z_cut <- quantile(ms[[mn]]$z_score, 0.20)
    flagged <- ms[[mn]]$z_score <= z_cut
    fl <- d_full[flagged, ]; cl <- d_full[!flagged, ]

    sd_f <- mean(apply(items_mat[flagged, ], 1, sd, na.rm=TRUE))
    sd_c <- mean(apply(items_mat[!flagged, ], 1, sd, na.rm=TRUE))

    # Agreement between methods
    agree <- NA
    for (mn2 in names(ms)) {
      if (mn2 == mn) next
      z_cut2 <- quantile(ms[[mn2]]$z_score, 0.20)
      flagged2 <- ms[[mn2]]$z_score <= z_cut2
      agree <- mean(flagged == flagged2)
    }

    cat(sprintf("%-10s  %6.3f %6.3f %6.3f %6.3f %6.2f %6.2f  %.1f%%\n",
        mn, mean(fl$CRT,na.rm=T), mean(cl$CRT,na.rm=T),
        mean(fl$Discernment,na.rm=T), mean(cl$Discernment,na.rm=T),
        sd_f, sd_c, 100*agree))
  }
}

# ============================================================================
# Method agreement (overlap of flagged respondents)
# ============================================================================

cat("\n\n============================================================\n")
cat("METHOD AGREEMENT (20%% FLAGGING)\n")
cat("============================================================\n")

for (si in 1:2) {
  ms <- if (si==1) methods_s1 else methods_s2
  N <- if (si==1) nrow(d1_c) else nrow(d2_c)

  cat(sprintf("\nStudy %d (N=%d, 20%% = %d flagged each):\n", si, N, round(N*0.2)))

  flags <- list()
  for (mn in names(ms)) {
    z_cut <- quantile(ms[[mn]]$z_score, 0.20)
    flags[[mn]] <- ms[[mn]]$z_score <= z_cut
  }

  # Pairwise overlap (Jaccard)
  mns <- names(flags)
  cat(sprintf("  %-10s %-10s  overlap  Jaccard  flagged_by_both\n", "Method1", "Method2"))
  for (i in 1:(length(mns)-1)) {
    for (j in (i+1):length(mns)) {
      both <- sum(flags[[mns[i]]] & flags[[mns[j]]])
      either <- sum(flags[[mns[i]]] | flags[[mns[j]]])
      jacc <- both / either
      cat(sprintf("  %-10s %-10s  %5d    %.3f    %.1f%%\n",
          mns[i], mns[j], both, jacc, 100*both/N))
    }
  }

  # All three agree
  all3 <- sum(flags$coupled & flags$weighted & flags$efa_d)
  any3 <- sum(flags$coupled | flags$weighted | flags$efa_d)
  cat(sprintf("  All 3 agree: %d (%.1f%% of N)\n", all3, 100*all3/N))
  cat(sprintf("  Any 1 flags: %d (%.1f%% of N)\n", any3, 100*any3/N))
}

# ============================================================================
# Consensus flagging: flag only if >=2 methods agree
# ============================================================================

cat("\n\n============================================================\n")
cat("CONSENSUS FLAGGING (>=2 methods agree, 20%% each)\n")
cat("============================================================\n")

for (si in 1:2) {
  d_full <- if (si==1) d1_c else d2_c
  ms <- if (si==1) methods_s1 else methods_s2

  flags <- list()
  for (mn in names(ms)) {
    z_cut <- quantile(ms[[mn]]$z_score, 0.20)
    flags[[mn]] <- ms[[mn]]$z_score <= z_cut
  }

  consensus <- (as.integer(flags$coupled) + as.integer(flags$weighted) + as.integer(flags$efa_d)) >= 2
  clean <- d_full[!consensus, ]

  r_all <- cor(d_full$CRT, d_full$Discernment, use="complete.obs")
  r_cln <- cor(clean$CRT, clean$Discernment, use="complete.obs")

  reg_all <- lm(Discernment ~ Age + Sex + Education + Conserv + CRT, data=d_full)
  reg_cln <- lm(Discernment ~ Age + Sex + Education + Conserv + CRT, data=clean)

  cat(sprintf("\nStudy %d: consensus flagged %d/%d (%.1f%%)\n",
      si, sum(consensus), nrow(d_full), 100*mean(consensus)))
  cat(sprintf("  r(CRT,Disc): %.3f → %.3f (%+.3f)\n", r_all, r_cln, r_cln-r_all))
  cat(sprintf("  R²:          %.3f → %.3f (%+.3f)\n",
      summary(reg_all)$r.squared, summary(reg_cln)$r.squared,
      summary(reg_cln)$r.squared - summary(reg_all)$r.squared))
  cat(sprintf("  beta(CRT):   %.3f → %.3f\n", coef(reg_all)["CRT"], coef(reg_cln)["CRT"]))

  # Table correlations
  disc_cols <- intersect(c("L_Discernment","C_Discernment","N_Discernment"), names(d_full))
  for (grp in c("Clinton","Trump")) {
    g <- clean[clean$ClintonTrump == ifelse(grp=="Clinton",1,2), ]
    cat(sprintf("  %s (N=%d):", grp, nrow(g)))
    for (dc in disc_cols) cat(sprintf("  %s=%+.3f", gsub("_Discernment","",dc), cor(g$CRT, g[[dc]], use="complete.obs")))
    cat("\n")
  }
}

cat("\n\n=== DONE ===\n")
