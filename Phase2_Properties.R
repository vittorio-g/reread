# Phase2_Properties.R
# Study the properties of the optimized index (RF + full ensemble + Scenario A
# at FPR=5%): MCC vs total_items, vs sample size, per-pattern, calibration, ROC.
#
# Uses existing sim_robust_scores.csv (24K respondents).

suppressPackageStartupMessages({
  library(MASS); library(dplyr); library(ggplot2); library(tidyr)
  library(randomForest); library(pROC)
})
select <- dplyr::select  # explicit shadow

df <- read.csv("sim_robust_scores.csv", stringsAsFactors = FALSE)
cat(sprintf("Loaded %d respondents\n", nrow(df)))

# Scenario A: clean vs corruption>80% only
mask_A <- (df$pattern == "clean") | (df$corruption > 0.80)
df_A <- df[mask_A, ]
df_A$gt <- df_A$pattern != "clean"
features <- c("z_rr", "z_rr_iter_t2", "irv", "longstring", "d2", "person_tot")

cat(sprintf("Scenario A: N = %d (positives=%d, negatives=%d)\n",
            nrow(df_A), sum(df_A$gt), sum(!df_A$gt)))

# 5-fold CV with RF
set.seed(42)
df_A$fold <- sample(seq_len(5), nrow(df_A), replace = TRUE)
df_A$prob <- NA_real_

for (f in 1:5) {
  tr <- df_A$fold != f; te <- df_A$fold == f
  X_tr <- as.matrix(df_A[tr, features, drop = FALSE])
  X_te <- as.matrix(df_A[te, features, drop = FALSE])
  y_tr <- factor(as.integer(df_A$gt[tr]), levels = c(0, 1))
  mod <- randomForest(X_tr, y_tr, ntree = 200,
                       mtry = max(1, floor(sqrt(ncol(X_tr)))))
  df_A$prob[te] <- predict(mod, X_te, type = "prob")[, "1"]
  cat(sprintf("Fold %d done\n", f))
}

# Calibrate threshold to FPR=5% on clean
thr <- quantile(df_A$prob[df_A$pattern == "clean"], 0.95)
df_A$flag <- df_A$prob > thr

cat(sprintf("\nThreshold @ FPR=5%%: %.3f\n", thr))

compute_mcc <- function(flag, gt) {
  TP <- sum(flag & gt); TN <- sum(!flag & !gt)
  FP <- sum(flag & !gt); FN <- sum(!flag & gt)
  den <- sqrt(as.numeric(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  list(mcc = if (den == 0) 0 else (TP*TN - FP*FN) / den,
       sens = TP/max(1,TP+FN), spec = TN/max(1,TN+FP),
       fpr = FP/max(1,FP+TN))
}

m_overall <- compute_mcc(df_A$flag, df_A$gt)
cat(sprintf("\nOverall: sens=%.3f spec=%.3f MCC=%.3f\n",
            m_overall$sens, m_overall$spec, m_overall$mcc))

# ============================================================
# Property 1 — MCC by total_items
# ============================================================
cat("\n=== Property 1: MCC by total_items ===\n")
prop_items <- df_A %>%
  group_by(total_items) %>%
  summarise(
    n = n(),
    n_pos = sum(gt),
    metrics = list(compute_mcc(flag, gt)),
    .groups = "drop"
  ) %>%
  mutate(
    sens = round(sapply(metrics, `[[`, "sens"), 3),
    spec = round(sapply(metrics, `[[`, "spec"), 3),
    mcc  = round(sapply(metrics, `[[`, "mcc"),  3)
  ) %>%
  select(total_items, n, n_pos, sens, spec, mcc)
print(as.data.frame(prop_items), row.names = FALSE)

# ============================================================
# Property 2 — MCC by nF (faceted)
# ============================================================
cat("\n=== Property 2: MCC by (nF, ipf) ===\n")
prop_nf <- df_A %>%
  group_by(nF, ipf, total_items) %>%
  summarise(
    n = n(),
    metrics = list(compute_mcc(flag, gt)),
    .groups = "drop"
  ) %>%
  mutate(
    sens = round(sapply(metrics, `[[`, "sens"), 3),
    mcc  = round(sapply(metrics, `[[`, "mcc"),  3)
  ) %>%
  select(nF, ipf, total_items, n, sens, mcc)
print(as.data.frame(prop_nf), row.names = FALSE)

# ============================================================
# Property 3 — Per-pattern detection at corruption ≥90%
# ============================================================
cat("\n=== Property 3: Per-pattern detection (corruption >80%) ===\n")
prop_pat <- df_A %>%
  filter(gt) %>%
  mutate(corr_bin = round(corruption * 10) / 10) %>%
  group_by(pattern, corr_bin) %>%
  summarise(n = n(), rate = round(mean(flag), 3), .groups = "drop")
print(as.data.frame(prop_pat), row.names = FALSE)

# ============================================================
# Property 4 — Stability across replications
# ============================================================
cat("\n=== Property 4: Stability across reps ===\n")
prop_rep <- df_A %>%
  group_by(nF, ipf, rep) %>%
  summarise(
    metrics = list(compute_mcc(flag, gt)),
    .groups = "drop"
  ) %>%
  mutate(mcc = sapply(metrics, `[[`, "mcc"))

stab <- prop_rep %>%
  group_by(nF, ipf) %>%
  summarise(
    mean_mcc = round(mean(mcc), 3),
    sd_mcc = round(sd(mcc), 3),
    se_mcc = round(sd(mcc)/sqrt(n()), 3),
    n_reps = n(),
    .groups = "drop"
  )
print(as.data.frame(stab), row.names = FALSE)

# ============================================================
# Property 5 — ROC curve
# ============================================================
roc_obj <- roc(df_A$gt, df_A$prob, quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))
cat(sprintf("\n=== AUC: %.4f ===\n", auc_val))

# ============================================================
# Save data + plots
# ============================================================
write.csv(df_A %>% select(-fold), "phase2_scored_data.csv", row.names = FALSE)
write.csv(prop_items, "phase2_mcc_by_items.csv", row.names = FALSE)
write.csv(prop_nf, "phase2_mcc_by_nF.csv", row.names = FALSE)
write.csv(prop_pat, "phase2_per_pattern.csv", row.names = FALSE)
write.csv(stab, "phase2_stability.csv", row.names = FALSE)

# Plot 1: MCC by total_items
p1 <- ggplot(prop_items, aes(x = total_items, y = mcc)) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "darkgreen",
             linewidth = 0.5) +
  annotate("text", x = 50, y = 0.72, label = "MCC = 0.7", size = 3,
           color = "darkgreen") +
  geom_line(linewidth = 1.2, color = "#2c5282") +
  geom_point(size = 4, color = "#2c5282") +
  geom_text(aes(label = sprintf("%.2f", mcc)), vjust = -1.0, size = 3.5) +
  scale_x_continuous(breaks = unique(prop_items$total_items)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(
    title = "MCC by questionnaire length (RF + full ensemble, Scenario A, FPR=5%)",
    subtitle = "Scenario A: clean vs corruption>80% only — middle (10-80%) excluded",
    x = "Total items (nF × ipf)",
    y = "MCC"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))
ggsave("plot_phase2_mcc_by_items.png", p1, width = 11, height = 6,
       dpi = 150, bg = "white")
cat("Saved: plot_phase2_mcc_by_items.png\n")

# Plot 2: per-pattern detection bar chart
prop_pat_filtered <- prop_pat %>% filter(corr_bin >= 0.9)
p2 <- ggplot(prop_pat_filtered, aes(x = pattern, y = rate, fill = factor(corr_bin))) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", rate)),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("0.9" = "#fee08b", "1" = "#d73027"),
                     name = "Corruption") +
  scale_y_continuous(limits = c(0, 1.1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Per-pattern detection rate at corruption >80%",
    subtitle = "RF + full ensemble + FPR=5% calibration",
    x = "Pattern", y = "Detection rate"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")
ggsave("plot_phase2_per_pattern.png", p2, width = 12, height = 6,
       dpi = 150, bg = "white")
cat("Saved: plot_phase2_per_pattern.png\n")

# Plot 3: ROC curve
roc_df <- data.frame(
  fpr = 1 - roc_obj$specificities,
  tpr = roc_obj$sensitivities
)
p3 <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
  geom_abline(slope = 1, linetype = "dashed", color = "grey60") +
  geom_path(linewidth = 1.2, color = "#e41a1c") +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "darkgreen") +
  annotate("text", x = 0.10, y = 0.4,
           label = sprintf("AUC = %.3f\nMCC at FPR=5%%: %.3f",
                            auc_val, m_overall$mcc),
           hjust = 0, size = 4, color = "#444444") +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  labs(title = "ROC curve — RF + full ensemble (Scenario A)",
       subtitle = "Trained on 80%, evaluated on held-out 20% (5-fold CV)",
       x = "False positive rate", y = "True positive rate") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))
ggsave("plot_phase2_roc.png", p3, width = 8, height = 8,
       dpi = 150, bg = "white")
cat("Saved: plot_phase2_roc.png\n")

# Plot 4: stability
p4 <- ggplot(prop_rep, aes(x = factor(rep), y = mcc)) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "darkgreen") +
  geom_boxplot(aes(group = paste(nF, ipf)), fill = "#fee08b", alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "#377eb8") +
  facet_wrap(~ paste0("nF=", nF, ", ipf=", ipf), ncol = 3) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(title = "Stability of MCC across 8 replications, by condition",
       subtitle = "Scenario A, FPR=5%, RF + full ensemble",
       x = "Replication", y = "MCC") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "#EEEEEE", color = NA),
        strip.text = element_text(face = "bold"))
ggsave("plot_phase2_stability.png", p4, width = 13, height = 7,
       dpi = 150, bg = "white")
cat("Saved: plot_phase2_stability.png\n")

cat("\n=== DONE Phase 2 ===\n")
