# Re-analysis: treat as "true careless" only respondents with corruption > 50%.
# Low-corruption respondents (10-50%) become part of the negative class.
# Recompute sensitivity, specificity, MCC under this labeling.

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(tidyr)
})
source("ReReReRe.R")

df <- read.csv("sim_detection_per_resp.csv", stringsAsFactors = FALSE)

# Compute auto-z threshold per cell
cond_thresholds <- df %>%
  distinct(nF, ipf, total_items) %>%
  rowwise() %>%
  mutate(auto_z = round(get_calibrated_z(total_items), 3)) %>%
  ungroup()

df <- df %>%
  left_join(cond_thresholds, by = c("nF","ipf","total_items")) %>%
  mutate(corruption_bin10 = ifelse(pattern == "clean", NA,
                                   round(corruption * 10) / 10),
         # NEW GT: corruption > 50% is "truly careless"; everything else (clean +
         # low-corruption) is non-careless.
         gt_careless = (pattern != "clean") & (corruption > 0.50),
         flagged_autoz    = z_standard <= auto_z,
         flagged_autoz_vp = z_varpen <= auto_z,
         flagged_z15      = z_standard <= 1.5,
         flagged_z15_vp   = z_varpen <= 1.5)

# ============================================================
# Aggregate metrics per variant
# ============================================================
metrics <- function(flag, truth) {
  TP <- sum( flag &  truth)
  TN <- sum(!flag & !truth)
  FP <- sum( flag & !truth)
  FN <- sum(!flag &  truth)
  den <- sqrt(as.numeric(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  mcc <- if (den == 0) 0 else (TP*TN - FP*FN) / den
  data.frame(
    sens = TP / max(1, TP+FN),
    spec = TN / max(1, TN+FP),
    fpr  = FP / max(1, FP+TN),
    mcc  = mcc,
    TP = TP, FP = FP, TN = TN, FN = FN
  )
}

m_z15    <- metrics(df$flagged_z15,      df$gt_careless) %>% mutate(variant = "z=1.5 fixed")
m_z15vp  <- metrics(df$flagged_z15_vp,   df$gt_careless) %>% mutate(variant = "z=1.5 + VP")
m_az     <- metrics(df$flagged_autoz,    df$gt_careless) %>% mutate(variant = "auto-z")
m_azvp   <- metrics(df$flagged_autoz_vp, df$gt_careless) %>% mutate(variant = "auto-z + VP")

summary_tbl <- bind_rows(m_z15, m_z15vp, m_az, m_azvp) %>%
  select(variant, sens, spec, fpr, mcc, TP, FP, TN, FN) %>%
  mutate(across(c(sens, spec, fpr, mcc), ~ round(.x, 3)))

cat("=== Overall metrics under GT = corruption > 50% ===\n")
print(as.data.frame(summary_tbl), row.names = FALSE)

# Per-pattern detection for "truly careless" (corruption > 50%)
pattern_metrics <- df %>%
  filter(gt_careless) %>%
  group_by(pattern) %>%
  summarise(
    n = n(),
    sens_z15    = round(mean(flagged_z15),      3),
    sens_z15_vp = round(mean(flagged_z15_vp),   3),
    sens_autoz  = round(mean(flagged_autoz),    3),
    sens_az_vp  = round(mean(flagged_autoz_vp), 3),
    .groups = "drop"
  )
cat("\n=== Per-pattern sensitivity (positive class only) ===\n")
print(as.data.frame(pattern_metrics), row.names = FALSE)

# Per-corruption detection curve (still plotted), with vertical line at 50%
det_curve <- df %>%
  filter(pattern != "clean") %>%
  group_by(pattern, corruption_bin10) %>%
  summarise(n = n(),
            det_z15   = mean(flagged_z15),
            det_autoz = mean(flagged_autoz),
            det_az_vp = mean(flagged_autoz_vp),
            .groups = "drop") %>%
  rename(corruption = corruption_bin10)

fpr_clean_autoz <- mean(df$flagged_autoz[df$pattern == "clean"])
# FPR on pooled negative class (clean + low corruption)
fpr_neg_autoz <- mean(df$flagged_autoz[!df$gt_careless])

pattern_colors <- c(
  random = "#e41a1c", longstring = "#377eb8", mixed = "#4daf4a",
  pure_straight = "#984ea3", acquiescent = "#ff7f00"
)

# ============================================================
# PLOT 1 — detection curve with GT boundary shaded
# ============================================================
p1 <- ggplot(det_curve, aes(x = corruption, y = det_autoz,
                             color = pattern, group = pattern)) +
  # Shade: left = false positive zone, right = true positive zone
  annotate("rect", xmin = 0.05, xmax = 0.55, ymin = 0, ymax = 1,
           fill = "#ffcccc", alpha = 0.3) +
  annotate("rect", xmin = 0.55, xmax = 1.05, ymin = 0, ymax = 1,
           fill = "#ccffcc", alpha = 0.3) +
  geom_vline(xintercept = 0.55, linetype = "solid",
             color = "grey30", linewidth = 0.6) +
  annotate("text", x = 0.30, y = 0.95,
           label = "NEGATIVE CLASS\n(low corruption = non-careless)",
           size = 3, color = "grey20", fontface = "bold") +
  annotate("text", x = 0.80, y = 0.95,
           label = "POSITIVE CLASS\n(corruption > 50% = true careless)",
           size = 3, color = "grey20", fontface = "bold") +
  geom_hline(yintercept = fpr_clean_autoz, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  annotate("text", x = 0.15, y = fpr_clean_autoz + 0.04,
           label = sprintf("clean-only FPR = %.1f%%",
                           fpr_clean_autoz * 100),
           size = 3, color = "grey30") +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3) +
  scale_color_manual(values = pattern_colors, name = "Pattern") +
  scale_x_continuous(breaks = seq(0.1, 1.0, 0.1),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Detection curve under GT = corruption > 50% (auto-z)",
    subtitle = sprintf("Left of line: flagging is a FP. Right: flagging is a TP. Overall MCC: %.3f | Sens: %.1f%% | Spec: %.1f%%",
                       m_az$mcc, m_az$sens * 100, m_az$spec * 100),
    x = "Corruption level",
    y = "Flagging rate at auto-z threshold",
    caption = "Red zone: respondents here are non-careless by definition, so flagging = false positive."
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

ggsave("plot_gt50_detection_curve.png", p1,
       width = 12, height = 6.5, dpi = 150, bg = "white")
cat("\nSaved: plot_gt50_detection_curve.png\n")

# ============================================================
# PLOT 2 — bar chart summary of metrics across variants
# ============================================================
bar_df <- summary_tbl %>%
  select(variant, sens, spec, mcc) %>%
  pivot_longer(-variant, names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         "sens" = "Sensitivity",
                         "spec" = "Specificity",
                         "mcc"  = "MCC"),
         metric = factor(metric, levels = c("Sensitivity", "Specificity", "MCC")),
         variant = factor(variant,
                          levels = c("z=1.5 fixed", "z=1.5 + VP",
                                     "auto-z", "auto-z + VP")))

p2 <- ggplot(bar_df, aes(x = variant, y = value, fill = metric)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", value)),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = c("Sensitivity" = "#377eb8",
                               "Specificity" = "#4daf4a",
                               "MCC" = "#e41a1c"), name = "") +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Overall metrics under GT = corruption > 50%",
    subtitle = "Comparison of flagging strategies (N=12,000 respondents; 6 conditions x 4 reps)",
    x = NULL,
    y = "Value"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"))

ggsave("plot_gt50_metrics_bars.png", p2,
       width = 10, height = 6, dpi = 150, bg = "white")
cat("Saved: plot_gt50_metrics_bars.png\n")

# ============================================================
# PLOT 3 — per-pattern sensitivity comparison
# ============================================================
pp_long <- pattern_metrics %>%
  select(pattern, sens_z15, sens_z15_vp, sens_autoz, sens_az_vp) %>%
  pivot_longer(-pattern, names_to = "variant", values_to = "sens") %>%
  mutate(variant = recode(variant,
                          "sens_z15"    = "z=1.5 fixed",
                          "sens_z15_vp" = "z=1.5 + VP",
                          "sens_autoz"  = "auto-z",
                          "sens_az_vp"  = "auto-z + VP"),
         variant = factor(variant,
                          levels = c("z=1.5 fixed", "z=1.5 + VP",
                                     "auto-z", "auto-z + VP")))

p3 <- ggplot(pp_long, aes(x = pattern, y = sens, fill = variant)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", sens)),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 2.8) +
  scale_fill_brewer(palette = "Set2", name = "Strategy") +
  scale_y_continuous(limits = c(0, 1.05),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Sensitivity per careless pattern (GT = corruption > 50%)",
    subtitle = "Restricted to true careless respondents only",
    x = "Pattern",
    y = "Sensitivity (true-positive rate)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0))

ggsave("plot_gt50_sens_by_pattern.png", p3,
       width = 12, height = 6, dpi = 150, bg = "white")
cat("Saved: plot_gt50_sens_by_pattern.png\n")

cat("=== DONE ===\n")
