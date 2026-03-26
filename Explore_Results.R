#### Explore Multiverse Results ####
# === REVISION 2026-03-10f ===
# Switched to z-score metric with fixed z_threshold.
# Columns: z_threshold, mean_z_good, mean_z_careless

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ============================================================
# LOAD DATA
# ============================================================

res <- read.csv("multiverse_results.csv")
res_type <- read.csv("multiverse_results_by_type.csv")

cat("Main results:", nrow(res), "rows\n")
cat("By-type results:", nrow(res_type), "rows\n\n")

# Items-per-factor pattern (fixed, cycling 10-6-3)
cat("Items per factor pattern: cycling [10, 6, 3]\n")
cat("Total items by nFactors:\n")
res %>% distinct(nFactors, total_items) %>% arrange(nFactors) %>% print()
cat("\n")

# Quick summary
cat("=== OVERALL PERFORMANCE ===\n")
cat("MCC:         ", round(mean(res$mcc, na.rm=T), 3),
    " (SD=", round(sd(res$mcc, na.rm=T), 3), ")\n")
cat("F1:          ", round(mean(res$f1, na.rm=T), 3),
    " (SD=", round(sd(res$f1, na.rm=T), 3), ")\n")
cat("Sensitivity: ", round(mean(res$sensitivity, na.rm=T), 3),
    " (SD=", round(sd(res$sensitivity, na.rm=T), 3), ")\n")
cat("Specificity: ", round(mean(res$specificity, na.rm=T), 3),
    " (SD=", round(sd(res$specificity, na.rm=T), 3), ")\n")
cat("Precision:   ", round(mean(res$precision, na.rm=T), 3),
    " (SD=", round(sd(res$precision, na.rm=T), 3), ")\n")
cat("Youden's J:  ", round(mean(res$youden, na.rm=T), 3),
    " (SD=", round(sd(res$youden, na.rm=T), 3), ")\n")
cat("Mean z good:     ", round(mean(res$mean_z_good, na.rm=T), 3), "\n")
cat("Mean z careless: ", round(mean(res$mean_z_careless, na.rm=T), 3), "\n\n")

# ============================================================
# 1. BEST PARAMETERS (by MCC — balanced across all 4 quadrants)
# ============================================================

cat("=== TOP 10 PARAMETER COMBINATIONS (MCC) ===\n")
res %>%
  arrange(desc(mcc)) %>%
  select(nFactors, total_items, n_respondents, pct_careless,
         corProp, z_threshold, sensitivity, specificity, precision, f1, mcc,
         mean_z_good, mean_z_careless) %>%
  head(10) %>%
  print()

cat("\n=== TOP 10 PARAMETER COMBINATIONS (F1) ===\n")
res %>%
  arrange(desc(f1)) %>%
  select(nFactors, total_items, n_respondents, pct_careless,
         corProp, z_threshold, sensitivity, specificity, precision, f1, mcc) %>%
  head(10) %>%
  print()

# ============================================================
# 2. HEATMAP: MCC by corProp x z_threshold
# ============================================================

p1 <- res %>%
  group_by(corProp, z_threshold) %>%
  summarise(mean_mcc = mean(mcc, na.rm=T), .groups="drop") %>%
  ggplot(aes(x = factor(corProp), y = factor(z_threshold), fill = mean_mcc)) +
  geom_tile() +
  geom_text(aes(label = round(mean_mcc, 2)), size = 2.5) +
  scale_fill_viridis_c(option = "plasma") +
  labs(title = "Mean MCC by corProp x z_threshold",
       x = "corProp", y = "z_threshold", fill = "MCC") +
  theme_minimal()

# ============================================================
# 3. SENSITIVITY, SPECIFICITY, PRECISION by z_threshold
# ============================================================

p2 <- res %>%
  group_by(z_threshold) %>%
  summarise(
    Sensitivity = mean(sensitivity, na.rm=T),
    Specificity = mean(specificity, na.rm=T),
    Precision = mean(precision, na.rm=T),
    F1 = mean(f1, na.rm=T),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(Sensitivity, Specificity, Precision, F1),
               names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = z_threshold, y = value, color = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Sensitivity, Specificity, Precision & F1 by z_threshold",
       subtitle = "z_threshold: flag if z-score <= threshold (higher threshold = more flagged)",
       x = "z_threshold", y = "Mean value", color = "") +
  theme_minimal()

# ============================================================
# 4. EFFECT OF SAMPLE SIZE
# ============================================================

p3 <- res %>%
  group_by(n_respondents, z_threshold) %>%
  summarise(mean_mcc = mean(mcc, na.rm=T), .groups="drop") %>%
  ggplot(aes(x = z_threshold, y = mean_mcc,
             color = factor(n_respondents), group = n_respondents)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "MCC by Sample Size and z_threshold",
       x = "z_threshold", y = "Mean MCC", color = "N") +
  theme_minimal()

# ============================================================
# 5. EFFECT OF QUESTIONNAIRE COMPLEXITY (nFactors)
# ============================================================

p4 <- res %>%
  group_by(nFactors, total_items, corProp) %>%
  summarise(mean_mcc = mean(mcc, na.rm=T), .groups="drop") %>%
  ggplot(aes(x = factor(nFactors), y = mean_mcc,
             color = factor(corProp), group = corProp)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set2") +
  labs(title = "MCC by nFactors and corProp",
       subtitle = "Items per factor cycle: 10, 6, 3",
       x = "nFactors", y = "Mean MCC", color = "corProp") +
  theme_minimal()

# ============================================================
# 6. EFFECT OF % CARELESS on metrics
# ============================================================

p5 <- res %>%
  group_by(pct_careless) %>%
  summarise(
    Sensitivity = mean(sensitivity, na.rm=T),
    Specificity = mean(specificity, na.rm=T),
    Precision = mean(precision, na.rm=T),
    F1 = mean(f1, na.rm=T),
    MCC = mean(mcc, na.rm=T),
    .groups="drop"
  ) %>%
  pivot_longer(cols = c(Sensitivity, Specificity, Precision, F1, MCC),
               names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = pct_careless, y = value, color = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Metrics by % Careless Respondents",
       subtitle = "Corruption levels: 50-100% | Z-score + fixed threshold",
       x = "% Careless respondents (injected)", y = "Value", color = "") +
  theme_minimal()

# ============================================================
# 7. DETECTION BY CARELESS TYPE
# ============================================================

p6 <- res_type %>%
  group_by(pattern, z_threshold) %>%
  summarise(mean_det = mean(detection_rate, na.rm=T), .groups="drop") %>%
  ggplot(aes(x = z_threshold, y = mean_det, color = pattern)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "Detection Rate by Careless Pattern Type",
       subtitle = "Corruption levels: 50-100% | Z-score + fixed threshold",
       x = "z_threshold", y = "Mean Detection Rate", color = "Pattern") +
  theme_minimal()

# ============================================================
# 7b. OPTIMAL z_threshold — where does MCC peak per condition?
# ============================================================

cat("\n=== OPTIMAL z_threshold BY nFactors ===\n")
res %>%
  group_by(nFactors, z_threshold) %>%
  summarise(mean_mcc = mean(mcc, na.rm=T), .groups="drop") %>%
  group_by(nFactors) %>%
  slice_max(mean_mcc, n=1) %>%
  print()

p6b <- res %>%
  group_by(nFactors, total_items, z_threshold) %>%
  summarise(mean_mcc = mean(mcc, na.rm=T), .groups="drop") %>%
  ggplot(aes(x = z_threshold, y = mean_mcc,
             color = factor(nFactors), group = nFactors)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_color_brewer(palette = "Set2") +
  labs(title = "MCC by z_threshold for each nFactors level",
       subtitle = "Identifies optimal z_threshold region per questionnaire size",
       x = "z_threshold", y = "Mean MCC", color = "nFactors") +
  theme_minimal()

# ============================================================
# 8. PRECISION vs SENSITIVITY TRADEOFF (by z_threshold)
# ============================================================

p7 <- res %>%
  group_by(z_threshold) %>%
  summarise(
    mean_prec = mean(precision, na.rm=T),
    mean_sens = mean(sensitivity, na.rm=T),
    .groups="drop"
  ) %>%
  ggplot(aes(x = mean_sens, y = mean_prec, color = z_threshold)) +
  geom_point(size = 3) +
  geom_path(linewidth = 0.8) +
  geom_text(aes(label = z_threshold), vjust = -0.8, size = 2.5) +
  scale_color_viridis_c() +
  labs(title = "Precision vs Sensitivity Tradeoff",
       subtitle = "Each point = one z_threshold value (averaged over all conditions)",
       x = "Mean Sensitivity", y = "Mean Precision", color = "z_threshold") +
  theme_minimal()

# ============================================================
# 9. Z-SCORE SEPARATION: mean_z_good vs mean_z_careless
# ============================================================

cat("\n=== Z-SCORE SEPARATION BY nFactors ===\n")
res %>%
  group_by(nFactors, total_items) %>%
  summarise(
    mean_z_good = mean(mean_z_good, na.rm=T),
    mean_z_careless = mean(mean_z_careless, na.rm=T),
    z_gap = mean(mean_z_good, na.rm=T) - mean(mean_z_careless, na.rm=T),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  print()

p8 <- res %>%
  group_by(nFactors, total_items, corProp) %>%
  summarise(
    z_good = mean(mean_z_good, na.rm=T),
    z_careless = mean(mean_z_careless, na.rm=T),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(z_good, z_careless),
               names_to = "group", values_to = "mean_z") %>%
  mutate(group = ifelse(group == "z_good", "Good", "Careless")) %>%
  ggplot(aes(x = factor(nFactors), y = mean_z,
             color = group, shape = factor(corProp))) +
  geom_point(size = 2.5, position = position_dodge(0.3)) +
  scale_color_manual(values = c("Good" = "steelblue", "Careless" = "firebrick")) +
  labs(title = "Mean Z-Score: Good vs Careless by nFactors",
       subtitle = "Larger gap = better discrimination",
       x = "nFactors", y = "Mean Z-Score", color = "", shape = "corProp") +
  theme_minimal()

# ============================================================
# 10. n < p FLAG: does it matter?
# ============================================================

cat("\n=== n < p COMPARISON ===\n")
res %>%
  group_by(n_lt_p) %>%
  summarise(
    mean_mcc = mean(mcc, na.rm=T),
    mean_f1 = mean(f1, na.rm=T),
    mean_sens = mean(sensitivity, na.rm=T),
    mean_spec = mean(specificity, na.rm=T),
    mean_prec = mean(precision, na.rm=T),
    mean_z_gap = mean(mean_z_good - mean_z_careless, na.rm=T),
    n_cells = n(),
    .groups="drop"
  ) %>%
  print()

# ============================================================
# DISPLAY PLOTS
# ============================================================

print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
print(p6b)
print(p8)

# Combined overview
combined <- (p1 | p2) / (p3 | p4) / (p5 | p6) / (p7 | p6b) / (p8 | plot_spacer()) +
  plot_annotation(title = "ReReReRe Multiverse Results Overview",
                  subtitle = "Z-Score + Fixed Threshold | Primary: MCC | Loadings U(0.4,0.8) | Corruption 50-100%",
                  theme = theme(plot.title = element_text(face = "bold", size = 16)))

ggsave("multiverse_overview.png", combined, width = 16, height = 28, dpi = 200)
cat("\nSaved: multiverse_overview.png\n")
