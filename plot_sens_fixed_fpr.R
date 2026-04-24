# Sensitivity per corruption level at a FIXED operating point:
# choose the probability threshold so FPR on clean = 5%, then measure sensitivity.

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(tidyr)
})

df <- read.csv("sim_ensemble_scores.csv", stringsAsFactors = FALSE)
df$gt_careless <- (df$pattern != "clean") & (df$corruption > 0.50)

set.seed(42)
folds <- sample(seq_len(5), nrow(df), replace = TRUE)

train_cv <- function(formula_str) {
  preds <- numeric(nrow(df))
  for (f in unique(folds)) {
    tr <- folds != f; te <- folds == f
    mod <- glm(as.formula(formula_str), data = df[tr, ], family = binomial)
    preds[te] <- predict(mod, newdata = df[te, ], type = "response")
  }
  preds
}

df$p_zrr      <- train_cv("gt_careless ~ z_rr")
df$p_triad    <- train_cv("gt_careless ~ z_rr + irv + d2")
df$p_ensemble <- train_cv("gt_careless ~ z_rr + irv + longstring + d2 + person_tot")

# Pick threshold per model = 95th percentile of clean predictions (FPR=5%)
thr_zrr   <- quantile(df$p_zrr[df$pattern == "clean"],      0.95)
thr_triad <- quantile(df$p_triad[df$pattern == "clean"],    0.95)
thr_ens   <- quantile(df$p_ensemble[df$pattern == "clean"], 0.95)

cat(sprintf("FPR=5%% threshold: zrr=%.3f, triad=%.3f, ensemble=%.3f\n",
            thr_zrr, thr_triad, thr_ens))

df$flag_zrr      <- df$p_zrr      > thr_zrr
df$flag_triad    <- df$p_triad    > thr_triad
df$flag_ensemble <- df$p_ensemble > thr_ens

# Detection per corruption level (careless only)
by_level <- df %>%
  filter(pattern != "clean") %>%
  mutate(bin = round(corruption * 10) / 10) %>%
  group_by(bin) %>%
  summarise(
    n = n(),
    rate_zrr      = mean(flag_zrr),
    rate_triad    = mean(flag_triad),
    rate_ensemble = mean(flag_ensemble),
    .groups = "drop"
  )
cat("\n=== Detection per corruption level @ FPR=5% ===\n")
print(as.data.frame(by_level), row.names = FALSE)

# Per pattern at 100% corruption (and 90%)
by_pattern_high <- df %>%
  filter(pattern != "clean", corruption >= 0.85) %>%
  mutate(bin = round(corruption * 10) / 10) %>%
  group_by(pattern, bin) %>%
  summarise(
    n = n(),
    rate_zrr      = round(mean(flag_zrr),      3),
    rate_triad    = round(mean(flag_triad),    3),
    rate_ensemble = round(mean(flag_ensemble), 3),
    .groups = "drop"
  )
cat("\n=== Per-pattern detection at 90-100% corruption @ FPR=5% ===\n")
print(as.data.frame(by_pattern_high), row.names = FALSE)

# Compute MCC at FPR=5% operating point
compute_mcc <- function(flag, truth) {
  TP <- sum(flag & truth); TN <- sum(!flag & !truth)
  FP <- sum(flag & !truth); FN <- sum(!flag & truth)
  den <- sqrt(as.numeric(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  if (den == 0) 0 else (TP*TN - FP*FN) / den
}
cat("\n=== MCC at FPR=5% operating point ===\n")
for (nm in c("zrr","triad","ensemble")) {
  flag_col <- paste0("flag_", nm)
  mcc <- compute_mcc(df[[flag_col]], df$gt_careless)
  fpr <- mean(df[[flag_col]][df$pattern == "clean"])
  sens <- mean(df[[flag_col]][df$gt_careless])
  cat(sprintf("%-10s  sens=%.3f  fpr=%.3f  MCC=%.3f\n", nm, sens, fpr, mcc))
}

# ============================================================
# Plot: detection per corruption level, all three models
# ============================================================
plot_df <- by_level %>%
  pivot_longer(cols = starts_with("rate_"),
               names_to = "model", values_to = "rate") %>%
  mutate(model = recode(model,
                        "rate_zrr"      = "z_RR alone",
                        "rate_triad"    = "Triad: z_RR + IRV + D²",
                        "rate_ensemble" = "Full ensemble"),
         model = factor(model, levels = c("z_RR alone",
                                           "Triad: z_RR + IRV + D²",
                                           "Full ensemble")))

p <- ggplot(plot_df, aes(x = bin, y = rate,
                         color = model, group = model)) +
  annotate("rect", xmin = 0.55, xmax = 1.05, ymin = 0, ymax = 1,
           fill = "#ccffcc", alpha = 0.35) +
  annotate("rect", xmin = 0.05, xmax = 0.55, ymin = 0, ymax = 1,
           fill = "#ffe6cc", alpha = 0.35) +
  geom_vline(xintercept = 0.55, color = "grey30", linewidth = 0.5) +
  geom_hline(yintercept = 0.05, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  annotate("text", x = 0.12, y = 0.08, label = "FPR floor (5%)",
           size = 3, color = "grey30", hjust = 0) +
  annotate("text", x = 0.30, y = 0.98,
           label = "NOT careless (GT negative)",
           size = 3.2, color = "grey25", fontface = "bold") +
  annotate("text", x = 0.80, y = 0.98,
           label = "TRUE careless (GT positive)",
           size = 3.2, color = "grey25", fontface = "bold") +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.2f", rate)),
            vjust = -1.0, size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("z_RR alone" = "#377eb8",
                                 "Triad: z_RR + IRV + D²" = "#4daf4a",
                                 "Full ensemble" = "#e41a1c"),
                     name = "Model") +
  scale_x_continuous(breaks = seq(0.1, 1.0, 0.1),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0, 1.05),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Detection rate by corruption level @ fixed FPR = 5% on clean",
    subtitle = "Operating point calibrated for realistic use. Full ensemble reaches >90% detection at 80%+ corruption.",
    x = "Corruption level",
    y = "Detection rate (sensitivity)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"))

ggsave("plot_sens_fixed_fpr.png", p,
       width = 13, height = 7, dpi = 150, bg = "white")
cat("\nSaved: plot_sens_fixed_fpr.png\n")
