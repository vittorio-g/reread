# Sensitivity (detection rate) per corruption level using the ensemble
# trained once with the default GT > 50% cutoff.

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(tidyr)
})

df <- read.csv("sim_ensemble_scores.csv", stringsAsFactors = FALSE)
df$gt_careless <- (df$pattern != "clean") & (df$corruption > 0.50)

set.seed(42)
folds <- sample(seq_len(5), nrow(df), replace = TRUE)

# Train three models on the default GT=50% labeling
train_cv <- function(formula_str) {
  preds <- numeric(nrow(df))
  for (f in unique(folds)) {
    tr <- folds != f; te <- folds == f
    mod <- glm(as.formula(formula_str), data = df[tr, ], family = binomial)
    preds[te] <- predict(mod, newdata = df[te, ], type = "response")
  }
  preds
}

df$pred_zrr      <- train_cv("gt_careless ~ z_rr")
df$pred_triad    <- train_cv("gt_careless ~ z_rr + irv + d2")
df$pred_ensemble <- train_cv("gt_careless ~ z_rr + irv + longstring + d2 + person_tot")

# Bin corruption into decimi, compute detection rate per bucket for each model
df <- df %>%
  mutate(corruption_bin = ifelse(pattern == "clean", "clean",
                                  sprintf("%d%%", round(corruption * 10) * 10)))

detect_by_bucket <- df %>%
  group_by(pattern, corruption_bin) %>%
  summarise(
    n = n(),
    rate_zrr      = mean(pred_zrr      > 0.5),
    rate_triad    = mean(pred_triad    > 0.5),
    rate_ensemble = mean(pred_ensemble > 0.5),
    .groups = "drop"
  )

cat("=== Detection rate per (pattern, corruption bucket) ===\n")
print(as.data.frame(detect_by_bucket), row.names = FALSE)

# Overall by corruption level (ignoring pattern), for careless-only respondents
by_level <- df %>%
  filter(pattern != "clean") %>%
  mutate(corruption_bin10 = round(corruption * 10) / 10) %>%
  group_by(corruption_bin10) %>%
  summarise(
    n = n(),
    rate_zrr      = mean(pred_zrr > 0.5),
    rate_triad    = mean(pred_triad > 0.5),
    rate_ensemble = mean(pred_ensemble > 0.5),
    .groups = "drop"
  )

fpr_clean <- df %>% filter(pattern == "clean") %>%
  summarise(
    fpr_zrr      = mean(pred_zrr > 0.5),
    fpr_triad    = mean(pred_triad > 0.5),
    fpr_ensemble = mean(pred_ensemble > 0.5)
  )

cat("\n=== Detection rate by overall corruption level ===\n")
print(as.data.frame(by_level), row.names = FALSE)

cat("\n=== FPR on clean ===\n")
print(fpr_clean)

# Plot
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

fpr_df <- data.frame(
  model = c("z_RR alone", "Triad: z_RR + IRV + D²", "Full ensemble"),
  fpr   = c(fpr_clean$fpr_zrr, fpr_clean$fpr_triad, fpr_clean$fpr_ensemble)
)

p <- ggplot(plot_df, aes(x = corruption_bin10, y = rate,
                         color = model, group = model)) +
  annotate("rect", xmin = 0.55, xmax = 1.05, ymin = 0, ymax = 1,
           fill = "#ccffcc", alpha = 0.35) +
  annotate("rect", xmin = 0.05, xmax = 0.55, ymin = 0, ymax = 1,
           fill = "#ffe6cc", alpha = 0.35) +
  geom_vline(xintercept = 0.55, color = "grey30", linewidth = 0.5) +
  geom_hline(data = fpr_df, aes(yintercept = fpr, color = model),
             linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
  annotate("text", x = 0.30, y = 0.97,
           label = "NOT careless (GT = non-positive)",
           size = 3.2, color = "grey25", fontface = "bold") +
  annotate("text", x = 0.80, y = 0.97,
           label = "TRUE careless (GT = positive)",
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
    title = "Detection rate by corruption level — ensemble is near-perfect on full careless",
    subtitle = "Models trained with GT > 50%. Dashed lines = false-positive rate on clean respondents.",
    x = "Corruption level",
    y = "Detection rate (p > 0.5)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"))

ggsave("plot_sens_by_corruption.png", p,
       width = 13, height = 7, dpi = 150, bg = "white")
cat("\nSaved: plot_sens_by_corruption.png\n")
