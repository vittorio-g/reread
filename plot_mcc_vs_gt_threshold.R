# Plot: MCC of the ensemble as a function of the GT cutoff
# (i.e., what corruption level we call "truly careless").

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(tidyr)
})

df <- read.csv("sim_ensemble_scores.csv", stringsAsFactors = FALSE)

# Define a range of GT thresholds
thresholds <- seq(0, 0.9, 0.1)

# 5-fold CV, each model refit per threshold (so the ensemble is optimized for
# each definition of "careless")
set.seed(42)
folds <- sample(seq_len(5), nrow(df), replace = TRUE)

compute_mcc <- function(flag, truth) {
  TP <- sum(flag & truth); TN <- sum(!flag & !truth)
  FP <- sum(flag & !truth); FN <- sum(!flag & truth)
  den <- sqrt(as.numeric(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  if (den == 0) 0 else (TP*TN - FP*FN) / den
}

results <- data.frame()

for (T in thresholds) {
  df$gt <- (df$pattern != "clean") & (df$corruption > T)
  n_pos <- sum(df$gt); n_neg <- sum(!df$gt)
  if (n_pos < 50 || n_neg < 50) next  # skip degenerate

  # Ensemble (all 5 detectors)
  preds_ens <- numeric(nrow(df))
  for (f in unique(folds)) {
    tr <- folds != f; te <- folds == f
    mod <- tryCatch(
      glm(gt ~ z_rr + irv + longstring + d2 + person_tot,
          data = df[tr, ], family = binomial),
      error = function(e) NULL
    )
    if (is.null(mod)) { preds_ens[te] <- 0.5; next }
    preds_ens[te] <- predict(mod, newdata = df[te, ], type = "response")
  }
  flag_ens <- preds_ens > 0.5

  # z_RR alone
  preds_zrr <- numeric(nrow(df))
  for (f in unique(folds)) {
    tr <- folds != f; te <- folds == f
    mod <- glm(gt ~ z_rr, data = df[tr, ], family = binomial)
    preds_zrr[te] <- predict(mod, newdata = df[te, ], type = "response")
  }
  flag_zrr <- preds_zrr > 0.5

  # Triad: z_rr + irv + d2
  preds_triad <- numeric(nrow(df))
  for (f in unique(folds)) {
    tr <- folds != f; te <- folds == f
    mod <- tryCatch(
      glm(gt ~ z_rr + irv + d2, data = df[tr, ], family = binomial),
      error = function(e) NULL
    )
    if (is.null(mod)) { preds_triad[te] <- 0.5; next }
    preds_triad[te] <- predict(mod, newdata = df[te, ], type = "response")
  }
  flag_triad <- preds_triad > 0.5

  results <- rbind(results, data.frame(
    gt_threshold = T,
    n_positive = n_pos,
    n_negative = n_neg,
    mcc_zrr      = compute_mcc(flag_zrr,   df$gt),
    mcc_triad    = compute_mcc(flag_triad, df$gt),
    mcc_ensemble = compute_mcc(flag_ens,   df$gt)
  ))

  cat(sprintf("T=%.2f: n_pos=%d  zrr=%.3f  triad=%.3f  ensemble=%.3f\n",
              T, n_pos, results$mcc_zrr[nrow(results)],
              results$mcc_triad[nrow(results)],
              results$mcc_ensemble[nrow(results)]))
}

write.csv(results, "mcc_vs_gt_threshold.csv", row.names = FALSE)

# Plot
plot_df <- results %>%
  pivot_longer(cols = c(mcc_zrr, mcc_triad, mcc_ensemble),
               names_to = "model", values_to = "mcc") %>%
  mutate(model = recode(model,
                        "mcc_zrr"      = "z_RR alone",
                        "mcc_triad"    = "Triad: z_RR + IRV + D²",
                        "mcc_ensemble" = "Full ensemble (5 detectors)"),
         model = factor(model, levels = c("z_RR alone",
                                           "Triad: z_RR + IRV + D²",
                                           "Full ensemble (5 detectors)")))

p <- ggplot(plot_df, aes(x = gt_threshold, y = mcc,
                         color = model, group = model)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
  annotate("text", x = 0.5, y = 0.05,
           label = "current default\n(GT > 50%)",
           size = 3, color = "grey30", hjust = -0.1) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3.5) +
  geom_text(aes(label = sprintf("%.3f", mcc)),
            vjust = -1.2, size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("z_RR alone" = "#377eb8",
                                 "Triad: z_RR + IRV + D²" = "#4daf4a",
                                 "Full ensemble (5 detectors)" = "#e41a1c"),
                     name = "Model") +
  scale_x_continuous(breaks = thresholds,
                     labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0, max(plot_df$mcc) + 0.08),
                     breaks = seq(0, 1, 0.1)) +
  labs(
    title = "MCC as a function of the 'careless' corruption threshold",
    subtitle = "X-axis: minimum corruption required to count as careless. Lower T = stricter task (more fuzzy respondents included).",
    x = "GT threshold (corruption > T → careless)",
    y = "MCC (5-fold CV)",
    caption = "12,000 respondents × 6 conditions × 4 reps, v3 injector. Logistic regression ensemble refit per threshold."
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"))

ggsave("plot_mcc_vs_gt_threshold.png", p,
       width = 12, height = 7, dpi = 150, bg = "white")
cat("\nSaved: plot_mcc_vs_gt_threshold.png\n")
cat("=== DONE ===\n")
