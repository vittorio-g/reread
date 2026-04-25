# Phase2b — Feature importance from Random Forest
# Also produces a calibration plot (predicted prob vs actual rate of careless)

suppressPackageStartupMessages({
  library(MASS); library(dplyr); library(ggplot2); library(tidyr)
  library(randomForest)
})
select <- dplyr::select

df <- read.csv("sim_robust_scores.csv", stringsAsFactors = FALSE)
mask_A <- (df$pattern == "clean") | (df$corruption > 0.80)
df_A <- df[mask_A, ]
df_A$gt <- df_A$pattern != "clean"
features <- c("z_rr", "z_rr_iter_t2", "irv", "longstring", "d2", "person_tot")

# Train one RF on full data for feature importance
set.seed(123)
X <- as.matrix(df_A[, features, drop = FALSE])
y <- factor(as.integer(df_A$gt), levels = c(0, 1))
rf_full <- randomForest(X, y, ntree = 500, importance = TRUE,
                         mtry = max(1, floor(sqrt(ncol(X)))))

# Importance metrics
imp <- importance(rf_full)
imp_df <- data.frame(
  feature = rownames(imp),
  mean_decrease_accuracy = round(imp[, "MeanDecreaseAccuracy"], 3),
  mean_decrease_gini     = round(imp[, "MeanDecreaseGini"],     1)
) %>% arrange(desc(mean_decrease_accuracy))

cat("=== Feature Importance (Random Forest) ===\n")
print(imp_df, row.names = FALSE)
write.csv(imp_df, "phase2b_feature_importance.csv", row.names = FALSE)

# Plot
imp_df$feature <- factor(imp_df$feature,
                          levels = imp_df$feature[order(imp_df$mean_decrease_accuracy)])
p <- ggplot(imp_df, aes(x = feature, y = mean_decrease_accuracy)) +
  geom_col(fill = "#2c5282", width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", mean_decrease_accuracy)),
            hjust = -0.1, size = 3.5) +
  coord_flip() +
  labs(title = "Random Forest feature importance",
       subtitle = "Mean decrease in accuracy when each feature is permuted (Scenario A)",
       x = NULL, y = "Mean decrease in accuracy") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))
ggsave("plot_phase2b_importance.png", p,
       width = 11, height = 6, dpi = 150, bg = "white")
cat("\nSaved: plot_phase2b_importance.png\n")

# Calibration plot
df_A$prob <- predict(rf_full, X, type = "prob")[, "1"]
calib <- df_A %>%
  mutate(prob_bin = cut(prob, breaks = seq(0, 1, 0.05), include.lowest = TRUE)) %>%
  group_by(prob_bin) %>%
  summarise(mean_prob = mean(prob), actual_rate = mean(as.integer(gt)),
            n = n(), .groups = "drop")

p2 <- ggplot(calib, aes(x = mean_prob, y = actual_rate)) +
  geom_abline(slope = 1, linetype = "dashed", color = "grey50") +
  geom_point(aes(size = n), color = "#e41a1c", alpha = 0.7) +
  geom_line(color = "#e41a1c") +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  scale_size_continuous(range = c(2, 8), guide = "none") +
  labs(title = "Calibration plot — RF predicted probability vs actual rate",
       subtitle = "Each point = bin of predicted probability. Diagonal = perfectly calibrated.",
       x = "Mean predicted probability", y = "Actual rate of careless") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))
ggsave("plot_phase2b_calibration.png", p2,
       width = 9, height = 8, dpi = 150, bg = "white")
cat("Saved: plot_phase2b_calibration.png\n")

cat("=== DONE Phase 2b ===\n")
