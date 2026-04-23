# Regenerate ensemble plots from saved CSVs (fixing the label bug)

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })

forward <- read.csv("sim_ensemble_forward.csv", stringsAsFactors = FALSE)
standalone <- read.csv("sim_ensemble_standalone.csv", stringsAsFactors = FALSE)

cat("=== FORWARD SELECTION ===\n")
print(forward, row.names = FALSE)
cat("\n=== STANDALONE ===\n")
print(standalone, row.names = FALSE)

# Per-pattern analysis on the full scores
df <- read.csv("sim_ensemble_scores.csv", stringsAsFactors = FALSE)
df$gt_careless <- (df$pattern != "clean") & (df$corruption > 0.50)

cat("\n=== PATTERN DISTRIBUTION ===\n")
print(df %>% count(pattern))

cat("\n=== GT DISTRIBUTION ===\n")
print(table(df$gt_careless))

cat("\n=== CARELESS PCT DISTRIBUTION PER PATTERN (mean, range) ===\n")
print(df %>% filter(pattern != "clean") %>%
  group_by(pattern) %>%
  summarise(n = n(),
            mean_pct = round(mean(corruption), 3),
            min_pct = min(corruption),
            max_pct = max(corruption),
            n_gt = sum(corruption > 0.50)))

# Per-pattern sensitivity of ensemble (refit full model)
set.seed(42)
folds <- sample(seq_len(5), nrow(df), replace = TRUE)

# Train ensemble full model on all features
feats_full <- c("z_rr", "irv", "longstring", "d2", "person_tot")
preds_full <- numeric(nrow(df))
for (f in unique(folds)) {
  tr <- folds != f; te <- folds == f
  mod <- glm(gt_careless ~ z_rr + irv + longstring + d2 + person_tot,
             data = df[tr, ], family = binomial)
  preds_full[te] <- predict(mod, newdata = df[te, ], type = "response")
}
df$pred_ensemble <- preds_full

# Standalone z_rr predictions
preds_zrr <- numeric(nrow(df))
for (f in unique(folds)) {
  tr <- folds != f; te <- folds == f
  mod <- glm(gt_careless ~ z_rr, data = df[tr, ], family = binomial)
  preds_zrr[te] <- predict(mod, newdata = df[te, ], type = "response")
}
df$pred_zrr <- preds_zrr

# Per-pattern sensitivity on gt_careless TRUE
pattern_sens <- df %>%
  filter(gt_careless) %>%
  group_by(pattern) %>%
  summarise(
    n = n(),
    sens_zrr = round(mean(pred_zrr > 0.5), 3),
    sens_ensemble = round(mean(pred_ensemble > 0.5), 3),
    gain = round(sens_ensemble - sens_zrr, 3)
  )
cat("\n=== PER-PATTERN SENSITIVITY (z_rr alone vs ensemble) ===\n")
print(as.data.frame(pattern_sens), row.names = FALSE)

# =====================
# PLOT 1 — forward selection
# =====================
forward$label_short <- sapply(strsplit(forward$added, " "), function(x) x[1])

p1 <- ggplot(forward, aes(x = step, y = mcc)) +
  geom_line(linewidth = 1.2, color = "#2c5282") +
  geom_point(size = 5, color = "#2c5282") +
  geom_text(aes(label = sprintf("%.3f", mcc)), vjust = -1.2, size = 4) +
  geom_text(aes(label = features), vjust = 2.2, size = 3, color = "grey30") +
  scale_x_continuous(breaks = forward$step) +
  scale_y_continuous(limits = c(min(forward$mcc) - 0.08,
                                 max(forward$mcc) + 0.08)) +
  labs(
    title = "Greedy forward selection — ensemble MCC per detector added",
    subtitle = "5-fold CV, logistic regression, GT = corruption > 50%",
    x = "Step",
    y = "MCC (cross-validated)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))
ggsave("plot_ensemble_forward.png", p1,
       width = 11, height = 6, dpi = 150, bg = "white")
cat("\nSaved: plot_ensemble_forward.png\n")

# =====================
# PLOT 2 — standalone bar chart with deltas
# =====================
standalone_ordered <- standalone %>% arrange(desc(mcc))
p2 <- ggplot(standalone_ordered,
              aes(x = reorder(detector, -mcc), y = mcc, fill = detector)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.3f", mcc)), vjust = -0.4, size = 4) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  scale_y_continuous(limits = c(0, max(standalone$mcc) + 0.1)) +
  labs(
    title = "Standalone detector strength (each alone, 5-fold CV MCC)",
    subtitle = "IRV alone beats z_RR under the v3 injector (equalized pattern semantics)",
    x = NULL, y = "MCC (CV)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))
ggsave("plot_ensemble_standalone.png", p2,
       width = 10, height = 6, dpi = 150, bg = "white")
cat("Saved: plot_ensemble_standalone.png\n")

# =====================
# PLOT 3 — per-pattern sensitivity comparison
# =====================
sens_long <- pattern_sens %>%
  pivot_longer(c(sens_zrr, sens_ensemble),
               names_to = "model", values_to = "sens") %>%
  mutate(model = recode(model, "sens_zrr" = "z_RR alone",
                        "sens_ensemble" = "Full ensemble"))

p3 <- ggplot(sens_long, aes(x = pattern, y = sens, fill = model)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", sens)),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = c("z_RR alone" = "#377eb8",
                                "Full ensemble" = "#e41a1c")) +
  scale_y_continuous(limits = c(0, 1.05),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Per-pattern sensitivity: z_RR alone vs full ensemble",
    subtitle = "Restricted to truly careless respondents (corruption > 50%)",
    x = "Pattern", y = "Sensitivity (at p > 0.5 threshold)",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")
ggsave("plot_ensemble_per_pattern.png", p3,
       width = 12, height = 6, dpi = 150, bg = "white")
cat("Saved: plot_ensemble_per_pattern.png\n")

cat("\n=== DONE ===\n")
