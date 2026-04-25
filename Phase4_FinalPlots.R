# Phase4 — Final summary plots integrating Phase 2 + Phase 3 results

suppressPackageStartupMessages({
  library(MASS); library(dplyr); library(ggplot2); library(tidyr)
})
select <- dplyr::select

# Load Phase 3 results
res_rate <- read.csv("phase3_results_by_rate.csv")
df_p3 <- read.csv("phase3_all_scored.csv")
df_p2 <- read.csv("phase2_scored_data.csv")

# Plot A — MCC by careless rate (validation)
pA <- ggplot(res_rate, aes(x = factor(careless_rate*100, levels=c("20","40","60")),
                            y = mcc)) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "darkgreen") +
  geom_col(fill = "#377eb8", width = 0.55) +
  geom_text(aes(label = sprintf("%.3f", mcc)), vjust = -0.5, size = 5,
            fontface = "bold") +
  geom_text(aes(label = sprintf("sens=%.3f", sens), y = 0.05),
            color = "white", size = 3.5) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(title = "Validation MCC by careless rate (fresh data, RF + full ensemble)",
       subtitle = "Scenario A: clean vs corruption>80% only · FPR=5% calibration",
       x = "Careless rate (%)", y = "MCC") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))
ggsave("plot_phase4_validation_mcc.png", pA,
       width = 10, height = 6, dpi = 150, bg = "white")

# Plot B — Per-pattern detection at >80% corruption, per rate
per_pat_rate <- df_p3 %>%
  filter(gt) %>%
  group_by(rate, pattern) %>%
  summarise(rate_det = mean(flag), n = n(), .groups = "drop")

pB <- ggplot(per_pat_rate, aes(x = pattern, y = rate_det,
                               fill = factor(rate*100))) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", rate_det)),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("20" = "#fee08b", "40" = "#fdae61", "60" = "#d73027"),
                     name = "Careless rate (%)") +
  scale_y_continuous(limits = c(0, 1.1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Per-pattern detection @ corruption >80%, by careless rate",
       subtitle = "RF + full ensemble validated on 54,000 fresh respondents",
       x = "Pattern", y = "Detection rate") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")
ggsave("plot_phase4_validation_patterns.png", pB,
       width = 13, height = 6, dpi = 150, bg = "white")

# Plot C — Phase2 + Phase3 combined: MCC by total_items (with shaded confidence)
# Phase 2: training = sim_robust (24K resp), Phase 3: validation = fresh (54K resp)
# Phase 2 MCC by total_items already saved, recompute Phase 3 by total_items
phase3_by_items <- df_p3 %>%
  group_by(total_items, rate) %>%
  do({
    tmp <- .
    TP <- sum(tmp$flag & tmp$gt); TN <- sum(!tmp$flag & !tmp$gt)
    FP <- sum(tmp$flag & !tmp$gt); FN <- sum(!tmp$flag & tmp$gt)
    den <- sqrt(as.numeric(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    mcc <- if (den == 0) 0 else (TP*TN - FP*FN) / den
    data.frame(mcc = mcc, sens = TP/max(1,TP+FN), n = nrow(tmp))
  }) %>% ungroup()

pC <- ggplot(phase3_by_items, aes(x = total_items, y = mcc,
                                    color = factor(rate*100), group = factor(rate*100))) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "darkgreen") +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5) +
  geom_text(aes(label = sprintf("%.2f", mcc)),
            vjust = -1.2, size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("20" = "#377eb8", "40" = "#4daf4a", "60" = "#e41a1c"),
                     name = "Careless rate (%)") +
  scale_x_continuous(breaks = c(48, 80, 96, 160, 180, 300)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(title = "Validation MCC vs questionnaire length, by careless rate",
       subtitle = "Fresh validation data · Scenario A · FPR=5% · MCC>=0.7 (green dashed)",
       x = "Total items (nF × ipf)", y = "MCC") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")
ggsave("plot_phase4_mcc_vs_items_validation.png", pC,
       width = 12, height = 7, dpi = 150, bg = "white")

# Plot D — Confusion-matrix style summary at FPR=5%
df_p3$class <- with(df_p3,
  ifelse(gt & flag, "TP",
  ifelse(!gt & !flag, "TN",
  ifelse(!gt & flag, "FP", "FN"))))

pD <- df_p3 %>%
  group_by(rate, class) %>%
  summarise(n = n(), .groups = "drop") %>%
  ggplot(aes(x = factor(rate*100), y = n, fill = class)) +
  geom_col(position = "fill") +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5), size = 3.5) +
  scale_fill_manual(values = c("TP" = "#1a9850", "TN" = "#a6d96a",
                                "FP" = "#fdae61", "FN" = "#d73027")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Confusion-matrix proportions at FPR=5%, by careless rate",
       subtitle = "Validation data, Scenario A",
       x = "Careless rate (%)", y = "Proportion") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")
ggsave("plot_phase4_confusion_proportions.png", pD,
       width = 9, height = 6, dpi = 150, bg = "white")

cat("=== ALL PLOTS SAVED ===\n")
cat("plot_phase4_validation_mcc.png\n")
cat("plot_phase4_validation_patterns.png\n")
cat("plot_phase4_mcc_vs_items_validation.png\n")
cat("plot_phase4_confusion_proportions.png\n")
