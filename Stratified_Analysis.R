#### Stratified Analysis: Detection by Corruption Level ####
# === REVISION 2026-03-10f ===
#
# Uses z-score metric (primary) with fixed z_threshold.
# Also reports percentile for comparison.
#
# The multiverse pools all corruption levels (50-100%).
# But a respondent with 50% corruption still has 50% genuine responses —
# probably harder to detect. This script answers:
#
# 1. Does detection improve at higher corruption levels?
# 2. Which patterns (random/longstring/mixed) are detected at each level?
# 3. What's the "operating range" of ReReReRe?
# 4. How does MCC change if we only count severe corruption as positive?
#
# This runs a few representative conditions (not the full multiverse)
# for quick diagnostics. ~5-10 minutes.

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(lavaan)
library(psych)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("Synthetic_Good_Responses_2.R")
source("Careless_machine_2.R")
source("ReReReRe.R")

# ============================================================
# PARAMETERS
# ============================================================

NFACTORS_LEVELS  <- c(4, 8, 15, 25)
ITEMS_CYCLE      <- c(10, 6, 3)
N_RESPONDENTS    <- 500        # large enough for stable rates
PCT_CARELESS     <- 0.25       # 25% careless → ~125 careless respondents
CORPROP          <- 0.10       # single corProp for speed
ITERATIONS       <- 100
MIN_PAIRS        <- 15
Z_THRESHOLDS     <- c(0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5)

CARELESS_TYPES   <- c(random = 1/3, longstring = 1/3, mixed = 1/3)
CARELESS_LEVELS  <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5)

make_nItems <- function(nf) rep_len(ITEMS_CYCLE, nf)

# ============================================================
# RUN DIAGNOSTICS
# ============================================================

all_respondent_data <- list()

for (nf in NFACTORS_LEVELS) {

  nItems_vec <- make_nItems(nf)
  total_items <- sum(nItems_vec)

  cat(sprintf("\n=== nFactors=%d (%d items) ===\n", nf, total_items))

  # Generate clean data
  clean <- tryCatch(
    simulated_good_responses(nConstructs = nf, nItems = nItems_vec, n = N_RESPONDENTS),
    error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(clean)) next

  # Inject careless
  inj <- inject_careless(
    data = clean,
    pct_careless = PCT_CARELESS,
    pct_types = CARELESS_TYPES,
    careless_levels = CARELESS_LEVELS
  )

  # Run ReReReRe
  rr <- ReReReRe(
    data = inj$data_corrupted,
    corProp = CORPROP,
    cutOff = 0,
    iterations = ITERATIONS,
    min_pairs = MIN_PAIRS,
    progress = TRUE
  )

  # Combine per-respondent data
  df <- data.frame(
    nFactors = nf,
    total_items = total_items,
    respondent = seq_len(N_RESPONDENTS),
    z_score = rr$z_score,
    pctile = rr$result,
    indCors = rr$indCors,
    rand_mean = rr$rand_mean,
    rand_sd = rr$rand_sd,
    careless = inj$labels$careless,
    careless_pct = inj$labels$careless_pct,
    pattern = inj$labels$pattern,
    stringsAsFactors = FALSE
  )

  # Bin corruption level (round to nearest 10%)
  df$corruption_bin <- ifelse(
    df$careless,
    paste0(round(df$careless_pct * 100 / 10) * 10, "%"),
    "clean"
  )

  all_respondent_data[[as.character(nf)]] <- df
}

# Combine all conditions
D <- do.call(rbind, all_respondent_data)
rownames(D) <- NULL

cat("\n\n========================================\n")
cat("STRATIFIED ANALYSIS RESULTS\n")
cat("========================================\n\n")

# ============================================================
# 1. MEAN Z-SCORE BY CORRUPTION LEVEL
# ============================================================

cat("=== 1. MEAN Z-SCORE BY CORRUPTION LEVEL ===\n\n")

corr_summary <- D %>%
  group_by(nFactors, corruption_bin) %>%
  summarise(
    n = n(),
    mean_z = round(mean(z_score, na.rm=T), 3),
    sd_z = round(sd(z_score, na.rm=T), 3),
    mean_pctile = round(mean(pctile, na.rm=T), 3),
    mean_indCors = round(mean(indCors, na.rm=T), 3),
    .groups = "drop"
  ) %>%
  arrange(nFactors, corruption_bin)

print(corr_summary, n = 100)

# Wide format for readability
cat("\n--- MEAN Z-SCORE (wide) ---\n")
corr_wide <- D %>%
  group_by(nFactors, corruption_bin) %>%
  summarise(mean_z = round(mean(z_score, na.rm=T), 3), .groups = "drop") %>%
  pivot_wider(names_from = corruption_bin, values_from = mean_z)
print(corr_wide, width = 200)

# ============================================================
# 2. DETECTION RATE BY CORRUPTION LEVEL (at z_threshold = 2)
# ============================================================

cat("\n\n=== 2. DETECTION RATE BY CORRUPTION LEVEL (z_threshold=2) ===\n\n")

DISPLAY_ZT <- 2

det_by_corr <- D %>%
  filter(careless) %>%
  mutate(flagged = z_score <= DISPLAY_ZT) %>%
  group_by(nFactors, corruption_bin) %>%
  summarise(
    n = n(),
    n_flagged = sum(flagged),
    detection_rate = round(mean(flagged), 3),
    mean_z = round(mean(z_score), 3),
    .groups = "drop"
  ) %>%
  arrange(nFactors, corruption_bin)

print(det_by_corr, n = 100)

# ============================================================
# 3. DETECTION RATE BY CORRUPTION x PATTERN
# ============================================================

cat("\n\n=== 3. DETECTION RATE BY CORRUPTION x PATTERN (z_threshold=2) ===\n\n")

det_by_corr_type <- D %>%
  filter(careless) %>%
  mutate(flagged = z_score <= DISPLAY_ZT) %>%
  group_by(nFactors, pattern, corruption_bin) %>%
  summarise(
    n = n(),
    detection_rate = round(mean(flagged), 3),
    mean_z = round(mean(z_score), 3),
    .groups = "drop"
  ) %>%
  arrange(nFactors, pattern, corruption_bin)

print(det_by_corr_type, n = 200)

# ============================================================
# 4. MCC AT VARYING EVAL THRESHOLDS
#    ("Only count corruption >= X% as truly careless")
# ============================================================

cat("\n\n=== 4. MCC BY EVAL THRESHOLD (varying what counts as 'truly careless') ===\n\n")

eval_thresholds <- c(0.01, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00)

mcc_by_eval <- expand.grid(
  nFactors = NFACTORS_LEVELS,
  z_threshold = Z_THRESHOLDS,
  eval_threshold = eval_thresholds,
  stringsAsFactors = FALSE
)

mcc_results <- lapply(seq_len(nrow(mcc_by_eval)), function(i) {
  row <- mcc_by_eval[i, ]
  dd <- D %>% filter(nFactors == row$nFactors)

  actual <- dd$careless_pct >= row$eval_threshold
  flagged <- dd$z_score <= row$z_threshold

  TP <- sum(flagged & actual)
  TN <- sum(!flagged & !actual)
  FP <- sum(flagged & !actual)
  FN <- sum(!flagged & actual)

  mcc_denom <- sqrt(as.numeric(TP+FP) * (TP+FN) * (TN+FP) * (TN+FN))
  mcc <- if (mcc_denom > 0) (TP*TN - FP*FN) / mcc_denom else NA_real_
  sens <- if ((TP+FN) > 0) TP / (TP+FN) else NA_real_
  spec <- if ((TN+FP) > 0) TN / (TN+FP) else NA_real_

  data.frame(
    nFactors = row$nFactors,
    z_threshold = row$z_threshold,
    eval_threshold = row$eval_threshold,
    n_positive = sum(actual),
    TP = TP, FP = FP, FN = FN, TN = TN,
    sensitivity = round(sens, 3),
    specificity = round(spec, 3),
    mcc = round(mcc, 3),
    stringsAsFactors = FALSE
  )
}) %>% do.call(rbind, .)

# Best MCC per nFactors x eval_threshold
cat("--- Best MCC per condition (optimizing over z_threshold) ---\n\n")
best_mcc <- mcc_results %>%
  group_by(nFactors, eval_threshold) %>%
  slice_max(mcc, n = 1, with_ties = FALSE) %>%
  select(nFactors, eval_threshold, n_positive, z_threshold, sensitivity, specificity, mcc) %>%
  arrange(nFactors, eval_threshold)

print(best_mcc, n = 100)

# Wide format: MCC by nFactors x eval_threshold (at best z_threshold)
cat("\n--- Best MCC (wide) ---\n")
best_wide <- best_mcc %>%
  select(nFactors, eval_threshold, mcc) %>%
  pivot_wider(names_from = eval_threshold, values_from = mcc,
              names_prefix = "eval_")
print(best_wide, width = 200)

# ============================================================
# 5. Z-SCORE vs PERCENTILE COHEN'S d COMPARISON
# ============================================================

cat("\n\n=== 5. Z-SCORE vs PERCENTILE COHEN'S d ===\n")
cat("(Confirms z-score outperforms percentile after longstring fix)\n\n")

cohens_d <- function(x_good, x_care) {
  n1 <- length(x_good); n2 <- length(x_care)
  m1 <- mean(x_good, na.rm=T); m2 <- mean(x_care, na.rm=T)
  s1 <- sd(x_good, na.rm=T); s2 <- sd(x_care, na.rm=T)
  pooled_sd <- sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2))
  (m1 - m2) / pooled_sd
}

d_comparison <- D %>%
  group_by(nFactors) %>%
  summarise(
    d_zscore = round(cohens_d(z_score[!careless], z_score[careless]), 3),
    d_pctile = round(cohens_d(pctile[!careless], pctile[careless]), 3),
    d_indCors = round(cohens_d(indCors[!careless], indCors[careless]), 3),
    z_advantage = paste0(round((cohens_d(z_score[!careless], z_score[careless]) /
                                cohens_d(pctile[!careless], pctile[careless]) - 1) * 100), "%"),
    mean_z_good = round(mean(z_score[!careless], na.rm=T), 2),
    mean_z_care = round(mean(z_score[careless], na.rm=T), 2),
    .groups = "drop"
  )
print(d_comparison)

# ============================================================
# 6. PLOTS
# ============================================================

# Plot 1: Mean z-score by corruption level and nFactors
p1 <- D %>%
  group_by(nFactors, corruption_bin) %>%
  summarise(mean_z = mean(z_score, na.rm=T), .groups="drop") %>%
  mutate(corruption_bin = factor(corruption_bin,
    levels = c("clean", "50%", "60%", "70%", "80%", "90%", "100%"))) %>%
  ggplot(aes(x = corruption_bin, y = mean_z,
             color = factor(nFactors), group = nFactors)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set2") +
  labs(title = "Mean Z-Score by Corruption Level",
       subtitle = "Higher corruption = lower z-score (more detectable)",
       x = "Corruption Level", y = "Mean Z-Score", color = "nFactors") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 2: Detection rate by corruption level (z_threshold=2)
p2 <- D %>%
  filter(careless) %>%
  mutate(flagged = z_score <= DISPLAY_ZT) %>%
  group_by(nFactors, corruption_bin) %>%
  summarise(det_rate = mean(flagged), .groups="drop") %>%
  mutate(corruption_bin = factor(corruption_bin,
    levels = c("50%", "60%", "70%", "80%", "90%", "100%"))) %>%
  ggplot(aes(x = corruption_bin, y = det_rate,
             color = factor(nFactors), group = nFactors)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set2") +
  labs(title = paste0("Detection Rate by Corruption Level (z_threshold=", DISPLAY_ZT, ")"),
       subtitle = "Method catches high corruption, misses partial",
       x = "Corruption Level", y = "Detection Rate", color = "nFactors") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 3: Detection by corruption x pattern
p3 <- D %>%
  filter(careless) %>%
  mutate(flagged = z_score <= DISPLAY_ZT) %>%
  group_by(nFactors, pattern, corruption_bin) %>%
  summarise(det_rate = mean(flagged), .groups="drop") %>%
  mutate(corruption_bin = factor(corruption_bin,
    levels = c("50%", "60%", "70%", "80%", "90%", "100%"))) %>%
  ggplot(aes(x = corruption_bin, y = det_rate,
             color = pattern, group = pattern)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  facet_wrap(~nFactors, labeller = label_both) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "Detection Rate by Corruption Level and Pattern Type",
       subtitle = paste0("z_threshold=", DISPLAY_ZT, " | Faceted by nFactors"),
       x = "Corruption Level", y = "Detection Rate", color = "Pattern") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 4: MCC by eval threshold
p4 <- best_mcc %>%
  ggplot(aes(x = eval_threshold, y = mcc,
             color = factor(nFactors), group = nFactors)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set2") +
  labs(title = "Best MCC by Evaluation Threshold",
       subtitle = "Higher threshold = only severe corruption counts as positive",
       x = "Eval Threshold (min corruption % to count as careless)",
       y = "Best MCC (optimizing over z_threshold)", color = "nFactors") +
  theme_minimal()

# Plot 5: Cohen's d comparison
p5 <- d_comparison %>%
  pivot_longer(cols = c(d_zscore, d_pctile, d_indCors),
               names_to = "metric", values_to = "cohens_d") %>%
  mutate(metric = case_when(
    metric == "d_zscore" ~ "Z-Score",
    metric == "d_pctile" ~ "Percentile",
    metric == "d_indCors" ~ "indCors"
  )) %>%
  ggplot(aes(x = factor(nFactors), y = cohens_d, fill = metric)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Z-Score" = "steelblue",
                                "Percentile" = "gray60",
                                "indCors" = "firebrick")) +
  labs(title = "Cohen's d: Z-Score vs Percentile vs indCors",
       subtitle = "Z-score advantage grows with questionnaire size",
       x = "nFactors", y = "Cohen's d (Good vs Careless)", fill = "") +
  theme_minimal()

print(p1)
print(p2)
print(p3)
print(p4)
print(p5)

# Combined
combined <- (p1 | p2) / (p3) / (p4 | p5) +
  plot_annotation(
    title = "ReReReRe Stratified Analysis (Z-Score)",
    subtitle = paste0("N=", N_RESPONDENTS, " | pct_careless=", PCT_CARELESS,
                      " | corProp=", CORPROP, " | iterations=", ITERATIONS),
    theme = theme(plot.title = element_text(face = "bold", size = 16))
  )

ggsave("stratified_analysis.png", combined, width = 16, height = 20, dpi = 200)
cat("\nSaved: stratified_analysis.png\n")
