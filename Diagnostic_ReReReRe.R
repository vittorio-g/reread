#### Diagnostic: Score Distributions for ReReReRe ####
#
# Purpose: generate datasets across conditions, run ReReReRe,
# and compare FOUR scoring approaches side-by-side:
#   1. Percentile (corComparedIndex) — original, rank-based
#   2. Z-score — (indCors - rand_mean) / rand_sd (per-respondent standardized)
#   3. Delta — indCors - rand_mean (per-respondent, unstandardized)
#   4. indCors — raw coupled absolute correlation
#
# The key question: which scoring approach best separates good vs careless?
#
# === REVISION 2026-03-10c ===
# Z-score approach failed: per-respondent rand_sd denominator is noisy and
# confounded with careless status, destroying Cohen's d (worse than percentile
# in 7/8 conditions). Additional issue: zero_val=1 for random pairs contaminates
# rand_mean and rand_sd for partial longstring respondents.
#
# This diagnostic compares all four approaches to find the best one empirically.
# ===

rm(list = ls())

library(dplyr)
library(ggplot2)
library(patchwork)
library(lavaan)
library(psych)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("Synthetic_Good_Responses_2.R")
source("Careless_machine_2.R")
source("ReReReRe.R")

# ============================================================
# DIAGNOSTIC CONDITIONS
# ============================================================

conditions <- data.frame(
  nFactors      = c(  4,   8,  10,  15,  20,  20,  20,  25),
  n_respondents = c(300, 300, 300, 300, 300, 300, 1000, 500),
  pct_careless  = c(.10, .10, .10, .10, .10, .25,  .10, .10),
  corProp       = c(.10, .10, .10, .05, .05, .05,  .05, .05),
  stringsAsFactors = FALSE
)

ITEMS_CYCLE     <- c(10, 6, 3)
CARELESS_TYPES  <- c(random = 1/3, longstring = 1/3, mixed = 1/3)
CARELESS_LEVELS <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5)
ITERATIONS      <- 100

make_nItems <- function(nf) rep_len(ITEMS_CYCLE, nf)

# ============================================================
# RUN DIAGNOSTICS
# ============================================================

all_scores <- list()

for (i in seq_len(nrow(conditions))) {
  cond <- conditions[i, ]
  nItems_vec <- make_nItems(cond$nFactors)
  total_items <- sum(nItems_vec)

  label <- sprintf("nF=%d (%d items), n=%d, pct=%.0f%%",
                   cond$nFactors, total_items, cond$n_respondents,
                   cond$pct_careless * 100)
  cat("Running:", label, "\n")

  clean <- simulated_good_responses(
    nConstructs = cond$nFactors,
    nItems = nItems_vec,
    n = cond$n_respondents
  )

  injection <- inject_careless(
    data = clean,
    pct_careless = cond$pct_careless,
    pct_types = CARELESS_TYPES,
    careless_levels = CARELESS_LEVELS
  )

  rr <- ReReReRe(
    data = injection$data_corrupted,
    corProp = cond$corProp,
    cutOff = 0,
    iterations = ITERATIONS,
    min_pairs = 15,
    progress = FALSE
  )

  # Compute delta (unstandardized per-respondent difference)
  delta <- rr$indCors - rr$rand_mean

  scores <- data.frame(
    condition = label,
    percentile = rr$result,       # P(random < coupled)
    z_score = rr$z_score,         # (indCors - rand_mean) / rand_sd
    delta = delta,                # indCors - rand_mean (no SD normalization)
    indCors = rr$indCors,         # raw coupled correlation
    rand_mean = rr$rand_mean,
    rand_sd = rr$rand_sd,
    is_careless = injection$labels$careless,
    careless_pct = injection$labels$careless_pct,
    pattern = injection$labels$pattern,
    stringsAsFactors = FALSE
  )

  # Classify for plotting
  scores$group <- ifelse(!scores$is_careless, "Good",
                    ifelse(scores$careless_pct >= 0.5, "Careless (>=50%)",
                           "Careless (<50%)"))

  all_scores[[i]] <- scores
}

df <- do.call(rbind, all_scores)

# ============================================================
# THE KEY TABLE: Cohen's d across all four metrics
# ============================================================

cat("\n=== COHEN'S d COMPARISON ACROSS ALL FOUR METRICS ===\n")
cat("(Higher = better discrimination between Good and Careless >=50%)\n\n")

cohens_d_fn <- function(good_vals, careless_vals) {
  mg <- mean(good_vals, na.rm=T); mc <- mean(careless_vals, na.rm=T)
  sg <- var(good_vals, na.rm=T); sc <- var(careless_vals, na.rm=T)
  sp <- sqrt((sg + sc) / 2)
  if (sp == 0) return(NA_real_)
  (mg - mc) / sp
}

comparison <- df %>%
  filter(group %in% c("Good", "Careless (>=50%)")) %>%
  group_by(condition) %>%
  summarise(
    d_percentile = cohens_d_fn(percentile[group == "Good"], percentile[group == "Careless (>=50%)"]),
    d_zscore     = cohens_d_fn(z_score[group == "Good"], z_score[group == "Careless (>=50%)"]),
    d_delta      = cohens_d_fn(delta[group == "Good"], delta[group == "Careless (>=50%)"]),
    d_indCors    = cohens_d_fn(indCors[group == "Good"], indCors[group == "Careless (>=50%)"]),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  mutate(best = case_when(
    d_percentile >= d_zscore & d_percentile >= d_delta & d_percentile >= d_indCors ~ "percentile",
    d_zscore >= d_percentile & d_zscore >= d_delta & d_zscore >= d_indCors ~ "z_score",
    d_delta >= d_percentile & d_delta >= d_zscore & d_delta >= d_indCors ~ "delta",
    TRUE ~ "indCors"
  ))

print(comparison, width = 120)

# ============================================================
# SUMMARY STATS BY GROUP FOR ALL METRICS
# ============================================================

cat("\n=== SUMMARY BY CONDITION AND GROUP ===\n")
df %>%
  group_by(condition, group) %>%
  summarise(
    n = n(),
    mean_pctile = round(mean(percentile, na.rm=T), 3),
    mean_z = round(mean(z_score, na.rm=T), 3),
    sd_z = round(sd(z_score, na.rm=T), 3),
    mean_delta = round(mean(delta, na.rm=T), 4),
    sd_delta = round(sd(delta, na.rm=T), 4),
    mean_indCors = round(mean(indCors, na.rm=T), 4),
    mean_rand = round(mean(rand_mean, na.rm=T), 4),
    mean_rand_sd = round(mean(rand_sd, na.rm=T), 5),
    .groups = "drop"
  ) %>%
  print(n = 50, width = 150)

# ============================================================
# PLOT 1: Delta distributions (the key test)
# ============================================================

p_delta <- ggplot(df, aes(x = delta, fill = group)) +
  geom_histogram(bins = 40, alpha = 0.6, position = "identity") +
  facet_wrap(~condition, scales = "free", ncol = 2) +
  scale_fill_manual(values = c("Good" = "steelblue",
                                "Careless (>=50%)" = "firebrick",
                                "Careless (<50%)" = "orange")) +
  labs(title = "Delta Distributions: Good vs Careless",
       subtitle = "delta = indCors - rand_mean | Preserves baseline calibration without noisy SD",
       x = "Delta (coupled - random baseline)", y = "Count", fill = "") +
  theme_minimal() +
  theme(legend.position = "bottom")

# ============================================================
# PLOT 2: Percentile distributions (for comparison)
# ============================================================

p_pctile <- ggplot(df, aes(x = percentile, fill = group)) +
  geom_histogram(bins = 40, alpha = 0.6, position = "identity") +
  facet_wrap(~condition, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Good" = "steelblue",
                                "Careless (>=50%)" = "firebrick",
                                "Careless (<50%)" = "orange")) +
  labs(title = "Percentile Distributions: Good vs Careless",
       subtitle = "corComparedIndex = P(random < coupled) | Note ceiling effect at 1.0",
       x = "Percentile", y = "Count", fill = "") +
  theme_minimal() +
  theme(legend.position = "bottom")

# ============================================================
# PLOT 3: Z-score distributions (for comparison)
# ============================================================

p_zscore <- ggplot(df, aes(x = z_score, fill = group)) +
  geom_histogram(bins = 40, alpha = 0.6, position = "identity") +
  facet_wrap(~condition, scales = "free", ncol = 2) +
  scale_fill_manual(values = c("Good" = "steelblue",
                                "Careless (>=50%)" = "firebrick",
                                "Careless (<50%)" = "orange")) +
  labs(title = "Z-Score Distributions: Good vs Careless",
       subtitle = "z = (indCors - rand_mean) / rand_sd | Noisy denominator inflates variance",
       x = "Z-Score", y = "Count", fill = "") +
  theme_minimal() +
  theme(legend.position = "bottom")

# ============================================================
# PLOT 4: Box plots comparing all metrics (one condition, the best one)
# ============================================================

# Pick the condition with best overall separation for the box plot comparison
best_cond <- comparison$condition[which.max(rowMeans(
  comparison[, c("d_percentile", "d_delta", "d_indCors")], na.rm = TRUE))]

df_best <- df %>% filter(condition == best_cond)

box_data <- df_best %>%
  select(group, percentile, z_score, delta, indCors) %>%
  tidyr::pivot_longer(cols = c(percentile, z_score, delta, indCors),
               names_to = "metric", values_to = "value")

p_box_compare <- ggplot(box_data, aes(x = group, y = value, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  facet_wrap(~metric, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Good" = "steelblue",
                                "Careless (>=50%)" = "firebrick",
                                "Careless (<50%)" = "orange")) +
  labs(title = paste("All Four Metrics Compared —", best_cond),
       x = "", y = "Score", fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")

# ============================================================
# PLOT 5: Cohen's d across conditions for each metric
# ============================================================

d_long <- comparison %>%
  select(condition, d_percentile, d_zscore, d_delta, d_indCors) %>%
  tidyr::pivot_longer(cols = starts_with("d_"),
               names_to = "metric", values_to = "cohens_d",
               names_prefix = "d_")

p_d_compare <- ggplot(d_long, aes(x = condition, y = cohens_d, fill = metric)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Cohen's d by Metric and Condition",
       subtitle = "Which scoring approach best separates Good from Careless?",
       x = "", y = "Cohen's d", fill = "Metric") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")

# ============================================================
# DETECTION BY PATTERN TYPE (key diagnostic for longstring fix)
# ============================================================

cat("\n=== MEAN PERCENTILE BY PATTERN TYPE ===\n")
cat("(After longstring fix: longstring should be MUCH lower than before)\n\n")

df %>%
  mutate(pattern_group = case_when(
    !is_careless ~ "clean",
    pattern == "random" ~ "random",
    pattern == "longstring" ~ "longstring",
    pattern == "mixed" ~ "mixed",
    TRUE ~ pattern
  )) %>%
  group_by(condition, pattern_group) %>%
  summarise(
    n = n(),
    mean_pctile = round(mean(percentile, na.rm=T), 3),
    sd_pctile = round(sd(percentile, na.rm=T), 3),
    mean_indCors = round(mean(indCors, na.rm=T), 3),
    mean_rand_mean = round(mean(rand_mean, na.rm=T), 3),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    id_cols = condition,
    names_from = pattern_group,
    values_from = mean_pctile,
    names_prefix = "pctile_"
  ) %>%
  print(width = 120)

cat("\n=== DETECTION RATE BY PATTERN (cutOff = 0.80) ===\n")
df %>%
  filter(is_careless, careless_pct >= 0.01) %>%
  mutate(flagged_80 = percentile <= 0.80) %>%
  group_by(condition, pattern) %>%
  summarise(
    n = n(),
    detection_rate = round(mean(flagged_80), 3),
    mean_pctile = round(mean(percentile, na.rm=T), 3),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    id_cols = condition,
    names_from = pattern,
    values_from = detection_rate,
    names_prefix = "det_"
  ) %>%
  print(width = 120)

# ============================================================
# NORMATIVE FLAGGING ON DELTA
# ============================================================

cat("\n=== NORMATIVE FLAGGING ON DELTA ===\n")
for (k_val in c(1.5, 2.0, 2.5, 3.0)) {
  cat(sprintf("\n--- k_mad = %.1f ---\n", k_val))
  df %>%
    group_by(condition) %>%
    summarise(
      d_med = round(median(delta, na.rm=T), 4),
      d_mad = round(mad(delta, na.rm=T), 4),
      threshold = round(median(delta, na.rm=T) - k_val * mad(delta, na.rm=T), 4),
      n_flagged = sum(delta < (median(delta, na.rm=T) - k_val * mad(delta, na.rm=T))),
      pct_flagged = round(100 * mean(delta < (median(delta, na.rm=T) - k_val * mad(delta, na.rm=T))), 1),
      sensitivity = {
        thresh <- median(delta, na.rm=T) - k_val * mad(delta, na.rm=T)
        flagged <- delta < thresh
        n_pos <- sum(careless_pct >= 0.01)
        if (n_pos > 0) round(sum(flagged & careless_pct >= 0.01) / n_pos, 3) else NA_real_
      },
      precision = {
        thresh <- median(delta, na.rm=T) - k_val * mad(delta, na.rm=T)
        flagged <- delta < thresh
        if (sum(flagged) > 0) round(sum(flagged & careless_pct >= 0.01) / sum(flagged), 3) else NA_real_
      },
      .groups = "drop"
    ) %>%
    print()
}

# ============================================================
# NORMATIVE FLAGGING ON PERCENTILE
# ============================================================

cat("\n=== NORMATIVE FLAGGING ON PERCENTILE ===\n")
for (k_val in c(1.5, 2.0, 2.5, 3.0)) {
  cat(sprintf("\n--- k_mad = %.1f ---\n", k_val))
  df %>%
    group_by(condition) %>%
    summarise(
      p_med = round(median(percentile, na.rm=T), 4),
      p_mad = round(mad(percentile, na.rm=T), 4),
      threshold = round(median(percentile, na.rm=T) - k_val * mad(percentile, na.rm=T), 4),
      n_flagged = sum(percentile < (median(percentile, na.rm=T) - k_val * mad(percentile, na.rm=T))),
      pct_flagged = round(100 * mean(percentile < (median(percentile, na.rm=T) - k_val * mad(percentile, na.rm=T))), 1),
      sensitivity = {
        thresh <- median(percentile, na.rm=T) - k_val * mad(percentile, na.rm=T)
        flagged <- percentile < thresh
        n_pos <- sum(careless_pct >= 0.01)
        if (n_pos > 0) round(sum(flagged & careless_pct >= 0.01) / n_pos, 3) else NA_real_
      },
      precision = {
        thresh <- median(percentile, na.rm=T) - k_val * mad(percentile, na.rm=T)
        flagged <- percentile < thresh
        if (sum(flagged) > 0) round(sum(flagged & careless_pct >= 0.01) / sum(flagged), 3) else NA_real_
      },
      .groups = "drop"
    ) %>%
    print()
}

# ============================================================
# DISPLAY AND SAVE
# ============================================================

print(p_delta)
print(p_pctile)
print(p_zscore)
print(p_box_compare)
print(p_d_compare)

combined <- (p_d_compare) / (p_delta | p_pctile) / (p_zscore | p_box_compare) +
  plot_annotation(title = "ReReReRe Diagnostic: Comparing Four Scoring Approaches",
                  subtitle = "Percentile vs Z-score vs Delta vs Raw indCors",
                  theme = theme(plot.title = element_text(face = "bold", size = 16)))

ggsave("diagnostic_distributions.png", combined, width = 18, height = 22, dpi = 200)
cat("\nSaved: diagnostic_distributions.png\n")
