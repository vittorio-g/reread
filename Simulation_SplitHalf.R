# Simulation_SplitHalf.R
# Compare standard RR vs split-half RR under GT=corruption>50%, auto-z threshold.
# Setup identical to Plot_Detection_vs_Corruption (same seeds): 6 conditions, 4 reps.

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(tidyr)
})
source("Synthetic_Good_Responses_2.R")
source("ReReReRe.R")
source("ReReReRe_SplitHalf.R")

# Load extended injector
inject_careless_extended <- local({
  src <- readLines("Simulation_VariancePenalty.R")
  start <- grep("^inject_careless_extended <- function", src)[1]
  end_idx <- start; brace_depth <- 0; started <- FALSE
  for (i in start:length(src)) {
    line <- src[i]
    opens <- lengths(regmatches(line, gregexpr("\\{", line)))
    closes <- lengths(regmatches(line, gregexpr("\\}", line)))
    brace_depth <- brace_depth + opens - closes
    if (opens > 0) started <- TRUE
    if (started && brace_depth == 0) { end_idx <- i; break }
  }
  eval(parse(text = paste(src[start:end_idx], collapse = "\n")))
  get("inject_careless_extended")
})

NF_LEVELS <- c(8, 16, 30)
IPF_LEVELS <- c(6, 10)
REPS <- 4
N_FIXED <- 500
PCT_CARELESS <- 0.50
ITER <- 100
CORPROP <- 0.03
SEED_BASE <- 2026
CARELESS_LEVELS <- seq(0.1, 1.0, 0.1)
GT_CUTOFF <- 0.50

conditions <- expand.grid(nF = NF_LEVELS, ipf = IPF_LEVELS, stringsAsFactors = FALSE)
conditions$total_items <- conditions$nF * conditions$ipf
TOTAL <- nrow(conditions) * REPS

cat(sprintf("%d conditions x %d reps = %d runs\n", nrow(conditions), REPS, TOTAL))
cat("Methods: std, std+VP, split-half (min), split-half (mean)\n\n")

rows <- list(); idx <- 0; t0 <- Sys.time()

for (ci in seq_len(nrow(conditions))) {
  nF <- conditions$nF[ci]; ipf <- conditions$ipf[ci]; tot <- conditions$total_items[ci]
  for (rep_id in seq_len(REPS)) {
    idx <- idx + 1
    set.seed(SEED_BASE + (rep_id - 1) * 31 + ci)

    clean <- simulated_good_responses(nConstructs = nF, nItems = rep(ipf, nF), n = N_FIXED)
    inj <- inject_careless_extended(clean, pct_careless = PCT_CARELESS,
                                    careless_levels = CARELESS_LEVELS)
    data_mat <- inj$data_corrupted

    rr_std <- ReReReRe(data_mat, corProp = CORPROP, iterations = ITER,
                       align_signs = TRUE, mode = "auto")
    rr_vp <- ReReReRe(data_mat, corProp = CORPROP, iterations = ITER,
                      align_signs = TRUE, mode = "auto",
                      variance_penalty = TRUE)
    sh <- score_split_half(data_mat, corProp = CORPROP, iterations = ITER,
                           align_signs = TRUE, n_splits = 5, aggregation = "min",
                           seed = SEED_BASE + (rep_id - 1) * 31 + ci)

    rows[[length(rows) + 1]] <- data.frame(
      nF = nF, ipf = ipf, total_items = tot, rep = rep_id,
      respondent = seq_len(nrow(data_mat)),
      pattern = inj$labels$pattern,
      corruption = round(inj$labels$careless_pct, 2),
      z_std = rr_std$z_score,
      z_std_vp = rr_vp$z_score,
      z_sh_min  = sh$z_min,
      z_sh_mean = sh$z_mean
    )

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    eta <- if (idx > 1) (elapsed / idx) * (TOTAL - idx) else NA
    cat(sprintf("[%d/%d] nF=%d ipf=%d rep=%d | %.1f min ETA %.0f\n",
                idx, TOTAL, nF, ipf, rep_id, elapsed, eta))
  }
}

df <- do.call(rbind, rows)
write.csv(df, "sim_splithalf_results.csv", row.names = FALSE)

# Auto-z per cell
cond_thresholds <- df %>%
  distinct(nF, ipf, total_items) %>%
  rowwise() %>%
  mutate(auto_z = round(get_calibrated_z(total_items), 3)) %>%
  ungroup()

df <- df %>%
  left_join(cond_thresholds, by = c("nF","ipf","total_items")) %>%
  mutate(gt_careless = (pattern != "clean") & (corruption > GT_CUTOFF),
         f_std      = z_std     <= auto_z,
         f_std_vp   = z_std_vp  <= auto_z,
         f_sh_min   = z_sh_min  <= auto_z,
         f_sh_mean  = z_sh_mean <= auto_z)

mcc_fun <- function(flag, truth) {
  TP <- sum(flag & truth); TN <- sum(!flag & !truth)
  FP <- sum(flag & !truth); FN <- sum(!flag & truth)
  den <- sqrt(as.numeric(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  c(sens = TP / max(1, TP+FN),
    spec = TN / max(1, TN+FP),
    mcc  = if (den == 0) 0 else (TP*TN - FP*FN) / den)
}

cat("\n=== Overall metrics (auto-z, GT=corruption>50%) ===\n")
for (col in c("f_std","f_std_vp","f_sh_min","f_sh_mean")) {
  m <- mcc_fun(df[[col]], df$gt_careless)
  cat(sprintf("%-12s  sens=%.3f  spec=%.3f  MCC=%.3f\n",
              col, m["sens"], m["spec"], m["mcc"]))
}

cat("\n=== Per-pattern sensitivity (truly careless, corruption>50%) ===\n")
per_pat <- df %>%
  filter(gt_careless) %>%
  group_by(pattern) %>%
  summarise(
    n = n(),
    std      = round(mean(f_std),     3),
    std_vp   = round(mean(f_std_vp),  3),
    sh_min   = round(mean(f_sh_min),  3),
    sh_mean  = round(mean(f_sh_mean), 3),
    .groups = "drop"
  )
print(as.data.frame(per_pat), row.names = FALSE)

cat("\n=== Detection rate by corruption bin (low-items: 10-50%) ===\n")
df_low <- df %>%
  filter(pattern != "clean") %>%
  mutate(corruption_bin10 = round(corruption * 10) / 10) %>%
  group_by(pattern, corruption_bin10) %>%
  summarise(
    n = n(),
    det_std     = round(mean(f_std),     3),
    det_sh_min  = round(mean(f_sh_min),  3),
    det_sh_mean = round(mean(f_sh_mean), 3),
    .groups = "drop"
  )
print(as.data.frame(df_low), row.names = FALSE)

# ============================================================
# PLOTS
# ============================================================
pattern_colors <- c(
  random = "#e41a1c", longstring = "#377eb8", mixed = "#4daf4a",
  pure_straight = "#984ea3", acquiescent = "#ff7f00"
)

# Plot 1 — MCC bars per strategy
mcc_all <- data.frame(
  strategy = c("std","std+VP","split-half min","split-half mean"),
  row.names = NULL,
  t(sapply(c("f_std","f_std_vp","f_sh_min","f_sh_mean"),
           function(col) mcc_fun(df[[col]], df$gt_careless)))
) %>% mutate(across(c(sens, spec, mcc), ~ round(.x, 3)))

bar_df <- mcc_all %>%
  pivot_longer(cols = c(sens, spec, mcc), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, "sens"="Sensitivity","spec"="Specificity","mcc"="MCC"),
         metric = factor(metric, levels = c("Sensitivity","Specificity","MCC")),
         strategy = factor(strategy,
                           levels = c("std","std+VP","split-half min","split-half mean")))

p1 <- ggplot(bar_df, aes(x = strategy, y = value, fill = metric)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", value)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = c("Sensitivity"="#377eb8","Specificity"="#4daf4a","MCC"="#e41a1c")) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Split-half vs standard under GT = corruption > 50% (auto-z)",
       subtitle = "n_splits=5, aggregation=min (or mean)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))
ggsave("plot_splithalf_metrics.png", p1, width = 11, height = 6, dpi = 150, bg = "white")

# Plot 2 — detection curve std vs split-half-min
cmp_df <- df %>%
  filter(pattern != "clean") %>%
  mutate(corruption_bin10 = round(corruption * 10) / 10) %>%
  group_by(pattern, corruption_bin10) %>%
  summarise(
    Standard          = mean(f_std),
    `Split-half min`  = mean(f_sh_min),
    .groups = "drop"
  ) %>%
  rename(corruption = corruption_bin10) %>%
  pivot_longer(cols = c(Standard, `Split-half min`),
               names_to = "method", values_to = "det_rate")

p2 <- ggplot(cmp_df, aes(x = corruption, y = det_rate,
                         color = pattern, group = pattern)) +
  annotate("rect", xmin = 0.05, xmax = 0.55, ymin = 0, ymax = 1,
           fill = "#ffcccc", alpha = 0.3) +
  annotate("rect", xmin = 0.55, xmax = 1.05, ymin = 0, ymax = 1,
           fill = "#ccffcc", alpha = 0.3) +
  geom_vline(xintercept = 0.55, color = "grey30", linewidth = 0.4) +
  geom_line(linewidth = 1.2) + geom_point(size = 2.5) +
  scale_color_manual(values = pattern_colors, name = "Pattern") +
  scale_x_continuous(breaks = seq(0.1, 1.0, 0.1),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  facet_wrap(~ method, ncol = 2) +
  labs(title = "Detection curve: standard vs split-half min",
       subtitle = "Left zone: flagging = FP. Right zone: flagging = TP.",
       x = "Corruption level", y = "Flagging rate at auto-z") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "#EEEEEE", color = NA),
        strip.text = element_text(face = "bold"))
ggsave("plot_splithalf_curves.png", p2, width = 13, height = 6, dpi = 150, bg = "white")

cat(sprintf("\nRuntime: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("Saved: sim_splithalf_results.csv, plot_splithalf_metrics.png, plot_splithalf_curves.png\n")
cat("=== DONE ===\n")
