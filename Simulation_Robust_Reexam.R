# Simulation_Robust_Reexam.R
# Full re-examination with focus on detecting partial carelessness (60-80%).
# Larger N (8 reps × 6 conditions = 48 datasets = 24,000 respondents).
# Includes hyperparameter tuning for iter_EFA and multi-threshold evaluation.

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(tidyr); library(MASS)
})
source("Synthetic_Good_Responses_2.R")
source("ReReReRe.R")
source("Careless_machine_v3.R")
source("ReReReRe_IterCoupled.R")

# --- Config ---
NF_LEVELS  <- c(8, 16, 30)
IPF_LEVELS <- c(6, 10)
REPS <- 8                 # doubled vs v2
N_FIXED <- 500
PCT_CARELESS <- 0.40
ITER <- 100
CORPROP <- 0.03
SEED_BASE <- 20260423     # fresh seeds, different from v2
CARELESS_LEVELS <- seq(0.1, 1.0, 0.1)
GT_THRESHOLDS <- c(0.40, 0.50, 0.60, 0.70, 0.80)
TRIM_PCTS <- c(0.10, 0.20, 0.30)
N_FOLDS <- 5

conditions <- expand.grid(nF = NF_LEVELS, ipf = IPF_LEVELS, stringsAsFactors = FALSE)
conditions$total_items <- conditions$nF * conditions$ipf
TOTAL <- nrow(conditions) * REPS

cat(sprintf("ROBUST RE-EXAMINATION: %d conditions x %d reps = %d datasets\n",
            nrow(conditions), REPS, TOTAL))
cat(sprintf("Total respondents: %d  |  GT thresholds tested: %s\n",
            TOTAL * N_FIXED, paste(GT_THRESHOLDS, collapse = ", ")))
cat(sprintf("iter_EFA trim_pct tuning: %s\n",
            paste(TRIM_PCTS, collapse = ", ")))

# --- Detector helpers ---
compute_irv <- function(d) apply(d, 1, sd, na.rm = TRUE)
compute_longstring <- function(d) apply(d, 1, function(r) {
  r <- r[!is.na(r)]; if (!length(r)) 0 else max(rle(r)$lengths)
})
compute_d2 <- function(d) {
  m <- as.matrix(d); m <- m[, apply(m, 2, function(x) sd(x, na.rm=TRUE) > 0), drop=FALSE]
  if (ncol(m) < 2) return(rep(0, nrow(d)))
  cov_mat <- cov(m, use = "pairwise.complete.obs")
  center <- colMeans(m, na.rm = TRUE)
  cov_inv <- tryCatch(ginv(cov_mat),
                      error = function(e) solve(cov_mat + diag(1e-6, ncol(cov_mat))))
  d2 <- apply(m, 1, function(x) {dlt <- x - center; as.numeric(t(dlt) %*% cov_inv %*% dlt)})
  d2[is.na(d2) | is.nan(d2)] <- 0; d2
}
compute_person_total <- function(d) {
  m <- as.matrix(d); im <- colMeans(m, na.rm = TRUE)
  apply(m, 1, function(r) if (sd(r, na.rm=TRUE)==0 || sd(im)==0) 0
        else suppressWarnings(cor(r, im, use = "pairwise.complete.obs")))
}

# --- Main loop ---
rows <- list(); idx <- 0; t0 <- Sys.time()
for (ci in seq_len(nrow(conditions))) {
  nF <- conditions$nF[ci]; ipf <- conditions$ipf[ci]; tot <- conditions$total_items[ci]
  for (rep_id in seq_len(REPS)) {
    idx <- idx + 1
    set.seed(SEED_BASE + (rep_id - 1) * 31 + ci)

    clean <- simulated_good_responses(nConstructs = nF, nItems = rep(ipf, nF), n = N_FIXED)
    inj <- inject_careless_v3(clean, pct_careless = PCT_CARELESS,
                              careless_levels = CARELESS_LEVELS)
    data_mat <- inj$data_corrupted

    # Standard RR with VP
    rr <- ReReReRe(data_mat, corProp = CORPROP, iterations = ITER,
                   align_signs = TRUE, mode = "auto",
                   auto_z = TRUE, variance_penalty = TRUE)
    z_rr <- rr$z_score

    # Iterative EFA with different trim_pct values (via initial_z_threshold proxy)
    # ReReReRe_F_iterative uses initial_z_threshold as the trim criterion.
    # Higher threshold = more aggressive trim. Map: trim_pct≈0.1→z=0.5, 0.2→1.0, 0.3→1.5
    iter_thrs <- c(0.5, 1.0, 1.5)  # rough mapping
    iter_z_by_trim <- list()
    for (it_i in seq_along(iter_thrs)) {
      iter_z_by_trim[[paste0("iter_trim", it_i)]] <- tryCatch({
        rr_it <- ReReReRe_F_iterative(data_mat, iterations = ITER,
                                      align_signs = TRUE,
                                      initial_z_threshold = iter_thrs[it_i])
        rr_it$z_score
      }, error = function(e) rep(NA_real_, nrow(data_mat)))
    }

    # Other detectors
    irv <- compute_irv(data_mat)
    ls  <- compute_longstring(data_mat)
    d2  <- compute_d2(data_mat)
    pt  <- compute_person_total(data_mat)

    scores <- data.frame(
      nF = nF, ipf = ipf, total_items = tot, rep = rep_id,
      respondent = seq_len(nrow(data_mat)),
      pattern = inj$labels$pattern,
      corruption = round(inj$labels$careless_pct, 2),
      z_rr             = -z_rr,
      z_rr_iter_t1     = -iter_z_by_trim[["iter_trim1"]],   # trim ~10%
      z_rr_iter_t2     = -iter_z_by_trim[["iter_trim2"]],   # trim ~20%
      z_rr_iter_t3     = -iter_z_by_trim[["iter_trim3"]],   # trim ~30%
      irv              = -irv,
      longstring       = ls,
      d2               = d2,
      person_tot       = -pt
    )
    rows[[length(rows) + 1]] <- scores

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    eta <- if (idx > 1) (elapsed / idx) * (TOTAL - idx) else NA
    cat(sprintf("[%d/%d] nF=%d ipf=%d rep=%d | %.1f min ETA %.0f\n",
                idx, TOTAL, nF, ipf, rep_id, elapsed, eta))
  }
}

df <- do.call(rbind, rows)

# Scale features
feature_cols <- c("z_rr", "z_rr_iter_t1", "z_rr_iter_t2", "z_rr_iter_t3",
                  "irv", "longstring", "d2", "person_tot")
for (col in feature_cols) {
  x <- df[[col]]
  mu <- mean(x, na.rm = TRUE); sd_ <- sd(x, na.rm = TRUE)
  if (!is.na(sd_) && sd_ > 0) df[[col]] <- (x - mu) / sd_
  df[[col]][is.na(df[[col]])] <- 0
}

write.csv(df, "sim_robust_scores.csv", row.names = FALSE)
cat(sprintf("\nSaved sim_robust_scores.csv (%d rows)\n", nrow(df)))

# ============================================================
# Analysis 1: Best iter_EFA trim_pct
# ============================================================
set.seed(42)
folds <- sample(seq_len(N_FOLDS), nrow(df), replace = TRUE)

cv_metrics <- function(df, feats, gt_col, folds) {
  preds <- numeric(nrow(df))
  for (f in unique(folds)) {
    tr <- folds != f; te <- folds == f
    fs <- paste0(gt_col, " ~ ", paste(feats, collapse = " + "))
    mod <- tryCatch(glm(as.formula(fs), data = df[tr, ], family = binomial),
                    error = function(e) NULL)
    if (is.null(mod)) { preds[te] <- 0.5; next }
    preds[te] <- predict(mod, newdata = df[te, ], type = "response")
  }
  # FPR=5% operating point
  thr <- quantile(preds[df$pattern == "clean"], 0.95)
  flag <- preds > thr
  truth <- df[[gt_col]]
  TP <- sum(flag & truth); TN <- sum(!flag & !truth)
  FP <- sum(flag & !truth); FN <- sum(!flag & truth)
  den <- sqrt(as.numeric(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  list(preds = preds, thr = thr, flag = flag,
       mcc  = if (den == 0) 0 else (TP*TN - FP*FN) / den,
       sens = TP/max(1,TP+FN), fpr = FP/max(1,FP+TN))
}

cat("\n=== ANALYSIS 1: iter_EFA trim_pct tuning ===\n")
df$gt50 <- (df$pattern != "clean") & (df$corruption > 0.50)
trim_tune <- data.frame()
for (iter_col in c("z_rr_iter_t1", "z_rr_iter_t2", "z_rr_iter_t3")) {
  feats <- c("z_rr", iter_col, "irv", "longstring", "d2", "person_tot")
  r <- cv_metrics(df, feats, "gt50", folds)
  trim_tune <- rbind(trim_tune, data.frame(
    trim_config = iter_col, sens = round(r$sens, 3), fpr = round(r$fpr, 3),
    mcc = round(r$mcc, 3)))
}
print(trim_tune, row.names = FALSE)

best_iter <- trim_tune$trim_config[which.max(trim_tune$mcc)]
cat(sprintf("\nBest iter_EFA config: %s\n", best_iter))

# ============================================================
# Analysis 2: Multi-threshold MCC for best ensemble
# ============================================================
cat("\n=== ANALYSIS 2: MCC across GT thresholds ===\n")
ensembles <- list(
  z_rr_only    = c("z_rr"),
  triad        = c("z_rr", "irv", "d2"),
  triad_plus_iter = c("z_rr", "irv", "d2", best_iter),
  all_5_no_iter = c("z_rr", "irv", "longstring", "d2", "person_tot"),
  all_6_with_iter = c("z_rr", best_iter, "irv", "longstring", "d2", "person_tot")
)

multi_thr <- data.frame()
for (T in GT_THRESHOLDS) {
  gt_col <- paste0("gt_", as.integer(T*100))
  df[[gt_col]] <- (df$pattern != "clean") & (df$corruption > T)
  n_pos <- sum(df[[gt_col]])
  if (n_pos < 100) next
  for (nm in names(ensembles)) {
    r <- cv_metrics(df, ensembles[[nm]], gt_col, folds)
    multi_thr <- rbind(multi_thr, data.frame(
      gt_threshold = T, ensemble = nm, n_positive = n_pos,
      sens = round(r$sens, 3), fpr = round(r$fpr, 3),
      mcc  = round(r$mcc, 3)
    ))
  }
}
print(multi_thr, row.names = FALSE)
write.csv(multi_thr, "robust_multi_threshold.csv", row.names = FALSE)

# ============================================================
# Analysis 3: Sensitivity per corruption bucket (focus on 60-80%)
# ============================================================
cat("\n=== ANALYSIS 3: Detection per corruption bucket (best ensemble) ===\n")
# Train best ensemble at GT>50%, evaluate per-bucket
df$gt_best <- (df$pattern != "clean") & (df$corruption > 0.50)
feats_best <- ensembles$all_6_with_iter
r_best <- cv_metrics(df, feats_best, "gt_best", folds)
df$flag_best <- r_best$flag

# Also the triad for comparison
feats_triad <- ensembles$triad
r_triad <- cv_metrics(df, feats_triad, "gt_best", folds)
df$flag_triad <- r_triad$flag

# z_rr alone
feats_zrr <- c("z_rr")
r_zrr <- cv_metrics(df, feats_zrr, "gt_best", folds)
df$flag_zrr <- r_zrr$flag

per_bucket <- df %>%
  filter(pattern != "clean") %>%
  mutate(bin = round(corruption * 10) / 10) %>%
  group_by(pattern, bin) %>%
  summarise(
    n = n(),
    rate_zrr   = round(mean(flag_zrr),   3),
    rate_triad = round(mean(flag_triad), 3),
    rate_best  = round(mean(flag_best),  3),
    .groups = "drop"
  )
cat("Per-pattern detection:\n")
print(as.data.frame(per_bucket %>% filter(bin >= 0.5)), row.names = FALSE)
write.csv(per_bucket, "robust_per_bucket_detection.csv", row.names = FALSE)

# Focus: 60-80% average across patterns
focus_60_80 <- per_bucket %>%
  filter(bin %in% c(0.6, 0.7, 0.8)) %>%
  group_by(bin) %>%
  summarise(
    mean_rate_zrr   = round(mean(rate_zrr),   3),
    mean_rate_triad = round(mean(rate_triad), 3),
    mean_rate_best  = round(mean(rate_best),  3),
    .groups = "drop"
  )
cat("\n=== Focus: partial-corruption range (60-80%, averaged across patterns) ===\n")
print(as.data.frame(focus_60_80), row.names = FALSE)

cat(sprintf("\nRuntime: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("=== SIMULATION DONE, NOW PLOTTING ===\n")

# ============================================================
# Plots
# ============================================================
pattern_colors <- c(
  random = "#e41a1c", longstring = "#377eb8", mixed = "#4daf4a",
  pure_straight = "#984ea3", acquiescent = "#ff7f00", fatigue = "#a65628"
)

# Plot 1: Sensitivity per corruption bucket, full ensemble vs z_rr alone, faceted by pattern
plot_df <- per_bucket %>%
  pivot_longer(cols = starts_with("rate_"),
               names_to = "model", values_to = "rate") %>%
  mutate(model = recode(model,
                        "rate_zrr"   = "z_RR alone",
                        "rate_triad" = "Triad (z_RR + IRV + D²)",
                        "rate_best"  = "Best ensemble (6 detectors)"),
         model = factor(model, levels = c("z_RR alone",
                                           "Triad (z_RR + IRV + D²)",
                                           "Best ensemble (6 detectors)")))

p1 <- ggplot(plot_df, aes(x = bin, y = rate, color = model, group = model)) +
  annotate("rect", xmin = 0.55, xmax = 0.85, ymin = 0, ymax = 1,
           fill = "#fff3cd", alpha = 0.5) +
  annotate("text", x = 0.70, y = 1.04,
           label = "PARTIAL CARELESSNESS (60-80%) — focus of re-examination",
           size = 3, color = "#664d03", fontface = "bold") +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  scale_color_manual(values = c("z_RR alone" = "#377eb8",
                                 "Triad (z_RR + IRV + D²)" = "#4daf4a",
                                 "Best ensemble (6 detectors)" = "#e41a1c"),
                     name = "Model") +
  scale_x_continuous(breaks = seq(0.1, 1.0, 0.1),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0, 1.10),
                     labels = scales::percent_format(accuracy = 1)) +
  facet_wrap(~ pattern, ncol = 3) +
  labs(
    title = "Sensitivity per corruption bucket — focus on partial carelessness (60-80%)",
    subtitle = sprintf("FPR=5%% operating point, trained at GT>50%%. Best ensemble uses %s.",
                       best_iter),
    x = "Corruption level", y = "Detection rate"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "#EEEEEE", color = NA),
        strip.text = element_text(face = "bold"))

ggsave("plot_robust_sens_per_pattern.png", p1,
       width = 13, height = 8, dpi = 150, bg = "white")
cat("Saved: plot_robust_sens_per_pattern.png\n")

# Plot 2: MCC vs GT threshold
p2 <- ggplot(multi_thr, aes(x = gt_threshold, y = mcc,
                             color = ensemble, group = ensemble)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.2f", mcc)),
            vjust = -1.0, size = 2.6, show.legend = FALSE) +
  scale_color_brewer(palette = "Set1", name = "Ensemble") +
  scale_x_continuous(breaks = GT_THRESHOLDS,
                     labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.1)) +
  labs(
    title = "MCC vs GT threshold — comparison of ensembles",
    subtitle = "FPR=5% operating point, 5-fold CV",
    x = "GT threshold (corruption > T → careless)",
    y = "MCC"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"))

ggsave("plot_robust_mcc_vs_threshold.png", p2,
       width = 12, height = 7, dpi = 150, bg = "white")
cat("Saved: plot_robust_mcc_vs_threshold.png\n")

# Plot 3: Focus plot — 60-80% detection summary
focus_plot_df <- per_bucket %>%
  filter(bin >= 0.5) %>%
  pivot_longer(cols = starts_with("rate_"),
               names_to = "model", values_to = "rate") %>%
  mutate(model = recode(model,
                        "rate_zrr"   = "z_RR alone",
                        "rate_triad" = "Triad",
                        "rate_best"  = "Best ensemble"))

p3 <- ggplot(focus_plot_df %>%
               group_by(bin, model) %>%
               summarise(rate = mean(rate), .groups = "drop"),
              aes(x = bin, y = rate, fill = model)) +
  geom_col(position = "dodge", width = 0.07) +
  geom_text(aes(label = sprintf("%.2f", rate)),
            position = position_dodge(width = 0.07),
            vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = c("z_RR alone" = "#377eb8",
                                "Triad" = "#4daf4a",
                                "Best ensemble" = "#e41a1c")) +
  scale_x_continuous(breaks = seq(0.5, 1.0, 0.1),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0, 1.05),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Detection rate @ corruption 50-100% — averaged across patterns",
    subtitle = "Best ensemble (iter_EFA + IRV + D² + LongString + PersonTotal + z_RR) at FPR=5%",
    x = "Corruption level", y = "Detection rate", fill = "Model"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"))

ggsave("plot_robust_focus_50_100.png", p3,
       width = 11, height = 7, dpi = 150, bg = "white")
cat("Saved: plot_robust_focus_50_100.png\n")

cat("\n=== DONE ===\n")
