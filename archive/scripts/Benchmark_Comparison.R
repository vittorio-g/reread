#### Benchmark Comparison: ReReReRe vs Standard Methods ####
# Compares ReReReRe against four standard careless detection methods:
#   1. Mahalanobis distance (multivariate outlier)
#   2. Person-total correlation (agreement with sample)
#   3. LongString index (max consecutive identical responses)
#   4. IRV — intra-individual response variability (SD per person)
#
# All methods are applied to the SAME simulated data per condition.
# Each method gets its own optimal threshold (post-hoc) to maximize MCC.
# AUC is also computed (threshold-free, via Mann-Whitney U).
#
# Design: representative subset of the multiverse grid.
# Seeds are set per condition for reproducibility.
#
# === REVISION 2026-03-10g ===

rm(list = ls())

library(dplyr)
library(tidyr)
library(lavaan)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("Synthetic_Good_Responses_2.R")
source("Careless_machine_2.R")
source("ReReReRe.R")

# ============================================================
# CONFIGURATION
# ============================================================

GRID <- expand.grid(
  nFactors      = c(4, 8, 10, 15, 20, 25),
  n_respondents = c(50, 100, 300, 1000),
  pct_careless  = c(0.05, 0.10, 0.25),
  stringsAsFactors = FALSE
)

# ReReReRe fixed parameters
CORPROP        <- 0.05          # best from multiverse
ITERATIONS     <- 100
MIN_PAIRS      <- 15
ITEMS_CYCLE    <- c(10, 6, 3)
CARELESS_TYPES <- c(random = 1/3, longstring = 1/3, mixed = 1/3)
CARELESS_LEVELS <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5)
EVAL_THRESHOLD <- 0.01
SEED_BASE      <- 2026

# z-threshold values for ReReReRe (post-hoc)
Z_THRESHOLDS <- c(0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5)

# ============================================================
# BENCHMARK METHOD IMPLEMENTATIONS
# ============================================================

# 1. Mahalanobis distance — flags multivariate outliers (HIGH = careless)
compute_mahad <- function(data) {
  x <- as.matrix(data)
  center <- colMeans(x)
  cv <- cov(x)
  tryCatch(
    mahalanobis(x, center, cv),
    error = function(e) {
      # Regularize singular/near-singular covariance (n < p cases)
      tryCatch(
        mahalanobis(x, center, cv + diag(1e-3, ncol(x))),
        error = function(e2) rep(NA_real_, nrow(x))
      )
    }
  )
}

# 2. Person-total correlation — flags low agreement with sample mean (LOW = careless)
compute_ptc <- function(data) {
  x <- as.matrix(data)
  col_means <- colMeans(x)
  apply(x, 1, function(row) cor(row, col_means))
}

# 3. LongString — max consecutive identical responses (HIGH = careless)
#    NOTE: Our post-fix longstrings are NON-contiguous (scattered positions),
#    so this index will underperform on our simulated longstring pattern.
#    This is a real-world limitation: if item order is scrambled,
#    classic LongString fails. This is fair to report.
compute_longstring <- function(data) {
  x <- as.matrix(data)
  apply(x, 1, function(row) {
    r <- rle(as.integer(row))
    max(r$lengths)
  })
}

# 4. IRV — intra-individual response variability (LOW = careless)
compute_irv <- function(data) {
  apply(as.matrix(data), 1, sd)
}

# ============================================================
# EVALUATION HELPERS
# ============================================================

# MCC
calc_mcc <- function(tp, tn, fp, fn) {
  denom <- sqrt(as.numeric(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  if (denom == 0) return(0)
  (tp * tn - fp * fn) / denom
}

# AUC via Mann-Whitney U
calc_auc <- function(scores, is_careless, higher_is_careless) {
  if (any(is.na(scores))) return(NA_real_)
  if (!higher_is_careless) scores <- -scores
  n_pos <- sum(is_careless)
  n_neg <- sum(!is_careless)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  ranks <- rank(scores, ties.method = "average")
  (sum(ranks[is_careless]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

# Find best threshold to maximize MCC (quantile-based search)
find_best_threshold <- function(scores, is_careless, higher_is_careless,
                                n_candidates = 200) {
  if (any(is.na(scores))) return(list(mcc = NA, threshold = NA,
                                       sensitivity = NA, specificity = NA))

  candidates <- quantile(scores, probs = seq(0.005, 0.995, length.out = n_candidates),
                          na.rm = TRUE)
  candidates <- unique(candidates)

  best_mcc    <- -Inf
  best_thresh <- NA
  best_sens   <- NA
  best_spec   <- NA

  for (thresh in candidates) {
    flagged <- if (higher_is_careless) scores >= thresh else scores <= thresh

    tp <- sum(flagged & is_careless)
    tn <- sum(!flagged & !is_careless)
    fp <- sum(flagged & !is_careless)
    fn <- sum(!flagged & is_careless)

    mcc <- calc_mcc(tp, tn, fp, fn)
    if (mcc > best_mcc) {
      best_mcc    <- mcc
      best_thresh <- thresh
      n_pos <- tp + fn
      best_sens <- if (n_pos > 0) tp / n_pos else NA
      best_spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA
    }
  }

  list(mcc = best_mcc, threshold = best_thresh,
       sensitivity = best_sens, specificity = best_spec)
}

# Find best z_threshold for ReReReRe
find_best_zt <- function(z_scores, is_careless) {
  best_mcc  <- -Inf
  best_zt   <- NA
  best_sens <- NA
  best_spec <- NA

  for (zt in Z_THRESHOLDS) {
    flagged <- z_scores <= zt

    tp <- sum(flagged & is_careless)
    tn <- sum(!flagged & !is_careless)
    fp <- sum(flagged & !is_careless)
    fn <- sum(!flagged & is_careless)

    mcc <- calc_mcc(tp, tn, fp, fn)
    if (mcc > best_mcc) {
      best_mcc  <- mcc
      best_zt   <- zt
      n_pos <- tp + fn
      best_sens <- if (n_pos > 0) tp / n_pos else NA
      best_spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA
    }
  }

  list(mcc = best_mcc, threshold = best_zt,
       sensitivity = best_sens, specificity = best_spec)
}

# ============================================================
# MAIN LOOP
# ============================================================

cat("=== BENCHMARK COMPARISON ===\n")
cat(sprintf("Grid: %d conditions | Methods: ReReReRe, Mahalanobis, PersonTotal, LongString, IRV\n",
            nrow(GRID)))
cat(sprintf("ReReReRe: corProp=%.2f, iterations=%d, min_pairs=%d\n\n", CORPROP, ITERATIONS, MIN_PAIRS))

results <- vector("list", nrow(GRID) * 5)
res_idx <- 0
t_start <- Sys.time()

for (i in seq_len(nrow(GRID))) {

  g <- GRID[i, ]
  nF  <- g$nFactors
  n   <- g$n_respondents
  pct <- g$pct_careless

  nItems_vec  <- rep_len(ITEMS_CYCLE, nF)
  total_items <- sum(nItems_vec)
  n_lt_p      <- n < total_items

  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  if (i > 1) {
    rate <- elapsed / (i - 1)
    eta <- sprintf("%.0fmin left", rate * (nrow(GRID) - i + 1))
  } else {
    eta <- "..."
  }

  cat(sprintf("[%d/%d] nF=%d n=%d pct=%.2f (items=%d, n<p=%s) | %.1fmin | %s  →  ",
              i, nrow(GRID), nF, n, pct, total_items, n_lt_p, elapsed, eta))

  # --- Reproducible seed ---
  set.seed(SEED_BASE + i)

  # --- Generate data ---
  clean_data <- tryCatch(
    simulated_good_responses(nConstructs = nF, nItems = nItems_vec, n = n),
    error = function(e) { cat("GEN FAIL\n"); return(NULL) }
  )
  if (is.null(clean_data)) next

  injection <- tryCatch(
    inject_careless(data = clean_data, pct_careless = pct,
                    pct_types = CARELESS_TYPES, careless_levels = CARELESS_LEVELS),
    error = function(e) { cat("INJ FAIL\n"); return(NULL) }
  )
  if (is.null(injection)) next

  data_mat   <- injection$data_corrupted
  labels     <- injection$labels
  is_careless <- labels$careless_pct >= EVAL_THRESHOLD

  # --- Method 1: ReReReRe ---
  rr <- tryCatch(
    ReReReRe(data = data_mat, corProp = CORPROP, cutOff = 0,
             iterations = ITERATIONS, min_pairs = MIN_PAIRS, progress = FALSE),
    error = function(e) NULL
  )
  if (!is.null(rr)) {
    rr_best <- find_best_zt(rr$z_score, is_careless)
    rr_auc  <- calc_auc(rr$z_score, is_careless, higher_is_careless = FALSE)
  } else {
    rr_best <- list(mcc = NA, threshold = NA, sensitivity = NA, specificity = NA)
    rr_auc  <- NA
  }

  # --- Method 2: Mahalanobis ---
  mahad <- compute_mahad(data_mat)
  mahad_best <- find_best_threshold(mahad, is_careless, higher_is_careless = TRUE)
  mahad_auc  <- calc_auc(mahad, is_careless, higher_is_careless = TRUE)

  # --- Method 3: Person-total correlation ---
  ptc <- compute_ptc(data_mat)
  ptc_best <- find_best_threshold(ptc, is_careless, higher_is_careless = FALSE)
  ptc_auc  <- calc_auc(ptc, is_careless, higher_is_careless = FALSE)

  # --- Method 4: LongString ---
  ls <- compute_longstring(data_mat)
  ls_best <- find_best_threshold(ls, is_careless, higher_is_careless = TRUE)
  ls_auc  <- calc_auc(ls, is_careless, higher_is_careless = TRUE)

  # --- Method 5: IRV ---
  irv <- compute_irv(data_mat)
  irv_best <- find_best_threshold(irv, is_careless, higher_is_careless = FALSE)
  irv_auc  <- calc_auc(irv, is_careless, higher_is_careless = FALSE)

  # --- Collect ---
  for (mi in list(
    list(name = "ReReReRe",    b = rr_best,    auc = rr_auc),
    list(name = "Mahalanobis", b = mahad_best,  auc = mahad_auc),
    list(name = "PersonTotal", b = ptc_best,    auc = ptc_auc),
    list(name = "LongString",  b = ls_best,     auc = ls_auc),
    list(name = "IRV",         b = irv_best,    auc = irv_auc)
  )) {
    res_idx <- res_idx + 1
    results[[res_idx]] <- data.frame(
      nFactors      = nF,
      n_respondents = n,
      total_items   = total_items,
      n_lt_p        = n_lt_p,
      pct_careless  = pct,
      method        = mi$name,
      best_mcc      = mi$b$mcc,
      best_threshold = mi$b$threshold,
      sensitivity   = mi$b$sensitivity,
      specificity   = mi$b$specificity,
      auc           = mi$auc,
      stringsAsFactors = FALSE
    )
  }

  cat(sprintf("RR=%.3f  Mah=%.3f  PTC=%.3f  LS=%.3f  IRV=%.3f\n",
              rr_best$mcc, mahad_best$mcc, ptc_best$mcc, ls_best$mcc, irv_best$mcc))
}

elapsed_total <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat(sprintf("\nDone in %.1f minutes.\n", elapsed_total))

# ============================================================
# COMBINE AND SAVE
# ============================================================

df <- bind_rows(results[seq_len(res_idx)])
write.csv(df, "benchmark_results.csv", row.names = FALSE)
cat("Saved: benchmark_results.csv\n\n")

# ============================================================
# ANALYSIS
# ============================================================

cat("=====================================================\n")
cat("  1. OVERALL METHOD COMPARISON\n")
cat("=====================================================\n")
overall <- df %>%
  group_by(method) %>%
  summarise(
    mean_mcc  = mean(best_mcc, na.rm = TRUE),
    median_mcc = median(best_mcc, na.rm = TRUE),
    max_mcc   = max(best_mcc, na.rm = TRUE),
    mean_auc  = mean(auc, na.rm = TRUE),
    mean_sens = mean(sensitivity, na.rm = TRUE),
    mean_spec = mean(specificity, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  arrange(desc(mean_mcc))
print(as.data.frame(overall), digits = 3)

cat("\n=====================================================\n")
cat("  2. MEAN BEST-MCC BY METHOD x nFactors\n")
cat("=====================================================\n")
mcc_nf <- df %>%
  group_by(nFactors, method) %>%
  summarise(mean_mcc = mean(best_mcc, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = mean_mcc)
print(as.data.frame(mcc_nf), digits = 3)

cat("\n=====================================================\n")
cat("  3. MEAN AUC BY METHOD x nFactors\n")
cat("=====================================================\n")
auc_nf <- df %>%
  group_by(nFactors, method) %>%
  summarise(mean_auc = mean(auc, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = mean_auc)
print(as.data.frame(auc_nf), digits = 3)

cat("\n=====================================================\n")
cat("  4. MEAN BEST-MCC BY METHOD x n_respondents\n")
cat("=====================================================\n")
mcc_n <- df %>%
  group_by(n_respondents, method) %>%
  summarise(mean_mcc = mean(best_mcc, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = mean_mcc)
print(as.data.frame(mcc_n), digits = 3)

cat("\n=====================================================\n")
cat("  5. MEAN BEST-MCC BY METHOD x pct_careless\n")
cat("=====================================================\n")
mcc_pct <- df %>%
  group_by(pct_careless, method) %>%
  summarise(mean_mcc = mean(best_mcc, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = mean_mcc)
print(as.data.frame(mcc_pct), digits = 3)

cat("\n=====================================================\n")
cat("  6. HEAD-TO-HEAD: ReReReRe vs each benchmark\n")
cat("=====================================================\n")
rr_df <- df %>% filter(method == "ReReReRe") %>%
  select(nFactors, n_respondents, pct_careless, rr_mcc = best_mcc, rr_auc = auc)

for (bench in c("Mahalanobis", "PersonTotal", "LongString", "IRV")) {
  b_df <- df %>% filter(method == bench) %>%
    select(nFactors, n_respondents, pct_careless, b_mcc = best_mcc, b_auc = auc)

  comp <- inner_join(rr_df, b_df, by = c("nFactors", "n_respondents", "pct_careless"))

  wins   <- sum(comp$rr_mcc > comp$b_mcc, na.rm = TRUE)
  ties   <- sum(comp$rr_mcc == comp$b_mcc, na.rm = TRUE)
  losses <- sum(comp$rr_mcc < comp$b_mcc, na.rm = TRUE)
  na_ct  <- sum(is.na(comp$rr_mcc) | is.na(comp$b_mcc))
  dmcc   <- mean(comp$rr_mcc - comp$b_mcc, na.rm = TRUE)
  dauc   <- mean(comp$rr_auc - comp$b_auc, na.rm = TRUE)

  cat(sprintf("  vs %-12s: W/T/L = %d/%d/%d (NA: %d) | mean dMCC = %+.3f | mean dAUC = %+.3f\n",
              bench, wins, ties, losses, na_ct, dmcc, dauc))
}

cat("\n=====================================================\n")
cat("  7. HEAD-TO-HEAD by nFactors (ReReReRe minus best benchmark)\n")
cat("=====================================================\n")
for (nf in sort(unique(df$nFactors))) {
  sub_rr <- df %>% filter(method == "ReReReRe", nFactors == nf)
  sub_bench <- df %>% filter(method != "ReReReRe", nFactors == nf) %>%
    group_by(n_respondents, pct_careless) %>%
    summarise(best_bench_mcc = max(best_mcc, na.rm = TRUE),
              best_bench_method = method[which.max(best_mcc)],
              .groups = "drop")

  comp <- inner_join(
    sub_rr %>% select(n_respondents, pct_careless, rr_mcc = best_mcc),
    sub_bench,
    by = c("n_respondents", "pct_careless")
  )

  delta <- comp$rr_mcc - comp$best_bench_mcc
  wins  <- sum(delta > 0, na.rm = TRUE)
  total <- sum(!is.na(delta))

  # What's the best benchmark for this nF?
  best_b <- df %>% filter(method != "ReReReRe", nFactors == nf) %>%
    group_by(method) %>%
    summarise(m = mean(best_mcc, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(m))

  cat(sprintf("  nF=%d: ReReReRe wins %d/%d | mean delta = %+.3f | best benchmark: %s (%.3f)\n",
              nf, wins, total, mean(delta, na.rm = TRUE),
              best_b$method[1], best_b$m[1]))
}

cat("\n=====================================================\n")
cat("  8. RECOMMENDED ZONE (nF>=10, n>=100)\n")
cat("=====================================================\n")
rec <- df %>% filter(nFactors >= 10, n_respondents >= 100)
rec_overall <- rec %>%
  group_by(method) %>%
  summarise(
    mean_mcc = mean(best_mcc, na.rm = TRUE),
    mean_auc = mean(auc, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  arrange(desc(mean_mcc))
print(as.data.frame(rec_overall), digits = 3)

# Head-to-head in recommended zone
rr_rec <- rec %>% filter(method == "ReReReRe") %>%
  select(nFactors, n_respondents, pct_careless, rr_mcc = best_mcc)

for (bench in c("Mahalanobis", "PersonTotal", "LongString", "IRV")) {
  b_rec <- rec %>% filter(method == bench) %>%
    select(nFactors, n_respondents, pct_careless, b_mcc = best_mcc)

  comp <- inner_join(rr_rec, b_rec, by = c("nFactors", "n_respondents", "pct_careless"))
  wins   <- sum(comp$rr_mcc > comp$b_mcc, na.rm = TRUE)
  total  <- sum(!is.na(comp$rr_mcc) & !is.na(comp$b_mcc))
  dmcc   <- mean(comp$rr_mcc - comp$b_mcc, na.rm = TRUE)

  cat(sprintf("    vs %-12s: wins %d/%d | mean dMCC = %+.3f\n",
              bench, wins, total, dmcc))
}

# ============================================================
# PLOT: Method comparison by nFactors
# ============================================================

cat("\nGenerating plots...\n")

# Plot 1: MCC by nFactors per method
plot_data <- df %>%
  group_by(nFactors, method) %>%
  summarise(mean_mcc = mean(best_mcc, na.rm = TRUE), .groups = "drop")

methods <- c("ReReReRe", "Mahalanobis", "PersonTotal", "LongString", "IRV")
cols <- c("ReReReRe" = "#E41A1C", "Mahalanobis" = "#377EB8",
          "PersonTotal" = "#4DAF4A", "LongString" = "#984EA3", "IRV" = "#FF7F00")
pchs <- c(16, 17, 15, 18, 8)
names(pchs) <- methods

nf_vals <- sort(unique(plot_data$nFactors))

png("benchmark_mcc_by_nFactors.png", width = 900, height = 600)
par(mar = c(5, 5, 3, 2))
plot(NULL, xlim = range(nf_vals), ylim = c(-0.05, 0.65),
     xlab = "Number of Factors", ylab = "Mean MCC (best threshold)",
     main = "ReReReRe vs Standard Methods: Mean MCC by Questionnaire Size",
     cex.lab = 1.3, cex.main = 1.3, xaxt = "n")
axis(1, at = nf_vals)
abline(h = seq(0, 0.6, 0.1), col = "grey90", lty = 2)
abline(h = 0, col = "grey50")

for (j in seq_along(methods)) {
  m <- methods[j]
  pd <- plot_data %>% filter(method == m)
  lines(pd$nFactors, pd$mean_mcc, col = cols[m], lwd = 2.5, type = "o",
        pch = pchs[m], cex = 1.5)
}
legend("topleft", legend = methods, col = cols[methods], lwd = 2.5,
       pch = pchs, pt.cex = 1.5, cex = 1.1, bg = "white")
dev.off()
cat("Saved: benchmark_mcc_by_nFactors.png\n")

# Plot 2: AUC by nFactors per method
auc_plot <- df %>%
  group_by(nFactors, method) %>%
  summarise(mean_auc = mean(auc, na.rm = TRUE), .groups = "drop")

png("benchmark_auc_by_nFactors.png", width = 900, height = 600)
par(mar = c(5, 5, 3, 2))
plot(NULL, xlim = range(nf_vals), ylim = c(0.4, 1.0),
     xlab = "Number of Factors", ylab = "Mean AUC",
     main = "ReReReRe vs Standard Methods: Mean AUC by Questionnaire Size",
     cex.lab = 1.3, cex.main = 1.3, xaxt = "n")
axis(1, at = nf_vals)
abline(h = seq(0.4, 1, 0.1), col = "grey90", lty = 2)
abline(h = 0.5, col = "grey50")

for (j in seq_along(methods)) {
  m <- methods[j]
  pd <- auc_plot %>% filter(method == m)
  lines(pd$nFactors, pd$mean_auc, col = cols[m], lwd = 2.5, type = "o",
        pch = pchs[m], cex = 1.5)
}
legend("bottomright", legend = methods, col = cols[methods], lwd = 2.5,
       pch = pchs, pt.cex = 1.5, cex = 1.1, bg = "white")
dev.off()
cat("Saved: benchmark_auc_by_nFactors.png\n")

# Plot 3: MCC by n_respondents per method
n_plot <- df %>%
  group_by(n_respondents, method) %>%
  summarise(mean_mcc = mean(best_mcc, na.rm = TRUE), .groups = "drop")

n_vals <- sort(unique(n_plot$n_respondents))

png("benchmark_mcc_by_n.png", width = 900, height = 600)
par(mar = c(5, 5, 3, 2))
plot(NULL, xlim = c(1, length(n_vals)), ylim = c(-0.05, 0.55),
     xlab = "Sample Size", ylab = "Mean MCC (best threshold)",
     main = "ReReReRe vs Standard Methods: Mean MCC by Sample Size",
     cex.lab = 1.3, cex.main = 1.3, xaxt = "n")
axis(1, at = seq_along(n_vals), labels = n_vals)
abline(h = seq(0, 0.5, 0.1), col = "grey90", lty = 2)
abline(h = 0, col = "grey50")

for (j in seq_along(methods)) {
  m <- methods[j]
  pd <- n_plot %>% filter(method == m) %>% arrange(n_respondents)
  x_pos <- match(pd$n_respondents, n_vals)
  lines(x_pos, pd$mean_mcc, col = cols[m], lwd = 2.5, type = "o",
        pch = pchs[m], cex = 1.5)
}
legend("topleft", legend = methods, col = cols[methods], lwd = 2.5,
       pch = pchs, pt.cex = 1.5, cex = 1.1, bg = "white")
dev.off()
cat("Saved: benchmark_mcc_by_n.png\n")

# Plot 4: Head-to-head delta (ReReReRe minus best benchmark) by nFactors x n
h2h <- df %>%
  filter(method != "ReReReRe") %>%
  group_by(nFactors, n_respondents, pct_careless) %>%
  summarise(best_bench = max(best_mcc, na.rm = TRUE), .groups = "drop") %>%
  inner_join(
    df %>% filter(method == "ReReReRe") %>%
      select(nFactors, n_respondents, pct_careless, rr_mcc = best_mcc),
    by = c("nFactors", "n_respondents", "pct_careless")
  ) %>%
  mutate(delta = rr_mcc - best_bench)

h2h_summary <- h2h %>%
  group_by(nFactors, n_respondents) %>%
  summarise(mean_delta = mean(delta, na.rm = TRUE), .groups = "drop")

png("benchmark_delta_heatmap.png", width = 800, height = 600)
par(mar = c(5, 5, 4, 7))
delta_mat <- h2h_summary %>%
  pivot_wider(names_from = n_respondents, values_from = mean_delta) %>%
  arrange(nFactors)
mat <- as.matrix(delta_mat[,-1])
rownames(mat) <- paste0("nF=", delta_mat$nFactors)

# Color: red = ReReReRe loses, white = tie, blue = ReReReRe wins
n_colors <- 21
pal <- colorRampPalette(c("#E41A1C", "white", "#377EB8"))(n_colors)
max_abs <- max(abs(mat), na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = n_colors + 1)

image(t(mat[nrow(mat):1, ]), col = pal, breaks = breaks,
      axes = FALSE,
      main = "ReReReRe minus Best Benchmark (Mean MCC)",
      xlab = "Sample Size", ylab = "nFactors")
axis(1, at = seq(0, 1, length.out = ncol(mat)),
     labels = colnames(mat))
axis(2, at = seq(0, 1, length.out = nrow(mat)),
     labels = rev(rownames(mat)), las = 1)
# Add text values
for (r in 1:nrow(mat)) {
  for (c in 1:ncol(mat)) {
    x <- (c - 1) / (ncol(mat) - 1)
    y <- 1 - (r - 1) / (nrow(mat) - 1)
    text(x, y, sprintf("%+.3f", mat[r, c]), cex = 1.0, font = 2)
  }
}
mtext("Blue = ReReReRe wins | Red = benchmark wins", side = 1, line = 3.5, cex = 0.9)
dev.off()
cat("Saved: benchmark_delta_heatmap.png\n")

cat("\n=== DONE ===\n")
