# Simulation_Ensemble.R
# Greedy forward selection: start with z_RR, add one detector at a time,
# report MCC improvement. Uses 5-fold CV logistic regression as combiner.
#
# Detectors tested (5 total, plus baseline z_RR):
#   1. z_RR      : ReReReRe z-score with auto-z + variance_penalty
#   2. IRV       : within-respondent SD (low = careless)
#   3. LongString: max consecutive identical run (high = careless)
#   4. D2        : Mahalanobis distance from sample centroid (high = careless)
#   5. PersonTot : correlation respondent vs item-means (low = careless)
#
# Ground truth: careless_pct > 0.50 (using v3 injector with unified semantics).

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(tidyr)
})
source("Synthetic_Good_Responses_2.R")
source("ReReReRe.R")
source("Careless_machine_v3.R")

# --- Config ---
NF_LEVELS  <- c(8, 16, 30)
IPF_LEVELS <- c(6, 10)
REPS <- 4
N_FIXED <- 500
PCT_CARELESS <- 0.40  # higher rate for more positives per cell
ITER <- 100
CORPROP <- 0.03
SEED_BASE <- 2026
CARELESS_LEVELS <- seq(0.1, 1.0, 0.1)
GT_CUTOFF <- 0.50
N_FOLDS <- 5

conditions <- expand.grid(nF = NF_LEVELS, ipf = IPF_LEVELS, stringsAsFactors = FALSE)
conditions$total_items <- conditions$nF * conditions$ipf
TOTAL <- nrow(conditions) * REPS

cat(sprintf("Ensemble simulation: %d conditions x %d reps = %d datasets\n",
            nrow(conditions), REPS, TOTAL))
cat("Detectors: z_RR, IRV, LongString, D2, PersonTotal\n\n")

# ============================================================
# Detector score helpers
# ============================================================
compute_irv <- function(data) apply(data, 1, sd, na.rm = TRUE)

compute_longstring <- function(data) {
  apply(data, 1, function(row) {
    row <- row[!is.na(row)]
    if (length(row) == 0) return(0)
    rl <- rle(row)
    max(rl$lengths)
  })
}

compute_d2 <- function(data) {
  mat <- as.matrix(data)
  mat <- mat[, apply(mat, 2, function(x) sd(x, na.rm = TRUE) > 0), drop = FALSE]
  if (ncol(mat) < 2) return(rep(0, nrow(data)))
  cov_mat <- cov(mat, use = "pairwise.complete.obs")
  center <- colMeans(mat, na.rm = TRUE)
  cov_inv <- tryCatch(MASS::ginv(cov_mat),
                      error = function(e) solve(cov_mat + diag(1e-6, ncol(cov_mat))))
  d2 <- apply(mat, 1, function(x) {
    delta <- x - center
    as.numeric(t(delta) %*% cov_inv %*% delta)
  })
  d2[is.na(d2) | is.nan(d2)] <- 0
  d2
}

compute_person_total <- function(data) {
  mat <- as.matrix(data)
  item_means <- colMeans(mat, na.rm = TRUE)
  apply(mat, 1, function(row) {
    if (sd(row, na.rm = TRUE) == 0 || sd(item_means) == 0) return(0)
    suppressWarnings(cor(row, item_means, use = "pairwise.complete.obs"))
  })
}

# ============================================================
# Main loop — compute all scores per respondent
# ============================================================
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

    rr <- ReReReRe(data_mat, corProp = CORPROP, iterations = ITER,
                   align_signs = TRUE, mode = "auto",
                   auto_z = TRUE, variance_penalty = TRUE)
    z_rr <- rr$z_score

    scores <- data.frame(
      nF = nF, ipf = ipf, total_items = tot, rep = rep_id,
      respondent = seq_len(nrow(data_mat)),
      pattern = inj$labels$pattern,
      corruption = round(inj$labels$careless_pct, 2),
      # Detectors (raw, sign fixed so "high = more careless" for ensemble input)
      z_rr        = -z_rr,                          # low z = careless → negate
      irv         = -compute_irv(data_mat),          # low irv = careless → negate
      longstring  = compute_longstring(data_mat),    # high = careless (already)
      d2          = compute_d2(data_mat),            # high = careless (already)
      person_tot  = -compute_person_total(data_mat)  # low = careless → negate
    )
    rows[[length(rows) + 1]] <- scores

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    eta <- if (idx > 1) (elapsed / idx) * (TOTAL - idx) else NA
    cat(sprintf("[%d/%d] nF=%d ipf=%d rep=%d | %.1f min ETA %.0f\n",
                idx, TOTAL, nF, ipf, rep_id, elapsed, eta))
  }
}

df <- do.call(rbind, rows)
df$gt_careless <- (df$pattern != "clean") & (df$corruption > GT_CUTOFF)

# Scale features to comparable ranges (z-score)
for (col in c("z_rr", "irv", "longstring", "d2", "person_tot")) {
  x <- df[[col]]
  df[[col]] <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

write.csv(df, "sim_ensemble_scores.csv", row.names = FALSE)
cat("\nSaved: sim_ensemble_scores.csv\n")

# ============================================================
# Greedy forward selection (5-fold CV)
# ============================================================
set.seed(42)
folds <- sample(seq_len(N_FOLDS), nrow(df), replace = TRUE)

cv_mcc <- function(df, feats, folds, gt_col = "gt_careless") {
  preds <- numeric(nrow(df))
  for (f in unique(folds)) {
    tr <- folds != f; te <- folds == f
    formula_str <- paste0("as.integer(", gt_col, ") ~ ",
                          paste(feats, collapse = " + "))
    mod <- tryCatch(
      glm(as.formula(formula_str), data = df[tr, ], family = binomial),
      error = function(e) NULL
    )
    if (is.null(mod)) { preds[te] <- 0.5; next }
    preds[te] <- predict(mod, newdata = df[te, ], type = "response")
  }
  flag <- preds > 0.5
  truth <- df[[gt_col]]
  TP <- sum(flag & truth); TN <- sum(!flag & !truth)
  FP <- sum(flag & !truth); FN <- sum(!flag & truth)
  sens <- TP / max(1, TP + FN)
  spec <- TN / max(1, TN + FP)
  den <- sqrt(as.numeric(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  mcc <- if (den == 0) 0 else (TP*TN - FP*FN) / den
  list(sens = sens, spec = spec, mcc = mcc, preds = preds)
}

detectors <- c("z_rr", "irv", "longstring", "d2", "person_tot")

# Step 0: baseline (z_RR alone)
res_base <- cv_mcc(df, "z_rr", folds)
results <- data.frame(step = 0,
                      added = "z_rr (baseline)",
                      features = "z_rr",
                      sens = round(res_base$sens, 3),
                      spec = round(res_base$spec, 3),
                      mcc  = round(res_base$mcc,  3))

current <- c("z_rr")
remaining <- setdiff(detectors, current)

cat("\n=== GREEDY FORWARD SELECTION ===\n")
cat(sprintf("Step 0: z_rr alone | sens=%.3f spec=%.3f MCC=%.3f\n",
            res_base$sens, res_base$spec, res_base$mcc))

while (length(remaining) > 0) {
  best_mcc <- -Inf; best_det <- NULL; best_res <- NULL
  for (det in remaining) {
    feats <- c(current, det)
    res <- cv_mcc(df, feats, folds)
    if (res$mcc > best_mcc) {
      best_mcc <- res$mcc; best_det <- det; best_res <- res
    }
  }
  current <- c(current, best_det)
  remaining <- setdiff(remaining, best_det)
  step <- nrow(results)
  results <- rbind(results, data.frame(
    step = step, added = best_det,
    features = paste(current, collapse = "+"),
    sens = round(best_res$sens, 3),
    spec = round(best_res$spec, 3),
    mcc  = round(best_res$mcc,  3)
  ))
  cat(sprintf("Step %d: +%s -> MCC=%.3f (+%.3f)  sens=%.3f spec=%.3f\n",
              step, best_det, best_res$mcc,
              best_res$mcc - results$mcc[step], best_res$sens, best_res$spec))
}

write.csv(results, "sim_ensemble_forward.csv", row.names = FALSE)

# ============================================================
# Also: each detector standalone (to see individual strength)
# ============================================================
standalone_results <- data.frame()
for (det in detectors) {
  res <- cv_mcc(df, det, folds)
  standalone_results <- rbind(standalone_results,
                               data.frame(detector = det,
                                          sens = round(res$sens, 3),
                                          spec = round(res$spec, 3),
                                          mcc  = round(res$mcc,  3)))
}
standalone_results <- standalone_results %>% arrange(desc(mcc))
write.csv(standalone_results, "sim_ensemble_standalone.csv", row.names = FALSE)

cat("\n=== STANDALONE DETECTOR STRENGTH (CV MCC) ===\n")
print(standalone_results, row.names = FALSE)

cat("\n=== FORWARD SELECTION RESULTS ===\n")
print(results, row.names = FALSE)

# ============================================================
# Plots
# ============================================================
p_forward <- ggplot(results, aes(x = step, y = mcc)) +
  geom_line(linewidth = 1.2, color = "#2c5282") +
  geom_point(size = 4, color = "#2c5282") +
  geom_text(aes(label = sprintf("%.3f", mcc)), vjust = -1.0, size = 3.8) +
  geom_text(aes(label = added), angle = 30, hjust = -0.1, vjust = 1.5,
            size = 3.5, color = "#444444") +
  scale_x_continuous(breaks = 0:5,
                     labels = c("baseline", paste0("+", results$added[-1]))) +
  scale_y_continuous(limits = c(min(results$mcc) - 0.05,
                                 max(results$mcc) + 0.06)) +
  labs(
    title = "Greedy forward selection — ensemble improvement per detector added",
    subtitle = sprintf("5-fold CV, logistic regression, GT=corruption>50%%, %d respondents",
                       nrow(df)),
    x = "Step (detector added)",
    y = "MCC (cross-validated)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 20, hjust = 1))
ggsave("plot_ensemble_forward.png", p_forward,
       width = 11, height = 6, dpi = 150, bg = "white")

# Standalone + combined bar chart
bar_df <- bind_rows(
  standalone_results %>% mutate(variant = paste0("solo: ", detector)),
  results %>% mutate(variant = features) %>% select(variant, sens, spec, mcc)
) %>%
  mutate(order = row_number()) %>%
  pivot_longer(c(sens, spec, mcc), names_to = "metric", values_to = "value")

p_bars <- ggplot(bar_df %>% filter(metric == "mcc"),
                  aes(x = reorder(variant, order), y = value)) +
  geom_col(fill = "#e41a1c", width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", value)), vjust = -0.3, size = 3.2) +
  scale_y_continuous(limits = c(0, max(bar_df$value[bar_df$metric=="mcc"]) + 0.05)) +
  labs(title = "Standalone detectors vs greedy ensembles (MCC)",
       x = NULL, y = "MCC (CV)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("plot_ensemble_bars.png", p_bars,
       width = 12, height = 6, dpi = 150, bg = "white")

cat(sprintf("\nRuntime: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("Saved: sim_ensemble_scores.csv, sim_ensemble_forward.csv, sim_ensemble_standalone.csv\n")
cat("Plots: plot_ensemble_forward.png, plot_ensemble_bars.png\n")
cat("=== DONE ===\n")
