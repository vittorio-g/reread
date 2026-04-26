# Phase 6 — Systematic threshold sweep
#
# Purpose: optimize the GT threshold tau_GT by SCANNING tau in {0.10, 0.20, ..., 0.90}
# and reporting all 13 classification metrics, in TWO scenarios:
#   - Scenario A: positives = corruption > tau, negatives = clean only
#   - Scenario B: positives = corruption > tau, negatives = corruption <= tau (incl clean)
#
# Output: phase6_threshold_sweep.csv with all metrics x scenarios x thresholds
#         plot_phase6_*.png
#
# Methodologically clean: NO threshold pre-chosen by the analyst.

suppressPackageStartupMessages({
  library(randomForest)
  library(pROC)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(PRROC)
})

set.seed(42)

cat("=== Phase 6: Systematic GT threshold sweep ===\n")
cat("Loading data...\n")
sc <- read.csv("sim_robust_scores.csv", stringsAsFactors = FALSE)
cat(sprintf("Loaded %d respondents.\n", nrow(sc)))
cat("Corruption distribution:\n")
print(table(round(sc$corruption, 2)))

# Use the SAME 6 features as Phase 5 / earlier work.
# z_rr_iter is built as the LAST iteration available per row (t3 if not NA, else t2, else t1)
sc$z_rr_iter <- with(sc, ifelse(!is.na(z_rr_iter_t3), z_rr_iter_t3,
                          ifelse(!is.na(z_rr_iter_t2), z_rr_iter_t2, z_rr_iter_t1)))

feature_cols <- c("z_rr", "z_rr_iter", "irv", "longstring", "d2", "person_tot")

# Standardize features
for (f in feature_cols) {
  sc[[f]] <- as.numeric(scale(sc[[f]]))
}

# Drop rows with any NA in features
sc_full <- sc[complete.cases(sc[, feature_cols]), ]
cat(sprintf("After NA filter: %d respondents.\n", nrow(sc_full)))

# ----------------------------------------------------------------------
# Helper: 13-metric set
# ----------------------------------------------------------------------
metrics_at_threshold <- function(prob, y, tau) {
  pred <- prob > tau
  TP <- sum(pred & y == 1)
  FP <- sum(pred & y == 0)
  TN <- sum(!pred & y == 0)
  FN <- sum(!pred & y == 1)
  sens <- if (TP + FN > 0) TP / (TP + FN) else NA
  spec <- if (TN + FP > 0) TN / (TN + FP) else NA
  ppv  <- if (TP + FP > 0) TP / (TP + FP) else NA
  npv  <- if (TN + FN > 0) TN / (TN + FN) else NA
  ba   <- (sens + spec) / 2
  yj   <- sens + spec - 1
  gm   <- sqrt(sens * spec)
  prec <- ppv
  rec  <- sens
  f1   <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0) 2 * prec * rec / (prec + rec) else NA
  f2   <- if (!is.na(prec) && !is.na(rec) && (4*prec + rec) > 0) 5 * prec * rec / (4 * prec + rec) else NA
  num <- as.double(TP)*as.double(TN) - as.double(FP)*as.double(FN)
  den <- sqrt(as.double(TP+FP) * as.double(TP+FN) * as.double(TN+FP) * as.double(TN+FN))
  mcc <- if (den > 0) num/den else 0
  po <- (TP+TN)/(TP+FP+TN+FN)
  pe <- ((TP+FP)*(TP+FN) + (TN+FP)*(TN+FN)) / (TP+FP+TN+FN)^2
  kap <- if (pe < 1) (po - pe)/(1 - pe) else 0
  list(sens=sens, spec=spec, ppv=ppv, npv=npv, bal_acc=ba, youden_j=yj,
       g_mean=gm, f1=f1, f2=f2, mcc=mcc, kappa=kap,
       TP=TP, FP=FP, TN=TN, FN=FN)
}

# ----------------------------------------------------------------------
# Cross-validated RF predictions
# ----------------------------------------------------------------------
cv_rf_predict <- function(X, y, K = 5, seed = 42) {
  set.seed(seed)
  n <- length(y)
  folds <- sample(rep(1:K, length.out = n))
  oof <- numeric(n)
  for (k in 1:K) {
    tr <- folds != k
    te <- folds == k
    rf <- randomForest(x = X[tr, , drop = FALSE], y = factor(y[tr], levels = c(0,1)),
                       ntree = 200, mtry = max(1, floor(sqrt(ncol(X)))))
    pp <- predict(rf, X[te, , drop = FALSE], type = "prob")
    oof[te] <- pp[, "1"]
  }
  oof
}

# ----------------------------------------------------------------------
# Run sweep
# ----------------------------------------------------------------------
tau_grid <- seq(0.10, 0.90, by = 0.10)
results <- list()

for (scenario in c("A", "B")) {
  cat(sprintf("\n=== Scenario %s ===\n", scenario))
  for (tau in tau_grid) {
    cat(sprintf("  tau_GT = %.2f ... ", tau))
    if (scenario == "A") {
      keep <- sc_full$corruption == 0 | sc_full$corruption > tau + 1e-9
    } else {
      keep <- rep(TRUE, nrow(sc_full))   # everyone is included
    }
    sub <- sc_full[keep, ]
    y <- as.integer(sub$corruption > tau + 1e-9)
    X <- as.matrix(sub[, feature_cols])

    n_pos <- sum(y == 1); n_neg <- sum(y == 0)
    if (n_pos < 100 || n_neg < 100) { cat("skip (too few in one class)\n"); next }

    prob <- cv_rf_predict(X, y, K = 5, seed = 42)

    # Calibrate threshold to FPR = 5% (95th percentile of negatives)
    tau_op <- quantile(prob[y == 0], 0.95)

    # AUC + AUPRC (threshold-free)
    auc_val <- as.numeric(pROC::auc(pROC::roc(y, prob, quiet = TRUE,
                                              direction = "<")))
    pr <- PRROC::pr.curve(scores.class0 = prob[y == 1],
                          scores.class1 = prob[y == 0],
                          curve = FALSE)
    auprc_val <- pr$auc.integral

    # All metrics at FPR=5% operating point
    m <- metrics_at_threshold(prob, y, tau_op)

    # ALSO: max-MCC operating point
    cands <- seq(0.05, 0.95, by = 0.01)
    mcc_grid <- sapply(cands, function(t) metrics_at_threshold(prob, y, t)$mcc)
    tau_mcc <- cands[which.max(mcc_grid)]
    m_mcc <- metrics_at_threshold(prob, y, tau_mcc)

    row <- data.frame(
      scenario = scenario,
      tau_GT = tau,
      n_pos = n_pos,
      n_neg = n_neg,
      n_total = n_pos + n_neg,
      tau_op_FPR5 = as.numeric(tau_op),
      auc = auc_val,
      auprc = auprc_val,
      # @ FPR=5%
      sens = m$sens, spec = m$spec, ppv = m$ppv, npv = m$npv,
      bal_acc = m$bal_acc, youden_j = m$youden_j, g_mean = m$g_mean,
      f1 = m$f1, f2 = m$f2, mcc = m$mcc, kappa = m$kappa,
      flagged = m$TP + m$FP,
      true_positives = m$TP,
      # @ MCC-optimal
      tau_op_MCC = tau_mcc,
      sens_mcc = m_mcc$sens, spec_mcc = m_mcc$spec, ppv_mcc = m_mcc$ppv,
      f1_mcc = m_mcc$f1, mcc_max = m_mcc$mcc,
      stringsAsFactors = FALSE
    )
    results[[length(results) + 1]] <- row
    cat(sprintf("MCC@FPR5=%.3f  F1=%.3f  AUPRC=%.3f  AUC=%.3f  | maxMCC=%.3f@tau=%.2f\n",
                m$mcc, m$f1, auprc_val, auc_val, m_mcc$mcc, tau_mcc))
  }
}

res <- do.call(rbind, results)
write.csv(res, "phase6_threshold_sweep.csv", row.names = FALSE)
cat("\nSaved phase6_threshold_sweep.csv\n")

# ----------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------
res_long <- res %>%
  select(scenario, tau_GT, n_pos, n_neg, mcc, f1, f2, auprc, auc, kappa,
         bal_acc, youden_j, g_mean, sens, spec, ppv, npv) %>%
  pivot_longer(c(mcc, f1, f2, auprc, auc, kappa, bal_acc, youden_j, g_mean,
                 sens, spec, ppv, npv),
               names_to = "metric", values_to = "value")

# Plot 1: all metrics vs tau_GT, faceted by scenario
p1 <- ggplot(res_long, aes(tau_GT, value, color = metric)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
  facet_wrap(~ scenario, labeller = labeller(scenario = c(A="Scenario A: clean vs > tau",
                                                          B="Scenario B: <= tau vs > tau"))) +
  scale_x_continuous(breaks = tau_grid) +
  labs(title = "Phase 6 — Metric values across GT threshold sweep",
       subtitle = "5-fold CV Random Forest, all 6 features. Operating point: FPR = 5%.",
       x = "tau_GT (GT corruption threshold for 'careless')", y = "Metric value") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave("plot_phase6_all_metrics.png", p1, width = 11, height = 6, dpi = 150)

# Plot 2: focus on the 4 paper metrics (MCC, F1, AUPRC, Kappa)
p2 <- res_long %>%
  filter(metric %in% c("mcc","f1","auprc","kappa")) %>%
  ggplot(aes(tau_GT, value, color = scenario)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_x_continuous(breaks = tau_grid) +
  labs(title = "Phase 6 — Headline metrics across GT threshold sweep",
       subtitle = "Random Forest @ FPR=5%. Both scenarios.",
       x = "tau_GT", y = "Metric value", color = "Scenario") +
  theme_minimal(base_size = 11)

ggsave("plot_phase6_headline_metrics.png", p2, width = 10, height = 7, dpi = 150)

# Plot 3: trade-off — F1 vs n_true_positives_at_FPR5 (coverage)
p3 <- ggplot(res, aes(tau_GT, color = scenario)) +
  geom_line(aes(y = mcc, linetype = "MCC"), linewidth = 1) +
  geom_line(aes(y = true_positives / n_pos, linetype = "Sens (true careless caught)"), linewidth = 1) +
  geom_line(aes(y = n_pos / n_total, linetype = "Base rate"), linewidth = 0.6, alpha = 0.6) +
  scale_x_continuous(breaks = tau_grid) +
  labs(title = "Phase 6 — Trade-off: classifier quality vs careless coverage",
       subtitle = "MCC = index quality. Sens = % of careless flagged. Base rate = positives in dataset.",
       x = "tau_GT (GT threshold)", y = "Value", linetype = "Quantity", color = "Scenario") +
  theme_minimal(base_size = 11)

ggsave("plot_phase6_tradeoff.png", p3, width = 10, height = 6, dpi = 150)

# ----------------------------------------------------------------------
# Print summary
# ----------------------------------------------------------------------
cat("\n=== Best tau_GT by metric ===\n")
for (sc in c("A","B")) {
  cat(sprintf("\n--- Scenario %s ---\n", sc))
  sub <- res[res$scenario == sc, ]
  for (m in c("mcc","f1","auprc","kappa","youden_j","bal_acc")) {
    if (any(!is.na(sub[[m]]))) {
      i <- which.max(sub[[m]])
      cat(sprintf("  argmax %-9s : tau_GT = %.2f  (value = %.3f)\n",
                  m, sub$tau_GT[i], sub[[m]][i]))
    }
  }
}

cat("\nDone. Files: phase6_threshold_sweep.csv, plot_phase6_*.png\n")
