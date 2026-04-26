# Sim_Article_Big_v2.R
# v2 of the article simulation, with a much larger grid:
#   - 8 questionnaire sizes (30, 50, 80, 100, 150, 200, 250, 300 items)
#   - 4 careless rates (0.10, 0.20, 0.40, 0.60)
#   - 12 replications per condition
#   - n = 500 respondents per dataset
#   - Full corruption grid (0, 0.1, ..., 1.0) injected via Careless_machine_v3
#   - 6 features per respondent: z_RR, z_RR_iter, IRV, LongString, D², PersonTotal
#   - Both Scenario A (tau=0.60) and Scenario B (tau=0.40) evaluated post-hoc
#   - 5-fold CV Random Forest + 13 metrics + ablation + per-pattern stats
#
# Total: 8 x 4 x 12 = 384 runs.

suppressPackageStartupMessages({
  library(dplyr); library(MASS); library(randomForest)
  library(pROC); library(PRROC)
})
source("Synthetic_Good_Responses_2.R")
source("ReReReRe.R")
source("Careless_machine_v3.R")
source("ReReReRe_IterCoupled.R")

# ------------------------------------------------------------
# Configuration (v2)
# ------------------------------------------------------------
SIZES <- list(
  list(nF =  6, ipf =  5, total =  30),
  list(nF = 10, ipf =  5, total =  50),
  list(nF = 10, ipf =  8, total =  80),
  list(nF = 10, ipf = 10, total = 100),
  list(nF = 15, ipf = 10, total = 150),
  list(nF = 20, ipf = 10, total = 200),
  list(nF = 25, ipf = 10, total = 250),
  list(nF = 30, ipf = 10, total = 300)
)
CARELESS_RATES <- c(0.10, 0.20, 0.40, 0.60)
REPS <- 12
N_FIXED <- 500
ITER <- 100
CORPROP <- 0.03
SEED_BASE <- 20260427
CARELESS_LEVELS <- seq(0.1, 1.0, 0.1)

TAU_A <- 0.60
TAU_B <- 0.40

OUT_PREFIX <- "sim_article_v2"

# ------------------------------------------------------------
# Detector helpers
# ------------------------------------------------------------
compute_irv <- function(d) apply(d, 1, sd, na.rm = TRUE)
compute_longstring <- function(d) apply(d, 1, function(r) {
  r <- r[!is.na(r)]; if (!length(r)) 0 else max(rle(r)$lengths)
})
compute_d2 <- function(d) {
  m <- as.matrix(d)
  m <- m[, apply(m, 2, function(x) sd(x, na.rm = TRUE) > 0), drop = FALSE]
  if (ncol(m) < 2) return(rep(0, nrow(d)))
  cov_mat <- cov(m, use = "pairwise.complete.obs")
  center <- colMeans(m, na.rm = TRUE)
  cov_inv <- tryCatch(ginv(cov_mat),
                      error = function(e) solve(cov_mat + diag(1e-6, ncol(cov_mat))))
  d2 <- apply(m, 1, function(x) { dlt <- x - center
    as.numeric(t(dlt) %*% cov_inv %*% dlt) })
  d2[is.na(d2) | is.nan(d2)] <- 0
  d2
}
compute_person_total <- function(d) {
  m <- as.matrix(d); im <- colMeans(m, na.rm = TRUE)
  apply(m, 1, function(r) {
    if (sd(r, na.rm = TRUE) == 0 || sd(im) == 0) 0
    else suppressWarnings(cor(r, im, use = "pairwise.complete.obs"))
  })
}

# ------------------------------------------------------------
# 13-metric calculator
# ------------------------------------------------------------
metrics_at_threshold <- function(prob, y, tau) {
  pred <- prob > tau
  TP <- sum(pred & y == 1); FP <- sum(pred & y == 0)
  TN <- sum(!pred & y == 0); FN <- sum(!pred & y == 1)
  sens <- if (TP + FN > 0) TP / (TP + FN) else NA
  spec <- if (TN + FP > 0) TN / (TN + FP) else NA
  ppv  <- if (TP + FP > 0) TP / (TP + FP) else NA
  npv  <- if (TN + FN > 0) TN / (TN + FN) else NA
  f1   <- if (!is.na(ppv) && !is.na(sens) && (ppv + sens) > 0) 2 * ppv * sens / (ppv + sens) else NA
  f2   <- if (!is.na(ppv) && !is.na(sens) && (4*ppv + sens) > 0) 5 * ppv * sens / (4*ppv + sens) else NA
  num  <- as.double(TP)*as.double(TN) - as.double(FP)*as.double(FN)
  den  <- sqrt(as.double(TP+FP) * as.double(TP+FN) *
               as.double(TN+FP) * as.double(TN+FN))
  mcc  <- if (den > 0) num/den else 0
  po   <- (TP+TN)/(TP+FP+TN+FN)
  pe   <- ((TP+FP)*(TP+FN) + (TN+FP)*(TN+FN)) / (TP+FP+TN+FN)^2
  kap  <- if (pe < 1) (po - pe)/(1 - pe) else 0
  ba   <- (sens + spec) / 2
  yj   <- sens + spec - 1
  gm   <- sqrt(sens * spec)
  list(sens = sens, spec = spec, ppv = ppv, npv = npv,
       bal_acc = ba, youden_j = yj, g_mean = gm,
       f1 = f1, f2 = f2, mcc = mcc, kappa = kap,
       TP = TP, FP = FP, TN = TN, FN = FN)
}

cv_rf_predict <- function(X, y, K = 5, seed = 42) {
  set.seed(seed)
  n <- length(y); folds <- sample(rep(1:K, length.out = n))
  oof <- numeric(n)
  for (k in 1:K) {
    tr <- folds != k; te <- folds == k
    rf <- randomForest(x = X[tr, , drop = FALSE],
                       y = factor(y[tr], levels = c(0, 1)),
                       ntree = 200, mtry = max(1, floor(sqrt(ncol(X)))))
    pp <- predict(rf, X[te, , drop = FALSE], type = "prob")
    oof[te] <- pp[, "1"]
  }
  oof
}

# ------------------------------------------------------------
# Main loop
# ------------------------------------------------------------
TOTAL <- length(SIZES) * length(CARELESS_RATES) * REPS
cat(sprintf("Total runs: %d (%d sizes x %d rates x %d reps)\n",
            TOTAL, length(SIZES), length(CARELESS_RATES), REPS))

scores_all <- list()
metrics_all <- list()
ablation_all <- list()
perpat_all <- list()

idx <- 0; t0 <- Sys.time()
for (sz in SIZES) {
  for (cr in CARELESS_RATES) {
    for (rep_id in 1:REPS) {
      idx <- idx + 1
      seed <- SEED_BASE + (rep_id - 1) * 10000 + sz$total * 10 + as.integer(cr * 100)
      set.seed(seed)

      cat(sprintf("\n[%d/%d] size=%d (nF=%d, ipf=%d) rate=%.2f rep=%d ...\n",
                  idx, TOTAL, sz$total, sz$nF, sz$ipf, cr, rep_id))
      step_t <- Sys.time()

      clean <- simulated_good_responses(nConstructs = sz$nF,
                                         nItems = rep(sz$ipf, sz$nF),
                                         n = N_FIXED)
      cat(sprintf("  clean done (%.1fs)\n",
                  as.numeric(difftime(Sys.time(), step_t, units = "secs"))))

      step_t <- Sys.time()
      inj <- inject_careless_v3(clean, pct_careless = cr,
                                 careless_levels = CARELESS_LEVELS)
      data_mat <- inj$data_corrupted
      cat(sprintf("  inject done (%.1fs)\n",
                  as.numeric(difftime(Sys.time(), step_t, units = "secs"))))

      step_t <- Sys.time()
      rr <- tryCatch(
        ReReReRe(data_mat, corProp = CORPROP, iterations = ITER,
                 align_signs = TRUE, mode = "auto",
                 auto_z = FALSE, variance_penalty = FALSE),
        error = function(e) NULL)
      z_rr <- if (is.null(rr)) rep(0, nrow(data_mat)) else rr$z_score
      cat(sprintf("  z_rr done (%.1fs)\n",
                  as.numeric(difftime(Sys.time(), step_t, units = "secs"))))

      step_t <- Sys.time()
      z_rr_iter <- tryCatch({
        rr_it <- score_iter_coupled(data_mat, corProp = CORPROP,
                                     iterations = ITER, align_signs = TRUE,
                                     trim_pct = 0.20)
        rr_it$z_score
      }, error = function(e) rep(0, nrow(data_mat)))
      cat(sprintf("  z_rr_iter done (%.1fs)\n",
                  as.numeric(difftime(Sys.time(), step_t, units = "secs"))))

      step_t <- Sys.time()
      irv <- compute_irv(data_mat)
      ls  <- compute_longstring(data_mat)
      d2  <- compute_d2(data_mat)
      pt  <- compute_person_total(data_mat)
      cat(sprintf("  aux done (%.1fs)\n",
                  as.numeric(difftime(Sys.time(), step_t, units = "secs"))))

      rec <- data.frame(
        size       = sz$total,
        nF         = sz$nF,
        ipf        = sz$ipf,
        rate       = cr,
        rep        = rep_id,
        respondent = seq_len(nrow(data_mat)),
        pattern    = inj$labels$pattern,
        corruption = round(inj$labels$careless_pct, 2),
        z_rr       = -z_rr,
        z_rr_iter  = -z_rr_iter,
        irv        = -irv,
        longstring =  ls,
        d2         =  d2,
        person_tot = -pt,
        stringsAsFactors = FALSE
      )
      scores_all[[length(scores_all) + 1]] <- rec

      features <- c("z_rr", "z_rr_iter", "irv", "longstring", "d2", "person_tot")

      for (sc in c("A", "B")) {
        tau <- if (sc == "A") TAU_A else TAU_B
        if (sc == "A") {
          keep <- rec$pattern == "clean" | rec$corruption > tau + 1e-9
        } else {
          keep <- rep(TRUE, nrow(rec))
        }
        sub <- rec[keep, ]
        y <- as.integer(sub$corruption > tau + 1e-9)
        if (sum(y == 1) < 30 || sum(y == 0) < 30) next
        Xfull <- as.matrix(sub[, features])

        prob_full <- cv_rf_predict(Xfull, y, K = 5, seed = seed)
        tau_op <- quantile(prob_full[y == 0], 0.95)
        m <- metrics_at_threshold(prob_full, y, tau_op)
        auc_v <- tryCatch(as.numeric(pROC::auc(pROC::roc(y, prob_full,
                                                           quiet = TRUE,
                                                           direction = "<"))),
                          error = function(e) NA)
        pr <- tryCatch(PRROC::pr.curve(scores.class0 = prob_full[y == 1],
                                        scores.class1 = prob_full[y == 0],
                                        curve = FALSE),
                        error = function(e) NULL)
        auprc_v <- if (!is.null(pr)) pr$auc.integral else NA

        metrics_all[[length(metrics_all) + 1]] <- data.frame(
          size = sz$total, rate = cr, rep = rep_id, scenario = sc, tau = tau,
          n_pos = sum(y == 1), n_neg = sum(y == 0),
          tau_op = as.numeric(tau_op),
          sens = m$sens, spec = m$spec, ppv = m$ppv, npv = m$npv,
          bal_acc = m$bal_acc, youden_j = m$youden_j, g_mean = m$g_mean,
          f1 = m$f1, f2 = m$f2, mcc = m$mcc, kappa = m$kappa,
          auc = auc_v, auprc = auprc_v,
          stringsAsFactors = FALSE)

        feats_no_rr <- c("irv", "longstring", "d2", "person_tot")
        Xnorr <- as.matrix(sub[, feats_no_rr])
        prob_norr <- cv_rf_predict(Xnorr, y, K = 5, seed = seed)
        tau_op_norr <- quantile(prob_norr[y == 0], 0.95)
        m_norr <- metrics_at_threshold(prob_norr, y, tau_op_norr)

        Xtriad <- as.matrix(sub[, c("z_rr_iter", "irv", "d2")])
        prob_triad <- cv_rf_predict(Xtriad, y, K = 5, seed = seed)
        tau_op_triad <- quantile(prob_triad[y == 0], 0.95)
        m_triad <- metrics_at_threshold(prob_triad, y, tau_op_triad)

        ablation_all[[length(ablation_all) + 1]] <- data.frame(
          size = sz$total, rate = cr, rep = rep_id, scenario = sc, tau = tau,
          mcc_full = m$mcc,    f1_full = m$f1,
          mcc_norr = m_norr$mcc, f1_norr = m_norr$f1,
          mcc_triad = m_triad$mcc, f1_triad = m_triad$f1,
          delta_mcc = m$mcc - m_norr$mcc,
          delta_f1  = m$f1  - m_norr$f1,
          stringsAsFactors = FALSE)

        sub$flag_full <- prob_full > tau_op
        sub$flag_norr <- prob_norr > tau_op_norr
        per_pat <- sub %>%
          filter(corruption > tau + 1e-9) %>%
          group_by(pattern, corruption) %>%
          summarise(n = n(),
                    rate_full  = mean(flag_full),
                    rate_norr  = mean(flag_norr),
                    delta = mean(flag_full) - mean(flag_norr),
                    .groups = "drop") %>%
          mutate(size = sz$total, rate_sample = cr,
                 rep = rep_id, scenario = sc, tau = tau)
        perpat_all[[length(perpat_all) + 1]] <- as.data.frame(per_pat)
      }

      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      eta <- if (idx > 0) (elapsed / idx) * (TOTAL - idx) else NA
      cat(sprintf("  -> done. elapsed %.1f min, ETA %.0f min\n", elapsed, eta))

      if (idx %% 10 == 0) {
        write.csv(do.call(rbind, scores_all),
                  paste0(OUT_PREFIX, "_scores.csv"),  row.names = FALSE)
        write.csv(do.call(rbind, metrics_all),
                  paste0(OUT_PREFIX, "_metrics.csv"), row.names = FALSE)
        write.csv(do.call(rbind, ablation_all),
                  paste0(OUT_PREFIX, "_ablation.csv"),row.names = FALSE)
        write.csv(do.call(rbind, perpat_all),
                  paste0(OUT_PREFIX, "_perpat.csv"),  row.names = FALSE)
        cat(sprintf("  [checkpoint saved at idx=%d]\n", idx))
      }
    }
  }
}

write.csv(do.call(rbind, scores_all),
          paste0(OUT_PREFIX, "_scores.csv"),  row.names = FALSE)
write.csv(do.call(rbind, metrics_all),
          paste0(OUT_PREFIX, "_metrics.csv"), row.names = FALSE)
write.csv(do.call(rbind, ablation_all),
          paste0(OUT_PREFIX, "_ablation.csv"),row.names = FALSE)
write.csv(do.call(rbind, perpat_all),
          paste0(OUT_PREFIX, "_perpat.csv"),  row.names = FALSE)

cat("\n=== DONE Sim_Article_Big_v2.R ===\n")
cat(sprintf("Total elapsed: %.1f minutes\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
