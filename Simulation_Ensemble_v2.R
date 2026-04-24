# Simulation_Ensemble_v2.R
# Same setup as v1, but adds two new features:
#   - z_rr_iter_coupled: iterative coupled RR (new)
#   - z_rr_iter_efa:     ReReReRe_F_iterative (existing)
#
# Goal: test if iterative variants improve detection of random-full-careless
# (the hardest case in the ensemble).

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(tidyr); library(MASS)
})
source("Synthetic_Good_Responses_2.R")
source("ReReReRe.R")
source("Careless_machine_v3.R")
source("ReReReRe_IterCoupled.R")

# --- Config (same as v1) ---
NF_LEVELS  <- c(8, 16, 30)
IPF_LEVELS <- c(6, 10)
REPS <- 4
N_FIXED <- 500
PCT_CARELESS <- 0.40
ITER <- 100
CORPROP <- 0.03
SEED_BASE <- 2026
CARELESS_LEVELS <- seq(0.1, 1.0, 0.1)
GT_CUTOFF <- 0.50
N_FOLDS <- 5

conditions <- expand.grid(nF = NF_LEVELS, ipf = IPF_LEVELS, stringsAsFactors = FALSE)
conditions$total_items <- conditions$nF * conditions$ipf
TOTAL <- nrow(conditions) * REPS

cat(sprintf("Ensemble v2: %d conditions x %d reps = %d datasets\n",
            nrow(conditions), REPS, TOTAL))
cat("Detectors: z_rr, z_rr_iter_coupled, z_rr_iter_efa, irv, longstring, d2, person_tot\n\n")

# Helpers
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
  d2 <- apply(m, 1, function(x) { dlt <- x - center; as.numeric(t(dlt) %*% cov_inv %*% dlt) })
  d2[is.na(d2) | is.nan(d2)] <- 0; d2
}
compute_person_total <- function(d) {
  m <- as.matrix(d); im <- colMeans(m, na.rm = TRUE)
  apply(m, 1, function(r) if (sd(r, na.rm=TRUE)==0 || sd(im)==0) 0
        else suppressWarnings(cor(r, im, use = "pairwise.complete.obs")))
}

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

    # Iterative coupled
    iter_c <- score_iter_coupled(data_mat, corProp = CORPROP, iterations = ITER,
                                 align_signs = TRUE, trim_pct = 0.20)
    z_rr_iter_coupled <- iter_c$z_score

    # Iterative EFA (existing)
    z_rr_iter_efa <- tryCatch({
      rr_efa <- ReReReRe_F_iterative(data_mat, iterations = ITER,
                                      align_signs = TRUE,
                                      initial_z_threshold = 1.0)
      rr_efa$z_score
    }, error = function(e) rep(NA_real_, nrow(data_mat)))

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
      z_rr              = -z_rr,
      z_rr_iter_coupled = -z_rr_iter_coupled,
      z_rr_iter_efa     = -z_rr_iter_efa,
      irv               = -irv,
      longstring        = ls,
      d2                = d2,
      person_tot        = -pt
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

# Scale features
for (col in c("z_rr","z_rr_iter_coupled","z_rr_iter_efa","irv","longstring","d2","person_tot")) {
  x <- df[[col]]
  mu <- mean(x, na.rm = TRUE); sd_ <- sd(x, na.rm = TRUE)
  if (!is.na(sd_) && sd_ > 0) df[[col]] <- (x - mu) / sd_
  df[[col]][is.na(df[[col]])] <- 0
}

write.csv(df, "sim_ensemble_v2_scores.csv", row.names = FALSE)

# Greedy forward selection
set.seed(42)
folds <- sample(seq_len(N_FOLDS), nrow(df), replace = TRUE)

cv_mcc <- function(df, feats, folds) {
  preds <- numeric(nrow(df))
  for (f in unique(folds)) {
    tr <- folds != f; te <- folds == f
    fs <- paste0("as.integer(gt_careless) ~ ", paste(feats, collapse = " + "))
    mod <- tryCatch(glm(as.formula(fs), data = df[tr, ], family = binomial),
                    error = function(e) NULL)
    if (is.null(mod)) { preds[te] <- 0.5; next }
    preds[te] <- predict(mod, newdata = df[te, ], type = "response")
  }
  flag <- preds > 0.5
  truth <- df$gt_careless
  TP <- sum(flag & truth); TN <- sum(!flag & !truth)
  FP <- sum(flag & !truth); FN <- sum(!flag & truth)
  den <- sqrt(as.numeric(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  list(mcc = if (den == 0) 0 else (TP*TN - FP*FN) / den,
       sens = TP/max(1,TP+FN), spec = TN/max(1,TN+FP),
       preds = preds)
}

detectors <- c("z_rr", "z_rr_iter_coupled", "z_rr_iter_efa",
               "irv", "longstring", "d2", "person_tot")

# Standalone
cat("\n=== STANDALONE ===\n")
standalone <- data.frame()
for (det in detectors) {
  r <- cv_mcc(df, det, folds)
  standalone <- rbind(standalone,
    data.frame(detector = det, sens = r$sens, spec = r$spec, mcc = r$mcc))
}
standalone <- standalone %>% arrange(desc(mcc)) %>%
  mutate(across(c(sens,spec,mcc), ~ round(.x, 3)))
print(standalone, row.names = FALSE)

# Forward selection
cat("\n=== FORWARD SELECTION ===\n")
res_base <- cv_mcc(df, "z_rr", folds)
results <- data.frame(step = 0, added = "z_rr (baseline)",
                      features = "z_rr",
                      sens = round(res_base$sens, 3),
                      spec = round(res_base$spec, 3),
                      mcc  = round(res_base$mcc,  3))
current <- c("z_rr"); remaining <- setdiff(detectors, current)
cat(sprintf("Step 0: z_rr alone | MCC=%.3f\n", res_base$mcc))
while (length(remaining) > 0) {
  best_mcc <- -Inf; best_det <- NULL; best_r <- NULL
  for (det in remaining) {
    r <- cv_mcc(df, c(current, det), folds)
    if (r$mcc > best_mcc) { best_mcc <- r$mcc; best_det <- det; best_r <- r }
  }
  current <- c(current, best_det); remaining <- setdiff(remaining, best_det)
  step <- nrow(results)
  results <- rbind(results, data.frame(
    step = step, added = best_det,
    features = paste(current, collapse = "+"),
    sens = round(best_r$sens, 3), spec = round(best_r$spec, 3),
    mcc = round(best_r$mcc, 3)
  ))
  cat(sprintf("Step %d: +%s -> MCC=%.3f (+%.3f)\n",
              step, best_det, best_r$mcc,
              best_r$mcc - results$mcc[step]))
}
write.csv(results, "sim_ensemble_v2_forward.csv", row.names = FALSE)

# Random@100% detection per model at FPR=5%
cat("\n=== Per-pattern detection @ FPR=5% operating point ===\n")

# Fit chosen ensembles and compute FPR-calibrated detection
formulas <- list(
  z_rr_only   = "gt_careless ~ z_rr",
  triad_v1    = "gt_careless ~ z_rr + irv + d2",
  all_non_iter= "gt_careless ~ z_rr + irv + longstring + d2 + person_tot",
  plus_iter_c = "gt_careless ~ z_rr + z_rr_iter_coupled + irv + longstring + d2 + person_tot",
  plus_iter_e = "gt_careless ~ z_rr + z_rr_iter_efa + irv + longstring + d2 + person_tot",
  all_features= paste("gt_careless ~", paste(detectors, collapse = " + "))
)

eval_at_fpr5 <- function(formula_str, df, folds) {
  preds <- numeric(nrow(df))
  for (f in unique(folds)) {
    tr <- folds != f; te <- folds == f
    mod <- glm(as.formula(formula_str), data = df[tr, ], family = binomial)
    preds[te] <- predict(mod, newdata = df[te, ], type = "response")
  }
  thr <- quantile(preds[df$pattern == "clean"], 0.95)
  flag <- preds > thr
  # Global
  TP <- sum(flag & df$gt_careless); TN <- sum(!flag & !df$gt_careless)
  FP <- sum(flag & !df$gt_careless); FN <- sum(!flag & df$gt_careless)
  den <- sqrt(as.numeric(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  mcc <- if (den == 0) 0 else (TP*TN - FP*FN) / den
  list(preds = preds, flag = flag, thr = thr, mcc = mcc,
       sens = TP/max(1,TP+FN), fpr = FP/max(1,FP+TN))
}

per_pattern_full <- list()
for (nm in names(formulas)) {
  e <- eval_at_fpr5(formulas[[nm]], df, folds)
  df[[paste0("flag_", nm)]] <- e$flag
  cat(sprintf("%-14s  sens=%.3f  fpr=%.3f  MCC=%.3f\n",
              nm, e$sens, e$fpr, e$mcc))
  per_pattern_full[[nm]] <- df %>%
    filter(pattern != "clean", corruption >= 0.8) %>%
    mutate(bin = round(corruption * 10) / 10) %>%
    group_by(pattern, bin) %>%
    summarise(rate = round(mean(e$flag[df$pattern != "clean" & df$corruption >= 0.8][seq_len(n())]), 3),
              .groups = "drop")
}

# Simpler per-pattern summary
cat("\n=== Detection @ 100% corruption per pattern, per model (FPR=5%) ===\n")
for (nm in names(formulas)) {
  flag_col <- paste0("flag_", nm)
  summ <- df %>%
    filter(pattern != "clean", corruption == 1.0) %>%
    group_by(pattern) %>%
    summarise(n = n(), rate = round(mean(.data[[flag_col]]), 3),
              .groups = "drop")
  cat(sprintf("\n[%s]\n", nm))
  print(as.data.frame(summ), row.names = FALSE)
}

cat(sprintf("\nRuntime: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("=== DONE ===\n")
