# Phase3_Validation.R
# Fresh validation simulation:
#   - Fresh seeds (different from prior runs)
#   - 3 careless rates (20%, 40%, 60%) to study robustness
#   - 6 reps × 6 conditions × 3 rates = 108 datasets
#   - All features computed with iter_EFA (best from Phase 1: trim_pct ~ 0.20)
# Then evaluate using RF + full ensemble (the optimized config from Phase 1).

suppressPackageStartupMessages({
  library(MASS); library(dplyr); library(ggplot2); library(tidyr)
  library(randomForest)
})
select <- dplyr::select  # explicit shadow
source("Synthetic_Good_Responses_2.R")
source("ReReReRe.R")
source("Careless_machine_v3.R")

NF_LEVELS  <- c(8, 16, 30)
IPF_LEVELS <- c(6, 10)
REPS <- 6
CARELESS_RATES <- c(0.20, 0.40, 0.60)
N_FIXED <- 500
ITER <- 100
CORPROP <- 0.03
SEED_BASE <- 42424242   # fresh, very different from prior
CARELESS_LEVELS <- seq(0.1, 1.0, 0.1)

# Detector helpers (same as Phase 2)
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
  d2 <- apply(m, 1, function(x) { dlt <- x - center; as.numeric(t(dlt) %*% cov_inv %*% dlt) })
  d2[is.na(d2) | is.nan(d2)] <- 0; d2
}
compute_person_total <- function(d) {
  m <- as.matrix(d); im <- colMeans(m, na.rm = TRUE)
  apply(m, 1, function(r) if (sd(r, na.rm=TRUE)==0 || sd(im)==0) 0
        else suppressWarnings(cor(r, im, use = "pairwise.complete.obs")))
}

conditions <- expand.grid(nF = NF_LEVELS, ipf = IPF_LEVELS,
                          careless_rate = CARELESS_RATES,
                          stringsAsFactors = FALSE)
conditions$total_items <- conditions$nF * conditions$ipf
TOTAL <- nrow(conditions) * REPS

cat(sprintf("Phase 3 validation: %d conditions x %d reps = %d datasets\n",
            nrow(conditions), REPS, TOTAL))
cat(sprintf("Total respondents: %d\n", TOTAL * N_FIXED))

# Compute scores per dataset
rows <- list(); idx <- 0; t0 <- Sys.time()
for (ci in seq_len(nrow(conditions))) {
  nF <- conditions$nF[ci]; ipf <- conditions$ipf[ci]
  rate <- conditions$careless_rate[ci]; tot <- conditions$total_items[ci]
  for (rep_id in seq_len(REPS)) {
    idx <- idx + 1
    set.seed(SEED_BASE + (rep_id - 1) * 31 + ci)

    clean <- simulated_good_responses(nConstructs = nF, nItems = rep(ipf, nF), n = N_FIXED)
    inj <- inject_careless_v3(clean, pct_careless = rate,
                              careless_levels = CARELESS_LEVELS)
    data_mat <- inj$data_corrupted

    rr <- ReReReRe(data_mat, corProp = CORPROP, iterations = ITER,
                   align_signs = TRUE, mode = "auto",
                   auto_z = TRUE, variance_penalty = TRUE)
    z_rr <- rr$z_score

    z_iter <- tryCatch({
      rr_it <- ReReReRe_F_iterative(data_mat, iterations = ITER,
                                     align_signs = TRUE,
                                     initial_z_threshold = 1.0)
      rr_it$z_score
    }, error = function(e) rep(NA_real_, nrow(data_mat)))

    irv <- compute_irv(data_mat)
    ls  <- compute_longstring(data_mat)
    d2  <- compute_d2(data_mat)
    pt  <- compute_person_total(data_mat)

    scores <- data.frame(
      nF = nF, ipf = ipf, total_items = tot, careless_rate = rate, rep = rep_id,
      respondent = seq_len(nrow(data_mat)),
      pattern = inj$labels$pattern,
      corruption = round(inj$labels$careless_pct, 2),
      z_rr         = -z_rr,
      z_rr_iter    = -z_iter,
      irv          = -irv,
      longstring   = ls,
      d2           = d2,
      person_tot   = -pt
    )
    rows[[length(rows) + 1]] <- scores

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    eta <- if (idx > 1) (elapsed / idx) * (TOTAL - idx) else NA
    cat(sprintf("[%d/%d] nF=%d ipf=%d rate=%.2f rep=%d | %.1f min ETA %.0f\n",
                idx, TOTAL, nF, ipf, rate, rep_id, elapsed, eta))
  }
}

df_val <- do.call(rbind, rows)
write.csv(df_val, "phase3_raw_scores.csv", row.names = FALSE)
cat(sprintf("\nSaved phase3_raw_scores.csv (%d rows)\n", nrow(df_val)))

# Standardize features per-rate (different rates have different baseline distributions)
features <- c("z_rr", "z_rr_iter", "irv", "longstring", "d2", "person_tot")

# Scenario A: clean vs corruption>80% only
df_val$gt <- (df_val$pattern != "clean") & (df_val$corruption > 0.80)
df_val_A <- df_val %>% filter(pattern == "clean" | corruption > 0.80)

# Per-rate evaluation
cat("\n=== Per-rate evaluation: RF + full ensemble + Scenario A + FPR=5% ===\n")
results_by_rate <- data.frame()
for (rate in CARELESS_RATES) {
  sub <- df_val_A %>% filter(careless_rate == rate)
  # Standardize features within this rate
  for (col in features) {
    x <- sub[[col]]
    mu <- mean(x, na.rm = TRUE); sd_ <- sd(x, na.rm = TRUE)
    if (!is.na(sd_) && sd_ > 0) sub[[col]] <- (x - mu) / sd_
    sub[[col]][is.na(sub[[col]])] <- 0
  }

  # 5-fold CV with RF
  set.seed(99)
  folds_v <- sample(seq_len(5), nrow(sub), replace = TRUE)
  preds <- numeric(nrow(sub))
  for (f in 1:5) {
    tr <- folds_v != f; te <- folds_v == f
    X_tr <- as.matrix(sub[tr, features, drop = FALSE])
    X_te <- as.matrix(sub[te, features, drop = FALSE])
    y_tr <- factor(as.integer(sub$gt[tr]), levels = c(0, 1))
    mod <- randomForest(X_tr, y_tr, ntree = 200,
                         mtry = max(1, floor(sqrt(ncol(X_tr)))))
    preds[te] <- predict(mod, X_te, type = "prob")[, "1"]
  }
  thr <- quantile(preds[sub$pattern == "clean"], 0.95)
  flag <- preds > thr
  TP <- sum(flag & sub$gt); TN <- sum(!flag & !sub$gt)
  FP <- sum(flag & !sub$gt); FN <- sum(!flag & sub$gt)
  den <- sqrt(as.numeric(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  mcc <- if (den == 0) 0 else (TP*TN - FP*FN) / den

  results_by_rate <- rbind(results_by_rate, data.frame(
    careless_rate = rate, n = nrow(sub),
    n_pos = sum(sub$gt), sens = round(TP/max(1,TP+FN), 3),
    fpr = round(FP/max(1,FP+TN), 3), mcc = round(mcc, 3)
  ))
  sub$flag <- flag
  sub$pred_prob <- preds
  assign(paste0("df_val_A_rate", rate*100), sub)
}
cat("\n")
print(results_by_rate, row.names = FALSE)
write.csv(results_by_rate, "phase3_results_by_rate.csv", row.names = FALSE)

# Combined per-pattern detection across rates
all_scored <- bind_rows(
  get("df_val_A_rate20") %>% mutate(rate = 0.20),
  get("df_val_A_rate40") %>% mutate(rate = 0.40),
  get("df_val_A_rate60") %>% mutate(rate = 0.60)
)
write.csv(all_scored, "phase3_all_scored.csv", row.names = FALSE)

cat("\n=== Per-pattern detection per rate (corruption >80%) ===\n")
per_pat <- all_scored %>%
  filter(gt) %>%
  group_by(careless_rate, pattern) %>%
  summarise(n = n(), rate_det = round(mean(flag), 3), .groups = "drop") %>%
  pivot_wider(names_from = careless_rate, values_from = c(n, rate_det),
              names_glue = "{.value}_{careless_rate}") %>%
  arrange(pattern)
print(as.data.frame(per_pat), row.names = FALSE)

cat(sprintf("\nRuntime: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("=== DONE Phase 3 ===\n")
