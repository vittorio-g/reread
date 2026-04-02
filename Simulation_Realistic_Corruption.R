###############################################################################
# Simulazione con corruzione realistica (10%-100%) + soglia GT al 50%
#
# Differenza dal vecchio setup:
#   PRIMA:  careless_levels = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0) — solo 50-100%
#   ORA:    careless_levels = seq(0.1, 1, 0.1) — 10%-100%
#           ground truth = careless_pct > 0.50 (solo quelli gravemente corrotti)
#           quelli con 10-40% sono "rumore" → conta come clean nella valutazione
#
# Confronto: ReReReRe (standard) vs ReReReRe_F (per-factor)
###############################################################################

setwd("C:/Users/vitto/Desktop/ReReReRe")
source("ReReReRe.R")
source("Synthetic_Good_Responses_2.R")
source("Careless_machine_2.R")

library(dplyr)

# ============================================================================
# PARAMETERS
# ============================================================================

NF_LEVELS   <- c(4, 8, 12, 20, 30)
IPF_LEVELS  <- c(3, 6, 10)
N_RESP      <- 300
PCT_CARELESS <- 0.15           # 15% careless
CARELESS_LEVELS <- seq(0.1, 1, 0.1)  # 10% to 100% corruption
GT_THRESHOLD <- 0.50           # ground truth: only >50% counts as careless
REPS        <- 5
ITERATIONS  <- 100
CORPROP     <- 0.03
Z_RANGE     <- seq(0.1, 3.0, 0.2)  # post-hoc z thresholds

SEED_BASE   <- 20260402

cat("============================================================\n")
cat("REALISTIC CORRUPTION SIMULATION\n")
cat(sprintf("  careless_levels: %s\n", paste(CARELESS_LEVELS*100, collapse=", ")))
cat(sprintf("  GT threshold: >%.0f%% corruption = careless\n", GT_THRESHOLD*100))
cat(sprintf("  Conditions: %d nF × %d ipf × %d reps = %d\n",
    length(NF_LEVELS), length(IPF_LEVELS), REPS,
    length(NF_LEVELS) * length(IPF_LEVELS) * REPS))
cat("============================================================\n\n")

# ============================================================================
# MCC function
# ============================================================================

compute_mcc <- function(pred, actual) {
  tp <- sum(pred & actual)
  tn <- sum(!pred & !actual)
  fp <- sum(pred & !actual)
  fn <- sum(!pred & actual)
  denom <- as.double(tp+fp) * as.double(tp+fn) * as.double(tn+fp) * as.double(tn+fn)
  if (denom == 0) return(0)
  (as.double(tp)*tn - as.double(fp)*fn) / sqrt(denom)
}

# ============================================================================
# CHECKPOINT
# ============================================================================

CHECKPOINT_FILE <- "sim_realistic_checkpoint.csv"
RESULTS_FILE    <- "sim_realistic_results.csv"

if (file.exists(RESULTS_FILE)) {
  results_all <- read.csv(RESULTS_FILE, stringsAsFactors=FALSE)
  cat(sprintf("Loaded %d existing results\n", nrow(results_all)))
} else {
  results_all <- data.frame()
}

done_key <- function(nf, ipf, rep) paste(nf, ipf, rep, sep="_")
done_keys <- if (nrow(results_all) > 0) {
  unique(paste(results_all$nF, results_all$ipf, results_all$rep, sep="_"))
} else character(0)

# ============================================================================
# MAIN LOOP
# ============================================================================

conditions <- expand.grid(nF=NF_LEVELS, ipf=IPF_LEVELS, rep=1:REPS)
total_cond <- nrow(conditions)
remaining <- sum(!done_key(conditions$nF, conditions$ipf, conditions$rep) %in% done_keys)

cat(sprintf("Total conditions: %d, remaining: %d\n\n", total_cond, remaining))

t_start <- Sys.time()
n_done <- 0

for (ci in seq_len(total_cond)) {
  nf  <- conditions$nF[ci]
  ipf <- conditions$ipf[ci]
  rep_i <- conditions$rep[ci]

  key <- done_key(nf, ipf, rep_i)
  if (key %in% done_keys) next

  seed <- SEED_BASE + (ci - 1) * 100 + rep_i
  total_items <- nf * ipf

  # Generate data
  set.seed(seed)
  items_per <- rep(ipf, nf)
  good_data <- tryCatch(
    simulated_good_responses(nConstructs=nf, nItems=items_per, n=N_RESP),
    error = function(e) NULL
  )
  if (is.null(good_data)) {
    cat(sprintf("  [%d/%d] nF=%d ipf=%d rep=%d — data generation FAILED\n", ci, total_cond, nf, ipf, rep_i))
    next
  }

  # Inject careless with FULL range 10-100%
  inj <- inject_careless(good_data, PCT_CARELESS,
                         careless_levels=CARELESS_LEVELS, seed=seed+500000)

  # Ground truth: only those with >50% corruption
  gt_careless <- inj$labels$careless & (inj$labels$careless_pct > GT_THRESHOLD)
  gt_noise    <- inj$labels$careless & (inj$labels$careless_pct <= GT_THRESHOLD)

  n_true_careless <- sum(gt_careless)
  n_noise <- sum(gt_noise)
  n_clean <- sum(!inj$labels$careless)

  # Run both methods
  rr_std <- tryCatch(
    ReReReRe(inj$data_corrupted, corProp=CORPROP, iterations=ITERATIONS,
             align_signs=TRUE, mode="auto"),
    error = function(e) NULL
  )

  rr_f <- tryCatch(
    ReReReRe_F(inj$data_corrupted, iterations=ITERATIONS, align_signs=TRUE),
    error = function(e) NULL
  )

  if (is.null(rr_std) || is.null(rr_f)) {
    cat(sprintf("  [%d/%d] nF=%d ipf=%d rep=%d — ReReReRe FAILED\n", ci, total_cond, nf, ipf, rep_i))
    next
  }

  # Evaluate at each z threshold
  batch <- data.frame()

  for (z_thr in Z_RANGE) {
    flag_std <- rr_std$z_score <= z_thr
    flag_f   <- rr_f$z_score <= z_thr

    mcc_std <- compute_mcc(flag_std, gt_careless)
    mcc_f   <- compute_mcc(flag_f, gt_careless)

    # Also compute sensitivity and specificity
    # (for GT: careless >50% = positive, everything else = negative)
    gt_neg <- !gt_careless  # includes clean + noise (10-40%)

    sens_std <- if (sum(gt_careless)>0) sum(flag_std & gt_careless)/sum(gt_careless) else NA
    spec_std <- if (sum(gt_neg)>0) sum(!flag_std & gt_neg)/sum(gt_neg) else NA
    sens_f   <- if (sum(gt_careless)>0) sum(flag_f & gt_careless)/sum(gt_careless) else NA
    spec_f   <- if (sum(gt_neg)>0) sum(!flag_f & gt_neg)/sum(gt_neg) else NA

    # How many noise (10-40%) respondents get wrongly flagged?
    noise_flagged_std <- if (n_noise>0) sum(flag_std[gt_noise])/n_noise else NA
    noise_flagged_f   <- if (n_noise>0) sum(flag_f[gt_noise])/n_noise else NA

    batch <- rbind(batch, data.frame(
      nF=nf, ipf=ipf, total_items=total_items, rep=rep_i,
      z_threshold=z_thr,
      n_true_careless=n_true_careless, n_noise=n_noise, n_clean=n_clean,
      mcc_std=mcc_std, mcc_f=mcc_f,
      sens_std=sens_std, spec_std=spec_std,
      sens_f=sens_f, spec_f=spec_f,
      noise_flagged_std=noise_flagged_std, noise_flagged_f=noise_flagged_f,
      mode_std=rr_std$mode_used[1], mode_f=rr_f$mode_used[1],
      n_pairs_std=rr_std$n_pairs[1], n_pairs_f=rr_f$n_pairs[1]
    ))
  }

  results_all <- rbind(results_all, batch)
  done_keys <- c(done_keys, key)

  # Save checkpoint
  write.csv(results_all, RESULTS_FILE, row.names=FALSE)

  n_done <- n_done + 1
  elapsed <- as.numeric(difftime(Sys.time(), t_start, units="mins"))
  eta <- if (n_done > 0) elapsed / n_done * (remaining - n_done) else NA

  cat(sprintf("  [%d/%d] nF=%d ipf=%d rep=%d — mcc_std=%.3f mcc_f=%.3f (z=1.5) | GT: %d careless, %d noise | ETA: %.0f min\n",
      ci, total_cond, nf, ipf, rep_i,
      batch$mcc_std[batch$z_threshold==1.5],
      batch$mcc_f[batch$z_threshold==1.5],
      n_true_careless, n_noise, eta))
}

cat(sprintf("\nDone. Total time: %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units="mins"))))

# ============================================================================
# ANALYSIS
# ============================================================================

cat("\n\n============================================================\n")
cat("ANALYSIS\n")
cat("============================================================\n")

res <- results_all

# --- Oracle best MCC per condition ---
best_std <- res %>% group_by(nF, ipf, total_items, rep) %>%
  summarize(mcc_std=max(mcc_std), .groups="drop")
best_f <- res %>% group_by(nF, ipf, total_items, rep) %>%
  summarize(mcc_f=max(mcc_f), .groups="drop")
best <- merge(best_std, best_f)

cat("\n--- OVERALL (oracle best z) ---\n")
cat(sprintf("  ReReReRe (std):  mean MCC = %.3f (SD=%.3f)\n", mean(best$mcc_std), sd(best$mcc_std)))
cat(sprintf("  ReReReRe_F:      mean MCC = %.3f (SD=%.3f)\n", mean(best$mcc_f), sd(best$mcc_f)))

# --- By item bin ---
best$item_bin <- cut(best$total_items, breaks=c(0, 30, 60, 100, 200, Inf),
                     labels=c("<30","30-60","60-100","100-200",">200"))

cat("\n--- BY ITEM BIN (oracle) ---\n")
cat(sprintf("%-10s  %8s  %8s  %8s\n", "Items", "Std", "F", "Winner"))
for (bin in levels(best$item_bin)) {
  sub <- best[best$item_bin==bin, ]
  if (nrow(sub)==0) next
  ms <- mean(sub$mcc_std); mf <- mean(sub$mcc_f)
  w <- if (ms > mf + 0.005) "Std" else if (mf > ms + 0.005) "F" else "~Tied"
  cat(sprintf("%-10s  %8.3f  %8.3f  %8s\n", bin, ms, mf, w))
}

# --- By nF ---
cat("\n--- BY nF (oracle) ---\n")
cat(sprintf("%4s  %8s  %8s  %8s\n", "nF", "Std", "F", "Winner"))
for (nf in sort(unique(best$nF))) {
  sub <- best[best$nF==nf, ]
  ms <- mean(sub$mcc_std); mf <- mean(sub$mcc_f)
  w <- if (ms > mf + 0.005) "Std" else if (mf > ms + 0.005) "F" else "~Tied"
  cat(sprintf("%4d  %8.3f  %8.3f  %8s\n", nf, ms, mf, w))
}

# --- By nF x ipf (oracle, heatmap data) ---
cat("\n--- WINNER MAP (nF x ipf, oracle) ---\n")
cat(sprintf("%4s %4s %6s  %8s  %8s  %8s  %8s\n", "nF", "ipf", "items", "Std", "F", "Winner", "Margin"))
for (nf in sort(unique(best$nF))) {
  for (ipf in sort(unique(best$ipf))) {
    sub <- best[best$nF==nf & best$ipf==ipf, ]
    if (nrow(sub)==0) next
    ms <- mean(sub$mcc_std); mf <- mean(sub$mcc_f)
    w <- if (ms > mf + 0.005) "Std" else if (mf > ms + 0.005) "F" else "~Tied"
    cat(sprintf("%4d %4d %6d  %8.3f  %8.3f  %8s  %8.3f\n", nf, ipf, nf*ipf, ms, mf, w, abs(ms-mf)))
  }
}

# --- At fixed z=1.5 ---
z15 <- res[res$z_threshold==1.5, ]

cat("\n--- AT z=1.5 (fixed) ---\n")
cat(sprintf("  Std:  mean MCC=%.3f, mean sens=%.3f, mean spec=%.3f\n",
    mean(z15$mcc_std), mean(z15$sens_std, na.rm=T), mean(z15$spec_std, na.rm=T)))
cat(sprintf("  F:    mean MCC=%.3f, mean sens=%.3f, mean spec=%.3f\n",
    mean(z15$mcc_f), mean(z15$sens_f, na.rm=T), mean(z15$spec_f, na.rm=T)))

# --- Noise flagging rate ---
cat("\n--- NOISE RESPONDENTS (10-40%% corruption) FLAGGED AT z=1.5 ---\n")
cat(sprintf("  Std: %.1f%% of noise flagged\n", 100*mean(z15$noise_flagged_std, na.rm=T)))
cat(sprintf("  F:   %.1f%% of noise flagged\n", 100*mean(z15$noise_flagged_f, na.rm=T)))

cat("\n--- NOISE FLAGGING BY ITEM BIN ---\n")
z15$item_bin <- cut(z15$total_items, breaks=c(0, 30, 60, 100, 200, Inf),
                    labels=c("<30","30-60","60-100","100-200",">200"))
for (bin in levels(z15$item_bin)) {
  sub <- z15[z15$item_bin==bin, ]
  if (nrow(sub)==0) next
  cat(sprintf("  %-10s  Std noise flagged: %.1f%%  F noise flagged: %.1f%%\n",
      bin, 100*mean(sub$noise_flagged_std, na.rm=T), 100*mean(sub$noise_flagged_f, na.rm=T)))
}

# --- Comparison with old setup (for context) ---
cat("\n\n============================================================\n")
cat("COMPARISON NOTE\n")
cat("============================================================\n")
cat("Old setup: careless_levels = 50-100%, all careless count as positive\n")
cat("New setup: careless_levels = 10-100%, only >50% counts as positive\n")
cat("The 10-40% corruption respondents act as realistic noise that\n")
cat("makes detection harder — they slightly disrupt the correlation\n")
cat("matrix but shouldn't be flagged.\n")

cat("\n=== DONE ===\n")
