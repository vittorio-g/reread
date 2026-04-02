###############################################################################
# Unified Simulation: All ReReReRe variants (OPTIMIZED)
#
# Key optimization: compute per-factor z ONCE, then evaluate all proplow
# thresholds and meanvar lambdas post-hoc from the stored factor_z_matrix.
#
# Methods:
#   1. std           — ReReReRe() standard (auto: weighted/coupled)
#   2. std_cf        — standard + cross-factor baseline
#   3. efa_d         — ReReReRe_F() global EFA-D
#   4. efa_d_cf      — EFA-D + cross-factor baseline
#   5. proplow_*     — per-factor, % factors below threshold (post-hoc)
#   6. meanvar_*     — per-factor, mean/combined score (post-hoc)
#   7. iterative     — iterative EFA (2 rounds)
#
# Realistic corruption: levels 10-100%, GT = careless_pct > 50%
###############################################################################

setwd("C:/Users/vitto/Desktop/ReReReRe")
source("ReReReRe.R")
source("Synthetic_Good_Responses_2.R")
source("Careless_machine_2.R")
library(dplyr)

# ============================================================================
# PARAMETERS
# ============================================================================

NF_LEVELS    <- c(4, 8, 12, 20, 30)
IPF_LEVELS   <- c(3, 6, 10)
N_RESP       <- 300
PCT_CARELESS <- 0.15
CARELESS_LEVELS <- seq(0.1, 1, 0.1)
GT_THRESHOLD <- 0.50
REPS         <- 5
ITERATIONS   <- 50
SEED_BASE    <- 20260402

Z_RANGE         <- seq(0.1, 3.0, 0.2)
FZ_THRESHOLDS   <- c(0.5, 1.0, 1.5, 2.0)
FLAG_THRESHOLDS <- seq(0.1, 0.8, 0.1)
LAMBDA_VALUES   <- c(0.5, 1.0, 1.5)

RESULTS_FILE    <- "sim_variants_results.csv"
CHECKPOINT_FILE <- "sim_variants_checkpoint.rds"

# ============================================================================
# HELPERS
# ============================================================================

compute_mcc <- function(pred, actual) {
  tp <- sum(pred & actual); tn <- sum(!pred & !actual)
  fp <- sum(pred & !actual); fn <- sum(!pred & actual)
  d <- as.double(tp+fp) * as.double(tp+fn) * as.double(tn+fp) * as.double(tn+fn)
  if (d == 0) return(0)
  (as.double(tp)*tn - as.double(fp)*fn) / sqrt(d)
}

find_oracle_mcc <- function(scores, gt, thresholds = seq(-2, 5, 0.1), lower_is_bad = TRUE) {
  best <- -1
  for (z in thresholds) {
    pred <- if (lower_is_bad) scores <= z else scores >= z
    m <- compute_mcc(pred, gt)
    if (m > best) best <- m
  }
  best
}

# ============================================================================
# CHECKPOINT — wipe old format (had redundant EFA issue)
# ============================================================================

# Start fresh — old checkpoint had duplicated EFA problem
results_all <- data.frame()
done_keys <- character(0)

if (file.exists(CHECKPOINT_FILE)) {
  cp <- tryCatch(readRDS(CHECKPOINT_FILE), error = function(e) NULL)
  if (!is.null(cp) && !is.null(cp$version) && cp$version == 2) {
    results_all <- cp$results
    done_keys <- cp$done_keys
    cat(sprintf("Resumed v2 checkpoint: %d conditions done\n", length(done_keys)))
  } else {
    cat("Old checkpoint found, starting fresh (optimized version)\n")
  }
}

# ============================================================================
# MAIN LOOP
# ============================================================================

conditions <- expand.grid(nF = NF_LEVELS, ipf = IPF_LEVELS, rep = 1:REPS)
total_cond <- nrow(conditions)
remaining <- sum(!paste(conditions$nF, conditions$ipf, conditions$rep, sep="_") %in% done_keys)

cat(sprintf("\n============================================================\n"))
cat(sprintf("ALL VARIANTS SIMULATION (optimized) — %d conditions, %d remaining\n", total_cond, remaining))
cat(sprintf("============================================================\n\n"))

t_start <- Sys.time()
n_done <- 0

for (ci in seq_len(total_cond)) {
  nf <- conditions$nF[ci]; ipf <- conditions$ipf[ci]; rep_i <- conditions$rep[ci]
  key <- paste(nf, ipf, rep_i, sep = "_")
  if (key %in% done_keys) next

  total_items <- nf * ipf
  seed <- SEED_BASE + (ci - 1) * 100 + rep_i

  set.seed(seed)
  good_data <- tryCatch(
    simulated_good_responses(nConstructs = nf, nItems = rep(ipf, nf), n = N_RESP),
    error = function(e) NULL
  )
  if (is.null(good_data)) {
    cat(sprintf("[%d/%d] nF=%d ipf=%d rep=%d SKIP\n", ci, total_cond, nf, ipf, rep_i))
    next
  }

  inj <- inject_careless(good_data, PCT_CARELESS, careless_levels = CARELESS_LEVELS, seed = seed + 500000)
  gt <- inj$labels$careless & (inj$labels$careless_pct > GT_THRESHOLD)
  gt_noise <- inj$labels$careless & (inj$labels$careless_pct <= GT_THRESHOLD)
  corrupted <- inj$data_corrupted
  mat <- as.matrix(corrupted[, sapply(corrupted, is.numeric), drop = FALSE])

  cat(sprintf("[%d/%d] nF=%d ipf=%d (%d items) rep=%d | GT:%d/%d ... ",
      ci, total_cond, nf, ipf, total_items, rep_i, sum(gt), sum(gt_noise)))

  batch <- data.frame()

  # === 1. Standard (auto) ===
  rr_std <- tryCatch(ReReReRe(corrupted, corProp=0.03, iterations=ITERATIONS, align_signs=TRUE, mode="auto"),
                     error = function(e) NULL)

  # === 2. Standard + cross-factor ===
  rr_std_cf <- tryCatch(ReReReRe(corrupted, corProp=0.03, iterations=ITERATIONS, align_signs=TRUE, mode="auto",
                                  cross_factor_baseline=TRUE),
                        error = function(e) NULL)

  # === 3. EFA-D global ===
  rr_f <- tryCatch(ReReReRe_F(corrupted, iterations=ITERATIONS, align_signs=TRUE),
                   error = function(e) NULL)

  # === 4. EFA-D + cross-factor ===
  rr_f_cf <- tryCatch(ReReReRe_F(corrupted, iterations=ITERATIONS, align_signs=TRUE, cross_factor_baseline=TRUE),
                      error = function(e) NULL)

  # === 5+6. Per-factor z (COMPUTED ONCE, evaluated post-hoc for proplow + meanvar) ===
  pf <- tryCatch(.compute_per_factor_z(mat, iterations=ITERATIONS, align_signs=TRUE,
                                        cross_factor_baseline=FALSE),
                 error = function(e) NULL)

  # === 7. Iterative ===
  rr_iter <- tryCatch(ReReReRe_F_iterative(corrupted, iterations=ITERATIONS, align_signs=TRUE),
                      error = function(e) NULL)

  # ------ Evaluate z-score methods ------
  z_methods <- list()
  if (!is.null(rr_std)) z_methods[["std"]] <- rr_std$z_score
  if (!is.null(rr_std_cf)) z_methods[["std_cf"]] <- rr_std_cf$z_score
  if (!is.null(rr_f)) z_methods[["efa_d"]] <- rr_f$z_score
  if (!is.null(rr_f_cf)) z_methods[["efa_d_cf"]] <- rr_f_cf$z_score
  if (!is.null(rr_iter)) z_methods[["iterative"]] <- rr_iter$z_score

  # Per-factor derived scores (computed once, evaluated at multiple params)
  if (!is.null(pf)) {
    fz <- pf$factor_z_matrix
    pf_mean <- rowMeans(fz, na.rm = TRUE)
    pf_var  <- apply(fz, 1, var, na.rm = TRUE)

    z_methods[["pf_mean"]] <- pf_mean

    for (lam in LAMBDA_VALUES) {
      z_methods[[sprintf("pf_comb_l%.1f", lam)]] <- pf_mean - lam * sqrt(pmax(pf_var, 0))
    }
  }

  for (method_name in names(z_methods)) {
    scores <- z_methods[[method_name]]
    oracle <- find_oracle_mcc(scores, gt)

    for (zt in Z_RANGE) {
      flagged <- scores <= zt
      mcc <- compute_mcc(flagged, gt)
      sens <- if (sum(gt)>0) sum(flagged & gt)/sum(gt) else NA
      spec <- if (sum(!gt)>0) sum(!flagged & !gt)/sum(!gt) else NA
      nf_rate <- if (sum(gt_noise)>0) sum(flagged[gt_noise])/sum(gt_noise) else NA

      batch <- rbind(batch, data.frame(
        nF=nf, ipf=ipf, total_items=total_items, rep=rep_i,
        method=method_name, z_threshold=zt,
        mcc=mcc, sens=sens, spec=spec, noise_flagged=nf_rate,
        oracle_mcc=oracle, stringsAsFactors=FALSE))
    }
  }

  # ------ PropLow (higher = worse, different threshold space) ------
  if (!is.null(pf)) {
    fz <- pf$factor_z_matrix
    for (fzt in FZ_THRESHOLDS) {
      prop_low <- rowMeans(fz < fzt, na.rm = TRUE)
      oracle_pl <- find_oracle_mcc(prop_low, gt, seq(0.05, 1, 0.05), lower_is_bad = FALSE)

      for (ft in FLAG_THRESHOLDS) {
        flagged <- prop_low >= ft
        mcc <- compute_mcc(flagged, gt)
        sens <- if (sum(gt)>0) sum(flagged & gt)/sum(gt) else NA
        spec <- if (sum(!gt)>0) sum(!flagged & !gt)/sum(!gt) else NA
        nf_rate <- if (sum(gt_noise)>0) sum(flagged[gt_noise])/sum(gt_noise) else NA

        batch <- rbind(batch, data.frame(
          nF=nf, ipf=ipf, total_items=total_items, rep=rep_i,
          method=sprintf("proplow_fz%.1f", fzt), z_threshold=ft,
          mcc=mcc, sens=sens, spec=spec, noise_flagged=nf_rate,
          oracle_mcc=oracle_pl, stringsAsFactors=FALSE))
      }
    }
  }

  results_all <- rbind(results_all, batch)
  done_keys <- c(done_keys, key)

  tryCatch({
    saveRDS(list(results=results_all, done_keys=done_keys, version=2), CHECKPOINT_FILE)
    write.csv(results_all, RESULTS_FILE, row.names=FALSE)
  }, error = function(e) cat(sprintf(" (save err: %s)", e$message)))

  n_done <- n_done + 1
  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  eta <- if (n_done > 0) elapsed / n_done * (remaining - n_done) else NA

  std_mcc <- if (!is.null(rr_std)) find_oracle_mcc(rr_std$z_score, gt) else NA
  f_mcc <- if (!is.null(rr_f)) find_oracle_mcc(rr_f$z_score, gt) else NA
  cat(sprintf("std=%.3f f=%.3f | %.1f min ETA %.0f\n", std_mcc, f_mcc, elapsed, eta))
}

total_time <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat(sprintf("\nDone. %.1f min\n", total_time))

# ============================================================================
# ANALYSIS
# ============================================================================

cat("\n\n============================================================\n")
cat("ANALYSIS\n")
cat("============================================================\n\n")

res <- results_all

# Core methods (z-score based)
core_names <- c("std","std_cf","efa_d","efa_d_cf","iterative","pf_mean",
                "pf_comb_l0.5","pf_comb_l1.0","pf_comb_l1.5")

# --- Oracle MCC all methods ---
oracle_all <- res %>%
  group_by(method) %>%
  summarize(oracle = mean(oracle_mcc, na.rm=T), .groups="drop") %>%
  arrange(desc(oracle))

cat("--- ORACLE MCC (all) ---\n")
print(as.data.frame(oracle_all), digits=3)

# --- Core methods by item bin ---
core_res <- res %>% filter(method %in% core_names)

core_bin <- core_res %>%
  mutate(bin = cut(total_items, c(0,30,60,100,200,Inf),
                   labels=c("<30","30-60","60-100","100-200",">200"))) %>%
  group_by(bin, method) %>%
  summarize(oracle = mean(oracle_mcc, na.rm=T), .groups="drop") %>%
  tidyr::pivot_wider(names_from=method, values_from=oracle)

cat("\n--- ORACLE BY ITEM BIN ---\n")
print(as.data.frame(core_bin), digits=3)

# --- By nF ---
core_nf <- core_res %>%
  group_by(nF, method) %>%
  summarize(oracle = mean(oracle_mcc, na.rm=T), .groups="drop") %>%
  tidyr::pivot_wider(names_from=method, values_from=oracle)

cat("\n--- ORACLE BY nF ---\n")
print(as.data.frame(core_nf), digits=3)

# --- Cross-factor effect ---
cat("\n--- CROSS-FACTOR BASELINE EFFECT ---\n")
cf_eff <- core_res %>%
  filter(method %in% c("std","std_cf","efa_d","efa_d_cf")) %>%
  group_by(method) %>%
  summarize(oracle = mean(oracle_mcc, na.rm=T), .groups="drop")
print(as.data.frame(cf_eff), digits=3)

# --- PropLow best ---
proplow_res <- res %>% filter(grepl("proplow", method))
if (nrow(proplow_res) > 0) {
  cat("\n--- PROPLOW BY THRESHOLD ---\n")
  pl_summary <- proplow_res %>%
    group_by(method) %>%
    summarize(oracle = mean(oracle_mcc, na.rm=T), .groups="drop") %>%
    arrange(desc(oracle))
  print(as.data.frame(pl_summary), digits=3)
}

# --- Winner map ---
cat("\n--- WINNER MAP (nF x ipf) ---\n")
top <- c("std","efa_d","efa_d_cf","iterative","pf_mean","pf_comb_l1.0")
wm <- core_res %>%
  filter(method %in% top) %>%
  group_by(nF, ipf, method) %>%
  summarize(oracle = mean(oracle_mcc, na.rm=T), .groups="drop") %>%
  group_by(nF, ipf) %>%
  slice_max(oracle, n=1, with_ties=FALSE) %>%
  ungroup()
cat(sprintf("%4s %4s %6s  %-20s %6s\n", "nF", "ipf", "items", "winner", "MCC"))
for (i in seq_len(nrow(wm))) {
  cat(sprintf("%4d %4d %6d  %-20s %.3f\n",
      wm$nF[i], wm$ipf[i], wm$nF[i]*wm$ipf[i], wm$method[i], wm$oracle[i]))
}

# --- Save report ---
report_dir <- "archive/variants_comparison"
dir.create(report_dir, recursive=TRUE, showWarnings=FALSE)

sink(file.path(report_dir, "report.txt"))
cat(sprintf("ALL VARIANTS | %s | %d cond | %.1f min\n\n", Sys.Date(), total_cond, total_time))
cat("=== ORACLE ALL ===\n"); print(as.data.frame(oracle_all), digits=3)
cat("\n=== BY BIN ===\n"); print(as.data.frame(core_bin), digits=3)
cat("\n=== BY nF ===\n"); print(as.data.frame(core_nf), digits=3)
cat("\n=== CF EFFECT ===\n"); print(as.data.frame(cf_eff), digits=3)
if (exists("pl_summary")) { cat("\n=== PROPLOW ===\n"); print(as.data.frame(pl_summary), digits=3) }
cat("\n=== WINNERS ===\n")
for (i in seq_len(nrow(wm))) cat(sprintf("  nF=%d ipf=%d: %s (%.3f)\n", wm$nF[i], wm$ipf[i], wm$method[i], wm$oracle[i]))
sink()

cat(sprintf("\nReport: %s/report.txt\n", report_dir))
cat("\n=== DONE ===\n")
