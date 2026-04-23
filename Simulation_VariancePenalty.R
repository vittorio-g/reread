# Simulation_VariancePenalty.R
# Test: idea A (z-score penalizzato per varianza within-respondent).
#
# Design:
#   - 5 tipi di careless: random, scattered-longstring (existing), pure-straight,
#     acquiescent, mixed. 20% ciascuno.
#   - Corruzione seq(0.1, 1.0, 0.1); GT = pct > 0.50.
#   - Evaluazione: MCC a z=1.5 per standard vs variance-penalized (post-hoc).
#   - Grid di (alpha, beta) esplorata post-hoc: una sola chiamata a ReReReRe per cella.
#
# Runtime atteso: ~25-40 minuti (meno della all-variants perché 1 sola chiamata).

suppressPackageStartupMessages({ library(dplyr) })
source("Synthetic_Good_Responses_2.R")
source("ReReReRe.R")

# ============================================================
# Extended careless injector (adds pure_straight + acquiescent)
# ============================================================
inject_careless_extended <- function(data, pct_careless = 0.15,
                                     pct_types = c(random = 0.20, longstring = 0.20,
                                                    pure_straight = 0.20,
                                                    acquiescent = 0.20,
                                                    mixed = 0.20),
                                     careless_levels = seq(0.1, 1.0, 0.1),
                                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- as.data.frame(data); N <- nrow(X); J <- ncol(X)
  for (jj in seq_len(J)) if (!is.numeric(X[[jj]])) X[[jj]] <- as.numeric(as.character(X[[jj]]))
  vals_all <- unlist(X); vals_all <- vals_all[!is.na(vals_all)]
  LEVELS <- sort(unique(as.integer(round(vals_all))))
  if (length(LEVELS) < 2 || length(LEVELS) > 15) LEVELS <- seq(min(LEVELS), max(LEVELS))
  high_LEVELS <- LEVELS[LEVELS >= quantile(LEVELS, 0.65)]  # top ~third for acquiescent

  pct_types <- pct_types / sum(pct_types)
  M <- round(N * pct_careless)
  type_names <- names(pct_types)
  counts <- setNames(floor(M * pct_types), type_names)
  rem <- M - sum(counts)
  if (rem > 0) counts[seq_len(rem)] <- counts[seq_len(rem)] + 1

  careless_ids <- sample(seq_len(N), M, replace = FALSE)

  # Assign type + corruption level to each careless id
  pat_vec <- character(0); pct_vec <- numeric(0)
  for (pat in type_names) {
    Mp <- counts[pat]
    if (Mp == 0) next
    L <- length(careless_levels)
    cntL <- rep(floor(Mp/L), L); rem2 <- Mp - sum(cntL)
    if (rem2 > 0) cntL[sample(L, rem2)] <- cntL[sample(L, rem2)] + 1
    pcts <- unlist(mapply(function(p, n) rep(p, n), careless_levels, cntL, SIMPLIFY = FALSE))
    if (length(pcts) > 1) pcts <- sample(pcts)
    pat_vec <- c(pat_vec, rep(pat, Mp)); pct_vec <- c(pct_vec, pcts)
  }
  perm <- sample(seq_along(pat_vec))
  pat_vec <- pat_vec[perm]; pct_vec <- pct_vec[perm]

  # Corruption functions
  do_random <- function(row, k) {
    idx <- sample(J, k); row[idx] <- sample(LEVELS, k, replace = TRUE); list(row = row, idx = idx)
  }
  do_longstring <- function(row, k) {
    idx <- sample(J, k)
    n_chunks <- if (k >= 4) sample(2:min(5, k), 1) else if (k >= 2) 2 else 1
    chunk <- sample(rep(seq_len(n_chunks), length.out = k))
    for (ch in seq_len(n_chunks)) {
      row[idx[chunk == ch]] <- sample(LEVELS, 1)
    }
    list(row = row, idx = sort(idx))
  }
  do_pure_straight <- function(row, k) {
    # Ignore k and pct level: replace ALL items with one constant value.
    val <- sample(LEVELS, 1)
    row[] <- val
    list(row = row, idx = seq_len(J))
  }
  do_acquiescent <- function(row, k) {
    # Replace ALL items with values from the upper-third of the Likert scale,
    # with small noise. Ignore k since acquiescence is by definition global.
    vals <- sample(high_LEVELS, J, replace = TRUE)
    list(row = vals, idx = seq_len(J))
  }
  do_mixed <- function(row, k) {
    if (k <= 1) return(do_random(row, k))
    k_long <- max(1, round(k/2)); k_rand <- k - k_long
    L <- do_longstring(row, k_long); row2 <- L$row
    avail <- setdiff(seq_len(J), L$idx)
    if (length(avail) > 0 && k_rand > 0) {
      ri <- sample(avail, min(k_rand, length(avail)))
      row2[ri] <- sample(LEVELS, length(ri), replace = TRUE)
      return(list(row = row2, idx = sort(unique(c(L$idx, ri)))))
    }
    list(row = row2, idx = L$idx)
  }

  Xc <- X
  labels <- data.frame(careless = rep(FALSE, N), careless_pct = 0,
                       pattern = "clean", stringsAsFactors = FALSE)
  for (m in seq_len(M)) {
    i <- careless_ids[m]; pct_i <- pct_vec[m]; pat_i <- pat_vec[m]
    k <- min(J, max(0, round(J * pct_i)))
    orig <- as.numeric(Xc[i, , drop = TRUE])
    res <- switch(pat_i,
      "random" = do_random(orig, k),
      "longstring" = do_longstring(orig, k),
      "pure_straight" = do_pure_straight(orig, k),
      "acquiescent" = do_acquiescent(orig, k),
      "mixed" = do_mixed(orig, k))
    Xc[i, ] <- res$row
    labels$careless[i] <- TRUE
    # For the global types (pure_straight, acquiescent), pct is always 1.0
    # since the whole row is replaced.
    labels$careless_pct[i] <- if (pat_i %in% c("pure_straight","acquiescent")) 1.0 else length(res$idx)/J
    labels$pattern[i] <- pat_i
  }
  list(data_corrupted = Xc, labels = labels)
}

# ============================================================
# Grid
# ============================================================
NF_LEVELS <- c(8, 16, 30)
IPF_LEVELS <- c(6, 10)
REPS <- 5
N_FIXED <- 300
PCT_CARELESS <- 0.20  # 20% careless, split across 5 types = 4% each
ITER <- 100
CORPROP <- 0.03
SEED_BASE <- 2026
CARELESS_LEVELS <- seq(0.1, 1.0, 0.1)
GT_CUTOFF <- 0.50

Z_GRID <- seq(0.1, 3.0, 0.1)
# Penalty grid (post-hoc)
ALPHA_GRID <- c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
BETA_GRID  <- c(0.3, 0.5, 0.8, 1.0, 1.5)

conditions <- expand.grid(nF = NF_LEVELS, ipf = IPF_LEVELS, stringsAsFactors = FALSE)
conditions$total_items <- conditions$nF * conditions$ipf
TOTAL_CONDITIONS <- nrow(conditions) * REPS

cat("============================================================\n")
cat(" VARIANCE-PENALTY SIMULATION\n")
cat("============================================================\n")
cat(sprintf("%d conditions x %d reps = %d runs\n",
            nrow(conditions), REPS, TOTAL_CONDITIONS))
cat(sprintf("Careless types: random, longstring, pure_straight, acquiescent, mixed (20%% each)\n"))
cat(sprintf("Penalty grid: %d alpha x %d beta = %d variants per condition\n\n",
            length(ALPHA_GRID), length(BETA_GRID),
            length(ALPHA_GRID) * length(BETA_GRID)))

# ============================================================
# Helpers
# ============================================================
compute_mcc <- function(z, lab, zt) {
  flag <- z <= zt
  TP <- sum(flag & lab); TN <- sum(!flag & !lab)
  FP <- sum(flag & !lab); FN <- sum(!flag & lab)
  den <- sqrt(as.numeric(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  if (den == 0) 0 else (TP*TN - FP*FN) / den
}
oracle_mcc <- function(z, lab) {
  max(sapply(Z_GRID, function(zt) compute_mcc(z, lab, zt)), na.rm = TRUE)
}

# ============================================================
# Main loop
# ============================================================
results <- list(); per_pattern_results <- list()
idx <- 0; t0 <- Sys.time()

for (ci in seq_len(nrow(conditions))) {
  nF <- conditions$nF[ci]; ipf <- conditions$ipf[ci]; tot <- conditions$total_items[ci]
  for (rep_id in seq_len(REPS)) {
    idx <- idx + 1
    set.seed(SEED_BASE + (rep_id - 1) * 17 + ci)

    clean <- simulated_good_responses(nConstructs = nF, nItems = rep(ipf, nF), n = N_FIXED)
    inj <- inject_careless_extended(clean, pct_careless = PCT_CARELESS,
                                    careless_levels = CARELESS_LEVELS)
    labels <- inj$labels$careless_pct > GT_CUTOFF
    pattern <- inj$labels$pattern
    data_mat <- inj$data_corrupted

    # One RR call
    rr <- ReReReRe(data_mat, corProp = CORPROP, iterations = ITER,
                   align_signs = TRUE, mode = "auto")
    z_base <- rr$z_score
    sd_resp <- apply(data_mat, 1, sd, na.rm = TRUE)

    # Baseline (no penalty)
    mcc_base_15  <- compute_mcc(z_base, labels, 1.5)
    mcc_base_or  <- oracle_mcc(z_base, labels)

    # Explore (alpha, beta) grid
    for (alpha in ALPHA_GRID) {
      for (beta in BETA_GRID) {
        z_adj <- z_base - alpha * exp(-sd_resp / beta)
        m15 <- compute_mcc(z_adj, labels, 1.5)
        m_or <- oracle_mcc(z_adj, labels)
        results[[length(results) + 1]] <- data.frame(
          nF = nF, ipf = ipf, total_items = tot, rep = rep_id,
          alpha = alpha, beta = beta,
          mcc_z15 = round(m15, 4),
          mcc_oracle = round(m_or, 4)
        )

        # Per-pattern MCC (at z=1.5 only, for compactness)
        for (pat in c("clean", "random", "longstring", "pure_straight",
                      "acquiescent", "mixed")) {
          # True positives for this pattern subset (plus all clean as negatives)
          mask <- (pattern == pat) | (pattern == "clean")
          lab_sub <- labels[mask]
          z_sub <- z_adj[mask]
          m_pat <- if (sum(lab_sub) > 0 && sum(!lab_sub) > 0)
            compute_mcc(z_sub, lab_sub, 1.5) else NA
          per_pattern_results[[length(per_pattern_results) + 1]] <- data.frame(
            nF = nF, ipf = ipf, total_items = tot, rep = rep_id,
            alpha = alpha, beta = beta, pattern = pat,
            mcc_pattern_z15 = round(m_pat, 4)
          )
        }
      }
    }

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    eta <- if (idx > 1) (elapsed / idx) * (TOTAL_CONDITIONS - idx) else NA
    cat(sprintf("[%d/%d] nF=%d ipf=%d rep=%d | base z=1.5: %.3f | oracle: %.3f | %.1f min ETA %.0f\n",
                idx, TOTAL_CONDITIONS, nF, ipf, rep_id,
                mcc_base_15, mcc_base_or, elapsed, eta))
  }
}

df <- do.call(rbind, results)
df_pat <- do.call(rbind, per_pattern_results)
write.csv(df, "sim_varpen_results.csv", row.names = FALSE)
write.csv(df_pat, "sim_varpen_per_pattern.csv", row.names = FALSE)

# ============================================================
# Analysis
# ============================================================
cat("\n============================================================\n")
cat(" RESULTS\n")
cat("============================================================\n\n")

cat("=== BASELINE (alpha=0) ===\n")
base <- df %>% filter(alpha == 0) %>%
  summarise(mcc_z15 = mean(mcc_z15), mcc_oracle = mean(mcc_oracle))
print(as.data.frame(base))

cat("\n=== BEST (alpha, beta) BY MCC z=1.5 (averaged over all cells/reps) ===\n")
best_grid <- df %>%
  group_by(alpha, beta) %>%
  summarise(mcc_z15 = mean(mcc_z15), mcc_oracle = mean(mcc_oracle), .groups = "drop") %>%
  arrange(desc(mcc_z15))
print(as.data.frame(best_grid), row.names = FALSE)

cat("\n=== BEST combo heatmap (alpha x beta, MCC z=1.5) ===\n")
hm <- best_grid %>%
  select(alpha, beta, mcc_z15) %>%
  tidyr::pivot_wider(names_from = beta, values_from = mcc_z15,
                     names_prefix = "beta=")
print(as.data.frame(hm), row.names = FALSE)

cat("\n=== PER-PATTERN MCC at z=1.5 (averaged) ===\n")
pat_summ <- df_pat %>%
  group_by(alpha, beta, pattern) %>%
  summarise(mcc = round(mean(mcc_pattern_z15, na.rm = TRUE), 3), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = pattern, values_from = mcc)
# Only print baseline and top 3 (alpha, beta) by overall mcc_z15
top3 <- head(best_grid, 3)
cat("Baseline (alpha=0):\n")
print(as.data.frame(pat_summ %>% filter(alpha == 0) %>% select(-beta) %>% distinct()),
      row.names = FALSE)
cat("\nTop 3 penalized variants:\n")
for (i in seq_len(nrow(top3))) {
  a <- top3$alpha[i]; b <- top3$beta[i]
  cat(sprintf("alpha=%.2f beta=%.2f (overall MCC=%.3f):\n",
              a, b, top3$mcc_z15[i]))
  print(as.data.frame(pat_summ %>% filter(alpha == a, beta == b)),
        row.names = FALSE)
  cat("\n")
}

cat(sprintf("Total runtime: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("Saved: sim_varpen_results.csv + sim_varpen_per_pattern.csv\n")
cat("=== DONE ===\n")
