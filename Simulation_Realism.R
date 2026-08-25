#' Simulation_Realism.R
#'
#' Confronto tra condizione "clean" (CFA pulito + 6 pattern legacy) e
#' "realistic" (CFA contaminato + popolazione mista + 9 pattern di careless,
#' inclusi anchor_noise / language_barrier / positional_fatigue).
#'
#' Output: sim_realism_scores.csv (1 riga per rispondente, 5 features +
#' separation_ratio + level del diagnostico).

suppressPackageStartupMessages({
  library(lavaan); library(dplyr)
})

ROOT <- "C:/Users/vitto/Desktop/ReReReRe"
setwd(ROOT)

source("Synthetic_Good_Responses_v3.R")
source("Careless_machine_realistic.R")
source("ReReReRe_Diagnostic.R")
source("ReReReRe.R")

# ----- Aux feature helpers -----
compute_irv <- function(X) {
  apply(X, 1, function(r) sd(as.numeric(r), na.rm = TRUE))
}
compute_longstring <- function(X) {
  apply(X, 1, function(r) {
    r <- as.numeric(r)
    rl <- rle(r)$lengths
    if (length(rl) == 0) return(0)
    max(rl)
  })
}
compute_d2 <- function(X) {
  X <- as.matrix(X)
  X <- X[, apply(X, 2, function(c) sd(c, na.rm = TRUE) > 0), drop = FALSE]
  if (ncol(X) < 2) return(rep(0, nrow(X)))
  cm <- colMeans(X, na.rm = TRUE)
  S <- cov(X, use = "pairwise.complete.obs")
  S_inv <- tryCatch(MASS::ginv(S), error = function(e) NULL)
  if (is.null(S_inv)) return(rep(0, nrow(X)))
  D <- numeric(nrow(X))
  for (i in seq_len(nrow(X))) {
    d <- as.numeric(X[i, ]) - cm
    D[i] <- as.numeric(d %*% S_inv %*% d)
  }
  D
}
compute_person_total <- function(X) {
  X <- as.matrix(X)
  col_means <- colMeans(X, na.rm = TRUE)
  apply(X, 1, function(r) suppressWarnings(cor(as.numeric(r),
                                                col_means,
                                                use = "pairwise.complete.obs")))
}

# ----- Grid -----
SIZES <- list(
  list(name = 30,  nF = 6,  ipf = 5),
  list(name = 60,  nF = 10, ipf = 6),
  list(name = 120, nF = 12, ipf = 10),
  list(name = 200, nF = 20, ipf = 10)
)
RATES <- c(0.20, 0.40)
COND  <- c("clean", "realistic")
N_REPS <- 4
N_RESP <- 300

# Realistic-condition contamination knobs
realistic_cfg <- list(
  cross_loading_prob  = 0.20,
  cross_loading_range = c(0.15, 0.35),
  method_strength     = 0.30,
  skew_prob           = 0.30,
  style_mix = list(normal = 0.55, extreme = 0.20,
                   central = 0.15, acquiescent_attentive = 0.10)
)

# Pattern sets per condizione
# Clean: 6 pattern legacy (random, longstring, pure_straight, acquiescent, mixed, fatigue)
# Realistic: tutti 9 pattern (legacy + anchor_noise + language_barrier + positional_fatigue)
clean_pct_types <- c(random = 1/6, longstring = 1/6, pure_straight = 1/6,
                     acquiescent = 1/6, mixed = 1/6, fatigue = 1/6,
                     anchor_noise = 0, language_barrier = 0,
                     positional_fatigue = 0)
realistic_pct_types <- rep(1/9, 9)
names(realistic_pct_types) <- c("random","longstring","pure_straight",
                                 "acquiescent","mixed","fatigue",
                                 "anchor_noise","language_barrier",
                                 "positional_fatigue")

cat(sprintf("[%s] Starting realism simulation\n",
            format(Sys.time(), "%H:%M:%S")))
cat(sprintf("  Grid: %d sizes x %d rates x %d conditions x %d reps = %d datasets\n",
            length(SIZES), length(RATES), length(COND), N_REPS,
            length(SIZES) * length(RATES) * length(COND) * N_REPS))

scores_rows <- list()
diag_rows   <- list()
done <- 0
total <- length(SIZES) * length(RATES) * length(COND) * N_REPS
t0 <- Sys.time()

for (size_cfg in SIZES) {
  size <- size_cfg$name; nF <- size_cfg$nF; ipf <- size_cfg$ipf
  for (rate in RATES) {
    for (cond in COND) {
      for (rep in seq_len(N_REPS)) {
        seed_v <- 20260428L + as.integer(size * 100 + rate * 10 + rep) +
                  if (cond == "realistic") 50000L else 0L

        # --- Generate good responses ---
        if (cond == "clean") {
          dat_good <- simulated_good_responses_v3(
            nConstructs = nF, nItems = rep(ipf, nF), n = N_RESP,
            cross_loading_prob = 0, method_strength = 0,
            skew_prob = 0,
            style_mix = list(normal = 1, extreme = 0, central = 0,
                             acquiescent_attentive = 0),
            seed = seed_v)
          pct_types_use <- clean_pct_types
        } else {
          dat_good <- simulated_good_responses_v3(
            nConstructs = nF, nItems = rep(ipf, nF), n = N_RESP,
            cross_loading_prob = realistic_cfg$cross_loading_prob,
            cross_loading_range = realistic_cfg$cross_loading_range,
            method_strength    = realistic_cfg$method_strength,
            skew_prob          = realistic_cfg$skew_prob,
            style_mix          = realistic_cfg$style_mix,
            seed = seed_v)
          pct_types_use <- realistic_pct_types
        }

        # --- Inject careless ---
        ic <- inject_careless_realistic(
          data = dat_good,
          pct_careless = rate,
          pct_types = pct_types_use,
          careless_levels = seq(0.1, 1.0, 0.1),
          seed = seed_v + 1L)
        Xc <- ic$data_corrupted
        labels <- ic$labels

        # --- Diagnostic ---
        diagn <- diagnose_structure(Xc, corProp = 0.03,
                                     align_signs = TRUE, warn = FALSE)
        diag_rows[[length(diag_rows) + 1]] <- data.frame(
          size = size, rate = rate, condition = cond, rep = rep,
          separation_ratio = diagn$separation_ratio,
          first_eig_ratio  = diagn$first_eig_ratio,
          n_top_k_items    = diagn$n_top_k_items,
          level            = diagn$level,
          stringsAsFactors = FALSE)

        # --- ReReReRe ---
        rr <- ReReReRe(Xc, corProp = 0.03, z_threshold = 1.5,
                       iterations = 100, align_signs = TRUE,
                       progress = FALSE)
        z_rr <- rr$z_score

        # --- Auxiliaries ---
        irv  <- compute_irv(Xc)
        long <- compute_longstring(Xc)
        d2   <- compute_d2(Xc)
        ptot <- compute_person_total(Xc)

        # --- Score table (one row per respondent) ---
        scores_rows[[length(scores_rows) + 1]] <- data.frame(
          size = size, rate = rate, condition = cond, rep = rep,
          respondent = seq_len(nrow(Xc)),
          pattern    = labels$pattern,
          corruption = labels$careless_pct,
          z_rr       = z_rr,
          irv        = irv,
          longstring = long,
          d2         = d2,
          person_tot = ptot,
          separation_ratio = diagn$separation_ratio,
          level      = diagn$level,
          stringsAsFactors = FALSE)

        done <- done + 1
        if (done %% 4 == 0 || done == total) {
          el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
          cat(sprintf("  %d/%d done (%.0fs)  size=%d rate=%.2f cond=%s rep=%d\n",
                      done, total, el, size, rate, cond, rep))
        }
      }
    }
  }
}

scores_df <- do.call(rbind, scores_rows)
diag_df   <- do.call(rbind, diag_rows)

write.csv(scores_df, "sim_realism_scores.csv", row.names = FALSE)
write.csv(diag_df,   "sim_realism_diagnostic.csv", row.names = FALSE)

cat(sprintf("\n[%s] DONE. Wrote sim_realism_scores.csv (%d rows) and sim_realism_diagnostic.csv (%d rows)\n",
            format(Sys.time(), "%H:%M:%S"), nrow(scores_df), nrow(diag_df)))

cat("\nDiagnostic summary by (size, condition):\n")
print(diag_df %>% group_by(size, condition) %>%
      summarise(mean_sep = mean(separation_ratio, na.rm = TRUE),
                mean_eig = mean(first_eig_ratio, na.rm = TRUE),
                n_ok     = sum(level == "ok"),
                n_marg   = sum(level == "marginal"),
                n_weak   = sum(level == "weak"),
                .groups = "drop"))
