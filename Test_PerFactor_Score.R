###############################################################################
# Test: Per-Factor Coherence Scores
#
# Instead of ONE global coherence score, compute a z-score PER FACTOR,
# then aggregate. This captures partial/local carelessness better.
#
# Aggregation strategies tested:
#   - mean_z:    mean of per-factor z-scores (≈ global, but cleaner)
#   - min_z:     minimum z across factors (worst factor = most suspicious)
#   - prop_low:  proportion of factors with z < 1.0 (how many factors fail)
#   - sd_z:      SD of z across factors (high = inconsistent profile)
#   - combined:  mean_z + penalty for low min_z
#
# Compared against: EFA-D global (current default)
###############################################################################

library(lavaan)
library(psych)
library(dplyr)
library(ggplot2)
library(tidyr)

source("Synthetic_Good_Responses_2.R")
source("Careless_machine_2.R")
source("ReReReRe.R")

# ============================================================================
# HELPERS
# ============================================================================

evaluate_mcc <- function(predicted, actual) {
  tp <- sum(predicted == 1 & actual == 1)
  tn <- sum(predicted == 0 & actual == 0)
  fp <- sum(predicted == 1 & actual == 0)
  fn <- sum(predicted == 0 & actual == 1)
  denom <- as.double(tp + fp) * as.double(tp + fn) *
           as.double(tn + fp) * as.double(tn + fn)
  if (denom == 0) return(0)
  return((as.double(tp) * tn - as.double(fp) * fn) / sqrt(denom))
}

find_best_z <- function(scores, labels, thresholds = seq(-3, 5, 0.1), lower_is_bad = TRUE) {
  best_mcc <- -1; best_z <- NA
  for (z in thresholds) {
    if (lower_is_bad) {
      pred <- as.integer(scores <= z)
    } else {
      pred <- as.integer(scores >= z)
    }
    mcc <- evaluate_mcc(pred, labels)
    if (mcc > best_mcc) { best_mcc <- mcc; best_z <- z }
  }
  return(list(z = best_z, mcc = best_mcc))
}

# ============================================================================
# PER-FACTOR COHERENCE
# ============================================================================

ReReReRe_perFactor <- function(data, iterations = 50, align_signs = TRUE,
                                min_items_per_factor = 3) {
  data <- data[, sapply(data, is.numeric), drop = FALSE]
  mat <- as.matrix(data)
  N <- nrow(mat); p <- ncol(mat)
  item_max <- apply(mat, 2, max, na.rm = TRUE)

  # Step 1: Estimate nF
  pa <- tryCatch({
    suppressMessages(suppressWarnings(
      fa.parallel(mat, fa = "fa", plot = FALSE, n.iter = 20)
    ))
  }, error = function(e) NULL)

  nF_est <- if (!is.null(pa) && !is.null(pa$nfact) && pa$nfact >= 1) pa$nfact
            else max(1, sum(eigen(cor(mat, use="pairwise.complete.obs"),
                                   symmetric=TRUE, only.values=TRUE)$values > 1))
  nF_est <- max(1, min(nF_est, floor(p / 2)))

  # Step 2: EFA
  efa_result <- tryCatch({
    suppressWarnings(fa(mat, nfactors = nF_est, rotate = "oblimin", fm = "minres",
                        scores = "none", warnings = FALSE))
  }, error = function(e) {
    tryCatch({
      suppressWarnings(fa(mat, nfactors = max(1, nF_est-1), rotate = "oblimin",
                          fm = "minres", scores = "none", warnings = FALSE))
    }, error = function(e2) NULL)
  })

  if (is.null(efa_result)) {
    return(list(mean_z = rep(NA, N), min_z = rep(NA, N),
                prop_low = rep(NA, N), sd_z = rep(NA, N),
                combined = rep(NA, N), global_z = rep(NA, N),
                nF_detected = nF_est, n_factors_used = 0))
  }

  loadings_mat <- as.matrix(efa_result$loadings[])
  nF_actual <- ncol(loadings_mat)
  primary_factor <- apply(abs(loadings_mat), 1, which.max)

  # Precompute sign matrix
  raw_cor_mat <- cor(mat, use = "pairwise.complete.obs")
  signMat <- sign(raw_cor_mat)
  signMat[upper.tri(signMat, diag = TRUE)] <- NA

  # Step 3: Per-factor computation
  factor_z_matrix <- matrix(NA_real_, nrow = N, ncol = nF_actual)
  global_pairs_A <- c(); global_pairs_B <- c(); global_weights <- c(); global_signs <- c()
  factors_used <- 0

  for (f in seq_len(nF_actual)) {
    items_f <- which(primary_factor == f)
    if (length(items_f) < min_items_per_factor) next

    # Within-factor pairs for this factor
    combos <- combn(items_f, 2)
    idx_A <- combos[1, ]
    idx_B <- combos[2, ]
    k_f <- length(idx_A)

    # Get observed |r| and sign
    pair_r <- numeric(k_f)
    for (j in seq_len(k_f)) pair_r[j] <- raw_cor_mat[idx_A[j], idx_B[j]]
    pair_abs_r <- abs(pair_r)
    pair_sign <- sign(pair_r)

    keep <- !is.na(pair_abs_r)
    if (sum(keep) < 2) next

    idx_A <- idx_A[keep]; idx_B <- idx_B[keep]
    pair_abs_r <- pair_abs_r[keep]; pair_sign <- pair_sign[keep]
    k_f <- length(idx_A)
    weights <- pair_abs_r; weights[weights < 1e-6] <- 1e-6

    # Store for global computation
    global_pairs_A <- c(global_pairs_A, idx_A)
    global_pairs_B <- c(global_pairs_B, idx_B)
    global_weights <- c(global_weights, weights)
    global_signs <- c(global_signs, pair_sign)

    # Build A/B matrices
    A_f <- mat[, idx_A, drop = FALSE]
    B_f <- mat[, idx_B, drop = FALSE]

    # Sign alignment
    if (align_signs) {
      nf <- which(pair_sign < 0)
      if (length(nf) > 0) {
        B_f[, nf] <- rep(item_max[idx_B[nf]] + 1, each = N) - B_f[, nf]
      }
    }

    # Coupled score for this factor
    coupled_f <- rowCor_weighted(A_f, B_f, weights)

    # Permutation baseline for this factor
    rand_f <- matrix(NA_real_, N, iterations)
    for (iter in seq_len(iterations)) {
      ri1 <- sample(p, k_f, replace = TRUE)
      ri2 <- sample(p, k_f, replace = TRUE)
      same <- ri1 == ri2
      while (any(same)) { ri2[same] <- sample(p, sum(same), replace=TRUE); same <- ri1==ri2 }
      Ar <- mat[, ri1, drop=FALSE]; Br <- mat[, ri2, drop=FALSE]
      if (align_signs) {
        ri <- pmax(ri1,ri2); ci <- pmin(ri1,ri2)
        rs <- signMat[cbind(ri,ci)]; rs[is.na(rs)] <- 1
        rn <- which(rs < 0)
        if (length(rn)>0) Br[,rn] <- rep(item_max[ri2[rn]]+1, each=N) - Br[,rn]
      }
      rand_f[, iter] <- rowCor_weighted(Ar, Br, weights)
    }

    rm_f <- rowMeans(rand_f, na.rm=TRUE)
    rsd_f <- apply(rand_f, 1, sd, na.rm=TRUE); rsd_f[rsd_f==0] <- 1e-10
    z_f <- (coupled_f - rm_f) / rsd_f

    factor_z_matrix[, f] <- z_f
    factors_used <- factors_used + 1
  }

  if (factors_used == 0) {
    return(list(mean_z = rep(NA, N), min_z = rep(NA, N),
                prop_low = rep(NA, N), sd_z = rep(NA, N),
                combined = rep(NA, N), global_z = rep(NA, N),
                nF_detected = nF_est, n_factors_used = 0))
  }

  # Step 4: Aggregate per-factor z-scores
  # Only use factors that were computed (non-NA columns)
  valid_cols <- which(colSums(!is.na(factor_z_matrix)) > 0)
  fz <- factor_z_matrix[, valid_cols, drop = FALSE]

  mean_z   <- rowMeans(fz, na.rm = TRUE)
  min_z    <- apply(fz, 1, min, na.rm = TRUE)
  prop_low <- rowMeans(fz < 1.0, na.rm = TRUE)  # proportion of factors below z=1
  sd_z     <- apply(fz, 1, sd, na.rm = TRUE)

  # Combined: mean_z with penalty for having any very low factor
  # If min_z is very low, it pulls the combined score down
  combined <- mean_z * 0.7 + min_z * 0.3

  # Step 5: Also compute global EFA-D for comparison
  k_global <- length(global_pairs_A)
  A_global <- mat[, global_pairs_A, drop=FALSE]
  B_global <- mat[, global_pairs_B, drop=FALSE]
  if (align_signs) {
    nf <- which(global_signs < 0)
    if (length(nf) > 0) {
      B_global[, nf] <- rep(item_max[global_pairs_B[nf]] + 1, each = N) - B_global[, nf]
    }
  }
  global_coupled <- rowCor_weighted(A_global, B_global, global_weights)

  global_rand <- matrix(NA_real_, N, iterations)
  for (iter in seq_len(iterations)) {
    ri1 <- sample(p, k_global, replace=TRUE); ri2 <- sample(p, k_global, replace=TRUE)
    same <- ri1==ri2; while(any(same)){ri2[same]<-sample(p,sum(same),replace=TRUE);same<-ri1==ri2}
    Ar <- mat[,ri1,drop=FALSE]; Br <- mat[,ri2,drop=FALSE]
    if (align_signs) {
      ri<-pmax(ri1,ri2);ci<-pmin(ri1,ri2);rs<-signMat[cbind(ri,ci)];rs[is.na(rs)]<-1
      rn<-which(rs<0); if(length(rn)>0) Br[,rn]<-rep(item_max[ri2[rn]]+1,each=N)-Br[,rn]
    }
    global_rand[,iter] <- rowCor_weighted(Ar, Br, global_weights)
  }
  grm <- rowMeans(global_rand,na.rm=TRUE); grsd <- apply(global_rand,1,sd,na.rm=TRUE)
  grsd[grsd==0] <- 1e-10
  global_z <- (global_coupled - grm) / grsd

  list(mean_z = mean_z, min_z = min_z, prop_low = prop_low, sd_z = sd_z,
       combined = combined, global_z = global_z,
       nF_detected = nF_est, n_factors_used = factors_used,
       factor_z_matrix = fz)
}

# ============================================================================
# SIMULATION
# ============================================================================

cat("\n============================================================\n")
cat("PER-FACTOR COHERENCE vs GLOBAL EFA-D\n")
cat("============================================================\n\n")

nF_levels  <- c(4, 6, 8, 10, 15, 20, 25, 30)
ipf_levels <- c(3, 6, 10)
N_val      <- 300
pct_levels <- c(0.05, 0.10, 0.25)
n_reps     <- 3
iterations <- 50

grid <- expand.grid(nF=nF_levels, ipf=ipf_levels, pct=pct_levels, rep=1:n_reps)
grid$total_items <- grid$nF * grid$ipf
total_cond <- nrow(grid)

cat(sprintf("Total conditions: %d\n\n", total_cond))

results <- data.frame()
t_start <- Sys.time()

for (i in seq_len(total_cond)) {
  nF <- grid$nF[i]; ipf <- grid$ipf[i]; pct <- grid$pct[i]; rep_i <- grid$rep[i]
  total_items <- grid$total_items[i]
  seed_val <- 60000 + i

  set.seed(seed_val)
  cat(sprintf("[%d/%d] nF=%d ipf=%d (%d items) pct=%.0f%% rep=%d ... ",
              i, total_cond, nF, ipf, total_items, pct*100, rep_i))

  good_data <- tryCatch(simulated_good_responses(nF, rep(ipf, nF), N_val),
                          error = function(e) { cat("SKIP\n"); NULL })
  if (is.null(good_data)) next

  inj <- inject_careless(good_data, pct, seed = seed_val + 10000)
  corrupted <- inj$data_corrupted
  labels_vec <- as.integer(inj$labels$careless)

  # Run per-factor
  pf <- tryCatch(ReReReRe_perFactor(corrupted, iterations=iterations),
                  error = function(e) { cat(sprintf("ERR: %s\n", e$message)); NULL })
  if (is.null(pf)) { cat("SKIP\n"); next }

  # Oracle best for each aggregation
  b_global   <- find_best_z(pf$global_z, labels_vec, seq(0.1, 3.0, 0.1))
  b_mean_z   <- find_best_z(pf$mean_z, labels_vec, seq(0.1, 3.0, 0.1))
  b_min_z    <- find_best_z(pf$min_z, labels_vec, seq(-2, 3.0, 0.1))
  b_combined <- find_best_z(pf$combined, labels_vec, seq(-1, 3.0, 0.1))
  b_prop     <- find_best_z(pf$prop_low, labels_vec, seq(0.1, 1.0, 0.05),
                              lower_is_bad = FALSE)  # HIGH prop = bad

  cat(sprintf("global=%.3f mean=%.3f min=%.3f comb=%.3f prop=%.3f (nF_det=%d, fUsed=%d)\n",
              b_global$mcc, b_mean_z$mcc, b_min_z$mcc, b_combined$mcc, b_prop$mcc,
              pf$nF_detected, pf$n_factors_used))

  results <- rbind(results, data.frame(
    nF=nF, ipf=ipf, total_items=total_items, pct=pct, rep=rep_i,
    mcc_global=b_global$mcc,
    mcc_mean_z=b_mean_z$mcc,
    mcc_min_z=b_min_z$mcc,
    mcc_combined=b_combined$mcc,
    mcc_prop=b_prop$mcc,
    nF_detected=pf$nF_detected,
    n_factors_used=pf$n_factors_used
  ))

  if (i %% 50 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), t_start, units="mins"))
    eta <- elapsed / i * (total_cond - i)
    cat(sprintf("  >> Checkpoint: %d/%d, %.1f min elapsed, ETA %.1f min\n",
                i, total_cond, elapsed, eta))
    write.csv(results, "test_perfactor_checkpoint.csv", row.names=FALSE)
  }
}

write.csv(results, "test_perfactor_results.csv", row.names=FALSE)
elapsed_total <- as.numeric(difftime(Sys.time(), t_start, units="mins"))

# ============================================================================
# ANALYSIS
# ============================================================================

results$item_bin <- cut(results$total_items, breaks=c(0,30,60,100,200,Inf),
                        labels=c("<30","30-60","60-100","100-200",">200"))

report_dir <- "archive/perfactor_comparison"
dir.create(report_dir, recursive=TRUE, showWarnings=FALSE)

res_long <- results %>%
  select(nF, ipf, total_items, pct, rep, item_bin,
         mcc_global, mcc_mean_z, mcc_min_z, mcc_combined, mcc_prop) %>%
  pivot_longer(cols=starts_with("mcc_"), names_to="method", values_to="mcc") %>%
  mutate(method = recode(method,
    mcc_global="Global EFA-D",
    mcc_mean_z="Per-factor: mean(z)",
    mcc_min_z="Per-factor: min(z)",
    mcc_combined="Per-factor: combined",
    mcc_prop="Per-factor: prop_low"))

colors5 <- c("Global EFA-D"="#E41A1C", "Per-factor: mean(z)"="#377EB8",
             "Per-factor: min(z)"="#4DAF4A", "Per-factor: combined"="#984EA3",
             "Per-factor: prop_low"="#FF7F00")

# Plot 1: Overall
p1 <- res_long %>%
  group_by(method) %>%
  summarise(mean_mcc=mean(mcc,na.rm=T), se=sd(mcc,na.rm=T)/sqrt(n()), .groups="drop") %>%
  ggplot(aes(x=reorder(method,-mean_mcc), y=mean_mcc, fill=method)) +
  geom_col(width=0.6) +
  geom_errorbar(aes(ymin=mean_mcc-se, ymax=mean_mcc+se), width=0.2) +
  geom_text(aes(label=sprintf("%.3f",mean_mcc)), vjust=-0.5, size=3.5) +
  labs(title="Overall Mean MCC: Per-Factor vs Global", x="", y="Mean MCC (oracle)") +
  theme_minimal(base_size=12) + theme(legend.position="none",
    axis.text.x=element_text(angle=20, hjust=1)) +
  scale_fill_manual(values=colors5) + ylim(0, NA)
ggsave(file.path(report_dir, "01_overall.png"), p1, width=10, height=5.5, dpi=150)

# Plot 2: By item bin
p2 <- res_long %>%
  group_by(item_bin, method) %>%
  summarise(mean_mcc=mean(mcc,na.rm=T), .groups="drop") %>%
  ggplot(aes(x=item_bin, y=mean_mcc, fill=method)) +
  geom_col(position=position_dodge(0.8), width=0.7) +
  labs(title="MCC by Item Range", x="Total items", y="Mean MCC (oracle)") +
  theme_minimal(base_size=12) + scale_fill_manual(values=colors5)
ggsave(file.path(report_dir, "02_by_item_bin.png"), p2, width=11, height=5.5, dpi=150)

# Plot 3: By nF
p3 <- res_long %>%
  group_by(nF, method) %>%
  summarise(mean_mcc=mean(mcc,na.rm=T), .groups="drop") %>%
  ggplot(aes(x=factor(nF), y=mean_mcc, color=method, group=method)) +
  geom_line(linewidth=1) + geom_point(size=2.5) +
  labs(title="MCC by nFactors", x="nFactors", y="Mean MCC (oracle)") +
  theme_minimal(base_size=12) + scale_color_manual(values=colors5)
ggsave(file.path(report_dir, "03_by_nF.png"), p3, width=10, height=5.5, dpi=150)

# Plot 4: By total items LOESS
p4 <- res_long %>%
  ggplot(aes(x=total_items, y=mcc, color=method)) +
  geom_point(alpha=0.1, size=0.8) +
  geom_smooth(method="loess", span=0.5, se=TRUE, linewidth=1) +
  labs(title="MCC by Total Items (LOESS)", x="Total items", y="MCC") +
  theme_minimal(base_size=12) + scale_color_manual(values=colors5)
ggsave(file.path(report_dir, "04_by_total_items.png"), p4, width=10, height=5.5, dpi=150)

# Plot 5: Heatmap nF x ipf for best per-factor method
best_pf <- results %>%
  group_by(nF, ipf) %>%
  summarise(global=mean(mcc_global,na.rm=T), mean_z=mean(mcc_mean_z,na.rm=T),
            min_z=mean(mcc_min_z,na.rm=T), combined=mean(mcc_combined,na.rm=T),
            prop=mean(mcc_prop,na.rm=T), .groups="drop") %>%
  rowwise() %>%
  mutate(best_pf = max(mean_z, min_z, combined, prop),
         gain = best_pf - global) %>%
  ungroup()
p5 <- ggplot(best_pf, aes(x=factor(ipf), y=factor(nF), fill=gain)) +
  geom_tile() +
  geom_text(aes(label=sprintf("%+.3f", gain)), size=3.5) +
  scale_fill_gradient2(low="#d73027", mid="white", high="#1a9850", midpoint=0) +
  labs(title="Best Per-Factor Method Gain over Global EFA-D",
       x="Items per factor", y="nFactors",
       subtitle="Green = per-factor wins") +
  theme_minimal(base_size=12)
ggsave(file.path(report_dir, "05_gain_heatmap.png"), p5, width=7, height=6, dpi=150)

# Plot 6: Advantage by total items
p6 <- results %>%
  mutate(best_pf = pmax(mcc_mean_z, mcc_min_z, mcc_combined, mcc_prop, na.rm=TRUE),
         gain = best_pf - mcc_global) %>%
  group_by(total_items) %>%
  summarise(mean_gain=mean(gain,na.rm=T), .groups="drop") %>%
  ggplot(aes(x=total_items, y=mean_gain)) +
  geom_col(aes(fill=mean_gain>0), width=8) +
  geom_hline(yintercept=0) +
  labs(title="Per-Factor Gain over Global by Total Items",
       x="Total items", y="Mean MCC gain") +
  theme_minimal(base_size=12) + theme(legend.position="none") +
  scale_fill_manual(values=c("TRUE"="#4DAF4A","FALSE"="#E41A1C"))
ggsave(file.path(report_dir, "06_gain_by_items.png"), p6, width=10, height=5, dpi=150)

# ============================================================================
# TEXT REPORT
# ============================================================================

sink(file.path(report_dir, "report.txt"))

cat("============================================================\n")
cat("PER-FACTOR COHERENCE vs GLOBAL EFA-D\n")
cat(sprintf("Date: %s | Conditions: %d | Runtime: %.1f min\n",
            Sys.Date(), total_cond, elapsed_total))
cat("============================================================\n\n")

cat("--- OVERALL MEANS ---\n\n")
overall <- res_long %>% group_by(method) %>%
  summarise(mean=mean(mcc,na.rm=T), sd=sd(mcc,na.rm=T), .groups="drop") %>%
  arrange(desc(mean))
print(as.data.frame(overall), digits=3)

cat("\n\n--- BY ITEM BIN ---\n\n")
bin_tab <- res_long %>%
  group_by(item_bin, method) %>%
  summarise(mcc=mean(mcc,na.rm=T), .groups="drop") %>%
  pivot_wider(names_from=method, values_from=mcc)
print(as.data.frame(bin_tab), digits=3)

cat("\n\n--- BY nF ---\n\n")
nf_tab <- res_long %>%
  group_by(nF, method) %>%
  summarise(mcc=mean(mcc,na.rm=T), .groups="drop") %>%
  pivot_wider(names_from=method, values_from=mcc)
print(as.data.frame(nf_tab), digits=3)

cat("\n\n--- GAIN: BEST PER-FACTOR vs GLOBAL (per nF x ipf) ---\n\n")
print(as.data.frame(best_pf %>% select(nF, ipf, global, best_pf, gain)), digits=3)

cat("\n\n--- OVERALL WINNER ---\n\n")
best_method <- overall$method[1]
best_mcc <- overall$mean[1]
global_mcc <- overall$mean[overall$method == "Global EFA-D"]
cat(sprintf("Best: %s (MCC=%.3f)\n", best_method, best_mcc))
cat(sprintf("Global EFA-D: %.3f\n", global_mcc))
cat(sprintf("Gain: %+.3f (%.1f%%)\n", best_mcc - global_mcc,
            (best_mcc - global_mcc)/global_mcc*100))

sink()

# Console summary
cat("\n\n============================================================\n")
cat("RESULTS SUMMARY\n")
cat("============================================================\n\n")

cat("--- OVERALL ---\n")
print(as.data.frame(overall), digits=3)

cat("\n--- BY ITEM BIN ---\n")
print(as.data.frame(bin_tab), digits=3)

cat(sprintf("\nPlots and report saved to %s/\n", report_dir))
cat("\n=== DONE ===\n")
