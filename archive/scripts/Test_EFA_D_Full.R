###############################################################################
# Comprehensive Comparison: EFA-D vs Standard vs Weighted
#
# Optimized design: fewer redundant levels, still wide coverage
#
# Design:
#   nF:  4, 6, 8, 10, 15, 20, 25, 30    (8 levels)
#   ipf: 3, 6, 10                         (3 levels)
#   N:   300                               (fixed — diminishing returns above)
#   pct: 5%, 10%, 25%                     (3 levels)
#   Reps: 5                                (more reps per cell for stability)
#
# Total: 8 x 3 x 3 x 5 = 360 conditions
# Estimated runtime: ~90-120 min
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

find_best_z <- function(z_scores, labels) {
  best_mcc <- -1; best_z <- NA
  for (z in seq(0.1, 3.0, 0.1)) {
    pred <- as.integer(z_scores <= z)
    mcc <- evaluate_mcc(pred, labels)
    if (mcc > best_mcc) { best_mcc <- mcc; best_z <- z }
  }
  return(list(z = best_z, mcc = best_mcc))
}

# ============================================================================
# EFA OPTION D: within-factor pairs + observed |r| weights
# ============================================================================

ReReReRe_efa_D <- function(data, iterations = 50, align_signs = TRUE,
                            nf_override = NULL) {
  data <- data[, sapply(data, is.numeric), drop = FALSE]
  mat <- as.matrix(data)
  N <- nrow(mat); p <- ncol(mat)
  item_max <- apply(mat, 2, max, na.rm = TRUE)

  # Step 1: Estimate nF
  nF_est <- nf_override
  if (is.null(nF_est)) {
    pa <- tryCatch({
      suppressMessages(suppressWarnings(
        fa.parallel(mat, fa = "fa", plot = FALSE, n.iter = 20)
      ))
    }, error = function(e) NULL)
    if (!is.null(pa) && !is.null(pa$nfact) && pa$nfact >= 1) {
      nF_est <- pa$nfact
    } else {
      ev <- eigen(cor(mat, use = "pairwise.complete.obs"),
                  symmetric = TRUE, only.values = TRUE)$values
      nF_est <- max(1, sum(ev > 1))
    }
  }
  nF_est <- max(1, min(nF_est, floor(p / 2)))

  # Step 2: EFA
  efa_result <- tryCatch({
    suppressWarnings(
      fa(mat, nfactors = nF_est, rotate = "oblimin", fm = "minres",
         scores = "none", warnings = FALSE)
    )
  }, error = function(e) {
    tryCatch({
      suppressWarnings(
        fa(mat, nfactors = max(1, nF_est - 1), rotate = "oblimin", fm = "minres",
           scores = "none", warnings = FALSE)
      )
    }, error = function(e2) NULL)
  })

  if (is.null(efa_result)) {
    return(data.frame(z_score = rep(NA_real_, N), indCors = rep(NA_real_, N),
                      rand_mean = rep(NA_real_, N), rand_sd = rep(NA_real_, N),
                      nF_detected = nF_est, n_pairs = NA_integer_))
  }

  # Step 3: Assign items to primary factor
  loadings_mat <- as.matrix(efa_result$loadings[])
  nF_actual <- ncol(loadings_mat)
  primary_factor <- apply(abs(loadings_mat), 1, which.max)

  # Step 4: Within-factor pairs
  pair_list <- list()
  for (f in seq_len(nF_actual)) {
    items_f <- which(primary_factor == f)
    if (length(items_f) < 2) next
    combos <- combn(items_f, 2)
    for (ci in seq_len(ncol(combos)))
      pair_list[[length(pair_list) + 1]] <- c(combos[1, ci], combos[2, ci])
  }
  if (length(pair_list) == 0) {
    return(data.frame(z_score = rep(NA_real_, N), indCors = rep(NA_real_, N),
                      rand_mean = rep(NA_real_, N), rand_sd = rep(NA_real_, N),
                      nF_detected = nF_est, n_pairs = 0L))
  }

  pair_mat <- do.call(rbind, pair_list)
  idx_A <- pair_mat[, 1]; idx_B <- pair_mat[, 2]

  # Step 5: Weight by observed |r|
  raw_cor_mat <- cor(mat, use = "pairwise.complete.obs")
  pair_r <- numeric(length(idx_A))
  for (j in seq_along(idx_A)) pair_r[j] <- raw_cor_mat[idx_A[j], idx_B[j]]
  pair_abs_r <- abs(pair_r)
  pair_sign <- sign(pair_r)

  keep <- !is.na(pair_abs_r)
  idx_A <- idx_A[keep]; idx_B <- idx_B[keep]
  pair_r <- pair_r[keep]; pair_abs_r <- pair_abs_r[keep]; pair_sign <- pair_sign[keep]
  k <- length(idx_A)
  weights <- pair_abs_r; weights[weights < 1e-6] <- 1e-6

  A_coupled <- mat[, idx_A, drop = FALSE]
  B_coupled <- mat[, idx_B, drop = FALSE]

  if (align_signs) {
    needs_flip <- which(pair_sign < 0)
    if (length(needs_flip) > 0) {
      flip_max <- item_max[idx_B[needs_flip]]
      B_coupled[, needs_flip] <- rep(flip_max + 1, each = N) - B_coupled[, needs_flip]
    }
  }

  indCors <- rowCor_weighted(A_coupled, B_coupled, weights)

  signMat <- sign(raw_cor_mat); signMat[upper.tri(signMat, diag = TRUE)] <- NA
  rand_scores <- matrix(NA_real_, N, iterations)
  for (iter in seq_len(iterations)) {
    ri1 <- sample(p, k, replace = TRUE); ri2 <- sample(p, k, replace = TRUE)
    same <- ri1 == ri2
    while (any(same)) { ri2[same] <- sample(p, sum(same), replace = TRUE); same <- ri1 == ri2 }
    Ar <- mat[, ri1, drop = FALSE]; Br <- mat[, ri2, drop = FALSE]
    if (align_signs) {
      ri <- pmax(ri1, ri2); ci <- pmin(ri1, ri2)
      rs <- signMat[cbind(ri, ci)]; rs[is.na(rs)] <- 1
      rn <- which(rs < 0)
      if (length(rn) > 0) Br[, rn] <- rep(item_max[ri2[rn]] + 1, each = N) - Br[, rn]
    }
    rand_scores[, iter] <- rowCor_weighted(Ar, Br, weights)
  }

  rm <- rowMeans(rand_scores, na.rm = TRUE)
  rsd <- apply(rand_scores, 1, sd, na.rm = TRUE); rsd[rsd == 0] <- 1e-10
  z_score <- (indCors - rm) / rsd

  data.frame(z_score = z_score, indCors = indCors,
             rand_mean = rm, rand_sd = rsd,
             nF_detected = nF_est, n_pairs = k)
}

# ============================================================================
# WEIGHTED (all pairs, |r| weights) — inline for comparison
# ============================================================================

ReReReRe_weighted_inline <- function(data, iterations = 50, align_signs = TRUE) {
  data <- data[, sapply(data, is.numeric), drop = FALSE]
  mat <- as.matrix(data); N <- nrow(mat); p <- ncol(mat)
  item_max <- apply(mat, 2, max, na.rm = TRUE)
  cor_matrix <- cor(mat, use = "pairwise.complete.obs")
  pair_idx <- which(lower.tri(cor_matrix), arr.ind = TRUE)
  pair_r <- cor_matrix[pair_idx]; pair_abs_r <- abs(pair_r)
  keep <- !is.na(pair_abs_r); pair_idx <- pair_idx[keep,]; pair_r <- pair_r[keep]
  pair_abs_r <- pair_abs_r[keep]; pair_sign <- sign(pair_r)
  k <- nrow(pair_idx); weights <- pair_abs_r
  A <- mat[, pair_idx[,1], drop=FALSE]; B <- mat[, pair_idx[,2], drop=FALSE]
  if (align_signs) {
    nf <- which(pair_sign < 0)
    if (length(nf) > 0) B[,nf] <- rep(item_max[pair_idx[nf,2]]+1, each=N) - B[,nf]
  }
  indCors <- rowCor_weighted(A, B, weights)
  signMat <- sign(cor_matrix); signMat[upper.tri(signMat, diag=TRUE)] <- NA
  rand_scores <- matrix(NA_real_, N, iterations)
  for (iter in seq_len(iterations)) {
    ri1 <- sample(p,k,replace=TRUE); ri2 <- sample(p,k,replace=TRUE)
    same <- ri1==ri2; while(any(same)){ri2[same]<-sample(p,sum(same),replace=TRUE);same<-ri1==ri2}
    Ar <- mat[,ri1,drop=FALSE]; Br <- mat[,ri2,drop=FALSE]
    if (align_signs) {
      ri<-pmax(ri1,ri2);ci<-pmin(ri1,ri2);rs<-signMat[cbind(ri,ci)];rs[is.na(rs)]<-1
      rn<-which(rs<0); if(length(rn)>0) Br[,rn]<-rep(item_max[ri2[rn]]+1,each=N)-Br[,rn]
    }
    rand_scores[,iter] <- rowCor_weighted(Ar, Br, weights)
  }
  rm <- rowMeans(rand_scores,na.rm=TRUE); rsd <- apply(rand_scores,1,sd,na.rm=TRUE)
  rsd[rsd==0]<-1e-10
  data.frame(z_score=(indCors-rm)/rsd, indCors=indCors, rand_mean=rm, rand_sd=rsd)
}

# ============================================================================
# SIMULATION
# ============================================================================

cat("\n============================================================\n")
cat("FULL COMPARISON: EFA-D vs Standard vs Weighted\n")
cat("============================================================\n\n")

nF_levels  <- c(4, 6, 8, 10, 15, 20, 25, 30)
ipf_levels <- c(3, 6, 10)
N_val      <- 300
pct_levels <- c(0.05, 0.10, 0.25)
n_reps     <- 5
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
  seed_val <- 50000 + i

  set.seed(seed_val)
  cat(sprintf("[%d/%d] nF=%d ipf=%d (%d items) pct=%.0f%% rep=%d ... ",
              i, total_cond, nF, ipf, total_items, pct*100, rep_i))

  good_data <- tryCatch(simulated_good_responses(nF, rep(ipf, nF), N_val),
                          error = function(e) { cat("SKIP\n"); NULL })
  if (is.null(good_data)) next

  inj <- inject_careless(good_data, pct, seed = seed_val + 10000)
  corrupted <- inj$data_corrupted
  labels_vec <- as.integer(inj$labels$careless)

  # Standard (coupled, corProp=0.03)
  r_std <- tryCatch(ReReReRe(corrupted, corProp=0.03, iterations=iterations,
                              align_signs=TRUE, mode="coupled"),
                     error=function(e) NULL)
  b_std <- if (!is.null(r_std)) find_best_z(r_std$z_score, labels_vec) else list(z=NA, mcc=NA)

  # Weighted (all pairs)
  r_wt <- tryCatch(ReReReRe_weighted_inline(corrupted, iterations=iterations),
                    error=function(e) NULL)
  b_wt <- if (!is.null(r_wt)) find_best_z(r_wt$z_score, labels_vec) else list(z=NA, mcc=NA)

  # EFA-D (within-factor + |r| weights)
  r_efaD <- tryCatch(ReReReRe_efa_D(corrupted, iterations=iterations),
                      error=function(e) NULL)
  b_efaD <- if (!is.null(r_efaD)) find_best_z(r_efaD$z_score, labels_vec) else list(z=NA, mcc=NA)

  nF_det <- if (!is.null(r_efaD)) r_efaD$nF_detected[1] else NA
  k_D    <- if (!is.null(r_efaD)) r_efaD$n_pairs[1] else NA

  cat(sprintf("std=%.3f wt=%.3f efaD=%.3f (nF_det=%s k=%s)\n",
              ifelse(is.na(b_std$mcc),0,b_std$mcc),
              ifelse(is.na(b_wt$mcc),0,b_wt$mcc),
              ifelse(is.na(b_efaD$mcc),0,b_efaD$mcc),
              ifelse(is.na(nF_det),"?",nF_det),
              ifelse(is.na(k_D),"?",k_D)))

  results <- rbind(results, data.frame(
    nF=nF, ipf=ipf, total_items=total_items, pct=pct, rep=rep_i,
    mcc_std=b_std$mcc, z_std=b_std$z,
    mcc_wt=b_wt$mcc, z_wt=b_wt$z,
    mcc_efaD=b_efaD$mcc, z_efaD=b_efaD$z,
    nF_detected=nF_det, n_pairs_D=k_D
  ))

  # Checkpoint every 50
  if (i %% 50 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), t_start, units="mins"))
    eta <- elapsed / i * (total_cond - i)
    cat(sprintf("  >> Checkpoint: %d/%d done, %.1f min elapsed, ETA %.1f min\n", i, total_cond, elapsed, eta))
    write.csv(results, "efa_d_full_comparison_results.csv", row.names=FALSE)
  }
}

write.csv(results, "efa_d_full_comparison_results.csv", row.names=FALSE)
elapsed_total <- as.numeric(difftime(Sys.time(), t_start, units="mins"))
cat(sprintf("\n=== Simulation complete: %d conditions in %.1f minutes ===\n\n", total_cond, elapsed_total))

# ============================================================================
# ANALYSIS & REPORT
# ============================================================================

results$item_bin <- cut(results$total_items, breaks=c(0,30,60,100,200,Inf),
                        labels=c("<30","30-60","60-100","100-200",">200"))

# Reshape for ggplot
res_long <- results %>%
  select(nF, ipf, total_items, pct, rep, item_bin, mcc_std, mcc_wt, mcc_efaD) %>%
  pivot_longer(cols=c(mcc_std, mcc_wt, mcc_efaD), names_to="method", values_to="mcc") %>%
  mutate(method = recode(method, mcc_std="Standard (coupled)",
                                  mcc_wt="Weighted (all pairs)",
                                  mcc_efaD="EFA-D (within-factor + |r|)"))

# Output directory
report_dir <- "archive/efa_d_comparison"
dir.create(report_dir, recursive=TRUE, showWarnings=FALSE)

# ========== PLOT 01: Overall by method ==========
p1 <- res_long %>%
  group_by(method) %>%
  summarise(mean_mcc = mean(mcc, na.rm=TRUE),
            se = sd(mcc, na.rm=TRUE)/sqrt(n()), .groups="drop") %>%
  ggplot(aes(x=reorder(method, -mean_mcc), y=mean_mcc, fill=method)) +
  geom_col(width=0.6) +
  geom_errorbar(aes(ymin=mean_mcc-se, ymax=mean_mcc+se), width=0.2) +
  geom_text(aes(label=sprintf("%.3f", mean_mcc)), vjust=-0.5, size=4) +
  labs(title="Overall Mean MCC by Method", x="", y="Mean MCC (oracle best z)") +
  theme_minimal(base_size=13) + theme(legend.position="none") +
  scale_fill_manual(values=c("Standard (coupled)"="#E41A1C",
                              "Weighted (all pairs)"="#377EB8",
                              "EFA-D (within-factor + |r|)"="#4DAF4A")) +
  ylim(0, NA)
ggsave(file.path(report_dir, "01_overall_by_method.png"), p1, width=8, height=5, dpi=150)

# ========== PLOT 02: By item bin ==========
p2 <- res_long %>%
  group_by(item_bin, method) %>%
  summarise(mean_mcc = mean(mcc, na.rm=TRUE), .groups="drop") %>%
  ggplot(aes(x=item_bin, y=mean_mcc, fill=method)) +
  geom_col(position=position_dodge(0.8), width=0.7) +
  geom_text(aes(label=sprintf("%.3f", mean_mcc)), position=position_dodge(0.8),
            vjust=-0.3, size=3) +
  labs(title="MCC by Item Range", x="Total items", y="Mean MCC (oracle)") +
  theme_minimal(base_size=13) +
  scale_fill_manual(values=c("Standard (coupled)"="#E41A1C",
                              "Weighted (all pairs)"="#377EB8",
                              "EFA-D (within-factor + |r|)"="#4DAF4A"))
ggsave(file.path(report_dir, "02_by_item_bin.png"), p2, width=10, height=5.5, dpi=150)

# ========== PLOT 03: By nF ==========
p3 <- res_long %>%
  group_by(nF, method) %>%
  summarise(mean_mcc = mean(mcc, na.rm=TRUE), .groups="drop") %>%
  ggplot(aes(x=factor(nF), y=mean_mcc, color=method, group=method)) +
  geom_line(linewidth=1.2) + geom_point(size=3) +
  labs(title="MCC by Number of Factors", x="nFactors", y="Mean MCC (oracle)") +
  theme_minimal(base_size=13) +
  scale_color_manual(values=c("Standard (coupled)"="#E41A1C",
                               "Weighted (all pairs)"="#377EB8",
                               "EFA-D (within-factor + |r|)"="#4DAF4A"))
ggsave(file.path(report_dir, "03_by_nF.png"), p3, width=10, height=5.5, dpi=150)

# ========== PLOT 04: By total_items (scatter + LOESS) ==========
p4 <- res_long %>%
  ggplot(aes(x=total_items, y=mcc, color=method)) +
  geom_point(alpha=0.15, size=1) +
  geom_smooth(method="loess", span=0.5, se=TRUE, linewidth=1.2) +
  labs(title="MCC by Total Items (scatter + LOESS)", x="Total items", y="MCC (oracle)") +
  theme_minimal(base_size=13) +
  scale_color_manual(values=c("Standard (coupled)"="#E41A1C",
                               "Weighted (all pairs)"="#377EB8",
                               "EFA-D (within-factor + |r|)"="#4DAF4A"))
ggsave(file.path(report_dir, "04_by_total_items_loess.png"), p4, width=10, height=5.5, dpi=150)

# ========== PLOT 05: Heatmap nF x ipf — EFA-D ==========
heat_efaD <- results %>%
  group_by(nF, ipf) %>%
  summarise(mcc = mean(mcc_efaD, na.rm=TRUE), .groups="drop")
p5 <- ggplot(heat_efaD, aes(x=factor(ipf), y=factor(nF), fill=mcc)) +
  geom_tile() +
  geom_text(aes(label=sprintf("%.3f", mcc)), size=3.5) +
  scale_fill_gradient2(low="#d73027", mid="#ffffbf", high="#1a9850", midpoint=0.3) +
  labs(title="EFA-D: MCC Heatmap (nF x items/factor)", x="Items per factor", y="nFactors") +
  theme_minimal(base_size=13)
ggsave(file.path(report_dir, "05_heatmap_nF_ipf_efaD.png"), p5, width=7, height=6, dpi=150)

# ========== PLOT 06: Heatmap nF x ipf — Difference (EFA-D minus Standard) ==========
heat_diff <- results %>%
  group_by(nF, ipf) %>%
  summarise(diff = mean(mcc_efaD - mcc_std, na.rm=TRUE), .groups="drop")
p6 <- ggplot(heat_diff, aes(x=factor(ipf), y=factor(nF), fill=diff)) +
  geom_tile() +
  geom_text(aes(label=sprintf("%+.3f", diff)), size=3.5) +
  scale_fill_gradient2(low="#d73027", mid="white", high="#1a9850", midpoint=0,
                       limits=c(-0.15, 0.15)) +
  labs(title="EFA-D minus Standard: MCC Difference", x="Items per factor", y="nFactors",
       subtitle="Green = EFA-D wins, Red = Standard wins") +
  theme_minimal(base_size=13)
ggsave(file.path(report_dir, "06_heatmap_diff_efaD_vs_std.png"), p6, width=7, height=6, dpi=150)

# ========== PLOT 07: Heatmap — Difference (EFA-D minus Weighted) ==========
heat_diff2 <- results %>%
  group_by(nF, ipf) %>%
  summarise(diff = mean(mcc_efaD - mcc_wt, na.rm=TRUE), .groups="drop")
p7 <- ggplot(heat_diff2, aes(x=factor(ipf), y=factor(nF), fill=diff)) +
  geom_tile() +
  geom_text(aes(label=sprintf("%+.3f", diff)), size=3.5) +
  scale_fill_gradient2(low="#d73027", mid="white", high="#1a9850", midpoint=0,
                       limits=c(-0.15, 0.15)) +
  labs(title="EFA-D minus Weighted: MCC Difference", x="Items per factor", y="nFactors",
       subtitle="Green = EFA-D wins, Red = Weighted wins") +
  theme_minimal(base_size=13)
ggsave(file.path(report_dir, "07_heatmap_diff_efaD_vs_wt.png"), p7, width=7, height=6, dpi=150)

# ========== PLOT 08: Winner map — which method wins in each cell ==========
winner_map <- results %>%
  group_by(nF, ipf) %>%
  summarise(std = mean(mcc_std, na.rm=TRUE),
            wt = mean(mcc_wt, na.rm=TRUE),
            efaD = mean(mcc_efaD, na.rm=TRUE), .groups="drop") %>%
  rowwise() %>%
  mutate(winner = c("Standard","Weighted","EFA-D")[which.max(c(std, wt, efaD))],
         margin = max(c(std, wt, efaD)) - sort(c(std, wt, efaD))[2]) %>%
  ungroup()
p8 <- ggplot(winner_map, aes(x=factor(ipf), y=factor(nF), fill=winner)) +
  geom_tile() +
  geom_text(aes(label=sprintf("%s\n(+%.3f)", winner, margin)), size=3) +
  scale_fill_manual(values=c("Standard"="#E41A1C", "Weighted"="#377EB8", "EFA-D"="#4DAF4A")) +
  labs(title="Winner per Cell (nF x ipf)", x="Items per factor", y="nFactors",
       subtitle="Margin = winner MCC minus runner-up") +
  theme_minimal(base_size=13)
ggsave(file.path(report_dir, "08_winner_map.png"), p8, width=8, height=6, dpi=150)

# ========== PLOT 09: By pct_careless ==========
p9 <- res_long %>%
  group_by(pct, method) %>%
  summarise(mean_mcc = mean(mcc, na.rm=TRUE), .groups="drop") %>%
  ggplot(aes(x=factor(pct*100), y=mean_mcc, fill=method)) +
  geom_col(position=position_dodge(0.8), width=0.7) +
  labs(title="MCC by Careless Rate", x="% Careless", y="Mean MCC (oracle)") +
  theme_minimal(base_size=13) +
  scale_fill_manual(values=c("Standard (coupled)"="#E41A1C",
                              "Weighted (all pairs)"="#377EB8",
                              "EFA-D (within-factor + |r|)"="#4DAF4A"))
ggsave(file.path(report_dir, "09_by_pct_careless.png"), p9, width=9, height=5, dpi=150)

# ========== PLOT 10: Boxplots by nF per method ==========
p10 <- res_long %>%
  ggplot(aes(x=factor(nF), y=mcc, fill=method)) +
  geom_boxplot(position=position_dodge(0.8), outlier.size=0.5) +
  labs(title="MCC Distribution by nF", x="nFactors", y="MCC (oracle)") +
  theme_minimal(base_size=13) +
  scale_fill_manual(values=c("Standard (coupled)"="#E41A1C",
                              "Weighted (all pairs)"="#377EB8",
                              "EFA-D (within-factor + |r|)"="#4DAF4A"))
ggsave(file.path(report_dir, "10_boxplots_by_nF.png"), p10, width=12, height=5.5, dpi=150)

# ========== PLOT 11: EFA-D advantage over BEST of other two ==========
results$best_other <- pmax(results$mcc_std, results$mcc_wt, na.rm=TRUE)
results$efaD_adv <- results$mcc_efaD - results$best_other

p11 <- results %>%
  group_by(total_items) %>%
  summarise(mean_adv = mean(efaD_adv, na.rm=TRUE),
            pct_wins = mean(efaD_adv > 0, na.rm=TRUE)*100, .groups="drop") %>%
  ggplot(aes(x=total_items, y=mean_adv)) +
  geom_col(aes(fill=mean_adv > 0), width=8) +
  geom_hline(yintercept=0, linewidth=0.5) +
  labs(title="EFA-D Advantage over Best of (Standard, Weighted)",
       x="Total items", y="Mean MCC advantage",
       subtitle="Positive = EFA-D wins, Negative = another method wins") +
  theme_minimal(base_size=13) + theme(legend.position="none") +
  scale_fill_manual(values=c("TRUE"="#4DAF4A", "FALSE"="#E41A1C"))
ggsave(file.path(report_dir, "11_efaD_advantage.png"), p11, width=10, height=5, dpi=150)

# ========== PLOT 12: nF detected vs true ==========
nf_accuracy <- results %>%
  group_by(nF) %>%
  summarise(mean_det = mean(nF_detected, na.rm=TRUE),
            sd_det = sd(nF_detected, na.rm=TRUE),
            ratio = mean(nF_detected/nF, na.rm=TRUE), .groups="drop")
p12 <- ggplot(nf_accuracy, aes(x=nF, y=mean_det)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
  geom_point(size=4, color="#4DAF4A") +
  geom_errorbar(aes(ymin=mean_det-sd_det, ymax=mean_det+sd_det), width=0.5, color="#4DAF4A") +
  labs(title="Parallel Analysis: Detected vs True nF",
       x="True nFactors", y="Detected nFactors",
       subtitle="Dashed = perfect detection") +
  theme_minimal(base_size=13) +
  coord_equal(xlim=c(0,35), ylim=c(0,35))
ggsave(file.path(report_dir, "12_nF_detected_vs_true.png"), p12, width=7, height=7, dpi=150)

# ============================================================================
# TEXT REPORT
# ============================================================================

sink(file.path(report_dir, "report.txt"))

cat("============================================================\n")
cat("COMPREHENSIVE REPORT: EFA-D vs Standard vs Weighted\n")
cat(sprintf("Date: %s\n", Sys.Date()))
cat(sprintf("Total conditions: %d | Runtime: %.1f minutes\n", total_cond, elapsed_total))
cat("============================================================\n\n")

cat("DESIGN:\n")
cat(sprintf("  nF: %s\n", paste(nF_levels, collapse=", ")))
cat(sprintf("  ipf: %s\n", paste(ipf_levels, collapse=", ")))
cat(sprintf("  N: %d (fixed)\n", N_val))
cat(sprintf("  pct_careless: %s\n", paste(pct_levels*100, collapse="%, ")))
cat(sprintf("  Reps: %d per cell\n", n_reps))
cat(sprintf("  Total items range: %d - %d\n", min(grid$total_items), max(grid$total_items)))
cat(sprintf("  iterations: %d per RR call\n", iterations))

cat("\n\n--- 1. OVERALL MEANS ---\n\n")
cat(sprintf("Standard (coupled, corProp=0.03): %.3f\n", mean(results$mcc_std, na.rm=TRUE)))
cat(sprintf("Weighted (all pairs, |r| wt):     %.3f\n", mean(results$mcc_wt, na.rm=TRUE)))
cat(sprintf("EFA-D (within-factor + |r| wt):   %.3f\n", mean(results$mcc_efaD, na.rm=TRUE)))

cat("\n\n--- 2. BY ITEM BIN ---\n\n")
bin_tab <- results %>%
  group_by(item_bin) %>%
  summarise(n=n(), std=mean(mcc_std,na.rm=T), wt=mean(mcc_wt,na.rm=T),
            efaD=mean(mcc_efaD,na.rm=T),
            best = c("Std","Wt","EFA-D")[which.max(c(mean(mcc_std,na.rm=T),
                                                      mean(mcc_wt,na.rm=T),
                                                      mean(mcc_efaD,na.rm=T)))],
            .groups="drop")
print(as.data.frame(bin_tab), digits=3)

cat("\n\n--- 3. BY nF ---\n\n")
nf_tab <- results %>%
  group_by(nF) %>%
  summarise(std=mean(mcc_std,na.rm=T), wt=mean(mcc_wt,na.rm=T),
            efaD=mean(mcc_efaD,na.rm=T), .groups="drop")
print(as.data.frame(nf_tab), digits=3)

cat("\n\n--- 4. BY nF x ipf (EFA-D MCC) ---\n\n")
cell_tab <- results %>%
  group_by(nF, ipf, total_items) %>%
  summarise(std=mean(mcc_std,na.rm=T), wt=mean(mcc_wt,na.rm=T),
            efaD=mean(mcc_efaD,na.rm=T),
            nF_det=mean(nF_detected,na.rm=T), .groups="drop") %>%
  arrange(total_items)
print(as.data.frame(cell_tab), digits=3)

cat("\n\n--- 5. EFA-D vs BEST OTHER ---\n\n")
cat(sprintf("EFA-D wins over best(Std,Wt): %.1f%% of cells\n",
            mean(results$efaD_adv > 0, na.rm=TRUE)*100))
cat(sprintf("Mean advantage when EFA-D wins: +%.3f\n",
            mean(results$efaD_adv[results$efaD_adv > 0], na.rm=TRUE)))
cat(sprintf("Mean disadvantage when EFA-D loses: %.3f\n",
            mean(results$efaD_adv[results$efaD_adv <= 0], na.rm=TRUE)))

by_bin_wins <- results %>%
  group_by(item_bin) %>%
  summarise(pct_efaD_wins = mean(efaD_adv > 0, na.rm=TRUE)*100,
            mean_adv = mean(efaD_adv, na.rm=TRUE), .groups="drop")
cat("\nBy item bin:\n")
print(as.data.frame(by_bin_wins), digits=3)

cat("\n\n--- 6. PARALLEL ANALYSIS ACCURACY ---\n\n")
cat("nF detected vs true:\n")
print(as.data.frame(nf_accuracy), digits=2)
cat(sprintf("\nOverall mean ratio (detected/true): %.2f\n",
            mean(results$nF_detected/results$nF, na.rm=TRUE)))

cat("\n\n--- 7. WINNER MAP (nF x ipf) ---\n\n")
print(as.data.frame(winner_map %>% select(nF, ipf, std, wt, efaD, winner, margin)), digits=3)

cat("\n\n--- 8. BY pct_careless ---\n\n")
pct_tab <- results %>%
  group_by(pct) %>%
  summarise(std=mean(mcc_std,na.rm=T), wt=mean(mcc_wt,na.rm=T),
            efaD=mean(mcc_efaD,na.rm=T), .groups="drop")
print(as.data.frame(pct_tab), digits=3)

cat("\n\n============================================================\n")
cat("KEY FINDING:\n")
cat("============================================================\n\n")

overall_std <- mean(results$mcc_std, na.rm=TRUE)
overall_wt <- mean(results$mcc_wt, na.rm=TRUE)
overall_efaD <- mean(results$mcc_efaD, na.rm=TRUE)
efaD_win_pct <- mean(results$efaD_adv > 0, na.rm=TRUE)*100

if (overall_efaD >= overall_std && overall_efaD >= overall_wt) {
  cat(sprintf("EFA-D is the best overall method (MCC=%.3f vs Std=%.3f, Wt=%.3f)\n",
              overall_efaD, overall_std, overall_wt))
} else {
  best_name <- c("Standard","Weighted","EFA-D")[which.max(c(overall_std, overall_wt, overall_efaD))]
  cat(sprintf("Best overall: %s\n", best_name))
}
cat(sprintf("EFA-D wins %.1f%% of individual cells vs best alternative\n", efaD_win_pct))

sink()

cat("\n============================================================\n")
cat("REPORT COMPLETE\n")
cat(sprintf("  Results: efa_d_full_comparison_results.csv\n"))
cat(sprintf("  Report:  %s/report.txt\n", report_dir))
cat(sprintf("  Plots:   %s/ (12 plots)\n", report_dir))
cat("============================================================\n")

# Print key results to console too
cat("\n--- OVERALL MEANS ---\n")
cat(sprintf("Standard:  %.3f\n", mean(results$mcc_std, na.rm=TRUE)))
cat(sprintf("Weighted:  %.3f\n", mean(results$mcc_wt, na.rm=TRUE)))
cat(sprintf("EFA-D:     %.3f\n", mean(results$mcc_efaD, na.rm=TRUE)))

cat("\n--- BY ITEM BIN ---\n")
print(as.data.frame(bin_tab), digits=3)

cat("\n--- WINNER MAP ---\n")
print(as.data.frame(winner_map %>% select(nF, ipf, winner, margin)), digits=3)

cat("\n=== DONE ===\n")
