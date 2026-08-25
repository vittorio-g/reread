# =====================================================================
# DECISIVE TEST (Reviewer point 3): does rr add INCREMENTAL validity
# over {LongString + Person-Total} on EXTERNAL data, stratified by length?
# External analog of Fig 5 (which lives only in simulation).
#
# Per dataset: repeated stratified k-fold CV of two logistic models
#   M0: y ~ longstring + person_total
#   M1: y ~ longstring + person_total + rr
# Compare OOF AUCs (DeLong, paired) + stratified bootstrap CI on dAUC.
# Also report standalone AUC(rr, PT, LS) and shipped ens_p AUC for context.
# =====================================================================

suppressWarnings(suppressMessages({
  ok <- requireNamespace("pROC", quietly = TRUE)
}))
if (!ok) { install.packages("pROC", repos="https://cloud.r-project.org") }
suppressWarnings(suppressMessages(library(pROC)))

set.seed(42)
BASE <- "C:/Users/vitto/Desktop/ReReReRe/Dataset/gt_benchmark_candidates"

# dataset registry: name, items J, scores path, labels path
reg <- list(
  list(name="Kay S2",        J=363, s="kay_idris_idria/kay_scores_s2.csv", l="kay_idris_idria/kay_labels_s2.csv"),
  list(name="warning IPIP",  J=300, s="warning_ipipneo300/_scores.csv",     l="warning_ipipneo300/_labels.csv"),
  list(name="Kay S1",        J=250, s="kay_idris_idria/kay_scores_s1.csv", l="kay_idris_idria/kay_labels_s1.csv"),
  list(name="opsy 16PF",     J=163, s="opsy_16pf_scores.csv",              l="opsy_16pf_labels.csv"),
  list(name="smarvus",       J=141, s="smarvus/smarvus_scores.csv",        l="smarvus/smarvus_labels.csv"),
  list(name="Kay S6",        J=67,  s="kay_idris_idria/kay_scores_s6.csv", l="kay_idris_idria/kay_labels_s6.csv"),
  list(name="Kay S5",        J=62,  s="kay_idris_idria/kay_scores_s5.csv", l="kay_idris_idria/kay_labels_s5.csv"),
  list(name="douglas 2023",  J=50,  s="douglas2023_dataquality/douglas_scores.csv", l="douglas2023_dataquality/douglas_labels.csv"),
  list(name="duckworth VCL", J=50,  s="duckworth_grit_vcl/duckworth_scores.csv",    l="duckworth_grit_vcl/duckworth_labels.csv"),
  list(name="krause youth",  J=49,  s="krause_ier_youth/krause_scores.csv",         l="krause_ier_youth/krause_labels.csv")
)

auc1 <- function(y, x) {  # orientation-free detection AUC (>= .5)
  a <- as.numeric(pROC::auc(pROC::roc(y, x, quiet=TRUE, direction="auto")))
  a
}

# repeated stratified k-fold OOF predictions for a glm formula
oof_pred <- function(dat, form, K=10, reps=20) {
  n <- nrow(dat); acc <- numeric(n)
  for (rp in 1:reps) {
    # stratified folds
    idx0 <- which(dat$y==0); idx1 <- which(dat$y==1)
    f0 <- sample(rep(1:K, length.out=length(idx0)))
    f1 <- sample(rep(1:K, length.out=length(idx1)))
    fold <- integer(n); fold[idx0] <- f0; fold[idx1] <- f1
    pr <- numeric(n)
    for (k in 1:K) {
      tr <- fold!=k; te <- fold==k
      m <- suppressWarnings(glm(form, data=dat[tr,], family=binomial()))
      pr[te] <- predict(m, newdata=dat[te,], type="response")
    }
    acc <- acc + pr
  }
  acc/reps
}

rows <- list()
for (d in reg) {
  sc <- read.csv(file.path(BASE, d$s))
  lb <- read.csv(file.path(BASE, d$l))
  n <- min(nrow(sc), nrow(lb))
  sc <- sc[1:n,]; y <- as.integer(lb[[1]][1:n])
  df <- data.frame(y=y, rr=sc$rr_z, longstring=sc$longstring, person_total=sc$person_total,
                   ens=sc$ens_p)
  n_na_rr <- sum(is.na(df$rr))
  # impute rr NA (straight-liners) with median so M1 is not structurally penalised
  if (n_na_rr>0) df$rr[is.na(df$rr)] <- median(df$rr, na.rm=TRUE)
  df <- df[complete.cases(df[,c("y","rr","longstring","person_total")]),]

  # standalone detection AUCs
  a_rr <- auc1(df$y, df$rr); a_pt <- auc1(df$y, df$person_total)
  a_ls <- auc1(df$y, df$longstring); a_ens <- auc1(df$y, df$ens)

  # incremental CV test
  p0 <- oof_pred(df, y ~ longstring + person_total)
  p1 <- oof_pred(df, y ~ longstring + person_total + rr)
  r0 <- pROC::roc(df$y, p0, quiet=TRUE, direction="auto")
  r1 <- pROC::roc(df$y, p1, quiet=TRUE, direction="auto")
  A0 <- as.numeric(pROC::auc(r0)); A1 <- as.numeric(pROC::auc(r1))
  dtest <- pROC::roc.test(r1, r0, method="delong", paired=TRUE)
  # stratified bootstrap CI on dAUC (OOF preds fixed)
  B <- 2000; db <- numeric(B)
  i0 <- which(df$y==0); i1 <- which(df$y==1)
  for (b in 1:B) {
    bs <- c(sample(i0, replace=TRUE), sample(i1, replace=TRUE))
    db[b] <- auc1(df$y[bs], p1[bs]) - auc1(df$y[bs], p0[bs])
  }
  ci <- quantile(db, c(.025,.975))

  rows[[length(rows)+1]] <- data.frame(
    dataset=d$name, J=d$J, n=nrow(df), prev=round(mean(df$y),3), na_rr=n_na_rr,
    AUC_rr=round(a_rr,3), AUC_PT=round(a_pt,3), AUC_LS=round(a_ls,3), AUC_ens=round(a_ens,3),
    AUC_M0_LSPT=round(A0,3), AUC_M1_plus_rr=round(A1,3),
    dAUC=round(A1-A0,3), ci_lo=round(ci[1],3), ci_hi=round(ci[2],3),
    delong_p=signif(dtest$p.value,3)
  )
  cat(sprintf("%-14s J=%3d n=%5d prev=%.2f | rr=%.3f PT=%.3f | M0=%.3f M1=%.3f dAUC=%+.3f [%.3f,%.3f] p=%.3g\n",
      d$name, d$J, nrow(df), mean(df$y), a_rr, a_pt, A0, A1, A1-A0, ci[1], ci[2], dtest$p.value))
}

res <- do.call(rbind, rows)
write.csv(res, "C:/Users/vitto/Desktop/ReReReRe/rr_incremental_external.csv", row.names=FALSE)

cat("\n================ STRATUM SUMMARY ================\n")
res$stratum <- ifelse(res$J>=150, "long (>=150)", ifelse(res$J>=60, "medium (60-149)", "short (<60)"))
for (st in c("long (>=150)","medium (60-149)","short (<60)")) {
  sub <- res[res$stratum==st,]
  cat(sprintf("%-16s: mean dAUC=%+.3f | sig(p<.05): %d/%d | datasets: %s\n",
      st, mean(sub$dAUC), sum(sub$delong_p<.05), nrow(sub), paste(sub$dataset, collapse=", ")))
}
cat("\nSaved: rr_incremental_external.csv\n")
