## infer_study.R — inferential stats for Study 1 main results (reviewer point 9).
## LOO-CV logistic combiner (gate g=1 on this study, so a plain logistic reproduces the shipped
## ensemble) for the ablation models; AUC + 95% DeLong CI, MCC + 95% bootstrap CI, and DeLong
## paired tests for the ablation contrasts (does rr's null pay on the real study?). n=157.
suppressMessages({library(pROC)}); set.seed(1)
sc  <- read.csv("_study_scores.csv")                 # rr, longstring, person_total, irv, d2, eta
lab <- read.csv("_study_labels.csv")[[1]]
n <- length(lab)
feats <- list(
  Full        = c("rr","irv","longstring","person_total","d2"),
  Shipped     = c("rr","longstring","person_total"),
  NoRR        = c("irv","longstring","person_total","d2"),
  rr_alone    = c("rr")
)
loo <- function(cols){                                # leave-one-out logistic predicted prob
  p <- numeric(n); df <- data.frame(y=lab, sc[,cols,drop=FALSE])
  for(i in seq_len(n)){
    fit <- suppressWarnings(glm(y ~ ., data=df[-i,,drop=FALSE], family=binomial))
    p[i] <- predict(fit, newdata=df[i,,drop=FALSE], type="response")
  }
  p
}
P <- sapply(names(feats), function(k) loo(feats[[k]]))
colnames(P) <- names(feats)
write.csv(round(P,6), "_study_loo.csv", row.names=FALSE)

mccf <- function(pred,lab,thr){f<-pred>=thr;tp<-sum(f&lab==1);fp<-sum(f&lab==0);fn<-sum(!f&lab==1);tn<-sum(!f&lab==0)
  d<-sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));if(d==0)0 else (tp*tn-fp*fn)/d}
best_mcc <- function(pred,lab){th<-sort(unique(pred));m<-sapply(th,function(t)mccf(pred,lab,t));list(mcc=max(m),thr=th[which.max(m)])}
boot_mcc_ci <- function(pred,lab,thr,B=2000){v<-numeric(B);for(b in 1:B){ix<-sample(n,n,replace=TRUE);v[b]<-mccf(pred[ix],lab[ix],thr)};quantile(v,c(.025,.975))}

fpr5_thr <- function(pred,lab) as.numeric(quantile(pred[lab==0], 0.95))   # deployable FPR=5% operating point
rocs <- lapply(names(feats), function(k) roc(lab, P[,k], direction="<", quiet=TRUE)); names(rocs)<-names(feats)
cat("=== Study 1 ablation: LOO AUC (95% DeLong CI) + MCC oracle & @FPR5% (95% bootstrap CI), n=",n," ===\n",sep="")
tab <- data.frame()
for(k in names(feats)){
  ci <- ci.auc(rocs[[k]], method="delong"); a <- as.numeric(auc(rocs[[k]]))
  bm <- best_mcc(P[,k], lab); mc <- boot_mcc_ci(P[,k], lab, bm$thr)
  t5 <- fpr5_thr(P[,k],lab); m5 <- mccf(P[,k],lab,t5); c5 <- boot_mcc_ci(P[,k],lab,t5)
  tab <- rbind(tab, data.frame(model=k, AUC=round(a,3), AUC_lo=round(ci[1],3), AUC_hi=round(ci[3],3),
    MCC_or=round(bm$mcc,3), MCC_or_lo=round(mc[1],3), MCC_or_hi=round(mc[2],3),
    MCC_fpr5=round(m5,3), MCC_fpr5_lo=round(c5[1],3), MCC_fpr5_hi=round(c5[2],3)))
  cat(sprintf("  %-9s AUC=%.3f [%.3f,%.3f]  MCC_oracle=%.3f [%.3f,%.3f]  MCC@FPR5=%.3f [%.3f,%.3f]\n",
    k,a,ci[1],ci[3],bm$mcc,mc[1],mc[2],m5,c5[1],c5[2]))
}
write.csv(tab, "infer_study_ablation.csv", row.names=FALSE)
cat("\n=== DeLong paired tests (does the rr null pay on the real study?) ===\n")
pr <- function(a,b){t<-roc.test(rocs[[a]],rocs[[b]],method="delong",paired=TRUE)
  cat(sprintf("  %-9s vs %-9s : dAUC=%+.3f  DeLong p=%.3g\n",a,b,as.numeric(auc(rocs[[a]]))-as.numeric(auc(rocs[[b]])),t$p.value))}
pr("Full","NoRR"); pr("Shipped","NoRR"); pr("Full","rr_alone"); pr("Shipped","rr_alone")
cat("\nwrote _study_loo.csv, infer_study_ablation.csv\n")
