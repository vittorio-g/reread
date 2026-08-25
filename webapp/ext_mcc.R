## ext_mcc.R — add MCC alongside AUC for the external benchmark (Table 4/5).
## Oracle (best-threshold) MCC: the maximum MCC over all cutpoints of the score.
## This is threshold-free in the same sense AUC is (no calibration confound), so
## ens / rc / PT are compared on equal footing. Uses the ORIENTED dump
## (higher = more careless for every column) so directions need no flipping.
suppressMessages(library(pROC))
OUT <- "ext_bench_oriented"
idx <- read.csv(file.path("ext_bench","index.csv"), stringsAsFactors=FALSE)

mcc_at <- function(lab, pred){
  tp <- as.double(sum(pred==1 & lab==1)); tn <- as.double(sum(pred==0 & lab==0))
  fp <- as.double(sum(pred==1 & lab==0)); fn <- as.double(sum(pred==0 & lab==1))
  den <- sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  if(!is.finite(den) || den==0) return(0)
  (tp*tn - fp*fn)/den
}
## oracle MCC: sweep unique score values as cutoffs (subsample cutoffs if huge)
oracle_mcc <- function(lab, s){
  ok <- is.finite(s) & is.finite(lab); lab <- lab[ok]; s <- s[ok]
  if(length(unique(lab))<2 || length(s)<20) return(NA_real_)
  cuts <- unique(quantile(s, probs=seq(0.001,0.999,length.out=400), na.rm=TRUE))
  best <- -1
  for(c in cuts){ m <- mcc_at(lab, as.integer(s >= c)); if(is.finite(m) && m>best) best <- m }
  best
}
aucv <- function(lab,s){ ok<-is.finite(s); if(sum(ok)<20||length(unique(lab[ok]))<2) return(NA_real_)
  as.numeric(auc(roc(lab[ok], s[ok], direction="<", quiet=TRUE))) }

rows <- list()
for(k in seq_len(nrow(idx))){
  nm <- idx$name[k]
  f <- file.path(OUT, paste0(nm,"_y.csv"))
  if(!file.exists(f)){ cat("skip (no oriented dump):", nm, "\n"); next }
  Y <- read.csv(f)
  lab <- Y$y
  S <- list(ens=Y$eta, rc=Y$rr, PT=Y$person_total)
  au <- sapply(S, function(s) aucv(lab,s))
  mc <- sapply(S, function(s) oracle_mcc(lab,s))
  rows[[length(rows)+1]] <- data.frame(dataset=nm, n=nrow(Y), J=idx$J[k], rate=mean(lab),
    AUC_ens=au["ens"], AUC_rc=au["rc"], AUC_PT=au["PT"],
    MCC_ens=mc["ens"], MCC_rc=mc["rc"], MCC_PT=mc["PT"])
  cat(sprintf("%-14s J=%-3d rate=%.3f | AUC ens %.3f rc %.3f PT %.3f | MCC ens %.3f rc %.3f PT %.3f\n",
    nm, idx$J[k], mean(lab), au["ens"],au["rc"],au["PT"], mc["ens"],mc["rc"],mc["PT"]))
}
D <- do.call(rbind, rows)
write.csv(D, "ext_mcc.csv", row.names=FALSE)
cat("\n=== means across datasets ===\n")
cat(sprintf("AUC: ens %.3f  rc %.3f  PT %.3f\n", mean(D$AUC_ens,na.rm=TRUE), mean(D$AUC_rc,na.rm=TRUE), mean(D$AUC_PT,na.rm=TRUE)))
cat(sprintf("MCC: ens %.3f  rc %.3f  PT %.3f\n", mean(D$MCC_ens,na.rm=TRUE), mean(D$MCC_rc,na.rm=TRUE), mean(D$MCC_PT,na.rm=TRUE)))
cat("\nwrote ext_mcc.csv\n")
