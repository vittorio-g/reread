## benchmark_study.R — head-to-head competitive benchmark on Study 1 (induced careless, n=157).
## rr (shipped index) and ReReRe (shipped ensemble) vs the standard careless indices,
## all on the SAME data, AUC + DeLong paired tests. Directions fixed a priori (high = careless).
suppressMessages({library(careless); library(pROC)})
set.seed(1)
DIR <- "."
mat0 <- as.matrix(read.csv(file.path(DIR,"_study_matrix.csv"), check.names=FALSE))
lab  <- read.csv(file.path(DIR,"_study_labels.csv"))[[1]]
it   <- read.csv(file.path(DIR,"_study_items.csv"), stringsAsFactors=FALSE)
sc   <- read.csv(file.path(DIR,"_study_scores.csv"))              # rr, eta from shipped engine
loo  <- read.csv(file.path(DIR,"_study_loo.csv"))                 # LOO ensemble predictions (honest, no fit optimism)
stopifnot(ncol(mat0)==nrow(it), nrow(mat0)==length(lab), nrow(sc)==length(lab))

## median-impute NA (a few) so every index is defined on all respondents
impute <- function(m){for(j in seq_len(ncol(m))){v<-m[,j];mm<-median(v,na.rm=TRUE);if(is.na(mm))mm<-0;m[is.na(v),j]<-mm};m}
mat <- impute(mat0)
## reverse-coded copy (for scale-structured indices): x_rev = (max+min) - x
matR <- mat
for(j in which(it$reverse==1)) matR[,j] <- (it$rmax[j]+it$rmin[j]) - mat[,j]

## coherent domains only (exclude the 10 Random decoy items) for even-odd / personal reliability
coh <- it$domain != "Random"
domains <- unique(it$domain[coh]); domains <- domains[domains!="Random"]

## ---- competitors (careless package) --------------------------------------
mahad_v <- tryCatch(as.numeric(careless::mahad(mat, plot=FALSE, flag=FALSE)), error=function(e){cat("mahad err:",conditionMessage(e),"\n");rep(NA,nrow(mat))})
irv_v   <- as.numeric(careless::irv(mat))                              # SD of responses
long_v  <- as.numeric(careless::longstring(mat))                       # max run
## even-odd: reverse-coded, columns ordered by coherent domain, factors = items per domain
ord <- unlist(lapply(domains, function(d) which(it$domain==d)))
fac <- as.numeric(table(factor(it$domain[ord], levels=domains)))
eo_v <- as.numeric(careless::evenodd(matR[,ord,drop=FALSE], factors=fac))
## psychometric synonyms / antonyms (raw data; package handles direction via correlation sign)
psyn_v <- tryCatch(as.numeric(careless::psychsyn(mat, critval=.60)),           error=function(e){rep(NA,nrow(mat))})
## antonyms: no pair reaches |r|>=.60 on these mixed-scale data; use -.50 (documented)
pant_v <- tryCatch(as.numeric(careless::psychsyn(mat, critval=-.50, anto=TRUE)),error=function(e){rep(NA,nrow(mat))})
cat("psychsyn pairs @.60; antonym pairs @-.50\n")

## ---- resampled personal reliability (Goldammer 2024 style) ---------------
## B random even/odd splits within each coherent domain; per split, per person,
## subscale mean on each half (z-standardized across persons within domain);
## within-person correlation across domains between the two halves; average over B.
B <- 200
prr_acc <- matrix(0, nrow(mat), B)
for(b in seq_len(B)){
  A <- matrix(NA, nrow(mat), length(domains)); Bm <- A
  for(di in seq_along(domains)){
    cols <- which(it$domain==domains[di]); cols <- sample(cols)
    h <- floor(length(cols)/2); if(h<1) next
    a <- cols[seq_len(h)]; bb <- cols[(h+1):length(cols)]
    ma <- rowMeans(matR[,a,drop=FALSE]); mb <- rowMeans(matR[,bb,drop=FALSE])
    A[,di]  <- (ma-mean(ma))/(sd(ma)+1e-9); Bm[,di] <- (mb-mean(mb))/(sd(mb)+1e-9)
  }
  prr_acc[,b] <- vapply(seq_len(nrow(mat)), function(i){
    x<-A[i,];y<-Bm[i,];ok<-is.finite(x)&is.finite(y)
    if(sum(ok)<3||sd(x[ok])==0||sd(y[ok])==0) return(NA_real_)
    suppressWarnings(cor(x[ok],y[ok]))}, numeric(1))
}
prr_v <- rowMeans(prr_acc, na.rm=TRUE)

## ---- assemble: careless-score = high means MORE careless (fixed a priori) --
## evenodd (careless>=1.2.0): higher already = more careless -> use as-is.
S <- list(
  "ReReRe (ensemble)"      = loo$Shipped,   # leave-one-out prediction, not full-sample fit
  "rr (index)"             = sc$rr,
  "Resampled pers. rel."   = -prr_v,
  "Psych synonyms"         = -psyn_v,
  "Person-Total corr."     = sc$person_total,
  "LongString"             = long_v,
  "Even-odd inconsistency" = eo_v,
  "IRV (low=careless)"     = -irv_v,
  "Mahalanobis D2"         = mahad_v,
  "Psych antonyms"         = -pant_v
)
cat("finite / careless-among-finite / attentive-among-finite per index:\n")
for(nmi in names(S)){ok<-is.finite(S[[nmi]]);cat(sprintf("  %-24s %3d  care=%2d att=%2d\n",nmi,sum(ok),sum(lab[ok]==1),sum(lab[ok]==0)))}
## AUC + 95% CI (direction fixed: cases[careless]=1 have HIGHER score => direction="<")
mkroc <- function(s){ok<-is.finite(s); if(sum(ok)<10||length(unique(lab[ok]))<2) return(NULL)
  roc(lab[ok], s[ok], direction="<", quiet=TRUE)}
rocs <- lapply(S, mkroc)
auc_ci <- t(sapply(rocs, function(r){ if(is.null(r)) return(c(NA,NA,NA))
  c(as.numeric(auc(r)), as.numeric(ci.auc(r, method="delong")))[c(1,2,4)]}))
colnames(auc_ci) <- c("AUC","lo95","hi95")

## DeLong paired tests vs rr and vs ensemble (paired on the cases where BOTH are defined)
pval <- function(a,b){
  ok <- is.finite(S[[a]]) & is.finite(S[[b]])
  if(sum(ok)<10 || length(unique(lab[ok]))<2) return(NA)
  r1<-roc(lab[ok],S[[a]][ok],direction="<",quiet=TRUE); r2<-roc(lab[ok],S[[b]][ok],direction="<",quiet=TRUE)
  tryCatch(roc.test(r1,r2,method="delong",paired=TRUE)$p.value, error=function(e)NA)
}
nm <- names(S)
p_vs_rr  <- sapply(nm, function(m) if(m=="rr (index)") NA else pval("rr (index)", m))
p_vs_ens <- sapply(nm, function(m) if(m=="ReReRe (ensemble)") NA else pval("ReReRe (ensemble)", m))

res <- data.frame(method=nm, AUC=round(auc_ci[,"AUC"],3), lo95=round(auc_ci[,"lo95"],3),
  hi95=round(auc_ci[,"hi95"],3), p_vs_rr=signif(p_vs_rr,3), p_vs_ensemble=signif(p_vs_ens,3),
  row.names=NULL)
res <- res[order(-res$AUC),]
cat("\n=== Study 1 competitive benchmark (n=",length(lab),", careless=",sum(lab),") ===\n",sep="")
print(res, row.names=FALSE)
write.csv(res, file.path(DIR,"benchmark_study.csv"), row.names=FALSE)
cat("\nwrote benchmark_study.csv\n")
