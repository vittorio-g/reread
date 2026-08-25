## exp_efa_strat.R — does forcing factor-diversified pair selection beat plain top-|r|?
## Compare on Study 1: (1) baseline top-|r|; (2) oracle round-robin within TRUE domains;
## (3) EFA round-robin (parallel-analysis factors, max-loading assignment). AUC vs careless label.
suppressMessages({library(pROC); library(psych)}); set.seed(1)
E <- readRDS("exp_efa_env.rds"); ut<-E$ut; rr_<-E$rr_; k<-E$k; P<-E$P; minp<-E$minp; it<-E$it; lab<-E$lab; n<-E$n
NP <- nrow(ut)
indcor <- function(rows){
  aa<-ut[rows,1];bb<-ut[rows,2];sg<-sign(rr_[rows]);A<-P[,aa,drop=FALSE];B<-P[,bb,drop=FALSE]
  neg<-which(sg<0);for(c in neg)B[,c]<-(1+minp[bb[c]])-B[,c]
  vapply(1:n,function(i){av<-A[i,];bv<-B[i,];if(sd(av)==0||sd(bv)==0)return(0);abs(cor(av,bv))},numeric(1))
}
z_of <- function(rows,B=200){co<-indcor(rows);rand<-matrix(NA,n,B);for(bb in 1:B)rand[,bb]<-indcor(sample(NP,length(rows)))
  z<-(co-rowMeans(rand))/apply(rand,1,sd);z[!is.finite(z)]<-min(z[is.finite(z)]);z}
aucz <- function(rows){as.numeric(auc(roc(lab,-z_of(rows),direction="<",quiet=TRUE)))}

absr <- abs(rr_); ord <- order(absr,decreasing=TRUE)
## (1) baseline
sel_base <- ord[1:k]
## round-robin within-group selection given an item->factor vector `fac`
roundrobin <- function(fac){
  di <- fac[ut[,1]]; dj <- fac[ut[,2]]; within <- !is.na(di)&!is.na(dj)&(di==dj)
  groups <- split(which(within), di[within])
  groups <- lapply(groups, function(ix) ix[order(absr[ix],decreasing=TRUE)])  # sort each factor's pairs by |r|
  picked <- integer(0); pos <- setNames(rep(1L,length(groups)),names(groups)); gn<-names(groups)
  repeat{ any<-FALSE
    for(g in gn){ if(pos[g]<=length(groups[[g]])){ picked<-c(picked,groups[[g]][pos[g]]); pos[g]<-pos[g]+1L; any<-TRUE; if(length(picked)>=k) break } }
    if(length(picked)>=k || !any) break }
  picked[1:min(k,length(picked))]
}
## (2) oracle: true domains (drop the Random decoy factor)
facT <- it$domain; facT[facT=="Random"] <- NA
sel_oracle <- roundrobin(facT)
## (3) EFA: parallel analysis for nfactors, then max-loading assignment
np <- tryCatch(fa.parallel(P, fa="fa", plot=FALSE)$nfact, error=function(e) 10)
np <- max(3, min(np, 15))
fit <- tryCatch(fa(P, nfactors=np, rotate="oblimin", fm="minres", warnings=FALSE), error=function(e) NULL)
if(!is.null(fit)){ L<-abs(fit$loadings[,,drop=FALSE]); facE<-apply(L,1,which.max) } else facE<-rep(NA,ncol(P))
sel_efa <- if(!is.null(fit)) roundrobin(facE) else sel_base

cat("nfactors (parallel analysis):", np, "\n\n")
cat(sprintf("%-34s k=%3d  AUC=%.3f\n","(1) baseline top-|r|",       length(sel_base),   aucz(sel_base)))
cat(sprintf("%-34s k=%3d  AUC=%.3f\n","(2) oracle round-robin (true) ",length(sel_oracle), aucz(sel_oracle)))
cat(sprintf("%-34s k=%3d  AUC=%.3f\n","(3) EFA round-robin",         length(sel_efa),    aucz(sel_efa)))
# how many distinct factors each selection touches
ndist <- function(sel,fac){length(unique(na.omit(c(fac[ut[sel,1]],fac[ut[sel,2]]))))}
cat("\nendpoints span (distinct true domains): base",ndist(sel_base,facT),
    " oracle",ndist(sel_oracle,facT)," efa",ndist(sel_efa,facT),"\n")
