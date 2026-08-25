## sim_bench.R — Study-3 competitive benchmark on simulated data (full factor structure known,
## so ALL competitors incl. even-odd & resampled personal reliability). AUC vs the injected
## careless label, aggregated by questionnaire length. Shows whether rr's per-respondent null
## pays (and how it scales) against the established indices. Directions fixed a priori.
suppressMessages({library(careless); library(pROC)})
set.seed(1)
OUT <- "sim_bench"
idx <- read.csv(file.path(OUT,"index.csv"), stringsAsFactors=FALSE)
aucv <- function(lab,s){ok<-is.finite(s); if(sum(ok)<20||length(unique(lab[ok]))<2) return(NA)
  as.numeric(auc(roc(lab[ok],s[ok],direction="<",quiet=TRUE)))}
persrel <- function(matR, fac, B=100){
  doms<-sort(unique(fac)); n<-nrow(matR); acc<-matrix(NA,n,B)
  for(b in seq_len(B)){A<-matrix(NA,n,length(doms));Bm<-A
    for(di in seq_along(doms)){cols<-which(fac==doms[di]);cols<-sample(cols);h<-floor(length(cols)/2);if(h<1)next
      a<-cols[seq_len(h)];bb<-cols[(h+1):length(cols)]
      ma<-rowMeans(matR[,a,drop=FALSE]);mb<-rowMeans(matR[,bb,drop=FALSE])
      A[,di]<-(ma-mean(ma))/(sd(ma)+1e-9);Bm[,di]<-(mb-mean(mb))/(sd(mb)+1e-9)}
    acc[,b]<-vapply(seq_len(n),function(i){x<-A[i,];y<-Bm[i,];ok<-is.finite(x)&is.finite(y)
      if(sum(ok)<3||sd(x[ok])==0||sd(y[ok])==0)return(NA_real_);suppressWarnings(cor(x[ok],y[ok]))},numeric(1))}
  rowMeans(acc,na.rm=TRUE)}

METH <- c("ReReRe","rr","PersRel","PsychSyn","EvenOdd","Mahalanobis","IRV","LongString")
rec <- list()
for(k in seq_len(nrow(idx))){
  tag<-idx$file[k]; J<-idx$items[k]
  mat <- as.matrix(read.csv(file.path(OUT,paste0(tag,"_M.csv")), header=FALSE))
  Y   <- read.csv(file.path(OUT,paste0(tag,"_y.csv")))
  meta<- read.csv(file.path(OUT,paste0(tag,"_meta.csv")))
  lab <- Y$y; fac<-meta$factor
  matR <- mat; rev<-which(meta$reverse==1); if(length(rev)) matR[,rev] <- 6 - mat[,rev]  # K=5 -> 6-x
  facN <- as.numeric(table(factor(fac, levels=sort(unique(fac)))))
  eo <- tryCatch(as.numeric(careless::evenodd(matR, factors=facN)), error=function(e)rep(NA,nrow(mat)))
  ps <- tryCatch(as.numeric(careless::psychsyn(mat, critval=.40)), error=function(e)rep(NA,nrow(mat)))  # .40: 5-pt sim r rarely reaches .60
  pr <- persrel(matR, fac)
  S <- list(ReReRe=Y$eta, rr=Y$rr, PersRel=-pr, PsychSyn=-ps, EvenOdd=eo,
            Mahalanobis=as.numeric(careless::mahad(mat,plot=FALSE,flag=FALSE)),
            IRV=-as.numeric(careless::irv(mat)), LongString=as.numeric(careless::longstring(mat)))
  a <- sapply(METH, function(m) aucv(lab, S[[m]]))
  rec[[k]] <- data.frame(items=J, rep=idx$rep[k], t(a))
  if(idx$rep[k]==0) cat(sprintf("J=%-3d | ReReRe %.3f rr %.3f PersRel %.3f PsychSyn %.3f EvenOdd %.3f Mah %.3f IRV %.3f Long %.3f\n",
    J,a["ReReRe"],a["rr"],a["PersRel"],a["PsychSyn"],a["EvenOdd"],a["Mahalanobis"],a["IRV"],a["LongString"]))
}
D <- do.call(rbind, rec)
write.csv(D, file.path(OUT,"sim_bench_raw.csv"), row.names=FALSE)
## aggregate: mean AUC per method per length
ag <- aggregate(D[,METH], by=list(items=D$items), FUN=function(x) mean(x,na.rm=TRUE))
cat("\n=== mean AUC by questionnaire length (",max(idx$rep)+1," reps/size) ===\n",sep="")
print(format(ag, digits=3), row.names=FALSE)
write.csv(ag, file.path(OUT,"sim_bench_meanauc.csv"), row.names=FALSE)
## key deltas: rr - best consistency competitor, and ensemble - each, per size with t-CI over reps
cat("\n=== delta AUC (rr - PersRel) and (rr - PsychSyn) by length, mean [95% t-CI over reps] ===\n")
for(J in sort(unique(D$items))){
  s<-D[D$items==J,]
  for(cmp in c("PersRel","PsychSyn","EvenOdd","Mahalanobis")){
    d<-s$rr - s[[cmp]]; m<-mean(d,na.rm=TRUE); se<-sd(d,na.rm=TRUE)/sqrt(sum(is.finite(d)))
    cat(sprintf("  J=%-3d rr-%-11s %+.3f [%+.3f,%+.3f]\n",J,cmp,m,m-1.96*se,m+1.96*se))
  }
}
cat("\nwrote sim_bench_raw.csv, sim_bench_meanauc.csv\n")
