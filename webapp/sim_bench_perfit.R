## sim_bench_perfit.R — person-fit (l_z, U3) head-to-head vs rc on Study-3 simulation.
## U3poly (fast, nonparametric) on ALL datasets; lzpoly (fits GRM, slow) on rep==0 per length.
## Consistency competitors (PersRel, PsychSyn, EvenOdd) + shipped eta/rr reused from sim_bench.R.
suppressMessages({library(careless); library(pROC); library(PerFit)})
setwd("C:/Users/vitto/Desktop/ReReReRe/webapp")
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

set.seed(1)
rec <- list()
for(k in seq_len(nrow(idx))){
  tag<-idx$file[k]; J<-idx$items[k]; rp<-idx$rep[k]
  mat <- as.matrix(read.csv(file.path(OUT,paste0(tag,"_M.csv")), header=FALSE))
  Y   <- read.csv(file.path(OUT,paste0(tag,"_y.csv")))
  meta<- read.csv(file.path(OUT,paste0(tag,"_meta.csv")))
  lab <- Y$y; fac<-meta$factor
  matR <- mat; rev<-which(meta$reverse==1); if(length(rev)) matR[,rev] <- 6 - mat[,rev]   # K=5 -> 6-x
  matP <- matR - 1                                                                        # PerFit wants 0..Ncat-1
  facN <- as.numeric(table(factor(fac, levels=sort(unique(fac)))))

  ## consistency competitors
  eo <- tryCatch(as.numeric(careless::evenodd(matR, factors=facN)), error=function(e)rep(NA,nrow(mat)))
  ps <- tryCatch(as.numeric(careless::psychsyn(mat, critval=.40)), error=function(e)rep(NA,nrow(mat)))
  pr <- persrel(matR, fac)

  ## person-fit: U3 (all datasets, higher = careless)
  u3v <- tryCatch(PerFit::U3poly(matP, Ncat=5)$PFscores[[1]], error=function(e){cat("U3 ERR",tag,conditionMessage(e),"\n");rep(NA,nrow(mat))})
  ## person-fit: l_z (rep==0 only; lower = careless -> negate)
  if(rp==0){
    t0<-Sys.time()
    lzraw <- tryCatch(PerFit::lzpoly(matP, Ncat=5)$PFscores[[1]], error=function(e){cat("lz ERR",tag,conditionMessage(e),"\n");rep(NA,nrow(mat))})
    lzv <- -lzraw
    cat(sprintf("lzpoly %s (J=%d) %.0fs\n", tag, J, as.numeric(difftime(Sys.time(),t0,units="secs"))))
  } else lzv <- rep(NA, nrow(mat))

  a <- c(ReReRe=aucv(lab,Y$eta), rr=aucv(lab,Y$rr), lz=aucv(lab,lzv), U3=aucv(lab,u3v),
         PersRel=aucv(lab,-pr), PsychSyn=aucv(lab,-ps), EvenOdd=aucv(lab,eo))
  rec[[k]] <- data.frame(items=J, rep=rp, t(a))
  if(rp==0) cat(sprintf("J=%-3d | ReReRe %.3f rr %.3f lz %.3f U3 %.3f PersRel %.3f PsychSyn %.3f EvenOdd %.3f\n",
    J,a["ReReRe"],a["rr"],a["lz"],a["U3"],a["PersRel"],a["PsychSyn"],a["EvenOdd"]))
}
D <- do.call(rbind, rec)
write.csv(D, file.path(OUT,"sim_bench_perfit_raw.csv"), row.names=FALSE)

METH <- c("ReReRe","rr","lz","U3","PersRel","PsychSyn","EvenOdd")
ag <- aggregate(D[,METH], by=list(items=D$items), FUN=function(x) mean(x,na.rm=TRUE))
cat("\n=== Study 3: mean AUC by questionnaire length (lz = rep 0 only; others = 10 reps) ===\n")
print(format(ag, digits=3), row.names=FALSE)
write.csv(ag, "sim_bench_perfit_meanauc.csv", row.names=FALSE)
cat("\nwrote sim_bench/sim_bench_perfit_raw.csv and sim_bench_perfit_meanauc.csv\n")
