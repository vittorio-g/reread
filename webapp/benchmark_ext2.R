## benchmark_ext2.R — Study-2 competitive benchmark on the 17 external datasets (same as tab:ext /
## tab:ltpa) with archived labels. Structure-free competitors (no factor map externally):
## Mahalanobis, IRV, LongString, psychometric synonyms. Per-dataset AUC + DeLong; cross-dataset
## paired Wilcoxon of ReReRe/rr vs each competitor. Directions fixed a priori (high=careless).
suppressMessages({library(careless); library(pROC)})
set.seed(1)
OUT <- "ext_bench2"
idx <- read.csv(file.path(OUT,"index.csv"), stringsAsFactors=FALSE)
impute <- function(m){for(j in seq_len(ncol(m))){v<-m[,j];mm<-median(v,na.rm=TRUE);if(is.na(mm))mm<-0;m[is.na(v),j]<-mm};m}
mkroc <- function(lab,s){ok<-is.finite(s); if(sum(ok)<20||length(unique(lab[ok]))<2) return(NULL); roc(lab[ok],s[ok],direction="<",quiet=TRUE)}
aucv  <- function(lab,s){r<-mkroc(lab,s); if(is.null(r))NA else as.numeric(auc(r))}
delong<- function(lab,a,b){ok<-is.finite(a)&is.finite(b); if(sum(ok)<20||length(unique(lab[ok]))<2)return(NA)
  r1<-roc(lab[ok],a[ok],direction="<",quiet=TRUE);r2<-roc(lab[ok],b[ok],direction="<",quiet=TRUE)
  tryCatch(roc.test(r1,r2,method="delong",paired=TRUE)$p.value,error=function(e)NA)}
COMP <- c("Mahalanobis","IRV","LongString","PsychSyn")
rows<-list(); pmat_ens<-list(); pmat_rr<-list()
for(k in seq_len(nrow(idx))){
  nm<-idx$name[k]
  M<-as.matrix(read.csv(file.path(OUT,paste0(nm,"_M.csv")),header=FALSE))
  Y<-read.csv(file.path(OUT,paste0(nm,"_y.csv"))); lab<-Y$y; mat<-impute(M)
  ps<-tryCatch(as.numeric(careless::psychsyn(mat,critval=.60)),error=function(e)rep(NA,nrow(mat)))
  if(sum(is.finite(ps))<0.5*nrow(mat)) ps<-tryCatch(as.numeric(careless::psychsyn(mat,critval=.40)),error=function(e)rep(NA,nrow(mat)))
  S<-list(ReReRe=Y$eta, rr=Y$rr, PT=Y$pt,
    Mahalanobis=tryCatch(as.numeric(careless::mahad(mat,plot=FALSE,flag=FALSE)),error=function(e)rep(NA,nrow(mat))),
    IRV=-as.numeric(careless::irv(mat)), LongString=as.numeric(careless::longstring(mat)), PsychSyn=-ps)
  a<-sapply(names(S),function(m) aucv(lab,S[[m]]))
  rows[[nm]]<-data.frame(dataset=nm,kind=idx$kind[k],n=nrow(mat),J=ncol(mat),rate=round(mean(lab),3),
    ReReRe=round(a["ReReRe"],3),rr=round(a["rr"],3),PT=round(a["PT"],3),Mah=round(a["Mahalanobis"],3),
    IRV=round(a["IRV"],3),Long=round(a["LongString"],3),Psyn=round(a["PsychSyn"],3),row.names=NULL)
  pmat_ens[[nm]]<-sapply(COMP,function(c) delong(lab,S[["ReReRe"]],S[[c]]))
  pmat_rr[[nm]] <-sapply(COMP,function(c) delong(lab,S[["rr"]],S[[c]]))
  cat(sprintf("%-14s J=%-3d rate=%.2f | ens %.3f rr %.3f PT %.3f | Mah %.3f IRV %.3f Long %.3f Psyn %s\n",
    nm,ncol(mat),mean(lab),a["ReReRe"],a["rr"],a["PT"],a["Mahalanobis"],a["IRV"],a["LongString"],
    ifelse(is.na(a["PsychSyn"]),"NA",sprintf("%.3f",a["PsychSyn"]))))
}
res<-do.call(rbind,rows); write.csv(res,file.path(OUT,"benchmark_ext2.csv"),row.names=FALSE)
cat("\n=== per-dataset AUC (sorted by J desc) ===\n"); print(res[order(-res$J),c("dataset","kind","J","rate","ReReRe","rr","PT","Mah","IRV","Long","Psyn")],row.names=FALSE)
cat("\nmean AUC:  ReReRe",round(mean(res$ReReRe,na.rm=TRUE),3)," rr",round(mean(res$rr,na.rm=TRUE),3),
  " Mah",round(mean(res$Mah,na.rm=TRUE),3)," IRV",round(mean(res$IRV,na.rm=TRUE),3),
  " Long",round(mean(res$Long,na.rm=TRUE),3)," Psyn",round(mean(res$Psyn,na.rm=TRUE),3),"\n")
## cross-dataset paired Wilcoxon (one-sided: ReReRe/rr > competitor), and win counts
cat("\n=== cross-dataset paired tests (n datasets) ===\n")
for(c in c("Mah","IRV","Long","Psyn")){
  de<-res$ReReRe-res[[c]]; dr<-res$rr-res[[c]]
  we<-tryCatch(wilcox.test(res$ReReRe,res[[c]],paired=TRUE,alternative="greater")$p.value,error=function(e)NA)
  wr<-tryCatch(wilcox.test(res$rr,res[[c]],paired=TRUE,alternative="greater")$p.value,error=function(e)NA)
  cat(sprintf("  vs %-11s ReReRe wins %2d/%2d (Wilcoxon p=%.3f) | rr wins %2d/%2d (p=%.3f)\n",
    c,sum(de>0,na.rm=TRUE),sum(is.finite(de)),we,sum(dr>0,na.rm=TRUE),sum(is.finite(dr)),wr))
}
## count datasets with a significant DeLong advantage (p<.05 & higher AUC)
cat("\n=== datasets with significant DeLong advantage (p<.05) ===\n")
for(ci in seq_along(COMP)){c<-COMP[ci]; cc<-c("Mah","IRV","Long","Psyn")[ci]
  se<-sum(sapply(names(pmat_ens),function(nm) {p<-pmat_ens[[nm]][ci]; a<-res[res$dataset==nm,]; isTRUE(p<.05 && a$ReReRe>a[[cc]])}))
  sr<-sum(sapply(names(pmat_rr), function(nm) {p<-pmat_rr[[nm]][ci];  a<-res[res$dataset==nm,]; isTRUE(p<.05 && a$rr>a[[cc]])}))
  cat(sprintf("  vs %-11s ReReRe sig-better in %2d, rr sig-better in %2d datasets\n",c,se,sr))
}
cat("\nwrote benchmark_ext2.csv\n")
