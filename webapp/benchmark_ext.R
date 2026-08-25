## benchmark_ext.R — head-to-head on external GT datasets (no scale structure available,
## so structure-free competitors: Mahalanobis, IRV, LongString, psychometric synonyms).
## rr / ReReRe (shipped, from Node) vs each, AUC + DeLong, per dataset. Directions fixed a priori.
suppressMessages({library(careless); library(pROC)})
set.seed(1)
OUT <- "ext_bench"
idx <- read.csv(file.path(OUT,"index.csv"), stringsAsFactors=FALSE)
impute <- function(m){for(j in seq_len(ncol(m))){v<-m[,j];mm<-median(v,na.rm=TRUE);if(is.na(mm))mm<-0;m[is.na(v),j]<-mm};m}
mkroc <- function(lab,s){ok<-is.finite(s); if(sum(ok)<20||length(unique(lab[ok]))<2) return(NULL)
  roc(lab[ok], s[ok], direction="<", quiet=TRUE)}
delong <- function(lab,a,b){ok<-is.finite(a)&is.finite(b); if(sum(ok)<20||length(unique(lab[ok]))<2) return(NA)
  r1<-roc(lab[ok],a[ok],direction="<",quiet=TRUE);r2<-roc(lab[ok],b[ok],direction="<",quiet=TRUE)
  tryCatch(roc.test(r1,r2,method="delong",paired=TRUE)$p.value,error=function(e)NA)}
aucv <- function(lab,s){r<-mkroc(lab,s); if(is.null(r))NA else as.numeric(auc(r))}

rows <- list()
for(k in seq_len(nrow(idx))){
  nm <- idx$name[k]
  M <- as.matrix(read.csv(file.path(OUT,paste0(nm,"_M.csv")), check.names=FALSE))
  Y <- read.csv(file.path(OUT,paste0(nm,"_y.csv")))
  lab <- Y$y; mat <- impute(M)
  mahad_v <- tryCatch(as.numeric(careless::mahad(mat,plot=FALSE,flag=FALSE)),error=function(e)rep(NA,nrow(mat)))
  irv_v   <- tryCatch(as.numeric(careless::irv(mat)),error=function(e)rep(NA,nrow(mat)))
  long_v  <- tryCatch(as.numeric(careless::longstring(mat)),error=function(e)rep(NA,nrow(mat)))
  psyn_v  <- tryCatch(as.numeric(careless::psychsyn(mat,critval=.60)),error=function(e)rep(NA,nrow(mat)))
  S <- list("ReReRe"=Y$eta,"rr"=Y$rr,"Mahalanobis"=mahad_v,"IRV"=-irv_v,"LongString"=long_v,"PsychSyn"=-psyn_v)
  au <- sapply(S, function(s) aucv(lab,s))
  # DeLong: rr and ensemble vs the 4 competitors
  comp <- c("Mahalanobis","IRV","LongString","PsychSyn")
  p_rr  <- sapply(comp, function(c) delong(lab,S[["rr"]],S[[c]]))
  p_ens <- sapply(comp, function(c) delong(lab,S[["ReReRe"]],S[[c]]))
  rows[[nm]] <- data.frame(dataset=nm,n=nrow(mat),J=ncol(mat),rate=round(mean(lab),3),
    AUC_ReReRe=round(au["ReReRe"],3),AUC_rr=round(au["rr"],3),AUC_Mah=round(au["Mahalanobis"],3),
    AUC_IRV=round(au["IRV"],3),AUC_Long=round(au["LongString"],3),AUC_Psyn=round(au["PsychSyn"],3),
    row.names=NULL)
  cat(sprintf("%-14s J=%-3d rate=%.2f | ReReRe %.3f  rr %.3f | Mah %.3f IRV %.3f Long %.3f Psyn %.3f\n",
    nm,ncol(mat),mean(lab),au["ReReRe"],au["rr"],au["Mahalanobis"],au["IRV"],au["LongString"],au["PsychSyn"]))
  # stash p-values
  rows[[nm]]$p_rr_Mah<-signif(p_rr["Mahalanobis"],3); rows[[nm]]$p_rr_IRV<-signif(p_rr["IRV"],3)
  rows[[nm]]$p_rr_Long<-signif(p_rr["LongString"],3); rows[[nm]]$p_rr_Psyn<-signif(p_rr["PsychSyn"],3)
  rows[[nm]]$p_ens_Mah<-signif(p_ens["Mahalanobis"],3); rows[[nm]]$p_ens_IRV<-signif(p_ens["IRV"],3)
  rows[[nm]]$p_ens_Long<-signif(p_ens["LongString"],3); rows[[nm]]$p_ens_Psyn<-signif(p_ens["PsychSyn"],3)
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUT,"benchmark_ext.csv"), row.names=FALSE)
cat("\n=== external AUC summary (sorted by J) ===\n")
o <- order(res$J); print(res[o, c("dataset","n","J","rate","AUC_ReReRe","AUC_rr","AUC_Mah","AUC_IRV","AUC_Long","AUC_Psyn")], row.names=FALSE)
# how often rr / ensemble beat each competitor (numerically and significantly at .05)
cat("\nrr AUC > competitor (of",nrow(res),"datasets):",
  " Mah",sum(res$AUC_rr>res$AUC_Mah,na.rm=TRUE)," IRV",sum(res$AUC_rr>res$AUC_IRV,na.rm=TRUE),
  " Long",sum(res$AUC_rr>res$AUC_Long,na.rm=TRUE)," Psyn",sum(res$AUC_rr>res$AUC_Psyn,na.rm=TRUE),"\n")
cat("ensemble AUC > competitor:",
  " Mah",sum(res$AUC_ReReRe>res$AUC_Mah,na.rm=TRUE)," IRV",sum(res$AUC_ReReRe>res$AUC_IRV,na.rm=TRUE),
  " Long",sum(res$AUC_ReReRe>res$AUC_Long,na.rm=TRUE)," Psyn",sum(res$AUC_ReReRe>res$AUC_Psyn,na.rm=TRUE),"\n")
cat("mean AUC:  ReReRe",round(mean(res$AUC_ReReRe,na.rm=TRUE),3)," rr",round(mean(res$AUC_rr,na.rm=TRUE),3),
  " Mah",round(mean(res$AUC_Mah,na.rm=TRUE),3)," IRV",round(mean(res$AUC_IRV,na.rm=TRUE),3),
  " Long",round(mean(res$AUC_Long,na.rm=TRUE),3)," Psyn",round(mean(res$AUC_Psyn,na.rm=TRUE),3),"\n")
cat("\nwrote benchmark_ext.csv\n")
