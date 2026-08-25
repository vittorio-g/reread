## rr_incremental_current.R  (NEW FILE — does not touch existing pipeline)
## External analog of the simulation ablation (Fig 5): incremental value of rr
## OVER the shipped auxiliaries {LongString + Person-Total}, on the CURRENT
## authoritative ext_bench data (same 4000-row subsample already frozen in *_M.csv/*_y.csv,
## so results are consistent with benchmark_ext.csv).
##   M0: y ~ longstring + person_total     M1: + rr
## OOF CV -> dAUC (DeLong) + oracle-MCC bootstrap. rr from _y.csv (shipped scorer);
## longstring/person_total via careless:: and profile-correlation.
suppressMessages({library(careless); library(pROC)})
set.seed(1)
OUT <- "ext_bench"
idx <- read.csv(file.path(OUT,"index.csv"), stringsAsFactors=FALSE)
impute <- function(m){for(j in seq_len(ncol(m))){v<-m[,j];mm<-median(v,na.rm=TRUE);if(is.na(mm))mm<-0;m[is.na(v),j]<-mm};m}
aucF <- function(y,s){ok<-is.finite(s)&is.finite(y); if(length(unique(y[ok]))<2)return(NA)
  as.numeric(pROC::auc(pROC::roc(y[ok],s[ok],levels=c(0,1),direction="<",quiet=TRUE)))}
mcc <- function(y,p){tp<-as.double(sum(p&y==1));tn<-as.double(sum(!p&y==0));fp<-as.double(sum(p&y==0));fn<-as.double(sum(!p&y==1))
  d<-sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)); if(d==0)0 else (tp*tn-fp*fn)/d}
omcc <- function(y,s){th<-quantile(s,seq(.02,.98,.02),names=FALSE,na.rm=TRUE);max(sapply(th,function(t)mcc(y,s>=t)))}
oof <- function(dat,form,K=10,reps=15){n<-nrow(dat);acc<-numeric(n)
  for(rp in 1:reps){i0<-which(dat$y==0);i1<-which(dat$y==1);fold<-integer(n)
    fold[i0]<-sample(rep(1:K,length.out=length(i0)));fold[i1]<-sample(rep(1:K,length.out=length(i1)))
    pr<-numeric(n);for(k in 1:K){tr<-fold!=k;te<-fold==k
      m<-suppressWarnings(glm(form,data=dat[tr,],family=binomial()));pr[te]<-predict(m,newdata=dat[te,],type="response")}
    acc<-acc+pr};acc/reps}

rows<-list()
for(k in seq_len(nrow(idx))){
  nm<-idx$name[k]
  M<-as.matrix(read.csv(file.path(OUT,paste0(nm,"_M.csv")),check.names=FALSE))
  Y<-read.csv(file.path(OUT,paste0(nm,"_y.csv")))
  lab<-Y$y; mat<-impute(M)
  # shipped-scorer rr (higher=careless, as used in benchmark_ext with direction "<")
  rr<-Y$rr
  # auxiliaries computed from the matrix (standard careless:: + profile correlation)
  long<-tryCatch(as.numeric(careless::longstring(mat)),error=function(e)rep(NA,nrow(mat)))
  mp<-colMeans(mat); pt<-apply(mat,1,function(r){s<-suppressWarnings(cor(r,mp)); if(is.na(s))0 else s})
  ptc<- -pt  # orient: lower profile-corr = more careless -> higher ptc = more careless
  # validation vs benchmark_ext.csv
  a_eta<-aucF(lab,Y$eta); a_rr<-aucF(lab,rr); a_long<-aucF(lab,long); a_pt<-aucF(lab,ptc)

  df<-data.frame(y=lab, rr=rr, longstring=long, person_total=ptc)
  df$rr[!is.finite(df$rr)]<-median(df$rr[is.finite(df$rr)]); df$longstring[!is.finite(df$longstring)]<-median(df$longstring,na.rm=TRUE)
  df<-df[is.finite(df$y),]
  p0<-oof(df,y~longstring+person_total); p1<-oof(df,y~longstring+person_total+rr)
  r0<-pROC::roc(df$y,p0,levels=c(0,1),direction="<",quiet=TRUE); r1<-pROC::roc(df$y,p1,levels=c(0,1),direction="<",quiet=TRUE)
  A0<-as.numeric(pROC::auc(r0)); A1<-as.numeric(pROC::auc(r1))
  dp<-tryCatch(pROC::roc.test(r1,r0,method="delong",paired=TRUE)$p.value,error=function(e)NA)
  m0<-omcc(df$y,p0); m1<-omcc(df$y,p1)
  B<-2000; i0<-which(df$y==0); i1<-which(df$y==1); db<-numeric(B)
  for(b in 1:B){bs<-c(sample(i0,replace=TRUE),sample(i1,replace=TRUE));db[b]<-omcc(df$y[bs],p1[bs])-omcc(df$y[bs],p0[bs])}
  ci<-quantile(db,c(.025,.975),names=FALSE)
  rows[[nm]]<-data.frame(dataset=nm,J=ncol(mat),n=nrow(df),npos=sum(df$y),rate=round(mean(df$y),3),
    AUC_eta=round(a_eta,3),AUC_rr=round(a_rr,3),AUC_Long=round(a_long,3),AUC_PT=round(a_pt,3),
    AUC_M0=round(A0,3),AUC_M1=round(A1,3),dAUC=round(A1-A0,3),delong_p=signif(dp,3),
    oMCC_M0=round(m0,3),oMCC_M1=round(m1,3),dMCC=round(m1-m0,3),
    dMCC_lo=round(ci[1],3),dMCC_hi=round(ci[2],3),P_dMCC_gt0=round(mean(db>0),3),row.names=NULL)
  cat(sprintf("%-14s J=%-3d n=%4d npos=%4d | eta=%.3f rr=%.3f Long=%.3f PT=%.3f | dAUC=%+.3f p=%.2g | dMCC=%+.3f[%+.3f,%+.3f] P=%.2f\n",
    nm,ncol(mat),nrow(df),sum(df$y),a_eta,a_rr,a_long,a_pt,A1-A0,dp,m1-m0,ci[1],ci[2],mean(db>0)))
}
res<-do.call(rbind,rows)
write.csv(res,file.path(OUT,"rr_incremental_current.csv"),row.names=FALSE)
cat("\n=== stratum summary (by J) ===\n")
res$stratum<-ifelse(res$J>=140,"long(>=140)",ifelse(res$J>=60,"med(60-139)","short(<60)"))
for(st in c("long(>=140)","med(60-139)","short(<60)")){s<-res[res$stratum==st,]
  if(nrow(s)>0) cat(sprintf("%-12s mean dMCC=%+.3f | sig dMCC(P>0.975): %d/%d | %s\n",st,mean(s$dMCC),sum(s$P_dMCC_gt0>=0.975),nrow(s),paste(s$dataset,collapse=", ")))}
cat("\nwrote ext_bench/rr_incremental_current.csv\n")
