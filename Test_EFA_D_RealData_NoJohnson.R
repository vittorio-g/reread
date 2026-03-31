###############################################################################
# EFA-D vs Legacy: Real Data Validation (excluding Johnson)
###############################################################################
setwd("C:/Users/vitto/Desktop/ReReReRe")
library(dplyr); library(pROC); library(haven)
source("ReReReRe.R")

compute_mcc <- function(tp, tn, fp, fn) {
  tp<-as.double(tp);tn<-as.double(tn);fp<-as.double(fp);fn<-as.double(fn)
  num<-(tp*tn)-(fp*fn); den<-sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  if(is.na(den)||den==0) return(0); num/den
}
evaluate <- function(pred, true) {
  tp<-sum(pred==1&true==1,na.rm=T);tn<-sum(pred==0&true==0,na.rm=T)
  fp<-sum(pred==1&true==0,na.rm=T);fn<-sum(pred==0&true==1,na.rm=T)
  list(MCC=compute_mcc(tp,tn,fp,fn),Sens=tp/(tp+fn),Spec=tn/(tn+fp))
}
compute_auc <- function(s,l,d) tryCatch(as.numeric(auc(roc(l,s,direction=d,quiet=T))),error=function(e)NA)
run_mah <- function(d,alpha=.001) {
  m<-as.matrix(d);p<-ncol(m);ctr<-colMeans(m,na.rm=T);cv<-cov(m,use="pairwise.complete.obs")
  if(any(is.na(cv))||det(cv)<1e-10) cv<-cv+diag(.01,p)
  d2<-mahalanobis(m,ctr,cv); list(d2=d2,flagged=d2>qchisq(1-alpha,df=p))
}
find_best_z <- function(z,labels) {
  bm<- -Inf;bz<-1.5
  for(zt in seq(0.1,3,.1)){ev<-evaluate(as.integer(z<=zt),labels);if(ev$MCC>bm){bm<-ev$MCC;bz<-zt}}
  list(z=bz,mcc=bm)
}
run3 <- function(items,labels,name) {
  cat(sprintf("  Running %s... ",name))
  rr<-tryCatch(ReReReRe(items,corProp=.03,iterations=100,align_signs=T,mode=name),
               error=function(e){cat(sprintf("ERROR: %s\n",e$message));NULL})
  if(is.null(rr)) return(list(auc=NA,mcc15=NA,oracle=NA,oz=NA))
  a<-compute_auc(rr$z_score,labels,">"); ev<-evaluate(as.integer(rr$z_score<=1.5),labels)
  b<-find_best_z(rr$z_score,labels)
  cat(sprintf("AUC=%.3f MCC15=%.3f Oracle=%.3f(z=%.1f)\n",
      ifelse(is.na(a),0,a),ev$MCC,b$mcc,b$z))
  list(auc=a,mcc15=ev$MCC,oracle=b$mcc,oz=b$z)
}

cat("\n============================================================\n")
cat("EFA-D vs COUPLED vs WEIGHTED: Real Data (no Johnson)\n")
cat("============================================================\n\n")

results <- list()

# 1. Schroeders
cat("=== 1. Schroeders HEXACO-60 ===\n")
d<-read.csv("external_datasets/01_Schroeders_2022/data_mod_resp.csv",sep=";")
lab<-d$Careless; items<-d%>%select(starts_with("HE"))
cc<-complete.cases(items); items<-items[cc,]; lab<-lab[cc]
cat(sprintf("  N=%d items=%d careless=%.1f%%\n",nrow(items),ncol(items),mean(lab)*100))
e<-run3(items,lab,"efa_d"); c_<-run3(items,lab,"coupled"); w<-run3(items,lab,"weighted")
mh<-run_mah(items); emh<-evaluate(mh$flagged,lab); amh<-compute_auc(mh$d2,lab,"<")
results[[1]]<-data.frame(DS="Schroeders",It=60,nF="~10",N=nrow(items),
  efaD_auc=e$auc,coup_auc=c_$auc,wt_auc=w$auc,mah_auc=amh,
  efaD_mcc15=e$mcc15,coup_mcc15=c_$mcc15,wt_mcc15=w$mcc15,mah_mcc=emh$MCC,
  efaD_or=e$oracle,coup_or=c_$oracle,wt_or=w$oracle)

# 2. Schneider
cat("\n=== 2. Schneider QoL-31 ===\n")
d<-read.csv("external_datasets/06_Schneider_QoL/carersp.csv")
lab<-d$c01; ic<-grep("^(dep|pain|cog|fat)\\d",names(d),value=T); items<-d[,ic]
cc<-complete.cases(items)&!is.na(lab); items<-items[cc,]; lab<-lab[cc]
cat(sprintf("  N=%d items=%d careless=%.1f%%\n",nrow(items),ncol(items),mean(lab)*100))
e<-run3(items,lab,"efa_d"); c_<-run3(items,lab,"coupled"); w<-run3(items,lab,"weighted")
mh<-run_mah(items); emh<-evaluate(mh$flagged,lab); amh<-compute_auc(mh$d2,lab,"<")
results[[2]]<-data.frame(DS="Schneider",It=31,nF="~5",N=nrow(items),
  efaD_auc=e$auc,coup_auc=c_$auc,wt_auc=w$auc,mah_auc=amh,
  efaD_mcc15=e$mcc15,coup_mcc15=c_$mcc15,wt_mcc15=w$mcc15,mah_mcc=emh$MCC,
  efaD_or=e$oracle,coup_or=c_$oracle,wt_or=w$oracle)

# 3. Niessen
cat("\n=== 3. Niessen IPIP-100 ===\n")
d<-read_sav("external_datasets/08_Niessen_2016/Raw_data.sav")%>%filter(Use_me==1)
lab<-as.integer(d$Conditon); ic<-grep("^[ECANO]\\d+$",names(d),value=T)
items<-data.frame(lapply(as.data.frame(d[,ic]),as.numeric))
cc<-complete.cases(items); items<-items[cc,]; lab<-lab[cc]
cat(sprintf("  N=%d items=%d careless=%.1f%%\n",nrow(items),ncol(items),mean(lab)*100))
e<-run3(items,lab,"efa_d"); c_<-run3(items,lab,"coupled"); w<-run3(items,lab,"weighted")
mh<-run_mah(items); emh<-evaluate(mh$flagged,lab); amh<-compute_auc(mh$d2,lab,"<")
results[[3]]<-data.frame(DS="Niessen",It=100,nF="~5",N=nrow(items),
  efaD_auc=e$auc,coup_auc=c_$auc,wt_auc=w$auc,mah_auc=amh,
  efaD_mcc15=e$mcc15,coup_mcc15=c_$mcc15,wt_mcc15=w$mcc15,mah_mcc=emh$MCC,
  efaD_or=e$oracle,coup_or=c_$oracle,wt_or=w$oracle)

# 4-6. Goldammer
for(si in 1:3) {
  cat(sprintf("\n=== %d. Goldammer Study %d ===\n", si+3, si))
  g<-read.csv(sprintf("external_datasets/13_Goldammer_2024/Study_%d.csv",si))
  if(si<=2){lab<-as.integer(g$careless_all)}else{lab<-as.integer(g$condition>0)}
  ic<-grep("^[ecano]_[a-z]+[0-9]+$",names(g),value=T)
  if(length(ic)==0) ic<-grep("^[ecano]_[a-z]+[0-9]+_t1$",names(g),value=T)
  items<-g[,ic]; cc<-complete.cases(items)&!is.na(lab); items<-items[cc,]; lab<-lab[cc]
  cat(sprintf("  N=%d items=%d careless=%.1f%%\n",nrow(items),ncol(items),mean(lab)*100))
  e<-run3(items,lab,"efa_d"); c_<-run3(items,lab,"coupled"); w<-run3(items,lab,"weighted")
  mh<-run_mah(items); emh<-evaluate(mh$flagged,lab); amh<-compute_auc(mh$d2,lab,"<")
  nm<-paste0("Goldammer_S",si)
  results[[si+3]]<-data.frame(DS=nm,It=ncol(items),nF=c("~8","~6","~7")[si],N=nrow(items),
    efaD_auc=e$auc,coup_auc=c_$auc,wt_auc=w$auc,mah_auc=amh,
    efaD_mcc15=e$mcc15,coup_mcc15=c_$mcc15,wt_mcc15=w$mcc15,mah_mcc=emh$MCC,
    efaD_or=e$oracle,coup_or=c_$oracle,wt_or=w$oracle)
}

# Summary
final<-do.call(rbind,results)
cat("\n\n============================================================\n")
cat("SUMMARY TABLE\n")
cat("============================================================\n\n")

cat("--- AUC ---\n")
cat(sprintf("%-15s %4s %4s  %6s %6s %6s  %6s  Winner\n","Dataset","It","nF","EFA-D","Coup","Wt","Mah"))
for(i in 1:nrow(final)){r<-final[i,]
  aucs<-c(r$efaD_auc,r$coup_auc,r$wt_auc)
  w<-c("EFA-D","Coupled","Weighted")[which.max(aucs)]
  cat(sprintf("%-15s %4d %4s  %.3f  %.3f  %.3f  %.3f  %s\n",
      r$DS,r$It,r$nF,r$efaD_auc,r$coup_auc,r$wt_auc,r$mah_auc,w))
}

cat("\n--- MCC z=1.5 ---\n")
cat(sprintf("%-15s %4s %4s  %6s %6s %6s  %6s  Winner\n","Dataset","It","nF","EFA-D","Coup","Wt","Mah"))
for(i in 1:nrow(final)){r<-final[i,]
  mccs<-c(r$efaD_mcc15,r$coup_mcc15,r$wt_mcc15)
  w<-c("EFA-D","Coupled","Weighted")[which.max(mccs)]
  cat(sprintf("%-15s %4d %4s  %.3f  %.3f  %.3f  %.3f  %s\n",
      r$DS,r$It,r$nF,r$efaD_mcc15,r$coup_mcc15,r$wt_mcc15,r$mah_mcc,w))
}

cat("\n--- Oracle MCC ---\n")
cat(sprintf("%-15s %4s %4s  %6s %6s %6s  Winner\n","Dataset","It","nF","EFA-D","Coup","Wt"))
for(i in 1:nrow(final)){r<-final[i,]
  mccs<-c(r$efaD_or,r$coup_or,r$wt_or)
  w<-c("EFA-D","Coupled","Weighted")[which.max(mccs)]
  cat(sprintf("%-15s %4d %4s  %.3f  %.3f  %.3f  %s\n",
      r$DS,r$It,r$nF,r$efaD_or,r$coup_or,r$wt_or,w))
}

cat("\n--- MEANS ---\n")
cat(sprintf("AUC:    EFA-D=%.3f Coupled=%.3f Weighted=%.3f Mah=%.3f\n",
    mean(final$efaD_auc,na.rm=T),mean(final$coup_auc,na.rm=T),mean(final$wt_auc,na.rm=T),mean(final$mah_auc,na.rm=T)))
cat(sprintf("MCC15:  EFA-D=%.3f Coupled=%.3f Weighted=%.3f Mah=%.3f\n",
    mean(final$efaD_mcc15,na.rm=T),mean(final$coup_mcc15,na.rm=T),mean(final$wt_mcc15,na.rm=T),mean(final$mah_mcc,na.rm=T)))
cat(sprintf("Oracle: EFA-D=%.3f Coupled=%.3f Weighted=%.3f\n",
    mean(final$efaD_or,na.rm=T),mean(final$coup_or,na.rm=T),mean(final$wt_or,na.rm=T)))

write.csv(final,"efa_d_realdata_no_johnson.csv",row.names=F)
cat("\n=== DONE ===\n")
