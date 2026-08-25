## Point 3 — weight robustness on the EXTERNAL datasets ("applied unchanged everywhere").
## Recompute eta with shipped / rounded / integer / unweighted / +-50% / per-dataset-refit
## weights and compare AUC vs GT, per dataset. Uses ext_bench_oriented/ (oriented features + gate g).
suppressMessages(library(pROC)); set.seed(7)
OUT<-"ext_bench_oriented"; idx<-read.csv(file.path(OUT,"index.csv"))
W<-c(b0=-0.2084, rr=2.6756, longstring=1.2425, person_total=1.4391)
aucd<-function(y,s) as.numeric(pROC::auc(pROC::roc(y,s,levels=c(0,1),direction="<",quiet=TRUE)))
eta<-function(w,d,g) w["b0"]+w["rr"]*d$rr + g*(w["longstring"]*d$longstring + w["person_total"]*d$person_total)
oof<-function(d,K=10,reps=10){n<-nrow(d);acc<-numeric(n)
  for(rp in 1:reps){i0<-which(d$y==0);i1<-which(d$y==1)
    fold<-integer(n);fold[i0]<-sample(rep(1:K,length.out=length(i0)));fold[i1]<-sample(rep(1:K,length.out=length(i1)))
    pr<-numeric(n);for(k in 1:K){tr<-fold!=k;te<-fold==k
      m<-suppressWarnings(glm(y~rr+longstring+person_total,d[tr,],family=binomial()));pr[te]<-predict(m,d[te,],type="response")}
    acc<-acc+pr};acc/reps}

variants<-list(shipped=W, rounded=c(b0=-0.21,rr=2.7,longstring=1.2,person_total=1.4),
  integers=c(b0=0,rr=3,longstring=1,person_total=1), unweighted=c(b0=0,rr=1,longstring=1,person_total=1))
rows<-list()
for(r in 1:nrow(idx)){ nm<-as.character(idx$name[r]); g<-idx$gate[r]
  d<-read.csv(file.path(OUT,paste0(nm,"_y.csv")))
  if(any(is.na(d$rr))) d$rr[is.na(d$rr)]<-median(d$rr,na.rm=TRUE)
  d<-d[complete.cases(d[,c("y","rr","longstring","person_total")]),]
  if(length(unique(d$y))<2) next
  # validation: shipped-weight eta must match the dumped eta
  chk<-max(abs(eta(W,d,g)-d$eta))
  A<-sapply(variants,function(w)aucd(d$y,eta(w,d,g)))
  # +-50% perturbation
  P<-numeric(300); for(j in 1:300){w<-W;w[2:4]<-W[2:4]*(1+runif(3,-.5,.5));P[j]<-aucd(d$y,eta(w,d,g))}
  # per-dataset refit (CV)
  Aref<-aucd(d$y,oof(d))
  rows[[length(rows)+1]]<-data.frame(dataset=nm,n=nrow(d),J=idx$J[r],gate=round(g,3),
    prev=round(mean(d$y),3), chk=signif(chk,2),
    AUC_shipped=round(A["shipped"],3), AUC_rounded=round(A["rounded"],3),
    AUC_int311=round(A["integers"],3), AUC_unweighted=round(A["unweighted"],3),
    AUC_pert_mean=round(mean(P),3), AUC_pert_min=round(min(P),3),
    AUC_refit=round(Aref,3),
    max_abs_dev=round(max(abs(A-A["shipped"])),3))
  cat(sprintf("%-14s chk=%.0e | shipped=%.3f round=%.3f int=%.3f unw=%.3f | pert[min %.3f] | refit=%.3f | maxdev=%.3f\n",
    nm,chk,A["shipped"],A["rounded"],A["integers"],A["unweighted"],min(P),Aref,max(abs(A-A["shipped"]))))
}
res<-do.call(rbind,rows); rownames(res)<-NULL
write.csv(res,file.path(OUT,"w_sensitivity_external.csv"),row.names=FALSE)
cat("\n================ SUMMARY across",nrow(res),"external datasets ================\n")
cat(sprintf("max validation error (shipped eta vs dumped): %.1e\n", max(res$chk)))
cat(sprintf("mean |AUC(rounded) - AUC(shipped)|   = %.4f  (max %.3f)\n", mean(abs(res$AUC_rounded-res$AUC_shipped)), max(abs(res$AUC_rounded-res$AUC_shipped))))
cat(sprintf("mean |AUC(int 3,1,1) - AUC(shipped)| = %.4f  (max %.3f)\n", mean(abs(res$AUC_int311-res$AUC_shipped)), max(abs(res$AUC_int311-res$AUC_shipped))))
cat(sprintf("mean |AUC(unweighted)- AUC(shipped)| = %.4f  (max %.3f)\n", mean(abs(res$AUC_unweighted-res$AUC_shipped)), max(abs(res$AUC_unweighted-res$AUC_shipped))))
cat(sprintf("mean AUC_pert_min across sets        = %.4f  (min %.3f)\n", mean(res$AUC_pert_min), min(res$AUC_pert_min)))
cat(sprintf("mean AUC(shipped)=%.3f  vs  mean AUC(per-dataset refit)=%.3f  (refit - shipped = %+.4f)\n",
  mean(res$AUC_shipped), mean(res$AUC_refit), mean(res$AUC_refit-res$AUC_shipped)))
cat("wrote",file.path(OUT,"w_sensitivity_external.csv"),"\n")
