# Bootstrap CI on incremental oracle-MCC (M1 = LS+PT+rr vs M0 = LS+PT)
# on the long, incoherence-relevant external batteries.
# OOF predictions fixed (10-fold x 20 reps CV); stratified paired eval-bootstrap,
# oracle threshold re-selected within each resample (symmetric across M0/M1).
suppressWarnings(suppressMessages(library(pROC)))
set.seed(2026)
BASE <- "C:/Users/vitto/Desktop/ReReReRe/Dataset/gt_benchmark_candidates"
reg <- list(
  list(name="Kay S2", J=363, gt="self-report exclusion", s="kay_idris_idria/kay_scores_s2.csv", l="kay_idris_idria/kay_labels_s2.csv"),
  list(name="smarvus",J=141, gt=">=2 instructed checks",  s="smarvus/smarvus_scores.csv",        l="smarvus/smarvus_labels.csv"),
  list(name="opsy 16PF",J=163,gt="extreme speeders",      s="opsy_16pf_scores.csv",              l="opsy_16pf_labels.csv")
)
mcc <- function(y,p){ tp<-as.double(sum(p&y==1));tn<-as.double(sum(!p&y==0))
  fp<-as.double(sum(p&y==0));fn<-as.double(sum(!p&y==1))
  d<-sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)); if(d==0) 0 else (tp*tn-fp*fn)/d }
oracle_mcc <- function(y,s){ th<-quantile(s,probs=seq(0.02,0.98,by=0.02),names=FALSE)
  max(sapply(th,function(t) mcc(y,s>=t))) }
oof <- function(dat,form,K=10,reps=20){ n<-nrow(dat); acc<-numeric(n)
  for(rp in 1:reps){ i0<-which(dat$y==0);i1<-which(dat$y==1)
    fold<-integer(n); fold[i0]<-sample(rep(1:K,length.out=length(i0))); fold[i1]<-sample(rep(1:K,length.out=length(i1)))
    pr<-numeric(n); for(k in 1:K){tr<-fold!=k;te<-fold==k
      m<-suppressWarnings(glm(form,data=dat[tr,],family=binomial())); pr[te]<-predict(m,newdata=dat[te,],type="response")}
    acc<-acc+pr }; acc/reps }

B <- 2000
out <- list()
for(d in reg){
  sc<-read.csv(file.path(BASE,d$s)); lb<-read.csv(file.path(BASE,d$l))
  n<-min(nrow(sc),nrow(lb)); sc<-sc[1:n,]; y<-as.integer(lb[[1]][1:n])
  df<-data.frame(y=y,rr=sc$rr_z,longstring=sc$longstring,person_total=sc$person_total)
  df$rr[is.na(df$rr)]<-median(df$rr,na.rm=TRUE); df<-df[complete.cases(df),]
  p0<-oof(df,y~longstring+person_total); p1<-oof(df,y~longstring+person_total+rr)
  m0<-oracle_mcc(df$y,p0); m1<-oracle_mcc(df$y,p1); pt<-oracle_mcc(df$y,df$person_total)
  i0<-which(df$y==0); i1<-which(df$y==1); db<-numeric(B)
  for(b in 1:B){ bs<-c(sample(i0,replace=TRUE),sample(i1,replace=TRUE))
    db[b]<-oracle_mcc(df$y[bs],p1[bs])-oracle_mcc(df$y[bs],p0[bs]) }
  ci<-quantile(db,c(.025,.975),names=FALSE); pgt<-mean(db>0)
  out[[length(out)+1]]<-data.frame(dataset=d$name,J=d$J,n=nrow(df),npos=sum(df$y),gt=d$gt,
    oMCC_PTonly=round(pt,3),oMCC_M0=round(m0,3),oMCC_M1=round(m1,3),
    dMCC=round(m1-m0,3),boot_mean=round(mean(db),3),ci_lo=round(ci[1],3),ci_hi=round(ci[2],3),P_gt0=round(pgt,3))
  cat(sprintf("%-10s J=%3d npos=%4d | PTonly=%.3f M0=%.3f M1=%.3f | dMCC=%+.3f  95%%CI[%+.3f,%+.3f]  P(d>0)=%.3f\n",
    d$name,d$J,sum(df$y),pt,m0,m1,m1-m0,ci[1],ci[2],pgt))
}
res<-do.call(rbind,out)
write.csv(res,"C:/Users/vitto/Desktop/ReReReRe/rr_bootstrap_dmcc.csv",row.names=FALSE)
cat("\nSaved: rr_bootstrap_dmcc.csv\n")
