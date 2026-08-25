# Why is rr's incremental AUC tiny? Test the redundancy-with-PT hypothesis,
# and check whether rr helps oracle-MCC more than AUC (metric fairness vs sim).
suppressWarnings(suppressMessages(library(pROC)))
set.seed(1)
BASE <- "C:/Users/vitto/Desktop/ReReReRe/Dataset/gt_benchmark_candidates"
reg <- list(
  list(name="Kay S2",J=363,s="kay_idris_idria/kay_scores_s2.csv",l="kay_idris_idria/kay_labels_s2.csv"),
  list(name="warning IPIP",J=300,s="warning_ipipneo300/_scores.csv",l="warning_ipipneo300/_labels.csv"),
  list(name="Kay S1",J=250,s="kay_idris_idria/kay_scores_s1.csv",l="kay_idris_idria/kay_labels_s1.csv"),
  list(name="opsy 16PF",J=163,s="opsy_16pf_scores.csv",l="opsy_16pf_labels.csv"),
  list(name="smarvus",J=141,s="smarvus/smarvus_scores.csv",l="smarvus/smarvus_labels.csv"),
  list(name="Kay S6",J=67,s="kay_idris_idria/kay_scores_s6.csv",l="kay_idris_idria/kay_labels_s6.csv"),
  list(name="Kay S5",J=62,s="kay_idris_idria/kay_scores_s5.csv",l="kay_idris_idria/kay_labels_s5.csv")
)
mcc <- function(y, pred, thr){ p<-pred>=thr
  tp<-as.double(sum(p&y==1));tn<-as.double(sum(!p&y==0));fp<-as.double(sum(p&y==0));fn<-as.double(sum(!p&y==1))
  d<-sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)); if(d==0) return(0); (tp*tn-fp*fn)/d }
oracle_mcc <- function(y, s){ th<-quantile(s, probs=seq(0.01,0.99,by=0.01)); max(sapply(th, function(t) mcc(y,s,t))) }
oof <- function(dat, form, K=10, reps=15){ n<-nrow(dat); acc<-numeric(n)
  for(rp in 1:reps){ i0<-which(dat$y==0);i1<-which(dat$y==1)
    fold<-integer(n); fold[i0]<-sample(rep(1:K,length.out=length(i0))); fold[i1]<-sample(rep(1:K,length.out=length(i1)))
    pr<-numeric(n); for(k in 1:K){tr<-fold!=k;te<-fold==k
      m<-suppressWarnings(glm(form,data=dat[tr,],family=binomial())); pr[te]<-predict(m,newdata=dat[te,],type="response")}
    acc<-acc+pr }; acc/reps }
cat(sprintf("%-14s %5s %6s %6s | %7s %7s %7s\n","dataset","r_rrPT","rho","VIF~","oMCC0","oMCC1","dMCC"))
for(d in reg){
  sc<-read.csv(file.path(BASE,d$s)); lb<-read.csv(file.path(BASE,d$l))
  n<-min(nrow(sc),nrow(lb)); sc<-sc[1:n,]; y<-as.integer(lb[[1]][1:n])
  df<-data.frame(y=y,rr=sc$rr_z,longstring=sc$longstring,person_total=sc$person_total)
  df$rr[is.na(df$rr)]<-median(df$rr,na.rm=TRUE); df<-df[complete.cases(df),]
  rP<-cor(df$rr,df$person_total); rho<-cor(df$rr,df$person_total,method="spearman")
  vif<-1/(1-summary(lm(rr~person_total+longstring,data=df))$r.squared)
  p0<-oof(df,y~longstring+person_total); p1<-oof(df,y~longstring+person_total+rr)
  m0<-oracle_mcc(df$y,p0); m1<-oracle_mcc(df$y,p1)
  cat(sprintf("%-14s %5.2f %6.2f %6.2f | %7.3f %7.3f %+7.3f\n",d$name,rP,rho,vif,m0,m1,m1-m0))
}
