# hexaco: GT = self-report V1/V2 disagree (<=3) = careless. Long battery (240),
# incoherence-relevant GT, but PT>ens and rr weak -> honest test of whether rr
# adds incremental MCC even here.
suppressWarnings(suppressMessages(library(pROC)))
set.seed(2026)
D <- "C:/Users/vitto/Desktop/ReReReRe/Dataset/gt_benchmark_candidates/opsy_hexaco"
raw <- read.csv(file.path(D,"data.csv"), sep="\t", check.names=FALSE)
sc  <- read.csv(file.path(D,"..","opsy_hexaco_scores.csv"))
# 240 item columns = all except last 4 (V1,V2,country,elapse)
itemcols <- setdiff(names(raw), c("V1","V2","country","elapse"))
stopifnot(length(itemcols)==240)
cc <- complete.cases(raw[,itemcols])
keep <- raw[cc,]
cat("raw=",nrow(raw)," complete=",nrow(keep)," scores=",nrow(sc),"\n")
n <- min(nrow(keep), nrow(sc)); keep <- keep[1:n,]; sc <- sc[1:n,]
y <- as.integer( (keep$V1<=3) | (keep$V2<=3) )
cat("n_pos=",sum(y)," rate=",round(mean(y),4),"  (readme: 435, 1.9%)\n")

auc1 <- function(y,x) as.numeric(pROC::auc(pROC::roc(y,x,quiet=TRUE,direction="auto")))
a_ens<-auc1(y,sc$ens_p); a_rr<-auc1(y,sc$rr_careless); a_pt<-auc1(y,sc$person_total)
cat(sprintf("VALIDATION AUC  ens=%.3f (0.697)  rr=%.3f (0.588)  PT=%.3f (0.709)\n",a_ens,a_rr,a_pt))

df <- data.frame(y=y, rr=sc$rr_z, longstring=sc$longstring, person_total=sc$person_total)
df$rr[is.na(df$rr)] <- median(df$rr,na.rm=TRUE); df <- df[complete.cases(df),]
mcc <- function(y,p){tp<-as.double(sum(p&y==1));tn<-as.double(sum(!p&y==0));fp<-as.double(sum(p&y==0));fn<-as.double(sum(!p&y==1))
  d<-sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)); if(d==0)0 else (tp*tn-fp*fn)/d}
omcc<-function(y,s){th<-quantile(s,seq(.02,.98,.02),names=FALSE);max(sapply(th,function(t)mcc(y,s>=t)))}
oof<-function(dat,form,K=10,reps=15){n<-nrow(dat);acc<-numeric(n)
  for(rp in 1:reps){i0<-which(dat$y==0);i1<-which(dat$y==1);fold<-integer(n)
    fold[i0]<-sample(rep(1:K,length.out=length(i0)));fold[i1]<-sample(rep(1:K,length.out=length(i1)))
    pr<-numeric(n);for(k in 1:K){tr<-fold!=k;te<-fold==k
      m<-suppressWarnings(glm(form,data=dat[tr,],family=binomial()));pr[te]<-predict(m,newdata=dat[te,],type="response")}
    acc<-acc+pr};acc/reps}
p0<-oof(df,y~longstring+person_total);p1<-oof(df,y~longstring+person_total+rr)
r0<-pROC::roc(df$y,p0,quiet=TRUE,direction="auto");r1<-pROC::roc(df$y,p1,quiet=TRUE,direction="auto")
A0<-as.numeric(pROC::auc(r0));A1<-as.numeric(pROC::auc(r1));dp<-pROC::roc.test(r1,r0,method="delong",paired=TRUE)$p.value
m0<-omcc(df$y,p0);m1<-omcc(df$y,p1)
B<-2000;i0<-which(df$y==0);i1<-which(df$y==1);db<-numeric(B)
for(b in 1:B){bs<-c(sample(i0,replace=TRUE),sample(i1,replace=TRUE));db[b]<-omcc(df$y[bs],p1[bs])-omcc(df$y[bs],p0[bs])}
ci<-quantile(db,c(.025,.975),names=FALSE)
cat(sprintf("\nINCREMENTAL  M0=%.3f M1=%.3f dAUC=%+.3f p=%.3g | oMCC M0=%.3f M1=%.3f dMCC=%+.3f CI[%+.3f,%+.3f] P(>0)=%.3f\n",
  A0,A1,A1-A0,dp,m0,m1,m1-m0,ci[1],ci[2],mean(db>0)))
res<-data.frame(dataset="opsy HEXACO",J=240,n=nrow(df),npos=sum(df$y),gt="self-report V1/V2",
  AUC_rr=round(a_rr,3),AUC_PT=round(a_pt,3),AUC_ens=round(a_ens,3),AUC_M0=round(A0,3),AUC_M1=round(A1,3),
  dAUC=round(A1-A0,3),delong_p=signif(dp,3),oMCC_M0=round(m0,3),oMCC_M1=round(m1,3),
  dMCC=round(m1-m0,3),dMCC_ci_lo=round(ci[1],3),dMCC_ci_hi=round(ci[2],3),P_dMCC_gt0=round(mean(db>0),3))
write.csv(res,"C:/Users/vitto/Desktop/ReReReRe/rr_incremental_hexaco.csv",row.names=FALSE)
cat("Saved rr_incremental_hexaco.csv\n")
