## Study 1 ablation INCLUDING resampled personal reliability (RPR, Goldammer 2024).
## Question: is rr redundant with its closest rival RPR, or complementary?
## LOO logistic per model; AUC (DeLong CI) + oracle-MCC; incremental DeLong tests.
suppressMessages({library(careless); library(pROC)})
set.seed(1)
mat0 <- as.matrix(read.csv("_study_matrix.csv", check.names=FALSE))
lab  <- read.csv("_study_labels.csv")[[1]]
it   <- read.csv("_study_items.csv", stringsAsFactors=FALSE)
sc   <- read.csv("_study_scores.csv")               # oriented shipped features rr/longstring/person_total
stopifnot(nrow(mat0)==length(lab), nrow(sc)==length(lab))

impute<-function(m){for(j in seq_len(ncol(m))){v<-m[,j];mm<-median(v,na.rm=TRUE);if(is.na(mm))mm<-0;m[is.na(v),j]<-mm};m}
mat<-impute(mat0); matR<-mat
for(j in which(it$reverse==1)) matR[,j]<-(it$rmax[j]+it$rmin[j])-mat[,j]
domains<-unique(it$domain[it$domain!="Random"])

## ---- resampled personal reliability (B random within-domain split-halves) ----
B<-200; prr_acc<-matrix(0,nrow(mat),B)
for(b in seq_len(B)){ A<-matrix(NA,nrow(mat),length(domains));Bm<-A
  for(di in seq_along(domains)){cols<-which(it$domain==domains[di]);cols<-sample(cols)
    h<-floor(length(cols)/2); if(h<1) next; a<-cols[seq_len(h)]; bb<-cols[(h+1):length(cols)]
    ma<-rowMeans(matR[,a,drop=FALSE]);mb<-rowMeans(matR[,bb,drop=FALSE])
    A[,di]<-(ma-mean(ma))/(sd(ma)+1e-9);Bm[,di]<-(mb-mean(mb))/(sd(mb)+1e-9)}
  prr_acc[,b]<-vapply(seq_len(nrow(mat)),function(i){x<-A[i,];y<-Bm[i,];ok<-is.finite(x)&is.finite(y)
    if(sum(ok)<3||sd(x[ok])==0||sd(y[ok])==0)return(NA_real_);suppressWarnings(cor(x[ok],y[ok]))},numeric(1))}
prr_v<-rowMeans(prr_acc,na.rm=TRUE)          # high = coherent (attentive)

## feature frame, all oriented so HIGH = MORE careless
D<-data.frame(y=as.integer(lab), rr=sc$rr, ls=sc$longstring, pt=sc$person_total, rpr=-prr_v)
D$rpr[is.na(D$rpr)]<-median(D$rpr,na.rm=TRUE)

aucd<-function(y,s) as.numeric(pROC::auc(pROC::roc(y,s,levels=c(0,1),direction="<",quiet=TRUE)))
loo<-function(form){n<-nrow(D);p<-numeric(n)
  for(i in seq_len(n)){m<-suppressWarnings(glm(form,D[-i,,drop=FALSE],family=binomial()))
    p[i]<-predict(m,D[i,,drop=FALSE],type="response")};p}
mcc<-function(y,p){tp<-sum(p&y==1);tn<-sum(!p&y==0);fp<-sum(p&y==0);fn<-sum(!p&y==1)
  d<-sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));if(d==0)0 else (tp*tn-fp*fn)/d}
omcc<-function(y,s){th<-quantile(s,seq(.02,.98,.02),names=FALSE);max(sapply(th,function(t)mcc(y,s>=t)))}
dl<-function(pa,pb){r1<-roc(D$y,pa,direction="<",quiet=TRUE);r2<-roc(D$y,pb,direction="<",quiet=TRUE)
  tryCatch(roc.test(r1,r2,method="delong",paired=TRUE)$p.value,error=function(e)NA)}

models<-list(
  "rr alone"                 = y~rr,
  "RPR alone"                = y~rpr,
  "LongString + PT"          = y~ls+pt,
  "RPR + LongString + PT"    = y~rpr+ls+pt,
  "rr + LongString + PT (shipped)" = y~rr+ls+pt,
  "rr + RPR + LongString + PT"     = y~rr+rpr+ls+pt
)
P<-lapply(models,loo)
ci<-function(p){r<-roc(D$y,p,direction="<",quiet=TRUE);a<-ci.auc(r,method="delong");sprintf("%.3f [%.3f, %.3f]",as.numeric(auc(r)),a[1],a[3])}
tab<-data.frame(model=names(models),
  AUC=sapply(P,function(p)round(aucd(D$y,p),3)),
  AUC_CI=sapply(P,ci),
  oracle_MCC=sapply(P,function(p)round(omcc(D$y,p),3)), row.names=NULL)
cat("=== Study 1 ablation with RPR (LOO, n=",nrow(D),") ===\n",sep=""); print(tab,row.names=FALSE)
write.csv(tab,"ablation_rpr_study1.csv",row.names=FALSE)

cat("\n--- Is rr redundant with RPR? (key incremental tests) ---\n")
cat(sprintf("rr alone AUC=%.3f  vs  RPR alone AUC=%.3f  (DeLong p=%.3g)\n",
  aucd(D$y,D$rr),aucd(D$y,D$rpr),dl(D$rr,D$rpr)))
cat(sprintf("add rr to {RPR+LS+PT}:  dAUC=%+.3f  (DeLong p=%.3g)  dMCC=%+.3f\n",
  aucd(D$y,P[["rr + RPR + LongString + PT"]])-aucd(D$y,P[["RPR + LongString + PT"]]),
  dl(P[["rr + RPR + LongString + PT"]],P[["RPR + LongString + PT"]]),
  omcc(D$y,P[["rr + RPR + LongString + PT"]])-omcc(D$y,P[["RPR + LongString + PT"]])))
cat(sprintf("add RPR to {rr+LS+PT} (shipped): dAUC=%+.3f  (DeLong p=%.3g)  dMCC=%+.3f\n",
  aucd(D$y,P[["rr + RPR + LongString + PT"]])-aucd(D$y,P[["rr + LongString + PT (shipped)"]]),
  dl(P[["rr + RPR + LongString + PT"]],P[["rr + LongString + PT (shipped)"]]),
  omcc(D$y,P[["rr + RPR + LongString + PT"]])-omcc(D$y,P[["rr + LongString + PT (shipped)"]])))
cat(sprintf("correlation(rr, RPR) = %.3f\n", cor(D$rr,D$rpr)))
cat("\nwrote ablation_rpr_study1.csv\n")
