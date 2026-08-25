## exp_efa_diag.R — is the concern real? How concentrated are rr's top-|r| coupled pairs
## across the study's TRUE factor domains? Plus baseline rr AUC. Study 1 (n=157, 109 items).
suppressMessages(library(pROC)); set.seed(1)
M   <- as.matrix(read.csv("_study_matrix.csv", check.names=FALSE))
it  <- read.csv("_study_items.csv", stringsAsFactors=FALSE)
lab <- read.csv("_study_labels.csv")[[1]]
imp <- function(m){for(j in 1:ncol(m)){v<-m[,j];mm<-median(v,na.rm=TRUE);if(is.na(mm))mm<-0;m[is.na(v),j]<-mm};m}
M <- imp(M); n<-nrow(M); J<-ncol(M)
imax <- apply(M,2,max); imax[imax==0]<-1; P <- sweep(M,2,imax,"/")
R <- cor(P); ut <- which(upper.tri(R),arr.ind=TRUE); rr_ <- R[upper.tri(R)]
k <- min(max(15,round(0.03*nrow(ut))), nrow(ut))
ord <- order(abs(rr_),decreasing=TRUE); sel <- ord[1:k]
a <- ut[sel,1]; b <- ut[sel,2]; da <- it$domain[a]; db <- it$domain[b]
cat("=== top-|r| coupled pairs: k =",k,"of",nrow(ut),"pairs, across", length(unique(it$domain)),"true domains ===\n\n")
within <- da==db
cat("within-domain pairs:",sum(within),"(",round(100*mean(within)),"%)  cross-domain:",sum(!within),"\n\n")
# how many distinct domains are represented; concentration
touch <- table(factor(c(da,db), levels=unique(it$domain)))
cat("pairs touching each domain (an endpoint in that domain):\n"); print(sort(touch,decreasing=TRUE))
cat("\ntop-domain share of endpoints:",round(100*max(touch)/sum(touch)),"%   domains with >=1 pair:",sum(touch>0),"/",length(touch),"\n")
# within-domain pair count per domain (the pairs a factor-stratified scheme would balance)
wd <- table(factor(da[within],levels=unique(it$domain)))
cat("\nwithin-domain pairs per domain:\n"); print(sort(wd,decreasing=TRUE))

## baseline rr (current selection) AUC, for reference
signv <- sign(rr_[sel]); minp <- (apply(M,2,min)/imax)
indcor <- function(rows){
  aa<-ut[rows,1];bb<-ut[rows,2];sg<-sign(rr_[rows]);A<-P[,aa,drop=FALSE];B<-P[,bb,drop=FALSE]
  neg<-which(sg<0);for(c in neg)B[,c]<-(1+minp[bb[c]])-B[,c]
  vapply(1:n,function(i){av<-A[i,];bv<-B[i,];if(sd(av)==0||sd(bv)==0)return(0);abs(cor(av,bv))},numeric(1))
}
z_of <- function(rows,B=200){co<-indcor(rows);rand<-matrix(NA,n,B);for(bb in 1:B)rand[,bb]<-indcor(sample(nrow(ut),length(rows)))
  (co-rowMeans(rand))/apply(rand,1,sd)}
z <- z_of(sel); z[!is.finite(z)] <- min(z[is.finite(z)])
au <- as.numeric(auc(roc(lab, -z, direction="<", quiet=TRUE)))   # low z = careless -> use -z so high=careless
cat("\nbaseline rr AUC (current top-|r| selection):", round(au,3), "\n")
saveRDS(list(ut=ut,rr_=rr_,k=k,P=P,minp=minp,it=it,lab=lab,n=n), "exp_efa_env.rds")
cat("(saved env for the stratified comparison)\n")
