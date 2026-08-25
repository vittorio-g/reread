## Point 8 — weight over-precision / overfitting, on Study 1 (home data, N=157).
## (c) bootstrap CI on coefficients; (a) rounding & perturbation robustness; (b) fixed vs refit (LOO).
suppressMessages(library(pROC))
set.seed(1)
S <- read.csv("_study_scores.csv"); L <- read.csv("_study_labels.csv")
n <- min(nrow(S), nrow(L)); S <- S[1:n,]; y <- as.integer(L$careless[1:n])
df <- data.frame(y=y, rr=S$rr, longstring=S$longstring, person_total=S$person_total)
W <- c(b0=-0.1057, rr=2.7645, longstring=1.2689, person_total=1.4452)
aucF <- function(y,s) as.numeric(pROC::auc(pROC::roc(y,s,levels=c(0,1),direction="<",quiet=TRUE)))
eta <- function(w,d) w["b0"] + w["rr"]*d$rr + w["longstring"]*d$longstring + w["person_total"]*d$person_total

cat("=== (validation) does a plain refit reproduce the shipped weights? ===\n")
fit <- glm(y ~ rr + longstring + person_total, df, family=binomial())
print(round(rbind(shipped=W, refit=coef(fit)),4))
cat(sprintf("AUC shipped weights = %.4f | AUC in-sample refit = %.4f\n\n", aucF(df$y,eta(W,df)), aucF(df$y,fitted(fit))))

cat("=== (c) bootstrap 95%% CI on coefficients (2000 resamples) ===\n")
B<-2000; C<-matrix(NA,B,4,dimnames=list(NULL,names(W)))
for(b in 1:B){ i<-sample(nrow(df),replace=TRUE)
  fb<-suppressWarnings(tryCatch(glm(y~rr+longstring+person_total,df[i,],family=binomial(),
        control=glm.control(maxit=50)),error=function(e)NULL))
  if(!is.null(fb)) C[b,]<-coef(fb) }
for(k in names(W)){ q<-quantile(C[,k],c(.025,.5,.975),na.rm=TRUE)
  cat(sprintf("  %-13s shipped=%+.4f | boot median=%+.3f  95%%CI [%+.3f, %+.3f]  (width %.2f)\n",
    k, W[k], q[2], q[1], q[3], q[3]-q[1])) }

cat("\n=== (a) rounding robustness: recompute ensemble with coarser weights ===\n")
variants <- list(
  shipped   = W,
  one_sigfig= c(b0=-0.2, rr=3, longstring=1, person_total=1),
  integers  = round(W),
  rr_plus_equal_aux = c(b0=0, rr=3, longstring=1.5, person_total=1.5),
  all_equal = c(b0=0, rr=1, longstring=1, person_total=1)   # unweighted sum of z's
)
for(nm in names(variants)){ w<-variants[[nm]]; cat(sprintf("  %-18s AUC=%.4f   weights=(%s)\n",
  nm, aucF(df$y, eta(w,df)), paste(sprintf("%.2g",w),collapse=", "))) }

cat("\n=== (a) perturbation: multiplicative jitter on the 3 slopes ===\n")
for(s in c(0.10,0.20,0.50)){ A<-numeric(500)
  for(j in 1:500){ w<-W; w[2:4]<-W[2:4]*(1+runif(3,-s,s)); A[j]<-aucF(df$y,eta(w,df)) }
  cat(sprintf("  +/-%2.0f%%: AUC mean=%.4f  sd=%.4f  min=%.4f  (shipped %.4f)\n",
    s*100, mean(A), sd(A), min(A), aucF(df$y,eta(W,df)))) }

cat("\n=== (b) freezing cost: LOO fixed weights vs LOO per-fold refit ===\n")
pf<-numeric(nrow(df))               # fixed
for(i in 1:nrow(df)) pf[i]<-eta(W,df[i,,drop=FALSE])   # fixed weights, no training
pr<-numeric(nrow(df))               # LOO refit
for(i in 1:nrow(df)){ m<-suppressWarnings(glm(y~rr+longstring+person_total,df[-i,],family=binomial()))
  pr[i]<-predict(m,newdata=df[i,,drop=FALSE],type="response") }
r_fix<-pROC::roc(df$y,pf,levels=c(0,1),direction="<",quiet=TRUE)
r_ref<-pROC::roc(df$y,pr,levels=c(0,1),direction="<",quiet=TRUE)
dp<-pROC::roc.test(r_ref,r_fix,method="delong",paired=TRUE)$p.value
cat(sprintf("  LOO AUC fixed weights = %.4f | LOO AUC per-fold refit = %.4f | dAUC=%+.4f (DeLong p=%.3f)\n",
  as.numeric(auc(r_fix)), as.numeric(auc(r_ref)), as.numeric(auc(r_ref))-as.numeric(auc(r_fix)), dp))
cat("\n(Study 1 N=157; eta is the shipped linear combiner, gate g=1 on this data.)\n")

## ---- dump results for plotting ----
coefdf <- data.frame(param=names(W), shipped=as.numeric(W),
  boot_lo=sapply(names(W),function(k)quantile(C[,k],.025,na.rm=TRUE)),
  boot_md=sapply(names(W),function(k)quantile(C[,k],.5,na.rm=TRUE)),
  boot_hi=sapply(names(W),function(k)quantile(C[,k],.975,na.rm=TRUE)), row.names=NULL)
write.csv(coefdf,"ext_bench/w8_coef_ci.csv",row.names=FALSE)
vardf <- data.frame(variant=names(variants),
  AUC=sapply(variants,function(w)aucF(df$y,eta(w,df))),
  label=c("shipped (2.68,1.24,1.44)","1 sig-fig (3,1,1)","integers (3,1,1)","rr + equal aux (3,1.5,1.5)","unweighted z-sum (1,1,1)"), row.names=NULL)
pert <- do.call(rbind, lapply(c(0.10,0.20,0.50), function(s){A<-numeric(1000)
  for(j in 1:1000){w<-W;w[2:4]<-W[2:4]*(1+runif(3,-s,s));A[j]<-aucF(df$y,eta(w,df))}
  data.frame(pct=s*100, mean=mean(A), lo=quantile(A,.025), hi=quantile(A,.975), min=min(A))}))
write.csv(vardf,"ext_bench/w8_variant_auc.csv",row.names=FALSE)
write.csv(pert,"ext_bench/w8_perturb.csv",row.names=FALSE)
write.csv(data.frame(kind=c("LOO fixed","LOO refit"),
  AUC=c(as.numeric(auc(r_fix)),as.numeric(auc(r_ref)))),"ext_bench/w8_loo.csv",row.names=FALSE)
cat("wrote ext_bench/w8_*.csv\n")
