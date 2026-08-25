## interactions_test.R — does adding interaction terms to the 3-index logistic combiner help?
## Two questions, in order of importance:
##   (1) IN-SAMPLE / LOO on Study 1 (n=157): does LOO AUC/MCC improve?
##   (2) TRANSFER: freeze coefficients on Study 1, apply to the 12 external datasets
##       (same oriented features), compare mean AUC. Interactions are the first thing to
##       overfit at n=157, so (2) is the decisive test.
suppressMessages(library(pROC))
set.seed(1)

S <- read.csv("_study_oriented.csv")
y <- S$y
X <- S[, c("rr","longstring","person_total")]
names(X) <- c("rc","ls","pt")

mcc <- function(lab, pred){
  tp<-as.double(sum(pred==1&lab==1)); tn<-as.double(sum(pred==0&lab==0))
  fp<-as.double(sum(pred==1&lab==0)); fn<-as.double(sum(pred==0&lab==1))
  d<-sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)); if(!is.finite(d)||d==0) return(0); (tp*tn-fp*fn)/d
}
oracle_mcc <- function(lab,s){
  ok<-is.finite(s); lab<-lab[ok]; s<-s[ok]
  cuts<-unique(quantile(s,probs=seq(.001,.999,length.out=300),na.rm=TRUE)); best<--1
  for(c in cuts){m<-mcc(lab,as.integer(s>=c)); if(is.finite(m)&&m>best)best<-m}; best
}
aucv <- function(lab,s){ok<-is.finite(s); if(length(unique(lab[ok]))<2) return(NA)
  as.numeric(auc(roc(lab[ok],s[ok],direction="<",quiet=TRUE)))}

MODELS <- list(
  additive   = y ~ rc + ls + pt,
  two_way    = y ~ rc + ls + pt + rc:ls + rc:pt + ls:pt,
  full       = y ~ rc * ls * pt,                       # 2-way + 3-way
  rc_x_pt    = y ~ rc + ls + pt + rc:pt                # single most plausible interaction
)

cat("=== (1) Study 1: leave-one-out (combiner refit on n-1 each fold) ===\n")
D <- data.frame(y=y, X)
loo_res <- list()
for(nm in names(MODELS)){
  p <- numeric(nrow(D))
  for(i in seq_len(nrow(D))){
    fit <- suppressWarnings(glm(MODELS[[nm]], data=D[-i,], family=binomial()))
    p[i] <- predict(fit, newdata=D[i,], type="response")
  }
  a <- aucv(y,p); m <- oracle_mcc(y,p)
  npar <- length(coef(suppressWarnings(glm(MODELS[[nm]], data=D, family=binomial()))))
  loo_res[[nm]] <- c(auc=a, mcc=m, npar=npar)
  cat(sprintf("  %-10s npar=%-2d  LOO AUC=%.4f  LOO oracle-MCC=%.3f\n", nm, npar, a, m))
}
## DeLong test: additive vs each interaction model (on LOO predictions)
cat("\n  DeLong (LOO predictions, additive vs each):\n")
p_add <- numeric(nrow(D))
for(i in seq_len(nrow(D))){ f<-suppressWarnings(glm(MODELS$additive,data=D[-i,],family=binomial())); p_add[i]<-predict(f,newdata=D[i,],type="response") }
for(nm in setdiff(names(MODELS),"additive")){
  p2 <- numeric(nrow(D))
  for(i in seq_len(nrow(D))){ f<-suppressWarnings(glm(MODELS[[nm]],data=D[-i,],family=binomial())); p2[i]<-predict(f,newdata=D[i,],type="response") }
  r1<-roc(y,p_add,direction="<",quiet=TRUE); r2<-roc(y,p2,direction="<",quiet=TRUE)
  pv<-tryCatch(roc.test(r1,r2,method="delong",paired=TRUE)$p.value,error=function(e)NA)
  cat(sprintf("    additive vs %-10s  p=%.3f\n", nm, pv))
}

cat("\n=== (2) TRANSFER: fit on ALL of Study 1, apply frozen to 12 external datasets ===\n")
fits <- lapply(MODELS, function(f) suppressWarnings(glm(f, data=D, family=binomial())))
idx <- read.csv(file.path("ext_bench","index.csv"), stringsAsFactors=FALSE)
rows <- list()
for(k in seq_len(nrow(idx))){
  nm <- idx$name[k]; f <- file.path("ext_bench_oriented", paste0(nm,"_y.csv"))
  if(!file.exists(f)) next
  E <- read.csv(f); names(E)[names(E)=="rr"]<-"rc"; names(E)[names(E)=="longstring"]<-"ls"; names(E)[names(E)=="person_total"]<-"pt"
  au <- sapply(names(MODELS), function(mn) aucv(E$y, predict(fits[[mn]], newdata=E, type="response")))
  rows[[length(rows)+1]] <- data.frame(dataset=nm, J=idx$J[k], t(au))
  cat(sprintf("  %-14s J=%-3d | %s\n", nm, idx$J[k],
      paste(sprintf("%s %.3f", names(MODELS), au), collapse="  ")))
}
T <- do.call(rbind, rows)
cat("\n  MEAN external AUC:\n")
for(nm in names(MODELS)) cat(sprintf("    %-10s %.4f\n", nm, mean(T[[nm]], na.rm=TRUE)))
## paired Wilcoxon: additive vs each, across datasets
cat("\n  Paired Wilcoxon across the 12 external datasets (additive vs each):\n")
for(nm in setdiff(names(MODELS),"additive")){
  w <- suppressWarnings(wilcox.test(T$additive, T[[nm]], paired=TRUE))
  cat(sprintf("    additive vs %-10s  wins=%d/%d  p=%.3f\n", nm,
      sum(T$additive > T[[nm]], na.rm=TRUE), sum(is.finite(T$additive)), w$p.value))
}
write.csv(T, "interactions_external.csv", row.names=FALSE)
cat("\n=== coefficients of the 2-way model (are interactions even significant?) ===\n")
print(summary(fits$two_way)$coefficients, digits=3)
cat("\nwrote interactions_external.csv\n")
