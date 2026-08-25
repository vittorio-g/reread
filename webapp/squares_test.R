## squares_test.R — does a SUPERLINEAR (squared) contribution of each index help?
## Idea: "the higher an index, the more it should weigh". The oriented features are signed
## robust-z (negative = attentive, positive = careless), so a plain x^2 would destroy direction
## (a very attentive respondent would look careless). We therefore test four forms:
##   additive   : rc + ls + pt                               (shipped baseline)
##   quadratic  : + I(x^2) for each                          (free curvature, keeps linear term)
##   signed_sq  : replace x with sign(x)*x^2                 (direction kept, magnitude amplified)
##   pos_sq     : + pmax(x,0)^2 for each                     (amplify ONLY the careless side)
## Same two-stage evaluation as interactions_test.R: (1) LOO on Study 1, (2) frozen-coefficient
## transfer to the 12 external datasets (the decisive test).
suppressMessages(library(pROC))
set.seed(1)

S <- read.csv("_study_oriented.csv"); y <- S$y
mk <- function(d){
  d$rc2 <- d$rc^2; d$ls2 <- d$ls^2; d$pt2 <- d$pt^2
  d$rcS <- sign(d$rc)*d$rc^2; d$lsS <- sign(d$ls)*d$ls^2; d$ptS <- sign(d$pt)*d$pt^2
  d$rcP <- pmax(d$rc,0)^2;    d$lsP <- pmax(d$ls,0)^2;    d$ptP <- pmax(d$pt,0)^2
  d
}
D <- mk(data.frame(y=y, rc=S$rr, ls=S$longstring, pt=S$person_total))

mcc <- function(lab,pred){tp<-as.double(sum(pred==1&lab==1));tn<-as.double(sum(pred==0&lab==0))
  fp<-as.double(sum(pred==1&lab==0));fn<-as.double(sum(pred==0&lab==1))
  d<-sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)); if(!is.finite(d)||d==0) return(0); (tp*tn-fp*fn)/d}
oracle_mcc <- function(lab,s){ok<-is.finite(s);lab<-lab[ok];s<-s[ok]
  cuts<-unique(quantile(s,probs=seq(.001,.999,length.out=300),na.rm=TRUE));best<--1
  for(c in cuts){m<-mcc(lab,as.integer(s>=c));if(is.finite(m)&&m>best)best<-m};best}
aucv <- function(lab,s){ok<-is.finite(s); if(length(unique(lab[ok]))<2) return(NA)
  as.numeric(auc(roc(lab[ok],s[ok],direction="<",quiet=TRUE)))}

MODELS <- list(
  additive  = y ~ rc + ls + pt,
  quadratic = y ~ rc + ls + pt + rc2 + ls2 + pt2,
  signed_sq = y ~ rcS + lsS + ptS,
  pos_sq    = y ~ rc + ls + pt + rcP + lsP + ptP
)

cat("=== (1) Study 1: leave-one-out ===\n")
preds <- list()
for(nm in names(MODELS)){
  p <- numeric(nrow(D))
  for(i in seq_len(nrow(D))){
    f <- suppressWarnings(glm(MODELS[[nm]], data=D[-i,], family=binomial()))
    p[i] <- predict(f, newdata=D[i,], type="response")
  }
  preds[[nm]] <- p
  npar <- length(coef(suppressWarnings(glm(MODELS[[nm]], data=D, family=binomial()))))
  cat(sprintf("  %-10s npar=%-2d  LOO AUC=%.4f  LOO oracle-MCC=%.3f\n", nm, npar, aucv(y,p), oracle_mcc(y,p)))
}
cat("\n  DeLong (LOO preds, additive vs each):\n")
r1 <- roc(y, preds$additive, direction="<", quiet=TRUE)
for(nm in setdiff(names(MODELS),"additive")){
  r2 <- roc(y, preds[[nm]], direction="<", quiet=TRUE)
  pv <- tryCatch(roc.test(r1,r2,method="delong",paired=TRUE)$p.value, error=function(e)NA)
  cat(sprintf("    additive vs %-10s p=%.3f\n", nm, pv))
}

cat("\n=== (2) TRANSFER: fit on all of Study 1, frozen -> 12 external datasets ===\n")
fits <- lapply(MODELS, function(f) suppressWarnings(glm(f, data=D, family=binomial())))
idx <- read.csv(file.path("ext_bench","index.csv"), stringsAsFactors=FALSE)
rows <- list()
for(k in seq_len(nrow(idx))){
  nm <- idx$name[k]; f <- file.path("ext_bench_oriented", paste0(nm,"_y.csv"))
  if(!file.exists(f)) next
  E <- read.csv(f)
  E <- mk(data.frame(y=E$y, rc=E$rr, ls=E$longstring, pt=E$person_total))
  au <- sapply(names(MODELS), function(mn) aucv(E$y, predict(fits[[mn]], newdata=E, type="response")))
  rows[[length(rows)+1]] <- data.frame(dataset=nm, J=idx$J[k], t(au))
  cat(sprintf("  %-14s J=%-3d | %s\n", nm, idx$J[k], paste(sprintf("%s %.3f", names(MODELS), au), collapse="  ")))
}
T <- do.call(rbind, rows)
cat("\n  MEAN external AUC:\n")
for(nm in names(MODELS)) cat(sprintf("    %-10s %.4f\n", nm, mean(T[[nm]], na.rm=TRUE)))
cat("\n  Paired Wilcoxon (additive vs each) across 12 external datasets:\n")
for(nm in setdiff(names(MODELS),"additive")){
  w <- suppressWarnings(wilcox.test(T$additive, T[[nm]], paired=TRUE))
  cat(sprintf("    additive vs %-10s wins=%d/%d p=%.3f\n", nm,
     sum(T$additive > T[[nm]], na.rm=TRUE), sum(is.finite(T$additive)), w$p.value))
}
write.csv(T, "squares_external.csv", row.names=FALSE)
cat("\n=== are the squared terms significant? (quadratic model) ===\n")
print(summary(fits$quadratic)$coefficients, digits=3)
cat("\n=== pos_sq model (careless-side amplification) ===\n")
print(summary(fits$pos_sq)$coefficients, digits=3)
cat("\nwrote squares_external.csv\n")
