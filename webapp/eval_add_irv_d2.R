## eval_add_irv_d2.R — do IRV and Mahalanobis D^2 earn a place in the ensemble?
## For every external dataset and every ground truth:
##   standalone AUC of IRV and D^2, and 10-fold CV AUC of
##      shipped triad {rc, longstring, person_total}   vs   triad + IRV + D^2
## CV is used so the 5-feature model is not rewarded merely for having more parameters.
suppressMessages(library(pROC)); set.seed(11)
B <- "../Dataset/gt_benchmark_candidates"
sets <- list(
  list("TISP", file.path(B,"tisp"), c("y1")),
  list("TUMI", file.path(B,"tumi"), c("y1","y2")),
  list("warnIPIP", file.path(B,"warning_ipipneo300"), c("y1","y2")),
  list("smarvus", file.path(B,"smarvus"), c("y1","y2")),
  list("KayS2", file.path(B,"kay_idris_idria/s2"), c("y1","y2")),
  list("KayS1", file.path(B,"kay_idris_idria/s1"), c("y1","y2")))
aucv <- function(y,s) as.numeric(pROC::auc(pROC::roc(y,s,quiet=TRUE,direction="auto")))
oof <- function(d, form, K=10, reps=5){
  n<-nrow(d); acc<-numeric(n)
  for(rp in 1:reps){
    i0<-which(d$y==0); i1<-which(d$y==1)
    fold<-integer(n); fold[i0]<-sample(rep(1:K,length.out=length(i0))); fold[i1]<-sample(rep(1:K,length.out=length(i1)))
    p<-numeric(n)
    for(k in 1:K){ tr<-fold!=k; te<-fold==k
      m<-suppressWarnings(glm(form,data=d[tr,],family=binomial()))
      p[te]<-predict(m,newdata=d[te,],type="response") }
    acc<-acc+p }
  acc/reps }
out<-list()
cat(sprintf("%-9s %-3s %6s %6s | %6s %6s | %7s %7s %7s\n",
            "dataset","GT","n","prev","AUCirv","AUCd2","triad","+irv+d2","delta"))
for(s in sets){ nm<-s[[1]]; d0<-read.csv(file.path(s[[2]],"_feats.csv"))
  for(yc in s[[3]]){
    d<-data.frame(y=d0[[yc]], rc=d0$rc, ls=d0$longstring, pt=d0$person_total, irv=d0$irv, d2=d0$d2)
    d<-d[complete.cases(d),]
    if(length(unique(d$y))<2) next
    a_irv<-aucv(d$y,d$irv); a_d2<-aucv(d$y,d$d2)
    p3<-oof(d, y~rc+ls+pt); p5<-oof(d, y~rc+ls+pt+irv+d2)
    A3<-aucv(d$y,p3); A5<-aucv(d$y,p5)
    r3<-roc(d$y,p3,quiet=TRUE,direction="auto"); r5<-roc(d$y,p5,quiet=TRUE,direction="auto")
    pv<-tryCatch(roc.test(r5,r3,method="delong",paired=TRUE)$p.value, error=function(e)NA)
    cat(sprintf("%-9s %-3s %6d %5.1f%% | %6.3f %6.3f | %7.3f %7.3f %+7.3f  p=%.3g\n",
                nm,yc,nrow(d),100*mean(d$y),a_irv,a_d2,A3,A5,A5-A3,pv))
    out[[length(out)+1]]<-data.frame(dataset=nm,gt=yc,n=nrow(d),prev=mean(d$y),
      auc_irv=a_irv,auc_d2=a_d2,auc_triad=A3,auc_plus=A5,delta=A5-A3,delong_p=pv)
  }}
res<-do.call(rbind,out); write.csv(res,"eval_add_irv_d2.csv",row.names=FALSE)
cat(sprintf("\nmean delta = %+.4f  | significant improvements: %d/%d\n",
            mean(res$delta), sum(res$delong_p<.05 & res$delta>0, na.rm=TRUE), nrow(res)))
