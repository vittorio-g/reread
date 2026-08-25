suppressWarnings(suppressMessages(library(pROC)))
B <- "C:/Users/vitto/Desktop/ReReReRe/Dataset/gt_benchmark_candidates"
# fixed-orientation detection AUC: higher score = careless (case=1)
aucF <- function(y,s){ ok<-!is.na(s); as.numeric(pROC::auc(pROC::roc(y[ok],s[ok],quiet=TRUE,
          levels=c(0,1), direction="<"))) }
chk <- function(name, sfile, y, tgt){
  sc<-read.csv(file.path(B,sfile)); n<-min(nrow(sc),length(y)); sc<-sc[1:n,]; y<-y[1:n]
  cat(sprintf("%-12s n=%5d npos=%4d | ens=%.3f/%.3f  rr(careless)=%.3f/%.3f  rr(-z)=%.3f  PT=%.3f/%.3f\n",
    name,n,sum(y), aucF(y,sc$ens_p),tgt[1], aucF(y,sc$rr_careless),tgt[2], aucF(y,-sc$rr_z), aucF(y,sc$person_total),tgt[3]))
}
# smarvus (>=2 instructed): target ens .864 rr .865 PT .886
ys<-read.csv(file.path(B,"smarvus/smarvus_labels.csv"))[[1]]
chk("smarvus","smarvus/smarvus_scores.csv",ys,c(.864,.865,.886))
# 16pf (extreme speeders): target .828 .793 .826
y16<-read.csv(file.path(B,"opsy_16pf_labels.csv"))[[1]]
chk("16pf","opsy_16pf_scores.csv",y16,c(.828,.793,.826))
# warning (self-rep diligence): .822 .670 .803
yw<-read.csv(file.path(B,"warning_ipipneo300/_labels.csv"))[[1]]
chk("warning","warning_ipipneo300/_scores.csv",yw,c(.822,.670,.803))
# Kay S1 dur_speed: .934 .852 .955
g1<-read.csv(file.path(B,"kay_idris_idria/kay_gt_s1.csv"))
chk("KayS1-speed","kay_idris_idria/kay_scores_s1.csv",g1$dur_speed,c(.934,.852,.955))
# Kay S1 fake_informants (what I wrongly used): .566 .548 .628
chk("KayS1-fake","kay_idris_idria/kay_scores_s1.csv",g1$fake_informants,c(.566,.548,.628))
