## Task 2: for each Study-2 dataset, can RPR / even-odd be computed from the COMMITTED matrix?
## Criterion: the item identifiers must expose an item->subscale map with >=2 subscales,
## each with >=4 items (needed for a within-subscale split-half correlated across subscales).
D<-"../Dataset/"
sets<-list(
 c("Kay S2","gt_benchmark_candidates/kay_idris_idria/kay_matrix_s2.csv"),
 c("warning IPIP-NEO","gt_benchmark_candidates/warning_ipipneo300/_matrix.csv"),
 c("Kay S1","gt_benchmark_candidates/kay_idris_idria/kay_matrix_s1.csv"),
 c("opsy 16PF","gt_benchmark_candidates/opsy_16pf_matrix.csv"),
 c("smarvus","gt_benchmark_candidates/smarvus/smarvus_matrix.csv"),
 c("Kay S6","gt_benchmark_candidates/kay_idris_idria/kay_matrix_s6.csv"),
 c("Kay S5","gt_benchmark_candidates/kay_idris_idria/kay_matrix_s5.csv"),
 c("Duckworth VCL","gt_benchmark_candidates/duckworth_grit_vcl/duckworth_matrix.csv"),
 c("Douglas 2023","gt_benchmark_candidates/douglas2023_dataquality/douglas_matrix.csv"),
 c("Krause youth","gt_benchmark_candidates/krause_ier_youth/krause_matrix.csv"),
 c("Mastroianni 2022","learning_to_pay_attention/mastroianni2022_attitude/_matrix.csv"),
 c("Ivanov 2021","learning_to_pay_attention/ivanov2021_racialresentment/_matrix.csv"),
 c("O'Grady 2019","learning_to_pay_attention/ogrady2019_moralfoundations/_matrix.csv"),
 c("Pennycook 2020","learning_to_pay_attention/pennycook2020_covid/_matrix.csv"),
 c("Moss 2023","learning_to_pay_attention/moss2023_ethical/_matrix.csv"),
 c("Alvarez 2019","learning_to_pay_attention/alvarez2019_inattentive/_matrix.csv")
)
stem<-function(x){ x<-gsub('"','',x)
  x<-sub("[._]?[0-9]+([._]?[rRxX])?$","",x)   # strip trailing item number (+ optional reverse/x marker)
  x<-sub("[._]+$","",x); x }
out<-list()
for(s in sets){ nm<-s[1]; h<-readLines(file.path(D,s[2]),n=1)
  cols<-strsplit(h,",")[[1]]; J<-length(cols)
  st<-stem(cols); tb<-sort(table(st),decreasing=TRUE)
  n_sub<-length(tb); n_sub4<-sum(tb>=4)
  top<-paste(head(paste0(names(tb),"(",as.integer(tb),")"),6),collapse=" ")
  out[[length(out)+1]]<-data.frame(dataset=nm,J=J,n_stems=n_sub,stems_ge4=n_sub4,example=top)
  cat(sprintf("%-17s J=%3d  stems=%2d  stems>=4items=%2d | %s\n",nm,J,n_sub,n_sub4,top))
}
res<-do.call(rbind,out); write.csv(res,"rpr_applicability_study2.csv",row.names=FALSE)
cat("\nHeuristic 'exposed' = stems_ge4 >= 2 (>=2 named subscales of >=4 items recoverable from headers).\n")
cat("Datasets meeting heuristic:", sum(res$stems_ge4>=2), "of", nrow(res),
    "-> RPR/even-odd NOT computable on", sum(res$stems_ge4<2), "of", nrow(res), "(mechanical).\n")
