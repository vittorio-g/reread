###############################################################################
# Per-Factor Coherence vs Standard Coupled — Real Data Validation
###############################################################################
setwd("C:/Users/vitto/Desktop/ReReReRe")
library(dplyr); library(pROC); library(haven); library(psych); library(lavaan)
source("ReReReRe.R")

compute_mcc <- function(tp,tn,fp,fn) {
  tp<-as.double(tp);tn<-as.double(tn);fp<-as.double(fp);fn<-as.double(fn)
  num<-(tp*tn)-(fp*fn); den<-sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  if(is.na(den)||den==0) return(0); num/den
}
evaluate <- function(pred,true) {
  tp<-sum(pred==1&true==1,na.rm=T);tn<-sum(pred==0&true==0,na.rm=T)
  fp<-sum(pred==1&true==0,na.rm=T);fn<-sum(pred==0&true==1,na.rm=T)
  list(MCC=compute_mcc(tp,tn,fp,fn))
}
compute_auc <- function(s,l,d) tryCatch(as.numeric(auc(roc(l,s,direction=d,quiet=T))),error=function(e)NA)
find_best_z <- function(scores,labels,thresholds=seq(0.1,3,.1),lower_is_bad=TRUE) {
  bm<- -Inf;bz<-1.5
  for(z in thresholds){
    pred<-if(lower_is_bad) as.integer(scores<=z) else as.integer(scores>=z)
    m<-evaluate(pred,labels)$MCC; if(m>bm){bm<-m;bz<-z}
  }
  list(z=bz,mcc=bm)
}

# Per-factor function (same as simulation script)
ReReReRe_perFactor <- function(data, iterations=100, align_signs=TRUE, min_items_per_factor=3) {
  data <- data[, sapply(data, is.numeric), drop=FALSE]
  mat <- as.matrix(data); N <- nrow(mat); p <- ncol(mat)
  item_max <- apply(mat, 2, max, na.rm=TRUE)

  pa <- tryCatch(suppressMessages(suppressWarnings(
    fa.parallel(mat, fa="fa", plot=FALSE, n.iter=20))), error=function(e) NULL)
  nF_est <- if(!is.null(pa)&&!is.null(pa$nfact)&&pa$nfact>=1) pa$nfact
            else max(1,sum(eigen(cor(mat,use="pairwise.complete.obs"),
                                  symmetric=TRUE,only.values=TRUE)$values>1))
  nF_est <- max(2, min(nF_est, floor(p/2)))

  efa <- tryCatch(suppressWarnings(
    fa(mat, nfactors=nF_est, rotate="oblimin", fm="minres", scores="none", warnings=FALSE)),
    error=function(e) tryCatch(suppressWarnings(
      fa(mat, nfactors=max(1,nF_est-1), rotate="oblimin", fm="minres",
         scores="none", warnings=FALSE)), error=function(e2) NULL))

  if(is.null(efa)) return(NULL)

  loadings_mat <- as.matrix(efa$loadings[])
  primary_factor <- apply(abs(loadings_mat), 1, which.max)
  raw_cor_mat <- cor(mat, use="pairwise.complete.obs")
  signMat <- sign(raw_cor_mat); signMat[upper.tri(signMat, diag=TRUE)] <- NA
  nF_actual <- ncol(loadings_mat)
  factor_z <- matrix(NA_real_, N, nF_actual)
  factors_used <- 0

  for(f in seq_len(nF_actual)) {
    items_f <- which(primary_factor == f)
    if(length(items_f) < min_items_per_factor) next
    combos <- combn(items_f, 2)
    idx_A <- combos[1,]; idx_B <- combos[2,]; k_f <- length(idx_A)
    pair_r <- sapply(seq_len(k_f), function(j) raw_cor_mat[idx_A[j], idx_B[j]])
    pair_abs_r <- abs(pair_r); pair_sign <- sign(pair_r)
    keep <- !is.na(pair_abs_r) & pair_abs_r > 0
    if(sum(keep)<2) next
    idx_A<-idx_A[keep]; idx_B<-idx_B[keep]; pair_abs_r<-pair_abs_r[keep]
    pair_sign<-pair_sign[keep]; k_f<-length(idx_A); weights<-pair_abs_r
    A_f<-mat[,idx_A,drop=FALSE]; B_f<-mat[,idx_B,drop=FALSE]
    if(align_signs){nf<-which(pair_sign<0)
      if(length(nf)>0) B_f[,nf]<-rep(item_max[idx_B[nf]]+1,each=N)-B_f[,nf]}
    coupled_f <- rowCor_weighted(A_f, B_f, weights)
    rand_f <- matrix(NA_real_, N, iterations)
    for(iter in seq_len(iterations)){
      ri1<-sample(p,k_f,replace=TRUE); ri2<-sample(p,k_f,replace=TRUE)
      same<-ri1==ri2; while(any(same)){ri2[same]<-sample(p,sum(same),replace=TRUE);same<-ri1==ri2}
      Ar<-mat[,ri1,drop=FALSE]; Br<-mat[,ri2,drop=FALSE]
      if(align_signs){ri<-pmax(ri1,ri2);ci<-pmin(ri1,ri2);rs<-signMat[cbind(ri,ci)];rs[is.na(rs)]<-1
        rn<-which(rs<0); if(length(rn)>0) Br[,rn]<-rep(item_max[ri2[rn]]+1,each=N)-Br[,rn]}
      rand_f[,iter] <- rowCor_weighted(Ar, Br, weights)
    }
    rm_f<-rowMeans(rand_f,na.rm=TRUE); rsd_f<-apply(rand_f,1,sd,na.rm=TRUE); rsd_f[rsd_f==0]<-1e-10
    factor_z[,f] <- (coupled_f - rm_f) / rsd_f
    factors_used <- factors_used + 1
  }

  if(factors_used==0) return(NULL)
  valid_cols <- which(colSums(!is.na(factor_z)) > 0)
  fz <- factor_z[, valid_cols, drop=FALSE]

  list(mean_z=rowMeans(fz,na.rm=T), min_z=apply(fz,1,min,na.rm=T),
       prop_low=rowMeans(fz<1.0,na.rm=T), combined=rowMeans(fz,na.rm=T)*0.7+apply(fz,1,min,na.rm=T)*0.3,
       nF_detected=nF_est, n_factors_used=factors_used)
}

# ============================================================================
# RUN ON ALL DATASETS
# ============================================================================

cat("\n============================================================\n")
cat("PER-FACTOR vs STANDARD COUPLED: Real Data\n")
cat("============================================================\n\n")

run_dataset <- function(name, items, labels) {
  cat(sprintf("\n=== %s ===\n", name))
  cat(sprintf("  N=%d, items=%d, careless=%.1f%%\n", nrow(items), ncol(items), mean(labels)*100))

  # Standard coupled
  cat("  Running standard coupled... ")
  rr <- tryCatch(ReReReRe(items, corProp=0.03, iterations=100, align_signs=TRUE, mode="coupled"),
                  error=function(e){cat(sprintf("ERROR: %s\n",e$message));NULL})
  if(!is.null(rr)){
    auc_std <- compute_auc(rr$z_score, labels, ">")
    mcc15_std <- evaluate(as.integer(rr$z_score<=1.5), labels)$MCC
    best_std <- find_best_z(rr$z_score, labels)
    cat(sprintf("AUC=%.3f MCC15=%.3f Oracle=%.3f\n", auc_std, mcc15_std, best_std$mcc))
  } else { auc_std<-NA; mcc15_std<-NA; best_std<-list(mcc=NA) }

  # Per-factor
  cat("  Running per-factor... ")
  pf <- tryCatch(ReReReRe_perFactor(items, iterations=100), error=function(e){cat(sprintf("ERROR\n"));NULL})
  if(!is.null(pf)){
    auc_mean <- compute_auc(pf$mean_z, labels, ">")
    auc_min <- compute_auc(pf$min_z, labels, ">")
    auc_comb <- compute_auc(pf$combined, labels, ">")
    b_mean <- find_best_z(pf$mean_z, labels)
    b_min <- find_best_z(pf$min_z, labels, seq(-2,3,.1))
    b_comb <- find_best_z(pf$combined, labels, seq(-1,3,.1))
    b_prop <- find_best_z(pf$prop_low, labels, seq(0.1,1,.05), lower_is_bad=FALSE)
    cat(sprintf("nF=%d fUsed=%d\n", pf$nF_detected, pf$n_factors_used))
    cat(sprintf("    mean: AUC=%.3f Oracle=%.3f\n", auc_mean, b_mean$mcc))
    cat(sprintf("    min:  AUC=%.3f Oracle=%.3f\n", auc_min, b_min$mcc))
    cat(sprintf("    comb: AUC=%.3f Oracle=%.3f\n", auc_comb, b_comb$mcc))
    cat(sprintf("    prop: Oracle=%.3f\n", b_prop$mcc))
  } else { auc_mean<-NA;auc_min<-NA;auc_comb<-NA
    b_mean<-list(mcc=NA);b_min<-list(mcc=NA);b_comb<-list(mcc=NA);b_prop<-list(mcc=NA) }

  data.frame(Dataset=name, Items=ncol(items), N=nrow(items),
    std_AUC=auc_std, std_MCC15=mcc15_std, std_oracle=best_std$mcc,
    pf_mean_AUC=auc_mean, pf_mean_oracle=b_mean$mcc,
    pf_min_AUC=auc_min, pf_min_oracle=b_min$mcc,
    pf_comb_AUC=auc_comb, pf_comb_oracle=b_comb$mcc,
    pf_prop_oracle=b_prop$mcc)
}

all_res <- list()

# 1. Schroeders
d<-read.csv("external_datasets/01_Schroeders_2022/data_mod_resp.csv",sep=";")
lab<-d$Careless; items<-d%>%select(starts_with("HE"))
cc<-complete.cases(items); all_res[[1]]<-run_dataset("Schroeders",items[cc,],lab[cc])

# 2. Schneider
d<-read.csv("external_datasets/06_Schneider_QoL/carersp.csv")
lab<-d$c01; ic<-grep("^(dep|pain|cog|fat)\\d",names(d),value=T); items<-d[,ic]
cc<-complete.cases(items)&!is.na(lab); all_res[[2]]<-run_dataset("Schneider",items[cc,],lab[cc])

# 3. Niessen
d<-read_sav("external_datasets/08_Niessen_2016/Raw_data.sav")%>%filter(Use_me==1)
lab<-as.integer(d$Conditon); ic<-grep("^[ECANO]\\d+$",names(d),value=T)
items<-data.frame(lapply(as.data.frame(d[,ic]),as.numeric))
cc<-complete.cases(items); all_res[[3]]<-run_dataset("Niessen",items[cc,],lab[cc])

# 4-6. Goldammer
for(si in 1:3){
  g<-read.csv(sprintf("external_datasets/13_Goldammer_2024/Study_%d.csv",si))
  if(si<=2){lab<-as.integer(g$careless_all)}else{lab<-as.integer(g$condition>0)}
  ic<-grep("^[ecano]_[a-z]+[0-9]+$",names(g),value=T)
  if(length(ic)==0) ic<-grep("^[ecano]_[a-z]+[0-9]+_t1$",names(g),value=T)
  items<-g[,ic]; cc<-complete.cases(items)&!is.na(lab)
  all_res[[si+3]]<-run_dataset(sprintf("Goldammer_S%d",si),items[cc,],lab[cc])
}

# Summary
final<-do.call(rbind,all_res)
write.csv(final,"perfactor_realdata_results.csv",row.names=FALSE)

cat("\n\n============================================================\n")
cat("SUMMARY\n")
cat("============================================================\n\n")

cat("--- AUC ---\n")
cat(sprintf("%-15s %4s  %6s  %6s  %6s  %6s\n","Dataset","It","Std","PF_mean","PF_min","PF_comb"))
for(i in 1:nrow(final)){r<-final[i,]
  cat(sprintf("%-15s %4d  %.3f   %.3f   %.3f   %.3f\n",
      r$Dataset,r$Items,r$std_AUC,r$pf_mean_AUC,r$pf_min_AUC,r$pf_comb_AUC))}

cat("\n--- Oracle MCC ---\n")
cat(sprintf("%-15s %4s  %6s  %6s  %6s  %6s  %6s\n","Dataset","It","Std","PF_mean","PF_min","PF_comb","PF_prop"))
for(i in 1:nrow(final)){r<-final[i,]
  cat(sprintf("%-15s %4d  %.3f   %.3f   %.3f   %.3f   %.3f\n",
      r$Dataset,r$Items,r$std_oracle,r$pf_mean_oracle,r$pf_min_oracle,r$pf_comb_oracle,r$pf_prop_oracle))}

cat("\n--- MEANS ---\n")
cat(sprintf("AUC:    Std=%.3f  PF_mean=%.3f  PF_min=%.3f  PF_comb=%.3f\n",
    mean(final$std_AUC,na.rm=T),mean(final$pf_mean_AUC,na.rm=T),
    mean(final$pf_min_AUC,na.rm=T),mean(final$pf_comb_AUC,na.rm=T)))
cat(sprintf("Oracle: Std=%.3f  PF_mean=%.3f  PF_min=%.3f  PF_comb=%.3f  PF_prop=%.3f\n",
    mean(final$std_oracle,na.rm=T),mean(final$pf_mean_oracle,na.rm=T),
    mean(final$pf_min_oracle,na.rm=T),mean(final$pf_comb_oracle,na.rm=T),
    mean(final$pf_prop_oracle,na.rm=T)))

cat("\n=== DONE ===\n")
