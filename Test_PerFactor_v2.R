###############################################################################
# Per-Factor Coherence vs Standard Coupled — Simulated Data
#
# Compares:
#   1. Standard coupled (current default, corProp=0.03)
#   2. Per-factor: mean(z)
#   3. Per-factor: min(z)
#   4. Per-factor: combined (0.7*mean + 0.3*min)
#   5. Per-factor: prop_low (% factors with z < 1)
#
# Design: 8 nF x 3 ipf x 3 pct x 5 reps, N=300
###############################################################################

library(lavaan); library(psych); library(dplyr); library(ggplot2); library(tidyr)
source("Synthetic_Good_Responses_2.R")
source("Careless_machine_2.R")
source("ReReReRe.R")

evaluate_mcc <- function(predicted, actual) {
  tp<-sum(predicted==1&actual==1); tn<-sum(predicted==0&actual==0)
  fp<-sum(predicted==1&actual==0); fn<-sum(predicted==0&actual==1)
  d<-as.double(tp+fp)*as.double(tp+fn)*as.double(tn+fp)*as.double(tn+fn)
  if(d==0) return(0); (as.double(tp)*tn-as.double(fp)*fn)/sqrt(d)
}

find_best_z <- function(scores, labels, thresholds=seq(0.1,3,.1), lower_is_bad=TRUE) {
  bm<- -1; bz<-NA
  for(z in thresholds){
    pred<-if(lower_is_bad) as.integer(scores<=z) else as.integer(scores>=z)
    m<-evaluate_mcc(pred,labels); if(m>bm){bm<-m;bz<-z}
  }
  list(z=bz,mcc=bm)
}

# ============================================================================
# PER-FACTOR COHERENCE (uses coupled-style top-k% for pair selection)
# ============================================================================

ReReReRe_perFactor <- function(data, corProp=0.03, iterations=50, align_signs=TRUE,
                                min_items_per_factor=3) {
  data <- data[, sapply(data, is.numeric), drop=FALSE]
  mat <- as.matrix(data); N <- nrow(mat); p <- ncol(mat)
  item_max <- apply(mat, 2, max, na.rm=TRUE)

  # EFA to get factor assignments
  pa <- tryCatch(suppressMessages(suppressWarnings(
    fa.parallel(mat, fa="fa", plot=FALSE, n.iter=20))), error=function(e) NULL)
  nF_est <- if(!is.null(pa)&&!is.null(pa$nfact)&&pa$nfact>=1) pa$nfact
            else max(1,sum(eigen(cor(mat,use="pairwise.complete.obs"),
                                  symmetric=TRUE,only.values=TRUE)$values>1))
  nF_est <- max(2, min(nF_est, floor(p/2)))

  efa <- tryCatch(suppressWarnings(
    fa(mat, nfactors=nF_est, rotate="oblimin", fm="minres",
       scores="none", warnings=FALSE)),
    error=function(e) tryCatch(suppressWarnings(
      fa(mat, nfactors=max(1,nF_est-1), rotate="oblimin", fm="minres",
         scores="none", warnings=FALSE)), error=function(e2) NULL))

  if(is.null(efa)) return(list(mean_z=rep(NA,N), min_z=rep(NA,N),
    prop_low=rep(NA,N), combined=rep(NA,N), nF_detected=nF_est, n_factors_used=0))

  loadings_mat <- as.matrix(efa$loadings[])
  primary_factor <- apply(abs(loadings_mat), 1, which.max)

  raw_cor_mat <- cor(mat, use="pairwise.complete.obs")
  signMat <- sign(raw_cor_mat); signMat[upper.tri(signMat, diag=TRUE)] <- NA

  # Per-factor z-scores
  nF_actual <- ncol(loadings_mat)
  factor_z <- matrix(NA_real_, N, nF_actual)
  factors_used <- 0

  for(f in seq_len(nF_actual)) {
    items_f <- which(primary_factor == f)
    if(length(items_f) < min_items_per_factor) next

    combos <- combn(items_f, 2)
    idx_A <- combos[1,]; idx_B <- combos[2,]
    k_f <- length(idx_A)

    pair_r <- sapply(seq_len(k_f), function(j) raw_cor_mat[idx_A[j], idx_B[j]])
    pair_abs_r <- abs(pair_r); pair_sign <- sign(pair_r)
    keep <- !is.na(pair_abs_r) & pair_abs_r > 0
    if(sum(keep) < 2) next
    idx_A<-idx_A[keep]; idx_B<-idx_B[keep]
    pair_abs_r<-pair_abs_r[keep]; pair_sign<-pair_sign[keep]
    k_f<-length(idx_A); weights<-pair_abs_r

    A_f<-mat[,idx_A,drop=FALSE]; B_f<-mat[,idx_B,drop=FALSE]
    if(align_signs){
      nf<-which(pair_sign<0)
      if(length(nf)>0) B_f[,nf]<-rep(item_max[idx_B[nf]]+1,each=N)-B_f[,nf]
    }
    coupled_f <- rowCor_weighted(A_f, B_f, weights)

    rand_f <- matrix(NA_real_, N, iterations)
    for(iter in seq_len(iterations)){
      ri1<-sample(p,k_f,replace=TRUE); ri2<-sample(p,k_f,replace=TRUE)
      same<-ri1==ri2; while(any(same)){ri2[same]<-sample(p,sum(same),replace=TRUE);same<-ri1==ri2}
      Ar<-mat[,ri1,drop=FALSE]; Br<-mat[,ri2,drop=FALSE]
      if(align_signs){
        ri<-pmax(ri1,ri2);ci<-pmin(ri1,ri2);rs<-signMat[cbind(ri,ci)];rs[is.na(rs)]<-1
        rn<-which(rs<0); if(length(rn)>0) Br[,rn]<-rep(item_max[ri2[rn]]+1,each=N)-Br[,rn]
      }
      rand_f[,iter] <- rowCor_weighted(Ar, Br, weights)
    }
    rm_f<-rowMeans(rand_f,na.rm=TRUE); rsd_f<-apply(rand_f,1,sd,na.rm=TRUE)
    rsd_f[rsd_f==0]<-1e-10
    factor_z[,f] <- (coupled_f - rm_f) / rsd_f
    factors_used <- factors_used + 1
  }

  if(factors_used == 0) return(list(mean_z=rep(NA,N), min_z=rep(NA,N),
    prop_low=rep(NA,N), combined=rep(NA,N), nF_detected=nF_est, n_factors_used=0))

  valid_cols <- which(colSums(!is.na(factor_z)) > 0)
  fz <- factor_z[, valid_cols, drop=FALSE]

  mean_z <- rowMeans(fz, na.rm=TRUE)
  min_z <- apply(fz, 1, min, na.rm=TRUE)
  prop_low <- rowMeans(fz < 1.0, na.rm=TRUE)
  combined <- mean_z * 0.7 + min_z * 0.3

  list(mean_z=mean_z, min_z=min_z, prop_low=prop_low, combined=combined,
       nF_detected=nF_est, n_factors_used=factors_used)
}

# ============================================================================
# SIMULATION
# ============================================================================

cat("\n============================================================\n")
cat("PER-FACTOR COHERENCE vs STANDARD COUPLED\n")
cat("============================================================\n\n")

nF_levels <- c(4,6,8,10,15,20,25,30)
ipf_levels <- c(3,6,10)
pct_levels <- c(0.05, 0.10, 0.25)
n_reps <- 5; N_val <- 300; iterations <- 50

grid <- expand.grid(nF=nF_levels, ipf=ipf_levels, pct=pct_levels, rep=1:n_reps)
grid$total_items <- grid$nF * grid$ipf
total_cond <- nrow(grid)
cat(sprintf("Total conditions: %d\n\n", total_cond))

results <- data.frame()
t_start <- Sys.time()

for(i in seq_len(total_cond)){
  nF<-grid$nF[i]; ipf<-grid$ipf[i]; pct<-grid$pct[i]; rep_i<-grid$rep[i]
  total_items<-grid$total_items[i]; seed_val<-70000+i

  set.seed(seed_val)
  cat(sprintf("[%d/%d] nF=%d ipf=%d (%d it) pct=%.0f%% rep=%d ... ",
              i,total_cond,nF,ipf,total_items,pct*100,rep_i))

  good_data<-tryCatch(simulated_good_responses(nF,rep(ipf,nF),N_val),error=function(e){cat("SKIP\n");NULL})
  if(is.null(good_data)) next
  inj<-inject_careless(good_data,pct,seed=seed_val+10000)
  corrupted<-inj$data_corrupted; labels_vec<-as.integer(inj$labels$careless)

  # Standard coupled
  r_std<-tryCatch(ReReReRe(corrupted,corProp=0.03,iterations=iterations,align_signs=TRUE,mode="coupled"),
                   error=function(e) NULL)
  b_std<-if(!is.null(r_std)) find_best_z(r_std$z_score,labels_vec) else list(z=NA,mcc=NA)

  # Per-factor
  pf<-tryCatch(ReReReRe_perFactor(corrupted,iterations=iterations), error=function(e) NULL)
  if(is.null(pf)){
    cat("PF_FAIL\n"); next
  }
  b_mean<-find_best_z(pf$mean_z,labels_vec,seq(0.1,3,.1))
  b_min<-find_best_z(pf$min_z,labels_vec,seq(-2,3,.1))
  b_comb<-find_best_z(pf$combined,labels_vec,seq(-1,3,.1))
  b_prop<-find_best_z(pf$prop_low,labels_vec,seq(0.1,1,.05),lower_is_bad=FALSE)

  cat(sprintf("std=%.3f mean=%.3f min=%.3f comb=%.3f prop=%.3f\n",
      ifelse(is.na(b_std$mcc),0,b_std$mcc), b_mean$mcc, b_min$mcc, b_comb$mcc, b_prop$mcc))

  results<-rbind(results, data.frame(
    nF=nF, ipf=ipf, total_items=total_items, pct=pct, rep=rep_i,
    mcc_std=b_std$mcc, mcc_mean=b_mean$mcc, mcc_min=b_min$mcc,
    mcc_comb=b_comb$mcc, mcc_prop=b_prop$mcc,
    nF_detected=pf$nF_detected, n_factors_used=pf$n_factors_used))

  if(i%%50==0){
    elapsed<-as.numeric(difftime(Sys.time(),t_start,units="mins"))
    cat(sprintf("  >> %d/%d, %.1f min, ETA %.1f min\n",i,total_cond,elapsed,elapsed/i*(total_cond-i)))
    write.csv(results,"test_perfactor_v2_checkpoint.csv",row.names=FALSE)
  }
}

write.csv(results, "test_perfactor_v2_results.csv", row.names=FALSE)
elapsed_total<-as.numeric(difftime(Sys.time(),t_start,units="mins"))

# ============================================================================
# REPORT
# ============================================================================

results$item_bin <- cut(results$total_items, breaks=c(0,30,60,100,200,Inf),
                        labels=c("<30","30-60","60-100","100-200",">200"))

report_dir <- "archive/perfactor_v2"
dir.create(report_dir, recursive=TRUE, showWarnings=FALSE)

res_long <- results %>%
  select(nF,ipf,total_items,pct,rep,item_bin,mcc_std,mcc_mean,mcc_min,mcc_comb,mcc_prop) %>%
  pivot_longer(cols=starts_with("mcc_"), names_to="method", values_to="mcc") %>%
  mutate(method=recode(method, mcc_std="Standard (coupled)",
    mcc_mean="PerFactor: mean(z)", mcc_min="PerFactor: min(z)",
    mcc_comb="PerFactor: combined", mcc_prop="PerFactor: prop_low"))

colors5 <- c("Standard (coupled)"="#E41A1C", "PerFactor: mean(z)"="#377EB8",
             "PerFactor: min(z)"="#4DAF4A", "PerFactor: combined"="#984EA3",
             "PerFactor: prop_low"="#FF7F00")

# Plot 1: Overall
p1<-res_long%>%group_by(method)%>%
  summarise(m=mean(mcc,na.rm=T),se=sd(mcc,na.rm=T)/sqrt(n()),.groups="drop")%>%
  ggplot(aes(x=reorder(method,-m),y=m,fill=method))+
  geom_col(width=.6)+geom_errorbar(aes(ymin=m-se,ymax=m+se),width=.2)+
  geom_text(aes(label=sprintf("%.3f",m)),vjust=-.5,size=3.5)+
  labs(title="Overall Mean MCC: Per-Factor vs Standard Coupled",x="",y="Mean MCC (oracle)")+
  theme_minimal(base_size=12)+theme(legend.position="none",axis.text.x=element_text(angle=20,hjust=1))+
  scale_fill_manual(values=colors5)+ylim(0,NA)
ggsave(file.path(report_dir,"01_overall.png"),p1,width=10,height=5.5,dpi=150)

# Plot 2: By item bin
p2<-res_long%>%group_by(item_bin,method)%>%summarise(m=mean(mcc,na.rm=T),.groups="drop")%>%
  ggplot(aes(x=item_bin,y=m,fill=method))+
  geom_col(position=position_dodge(.8),width=.7)+
  labs(title="MCC by Item Range",x="Total items",y="Mean MCC")+
  theme_minimal(base_size=12)+scale_fill_manual(values=colors5)
ggsave(file.path(report_dir,"02_by_item_bin.png"),p2,width=11,height=5.5,dpi=150)

# Plot 3: By nF
p3<-res_long%>%group_by(nF,method)%>%summarise(m=mean(mcc,na.rm=T),.groups="drop")%>%
  ggplot(aes(x=factor(nF),y=m,color=method,group=method))+
  geom_line(linewidth=1)+geom_point(size=2.5)+
  labs(title="MCC by nFactors",x="nFactors",y="Mean MCC")+
  theme_minimal(base_size=12)+scale_color_manual(values=colors5)
ggsave(file.path(report_dir,"03_by_nF.png"),p3,width=10,height=5.5,dpi=150)

# Plot 4: Scatter + LOESS
p4<-res_long%>%ggplot(aes(x=total_items,y=mcc,color=method))+
  geom_point(alpha=.1,size=.8)+geom_smooth(method="loess",span=.5,se=TRUE,linewidth=1)+
  labs(title="MCC by Total Items (LOESS)",x="Total items",y="MCC")+
  theme_minimal(base_size=12)+scale_color_manual(values=colors5)
ggsave(file.path(report_dir,"04_by_total_items.png"),p4,width=10,height=5.5,dpi=150)

# Plot 5: Heatmap gain best-PF vs standard
heat<-results%>%group_by(nF,ipf)%>%
  summarise(std=mean(mcc_std,na.rm=T),
            best_pf=mean(pmax(mcc_mean,mcc_min,mcc_comb,mcc_prop,na.rm=TRUE),na.rm=T),
            .groups="drop")%>%mutate(gain=best_pf-std)
p5<-ggplot(heat,aes(x=factor(ipf),y=factor(nF),fill=gain))+geom_tile()+
  geom_text(aes(label=sprintf("%+.3f",gain)),size=3.5)+
  scale_fill_gradient2(low="#d73027",mid="white",high="#1a9850",midpoint=0)+
  labs(title="Best Per-Factor Gain over Standard Coupled",x="Items/factor",y="nFactors")+
  theme_minimal(base_size=12)
ggsave(file.path(report_dir,"05_gain_heatmap.png"),p5,width=7,height=6,dpi=150)

# Plot 6: Winner map
wm<-results%>%group_by(nF,ipf)%>%
  summarise(std=mean(mcc_std,na.rm=T),mean_z=mean(mcc_mean,na.rm=T),
            min_z=mean(mcc_min,na.rm=T),comb=mean(mcc_comb,na.rm=T),
            prop=mean(mcc_prop,na.rm=T),.groups="drop")%>%
  rowwise()%>%
  mutate(winner=c("Standard","mean(z)","min(z)","combined","prop_low")[
    which.max(c(std,mean_z,min_z,comb,prop))],
    margin=max(c(std,mean_z,min_z,comb,prop))-sort(c(std,mean_z,min_z,comb,prop))[4])%>%
  ungroup()
p6<-ggplot(wm,aes(x=factor(ipf),y=factor(nF),fill=winner))+geom_tile()+
  geom_text(aes(label=sprintf("%s\n(+%.3f)",winner,margin)),size=2.5)+
  labs(title="Winner per Cell",x="Items/factor",y="nFactors")+
  theme_minimal(base_size=12)
ggsave(file.path(report_dir,"06_winner_map.png"),p6,width=8,height=6,dpi=150)

# Text report
sink(file.path(report_dir,"report.txt"))
cat(sprintf("PER-FACTOR vs STANDARD COUPLED | %s | %d cond | %.1f min\n\n",
            Sys.Date(),total_cond,elapsed_total))

cat("--- OVERALL ---\n")
overall<-res_long%>%group_by(method)%>%summarise(m=mean(mcc,na.rm=T),sd=sd(mcc,na.rm=T),.groups="drop")%>%arrange(desc(m))
print(as.data.frame(overall),digits=3)

cat("\n\n--- BY ITEM BIN ---\n")
bt<-res_long%>%group_by(item_bin,method)%>%summarise(m=mean(mcc,na.rm=T),.groups="drop")%>%
  pivot_wider(names_from=method,values_from=m)
print(as.data.frame(bt),digits=3)

cat("\n\n--- BY nF ---\n")
nt<-res_long%>%group_by(nF,method)%>%summarise(m=mean(mcc,na.rm=T),.groups="drop")%>%
  pivot_wider(names_from=method,values_from=m)
print(as.data.frame(nt),digits=3)

cat("\n\n--- WINNER MAP ---\n")
print(as.data.frame(wm%>%select(nF,ipf,std,mean_z,min_z,comb,prop,winner,margin)),digits=3)
sink()

# Console
cat("\n\n============================================================\n")
cat("RESULTS\n============================================================\n\n")
print(as.data.frame(overall),digits=3)
cat("\n--- BY ITEM BIN ---\n")
print(as.data.frame(bt),digits=3)
cat(sprintf("\nReport: %s/\n",report_dir))
cat("\n=== DONE ===\n")
