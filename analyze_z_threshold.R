#### Analyze optimal z threshold from final multiverse results ####
library(data.table)

cat("Loading final multiverse results...\n")
dt <- fread("final_multiverse_results.csv")
cat(sprintf("Loaded: %d rows\n", nrow(dt)))
cat(sprintf("Columns: %s\n", paste(names(dt), collapse=", ")))

# Check z_threshold values
cat(sprintf("\nz_threshold levels: %s\n", paste(sort(unique(dt$z_threshold)), collapse=", ")))
cat(sprintf("corProp levels: %s\n", paste(sort(unique(dt$corProp)), collapse=", ")))
cat(sprintf("nFactors levels: %s\n", paste(sort(unique(dt$nFactors)), collapse=", ")))
cat(sprintf("n_respondents levels: %s\n", paste(sort(unique(dt$n_respondents)), collapse=", ")))

# === 1. Overall best z_threshold ===
cat("\n=== 1. OVERALL: Mean MCC by z_threshold ===\n")
by_z <- dt[, .(mean_mcc = mean(mcc, na.rm=TRUE),
               median_mcc = median(mcc, na.rm=TRUE),
               mean_auc = mean(auc, na.rm=TRUE),
               mean_sens = mean(sensitivity, na.rm=TRUE),
               mean_spec = mean(specificity, na.rm=TRUE)),
           by = z_threshold][order(z_threshold)]
print(by_z)

cat(sprintf("\nBest overall z: %.1f (MCC=%.3f)\n",
    by_z$z_threshold[which.max(by_z$mean_mcc)],
    max(by_z$mean_mcc)))

# === 2. Best z_threshold by nFactors ===
cat("\n=== 2. Best z_threshold by nFactors ===\n")
by_nf_z <- dt[, .(mean_mcc = mean(mcc, na.rm=TRUE)),
              by = .(nFactors, z_threshold)]
best_z_by_nf <- by_nf_z[, .SD[which.max(mean_mcc)], by = nFactors][order(nFactors)]
cat(sprintf("%-6s %-8s %-10s\n", "nF", "Best_z", "MCC"))
for (i in 1:nrow(best_z_by_nf)) {
  cat(sprintf("%-6d %-8.1f %-10.3f\n",
      best_z_by_nf$nFactors[i], best_z_by_nf$z_threshold[i], best_z_by_nf$mean_mcc[i]))
}

# === 3. Best z for "recommended" zone (nF>=10, n>=100) ===
cat("\n=== 3. Recommended zone (nF>=10, n>=100) ===\n")
rec <- dt[nFactors >= 10 & n_respondents >= 100]
by_z_rec <- rec[, .(mean_mcc = mean(mcc, na.rm=TRUE),
                     mean_sens = mean(sensitivity, na.rm=TRUE),
                     mean_spec = mean(specificity, na.rm=TRUE)),
                by = z_threshold][order(z_threshold)]
print(by_z_rec)
cat(sprintf("\nBest recommended z: %.1f (MCC=%.3f)\n",
    by_z_rec$z_threshold[which.max(by_z_rec$mean_mcc)],
    max(by_z_rec$mean_mcc)))

# === 4. Best z for "sweet spot" (nF>=15, n>=300) ===
cat("\n=== 4. Sweet spot (nF>=15, n>=300) ===\n")
sweet <- dt[nFactors >= 15 & n_respondents >= 300]
by_z_sweet <- sweet[, .(mean_mcc = mean(mcc, na.rm=TRUE),
                         mean_sens = mean(sensitivity, na.rm=TRUE),
                         mean_spec = mean(specificity, na.rm=TRUE)),
                    by = z_threshold][order(z_threshold)]
print(by_z_sweet)
cat(sprintf("\nBest sweet-spot z: %.1f (MCC=%.3f)\n",
    by_z_sweet$z_threshold[which.max(by_z_sweet$mean_mcc)],
    max(by_z_sweet$mean_mcc)))

# === 5. Fixed corProp=0.05, vary z only ===
cat("\n=== 5. Fixed corProp=0.05: MCC by z_threshold ===\n")
fixed_cp <- dt[corProp == 0.05]
by_z_fixed <- fixed_cp[, .(mean_mcc = mean(mcc, na.rm=TRUE),
                            mean_sens = mean(sensitivity, na.rm=TRUE),
                            mean_spec = mean(specificity, na.rm=TRUE)),
                       by = z_threshold][order(z_threshold)]
print(by_z_fixed)
cat(sprintf("\nBest z (corProp=0.05): %.1f (MCC=%.3f)\n",
    by_z_fixed$z_threshold[which.max(by_z_fixed$mean_mcc)],
    max(by_z_fixed$mean_mcc)))

# === 6. Best z by nFactors, fixed corProp=0.05 ===
cat("\n=== 6. Best z by nFactors (corProp=0.05) ===\n")
by_nf_z_fixed <- fixed_cp[, .(mean_mcc = mean(mcc, na.rm=TRUE)),
                           by = .(nFactors, z_threshold)]
best_z_nf_fixed <- by_nf_z_fixed[, .SD[which.max(mean_mcc)], by = nFactors][order(nFactors)]
cat(sprintf("%-6s %-8s %-10s\n", "nF", "Best_z", "MCC"))
for (i in 1:nrow(best_z_nf_fixed)) {
  cat(sprintf("%-6d %-8.1f %-10.3f\n",
      best_z_nf_fixed$nFactors[i], best_z_nf_fixed$z_threshold[i], best_z_nf_fixed$mean_mcc[i]))
}

# === 7. Comparison: z=1.0 vs z=1.5 vs z=2.0 (corProp=0.05) ===
cat("\n=== 7. Head-to-head: z=1.0 vs z=1.5 vs z=2.0 (corProp=0.05) ===\n")
for (z in c(1.0, 1.5, 2.0)) {
  sub <- fixed_cp[z_threshold == z]
  all_mcc <- mean(sub$mcc, na.rm=TRUE)
  rec_mcc <- mean(sub[nFactors >= 10 & n_respondents >= 100]$mcc, na.rm=TRUE)
  sweet_mcc <- mean(sub[nFactors >= 15 & n_respondents >= 300]$mcc, na.rm=TRUE)
  cat(sprintf("z=%.1f: All=%.3f  Recommended=%.3f  SweetSpot=%.3f\n",
      z, all_mcc, rec_mcc, sweet_mcc))
}

# === 8. MCC heatmap data: nFactors x z_threshold (corProp=0.05) ===
cat("\n=== 8. MCC heatmap: nFactors x z_threshold (corProp=0.05) ===\n")
heatmap <- dcast(by_nf_z_fixed, nFactors ~ z_threshold, value.var = "mean_mcc")
print(heatmap, digits=3)

cat("\nDone.\n")
