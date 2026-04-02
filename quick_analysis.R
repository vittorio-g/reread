res <- read.csv("sim_realistic_results.csv")
cat("Z values:", sort(unique(res$z_threshold)), "\n")
cat("Conditions:", nrow(res)/length(unique(res$z_threshold)), "\n\n")

# z=1.5 is in the range (0.1, 0.3, ..., 1.5, 1.7, ...)
# seq(0.1, 3.0, 0.2) = 0.1 0.3 0.5 0.7 0.9 1.1 1.3 1.5 1.7 ...
z15 <- res[abs(res$z_threshold - 1.5) < 0.05, ]
cat("Rows at z=1.5:", nrow(z15), "\n")

if (nrow(z15) > 0) {
  cat(sprintf("\n--- AT z=1.5 ---\n"))
  cat(sprintf("  Std: MCC=%.3f, sens=%.3f, spec=%.3f\n",
      mean(z15$mcc_std, na.rm=T), mean(z15$sens_std, na.rm=T), mean(z15$spec_std, na.rm=T)))
  cat(sprintf("  F:   MCC=%.3f, sens=%.3f, spec=%.3f\n",
      mean(z15$mcc_f, na.rm=T), mean(z15$sens_f, na.rm=T), mean(z15$spec_f, na.rm=T)))
  cat(sprintf("  Noise flagged: Std=%.1f%%, F=%.1f%%\n",
      100*mean(z15$noise_flagged_std, na.rm=T), 100*mean(z15$noise_flagged_f, na.rm=T)))
}

# Try other z values
for (zt in c(0.9, 1.1, 1.3, 1.7)) {
  zs <- res[abs(res$z_threshold - zt) < 0.05, ]
  if (nrow(zs) > 0) {
    cat(sprintf("\n--- AT z=%.1f ---\n", zt))
    cat(sprintf("  Std: MCC=%.3f  F: MCC=%.3f  Noise: Std=%.1f%% F=%.1f%%\n",
        mean(zs$mcc_std, na.rm=T), mean(zs$mcc_f, na.rm=T),
        100*mean(zs$noise_flagged_std, na.rm=T), 100*mean(zs$noise_flagged_f, na.rm=T)))
  }
}

# By item bin at z=1.3
z13 <- res[abs(res$z_threshold - 1.3) < 0.05, ]
if (nrow(z13) > 0) {
  z13$item_bin <- cut(z13$total_items, breaks=c(0,30,60,100,200,Inf),
                      labels=c("<30","30-60","60-100","100-200",">200"))
  cat("\n--- BY ITEM BIN at z=1.3 ---\n")
  for (b in levels(z13$item_bin)) {
    s <- z13[z13$item_bin==b, ]
    if (nrow(s)==0) next
    cat(sprintf("  %-10s Std=%.3f F=%.3f  Noise_Std=%.1f%% Noise_F=%.1f%%\n",
        b, mean(s$mcc_std, na.rm=T), mean(s$mcc_f, na.rm=T),
        100*mean(s$noise_flagged_std, na.rm=T), 100*mean(s$noise_flagged_f, na.rm=T)))
  }
}
