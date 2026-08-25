## core_only_headtohead.R — head-to-head restricted to the CORE datasets of Table 4
## (boundary/out-of-envelope sets removed from the paper), mean AUC + paired Wilcoxon.
D <- read.csv("ext_bench2/benchmark_ext2.csv", stringsAsFactors=FALSE)
C <- D[D$kind=="core",]
cat("core datasets:", nrow(C), "->", paste(C$dataset, collapse=", "), "\n\n")
M <- c("ReReRe","rr","PT","Mah","IRV","Long","Psyn")
cat("=== mean AUC (core only) ===\n")
for(m in M) cat(sprintf("  %-8s %.3f\n", m, mean(C[[m]], na.rm=TRUE)))
cat("\n=== paired Wilcoxon: ensemble vs each competitor (core only) ===\n")
for(m in setdiff(M,c("ReReRe","rr"))){
  ok <- is.finite(C$ReReRe) & is.finite(C[[m]])
  w <- suppressWarnings(wilcox.test(C$ReReRe[ok], C[[m]][ok], paired=TRUE))
  cat(sprintf("  ens vs %-6s wins=%d/%d  p=%.4f\n", m, sum(C$ReReRe[ok]>C[[m]][ok]), sum(ok), w$p.value))
}
cat("\n=== paired Wilcoxon: rr vs each competitor (core only) ===\n")
for(m in setdiff(M,c("ReReRe","rr"))){
  ok <- is.finite(C$rr) & is.finite(C[[m]])
  w <- suppressWarnings(wilcox.test(C$rr[ok], C[[m]][ok], paired=TRUE))
  cat(sprintf("  rr  vs %-6s wins=%d/%d  p=%.4f\n", m, sum(C$rr[ok]>C[[m]][ok]), sum(ok), w$p.value))
}
cat("\n=== ensemble vs rc: how often does the gate keep ens >= rc? (core) ===\n")
cat(sprintf("  ens >= rr on %d/%d datasets\n", sum(C$ReReRe >= C$rr, na.rm=TRUE), sum(is.finite(C$ReReRe))))
