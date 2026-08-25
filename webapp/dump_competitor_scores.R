## dump_competitor_scores.R — per-respondent scores of the structure-free competitor
## indices on the Study-2 external datasets, so that every index can be put through the
## SAME label-free cut as the ensemble (see ext_deployable.js).
##
## benchmark_ext2.R computed these only to take an AUC and threw the scores away; a
## deployable comparison needs the scores themselves. Directions are the ones fixed a
## priori there: high = careless.
suppressMessages({library(careless)})
set.seed(1)
OUT <- "ext_bench2"
idx <- read.csv(file.path(OUT, "index.csv"), stringsAsFactors = FALSE)
idx <- idx[idx$kind == "core", ]           # the ten datasets of tab:ext
impute <- function(m) {
  for (j in seq_len(ncol(m))) { v <- m[, j]; mm <- median(v, na.rm = TRUE); if (is.na(mm)) mm <- 0; m[is.na(v), j] <- mm }
  m
}
for (k in seq_len(nrow(idx))) {
  nm <- idx$name[k]
  t0 <- Sys.time()
  M <- as.matrix(read.csv(file.path(OUT, paste0(nm, "_M.csv")), header = FALSE))
  mat <- impute(M)
  ps <- tryCatch(as.numeric(careless::psychsyn(mat, critval = .60)), error = function(e) rep(NA, nrow(mat)))
  if (sum(is.finite(ps)) < 0.5 * nrow(mat))
    ps <- tryCatch(as.numeric(careless::psychsyn(mat, critval = .40)), error = function(e) rep(NA, nrow(mat)))
  D <- data.frame(
    Mahalanobis = tryCatch(as.numeric(careless::mahad(mat, plot = FALSE, flag = FALSE)), error = function(e) rep(NA, nrow(mat))),
    IRV         = -as.numeric(careless::irv(mat)),
    LongString  =  as.numeric(careless::longstring(mat)),
    PsychSyn    = -ps
  )
  write.csv(D, file.path(OUT, paste0(nm, "_comp.csv")), row.names = FALSE)
  cat(sprintf("%-14s n=%-6d J=%-4d  %.1fs\n", nm, nrow(mat), ncol(mat),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  flush.console()
}
cat("done\n")
