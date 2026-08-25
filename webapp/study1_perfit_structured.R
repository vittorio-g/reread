## study1_perfit_structured.R — structure-based person-fit (l_z) on Study 1, giving l_z its
## fairest shot by using the KNOWN item->subscale map (same grouping even-odd/RPR use: the
## `domain` column of _study_items.csv, Random decoy excluded). Fit a unidimensional GRM
## WITHIN each subscale, get per-subscale Zh, combine as the item-count-weighted mean of the
## per-subscale Zh. Careless = worse fit = lower combined Zh -> predictor = -combined.
suppressMessages({library(mirt); library(pROC)})
setwd("C:/Users/vitto/Desktop/ReReReRe/webapp")
M   <- as.matrix(read.csv("_study_matrix.csv", check.names=FALSE))
lab <- read.csv("_study_labels.csv")$careless
items <- read.csv("_study_items.csv", stringsAsFactors=FALSE, na.strings="")  # keep literal "NA" domain (PANAS Negative Affect)
stopifnot(all(colnames(M) == items$id))          # column order matches the map
cat("matrix:", nrow(M), "x", ncol(M), " careless:", sum(lab), "\n")

## item -> subscale map = domain column; exclude the Random decoy (as even-odd/RPR do)
items$sub <- items$domain
use_items <- which(items$sub != "Random")
subs <- sort(unique(items$sub[use_items]))
cat("subscales (Random excluded):", length(subs), "->", paste(subs, collapse=", "), "\n")

n <- nrow(M)
ZH <- matrix(NA_real_, n, length(subs)); colnames(ZH) <- subs
wts <- numeric(length(subs)); dropped <- character(0)
for(si in seq_along(subs)){
  cols <- which(items$sub == subs[si])
  X <- M[, cols, drop=FALSE]
  nun <- apply(X, 2, function(x) length(unique(x[is.finite(x)])))
  X <- X[, nun >= 2, drop=FALSE]               # drop constant items within subscale
  J <- ncol(X)
  if(J < 4){ dropped <- c(dropped, sprintf("%s(J=%d,<4)", subs[si], J)); next }
  mod <- tryCatch(mirt(as.data.frame(X), 1, itemtype="graded", verbose=FALSE,
                       technical=list(NCYCLES=3000)),
                  error=function(e){cat("fit ERR", subs[si], conditionMessage(e), "\n"); NULL})
  if(is.null(mod) || !extract.mirt(mod, "converged")){
    dropped <- c(dropped, sprintf("%s(no-converge)", subs[si])); next }
  ZH[, si] <- personfit(mod)$Zh
  wts[si]  <- J
  cat(sprintf("  %-18s J=%2d converged Zh mean careful=%+.3f careless=%+.3f\n",
      subs[si], J, mean(ZH[lab==0,si]), mean(ZH[lab==1,si])))
}
keep <- which(wts > 0)
cat("used", length(keep), "subscales,", sum(wts), "items; dropped:",
    if(length(dropped)) paste(dropped, collapse=", ") else "none", "\n")

## combine: item-count-weighted mean of per-subscale Zh
W <- wts[keep]
comb <- as.numeric(ZH[, keep, drop=FALSE] %*% W) / sum(W)
cat(sprintf("combined Zh mean careful=%+.3f careless=%+.3f (careless should be LOWER)\n",
    mean(comb[lab==0]), mean(comb[lab==1])))

aucv <- function(lab, s){ ok <- is.finite(s)
  as.numeric(auc(roc(lab[ok], s[ok], direction="<", quiet=TRUE))) }
auc_struct <- aucv(lab, -comb)                    # careless = higher predictor = -Zh
cat(sprintf("\n=== Study 1: STRUCTURE-BASED l_z AUC vs careless = %.3f ===\n", auc_struct))

## append to the Study 1 AUC CSV
csv <- "study1_perfit_auc.csv"
res <- read.csv(csv, stringsAsFactors=FALSE)
res <- res[res$index != "person_fit_lz_structured", ]     # idempotent
res <- rbind(res, data.frame(index="person_fit_lz_structured", auc=round(auc_struct,3)))
write.csv(res, csv, row.names=FALSE)
print(res)
