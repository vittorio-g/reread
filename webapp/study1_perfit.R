## study1_perfit.R — person-fit (l_z via mirt GRM) head-to-head vs rc on Study 1 real data.
suppressMessages({library(mirt); library(pROC)})
setwd("C:/Users/vitto/Desktop/ReReReRe/webapp")
M   <- as.matrix(read.csv("_study_matrix.csv", check.names=FALSE))
lab <- read.csv("_study_labels.csv")$careless
cat("matrix:", nrow(M), "x", ncol(M), " labels:", length(lab), " careless:", sum(lab), "\n")

## Drop degenerate (constant) items
nun <- apply(M, 2, function(x) length(unique(x[is.finite(x)])))
keep <- which(nun >= 2)
cat("dropped", ncol(M)-length(keep), "constant/degenerate items; kept", length(keep), "\n")
Mk <- M[, keep, drop=FALSE]

## Fit unidimensional graded response model (handles per-item mixed #categories)
mod <- mirt(as.data.frame(Mk), 1, itemtype="graded", verbose=FALSE,
            technical=list(NCYCLES=2000))
cat("converged:", extract.mirt(mod,"converged"), " (", extract.mirt(mod,"iterations"), "cycles )\n")

pf <- personfit(mod)
cat("personfit columns:", paste(colnames(pf), collapse=", "), "\n")

Zh <- pf$Zh
cat(sprintf("mean Zh careful=%.3f careless=%.3f (careless should be LOWER)\n",
            mean(Zh[lab==0], na.rm=TRUE), mean(Zh[lab==1], na.rm=TRUE)))

aucv <- function(lab, s){ ok <- is.finite(s)
  as.numeric(auc(roc(lab[ok], s[ok], direction="<", quiet=TRUE))) }
## predictor for careless = -Zh (low l_z = misfit = careless)
auc_lz <- aucv(lab, -Zh)
cat(sprintf("\n=== Study 1: l_z (Zh) AUC vs careless = %.3f ===\n", auc_lz))

if("Zh_theta_est" %in% colnames(pf) || "lzstar" %in% tolower(colnames(pf))){
  # not standard; check below
}
## lzstar if present
lzs <- NULL
for(cn in colnames(pf)) if(grepl("star", cn, ignore.case=TRUE) || cn=="Zh_star") lzs <- pf[[cn]]
if(!is.null(lzs)){
  cat(sprintf("=== Study 1: l_z* AUC = %.3f ===\n", aucv(lab, -lzs)))
} else cat("l_z* not returned by personfit (only Zh available)\n")

res <- data.frame(index=c("person_fit_lz","rc_paper","ensemble_paper"),
                  auc=c(round(auc_lz,3), 0.898, 0.984))
write.csv(res, "study1_perfit_auc.csv", row.names=FALSE)
print(res)
