setwd("C:/Users/vitto/Desktop/ReReReRe")
df <- suppressWarnings(foreign::read.spss("external_datasets/12_Johnson_IPIP300/ipip20993.sav",
       to.data.frame=TRUE, use.value.labels=FALSE))
nm <- names(df)
items <- grep("^[iI][0-9]+$", nm, value=TRUE)
if (length(items) < 50)
  items <- nm[sapply(df, function(x){x<-suppressWarnings(as.numeric(x));v<-x[!is.na(x)];length(v)>0&&all(v>=0&v<=5)&&length(unique(v[v>0]))>=3})]
cat("item cols:", length(items), "| first:", paste(head(items,6),collapse=","), "\n")
X <- df[, items]; for (j in seq_along(X)) X[[j]] <- suppressWarnings(as.numeric(X[[j]]))
X[X==0] <- NA
Xc <- X[complete.cases(X),]; cat("complete cases:", nrow(Xc), "\n")
set.seed(20260706)
NI <- min(100, ncol(Xc)); NR <- 300
Xs <- Xc[sample(nrow(Xc), NR), seq_len(NI)]; colnames(Xs) <- paste0("item_", seq_len(NI))
write.csv(cbind(id=paste0("R", sprintf("%03d", seq_len(NR))), Xs), "webapp/demo_real_base.csv", row.names=FALSE)
cat("WROTE", NR, "x", NI, "| range", min(unlist(Xs),na.rm=T), "-", max(unlist(Xs),na.rm=T), "\n")
