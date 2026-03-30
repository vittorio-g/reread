source("Synthetic_Good_Responses_2.R")
source("Careless_machine_2.R")
source("ReReReRe.R")

# Quick timing for nF=2, 10, 20, 40
for (nf in c(2, 10, 20, 40)) {
  t0 <- Sys.time()
  nItems_vec <- rep(6, nf)
  dat <- simulated_good_responses(nConstructs = nf, nItems = nItems_vec, n = 300, seed = 1)
  inj <- inject_careless(dat, pct_careless = 0.10, seed = 1)
  rr <- ReReReRe(inj$data_corrupted, corProp = 0.05, iterations = 100)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("nF=%2d (%3d items): %.1f sec\n", nf, sum(nItems_vec), elapsed))
}
