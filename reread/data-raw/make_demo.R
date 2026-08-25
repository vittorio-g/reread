# Generates the bundled demo dataset (simulated; no licensing constraints).
library(reread)
set.seed(20260716)
clean <- simulate_clean(n_factors = 15, items_per_factor = 6, n = 400, seed = 20260716)
inj   <- inject_careless(clean$data, prevalence = 0.2, seed = 99)
demo_careless <- list(
  responses = inj$data,             # 400 x 90 integer matrix of 1-5 Likert responses
  careless  = inj$labels,           # 0/1 ground-truth label
  severity  = inj$severity,         # corrupted fraction (0 = clean)
  pattern   = inj$pattern)          # careless pattern per respondent (NA = clean)
save(demo_careless, file = "data/demo_careless.rda", compress = "xz")
cat("wrote data/demo_careless.rda:", nrow(demo_careless$responses), "x",
    ncol(demo_careless$responses), " careless =", sum(demo_careless$careless), "\n")
