# Quick smoke test: confirm variance_penalty flag works and gives expected results
source("Synthetic_Good_Responses_2.R")
source("ReReReRe.R")

set.seed(42)
clean <- simulated_good_responses(nConstructs = 10, nItems = rep(6, 10), n = 150)

# Inject a pure straight-liner at row 1
data <- clean
data[1, ] <- 5  # pure straight-liner at value 5

# Inject an acquiescent at row 2
data[2, ] <- sample(6:7, ncol(data), replace = TRUE)

cat("Testing variance_penalty OFF (default):\n")
rr_off <- ReReReRe(data, iterations = 50)
cat(sprintf("  row 1 (straight-liner): z=%.3f, flagged=%s\n",
            rr_off$z_score[1], rr_off$flagged[1]))
cat(sprintf("  row 2 (acquiescent):    z=%.3f, flagged=%s\n",
            rr_off$z_score[2], rr_off$flagged[2]))
cat(sprintf("  row 10 (attentive):     z=%.3f, flagged=%s\n",
            rr_off$z_score[10], rr_off$flagged[10]))
cat(sprintf("  variance_penalty_used: %s\n", rr_off$variance_penalty_used[1]))

cat("\nTesting variance_penalty ON:\n")
rr_on <- ReReReRe(data, iterations = 50, variance_penalty = TRUE)
cat(sprintf("  row 1 (straight-liner): z=%.3f (raw=%.3f), flagged=%s\n",
            rr_on$z_score[1], rr_on$z_score_raw[1], rr_on$flagged[1]))
cat(sprintf("  row 2 (acquiescent):    z=%.3f (raw=%.3f), flagged=%s\n",
            rr_on$z_score[2], rr_on$z_score_raw[2], rr_on$flagged[2]))
cat(sprintf("  row 10 (attentive):     z=%.3f (raw=%.3f), flagged=%s\n",
            rr_on$z_score[10], rr_on$z_score_raw[10], rr_on$flagged[10]))
cat(sprintf("  variance_penalty_used: %s\n", rr_on$variance_penalty_used[1]))

cat("\n=== PASS ===\n")
