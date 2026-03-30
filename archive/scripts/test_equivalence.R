## Test: sign flip vs reverse coding equivalence for uniform Likert scale
set.seed(42)

# Simulate one respondent, k=20 coupled pairs, all 1-7 scale
k <- 20
A <- sample(1:7, k, replace=TRUE)
B <- sample(1:7, k, replace=TRUE)
neg <- sample(k, 8)  # 8 "negative correlation" pairs

# Sign flip: B * -1
B_flip <- B; B_flip[neg] <- -B[neg]

# Reverse code: (7+1) - B = 8 - B
B_rev <- B; B_rev[neg] <- 8 - B[neg]

cat("=== Uniform scale (all 1-7) ===\n")
cat(sprintf("|cor(A, B_flip)| = %.6f\n", abs(cor(A, B_flip))))
cat(sprintf("|cor(A, B_rev)|  = %.6f\n", abs(cor(A, B_rev))))
cat(sprintf("Difference: %.8f\n\n", abs(cor(A, B_flip)) - abs(cor(A, B_rev))))

# Mathematical argument:
# For uniform scale, B_rev = (max+1) - B = -B + C (constant)
# In the k-vector: B_rev[j] = -B[j] + C for negative pairs, B[j] for positive
# B_flip[j] = -B[j] for negative pairs, B[j] for positive
# B_rev = B_flip + C_vector where C_vector[j] = C for neg pairs, 0 for pos pairs
# cor(A, B_rev) != cor(A, B_flip) because adding a non-constant vector changes correlation

# But for UNIFORM scale where all items have the same max:
# The constant C = max+1 is the same for all negative pairs
# Still, it's added only to some elements, making it a non-uniform shift

# Run 1000 simulations to check the magnitude of the difference
diffs <- replicate(1000, {
  A_sim <- sample(1:7, k, replace=TRUE)
  B_sim <- sample(1:7, k, replace=TRUE)
  neg_sim <- sample(k, 8)
  B_f <- B_sim; B_f[neg_sim] <- -B_sim[neg_sim]
  B_r <- B_sim; B_r[neg_sim] <- 8 - B_sim[neg_sim]
  abs(cor(A_sim, B_f)) - abs(cor(A_sim, B_r))
})

cat("=== 1000 simulations (uniform 1-7 scale) ===\n")
cat(sprintf("Mean |cor| difference: %.6f\n", mean(diffs)))
cat(sprintf("SD: %.6f\n", sd(diffs)))
cat(sprintf("Max absolute diff: %.6f\n", max(abs(diffs))))
cat(sprintf("Mean absolute diff: %.6f\n", mean(abs(diffs))))

cat("\nConclusion: ")
if (mean(abs(diffs)) < 0.01) {
  cat("Differences are negligible for uniform scales.\n")
  cat("Existing simulation results (with sign flip) are valid.\n")
} else {
  cat("Differences are non-negligible. Simulations should be re-run.\n")
}
