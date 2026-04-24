suppressPackageStartupMessages(library(dplyr))
df <- read.csv("sim_ensemble_scores.csv")
df$gt <- (df$pattern != "clean") & (df$corruption > 0.50)
set.seed(42)
folds <- sample(1:5, nrow(df), replace = TRUE)
pred <- numeric(nrow(df))
for (f in 1:5) {
  tr <- folds != f; te <- folds == f
  m <- glm(gt ~ z_rr + irv + longstring + d2 + person_tot,
           data = df[tr, ], family = binomial)
  pred[te] <- predict(m, newdata = df[te, ], type = "response")
}
df$pred_ens <- pred

summ <- df %>%
  filter(pattern != "clean") %>%
  mutate(bin = round(corruption * 10) / 10) %>%
  group_by(pattern, bin) %>%
  summarise(n = n(),
            det = round(mean(pred_ens > 0.5), 3),
            .groups = "drop") %>%
  arrange(bin, pattern)

# Print only high-corruption buckets
cat("=== Detection by pattern x corruption (ensemble) ===\n")
print(as.data.frame(summ %>% filter(bin >= 0.8)), row.names = FALSE)
cat("\n=== All buckets ===\n")
print(as.data.frame(summ), row.names = FALSE)
