suppressPackageStartupMessages(library(dplyr))
df <- read.csv("sim_ensemble_scores.csv")

# What does z_rr alone do on random-100%?
df_r100 <- df %>% filter(pattern == "random", corruption == 1.0)
df_clean <- df %>% filter(pattern == "clean")

cat("=== z_RR score distribution (standardized, negated: high = careless) ===\n")
cat(sprintf("random @ 100%% corruption: mean=%.2f, median=%.2f, sd=%.2f (n=%d)\n",
            mean(df_r100$z_rr), median(df_r100$z_rr), sd(df_r100$z_rr), nrow(df_r100)))
cat(sprintf("clean respondents:         mean=%.2f, median=%.2f, sd=%.2f (n=%d)\n",
            mean(df_clean$z_rr), median(df_clean$z_rr), sd(df_clean$z_rr), nrow(df_clean)))

cat("\n=== IRV distribution ===\n")
cat(sprintf("random @ 100%%: mean=%.2f, sd=%.2f\n",
            mean(df_r100$irv), sd(df_r100$irv)))
cat(sprintf("clean:         mean=%.2f, sd=%.2f\n",
            mean(df_clean$irv), sd(df_clean$irv)))

# Simple z_rr threshold
cat("\n=== z_RR alone threshold analysis ===\n")
for (tau in c(0, 0.5, 1.0, 1.5, 2.0)) {
  det_r100 <- mean(df_r100$z_rr > tau)
  fpr_clean <- mean(df_clean$z_rr > tau)
  cat(sprintf("threshold z_rr > %.1f (std): det on random@100 = %.3f, FPR clean = %.3f\n",
              tau, det_r100, fpr_clean))
}

# Train glm with ONLY z_rr on just random-100% vs clean
cat("\n=== Logistic regression: z_rr alone, random-100% vs clean only ===\n")
df2 <- rbind(df_r100 %>% mutate(y = 1), df_clean %>% mutate(y = 0))
mod <- glm(y ~ z_rr, data = df2, family = binomial)
preds <- predict(mod, newdata = df2, type = "response")
flag <- preds > 0.5
cat(sprintf("detection on random@100 = %.3f (in-sample)\n",
            mean(flag[df2$y == 1])))
cat(sprintf("FPR on clean = %.3f\n", mean(flag[df2$y == 0])))
