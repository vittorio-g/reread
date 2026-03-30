#### Johnson IPIP-NEO-300: Inject-and-Detect Validation ####
# 20,993 respondents, 300 items, 30 facets (10 items each), ~148 reverse-coded
# No ground truth — we INJECT careless respondents and test detection
# This stress-tests align_signs on heavy reverse coding

setwd("C:/Users/vitto/Desktop/ReReReRe")
library(haven)
source("ReReReRe.R")
source("Careless_machine_2.R")

# ============================================================
# Load and prepare data
# ============================================================
cat("=== Loading Johnson IPIP-NEO-300 ===\n")
d <- read_sav("external_datasets/12_Johnson_IPIP300/ipip20993.sav")
cat(sprintf("Raw: N=%d, cols=%d\n", nrow(d), ncol(d)))

# Extract 300 items (i1 to i300)
item_cols <- paste0("i", 1:300)
items <- as.data.frame(d[, item_cols])
# Force to plain numeric (strip haven labels)
for (col in names(items)) items[[col]] <- as.numeric(items[[col]])

# 0 = missing/unanswered (scale is 1-5). Recode 0 → NA
items[items == 0] <- NA

# Remove rows with any NA
complete <- complete.cases(items)
items <- items[complete, ]
cat(sprintf("Complete cases (no 0s): %d\n", nrow(items)))

# Check scale range
cat(sprintf("Item range: %d-%d\n", min(items, na.rm=TRUE), max(items, na.rm=TRUE)))

# Sample down to manageable size (5000 for speed)
set.seed(42)
if (nrow(items) > 5000) {
  idx <- sample(nrow(items), 5000)
  items <- items[idx, ]
  cat(sprintf("Sampled to: %d respondents\n", nrow(items)))
}

# ============================================================
# MCC helper
# ============================================================
compute_mcc <- function(tp, tn, fp, fn) {
  tp <- as.double(tp); tn <- as.double(tn)
  fp <- as.double(fp); fn <- as.double(fn)
  num <- (tp * tn) - (fp * fn)
  den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  if (is.na(den) || den == 0) return(0)
  num / den
}

# ============================================================
# Inject at multiple rates and test detection
# ============================================================
pct_levels <- c(0.05, 0.10, 0.20)
z_range <- seq(0.5, 3.0, 0.1)
reps <- 3

results <- list()
idx <- 0

for (pct in pct_levels) {
  for (rep in 1:reps) {
    seed <- 42 + rep * 100 + round(pct * 1000)

    cat(sprintf("\n--- pct=%.0f%%, rep=%d ---\n", pct*100, rep))

    # Inject careless
    inj <- inject_careless(items, pct_careless = pct, seed = seed)
    labels <- as.integer(inj$labels$careless)

    cat(sprintf("  Injected: %d careless / %d total\n", sum(labels), length(labels)))

    # Run ReReReRe with defaults (corProp=0.03)
    rr <- ReReReRe(inj$data_corrupted, corProp = 0.03, iterations = 100,
                   align_signs = TRUE, progress = FALSE)

    # Also run auto_z
    rr_auto <- ReReReRe(inj$data_corrupted, corProp = 0.03, iterations = 100,
                        align_signs = TRUE, auto_z = TRUE, progress = FALSE)

    cat(sprintf("  Auto-z chose: z=%.2f\n", rr_auto$z_threshold_used[1]))

    # Evaluate at all z thresholds
    for (z in z_range) {
      flagged <- as.integer(rr$z_score <= z)
      tp <- sum(flagged == 1 & labels == 1)
      tn <- sum(flagged == 0 & labels == 0)
      fp <- sum(flagged == 1 & labels == 0)
      fn <- sum(flagged == 0 & labels == 1)
      mcc <- compute_mcc(tp, tn, fp, fn)

      idx <- idx + 1
      results[[idx]] <- data.frame(
        pct_careless = pct,
        rep_id = rep,
        z_threshold = z,
        mcc = mcc,
        sensitivity = if(tp+fn>0) tp/(tp+fn) else 0,
        specificity = if(tn+fp>0) tn/(tn+fp) else 0,
        stringsAsFactors = FALSE
      )
    }

    # Auto-z evaluation
    auto_flagged <- as.integer(rr_auto$flagged)
    tp_a <- sum(auto_flagged == 1 & labels == 1)
    tn_a <- sum(auto_flagged == 0 & labels == 0)
    fp_a <- sum(auto_flagged == 1 & labels == 0)
    fn_a <- sum(auto_flagged == 0 & labels == 1)
    mcc_a <- compute_mcc(tp_a, tn_a, fp_a, fn_a)

    cat(sprintf("  Auto-z MCC: %.3f (sens=%.3f, spec=%.3f)\n",
                mcc_a,
                if(tp_a+fn_a>0) tp_a/(tp_a+fn_a) else 0,
                if(tn_a+fp_a>0) tn_a/(tn_a+fp_a) else 0))

    # Best oracle z
    all_mcc <- sapply(z_range, function(z) {
      fl <- as.integer(rr$z_score <= z)
      compute_mcc(sum(fl==1 & labels==1), sum(fl==0 & labels==0),
                  sum(fl==1 & labels==0), sum(fl==0 & labels==1))
    })
    best_z <- z_range[which.max(all_mcc)]
    cat(sprintf("  Oracle best: z=%.1f, MCC=%.3f\n", best_z, max(all_mcc)))
  }
}

# ============================================================
# Save and summarize
# ============================================================
raw_df <- do.call(rbind, results)
write.csv(raw_df, "archive/johnson_inject_detect.csv", row.names=FALSE)

cat("\n\n====================================================\n")
cat("   JOHNSON IPIP-NEO-300: INJECT-AND-DETECT RESULTS\n")
cat("   300 items, 30 facets, ~148 reverse-coded\n")
cat("   N=5000 (sampled), corProp=0.03\n")
cat("====================================================\n\n")

# Best MCC per pct
for (pct in pct_levels) {
  sub <- raw_df[raw_df$pct_careless == pct, ]
  agg <- aggregate(mcc ~ z_threshold, data=sub, FUN=mean)
  best_idx <- which.max(agg$mcc)
  cat(sprintf("pct=%.0f%%: best z=%.1f, MCC=%.3f (mean over %d reps)\n",
      pct*100, agg$z_threshold[best_idx], agg$mcc[best_idx], reps))
}

cat("\n=== DONE ===\n")
