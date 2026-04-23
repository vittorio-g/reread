# Plot_Detection_vs_Corruption.R
# Analisi: quanto bene vengono detectati i careless in funzione del grado di corruzione?
#
# Per ogni rispondente salviamo (pattern, corruption_level, z_score).
# Poi graficiamo detection rate e distribuzione degli z per livello di corruzione.

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(tidyr)
})
source("Synthetic_Good_Responses_2.R")
source("ReReReRe.R")
# Extract just the inject_careless_extended function body from Simulation_VariancePenalty.R
# WITHOUT running the main loop (source would execute the whole file).
inject_careless_extended <- local({
  src <- readLines("Simulation_VariancePenalty.R")
  start <- grep("^inject_careless_extended <- function", src)[1]
  end_idx <- start
  brace_depth <- 0; started <- FALSE
  for (i in start:length(src)) {
    line <- src[i]
    opens <- lengths(regmatches(line, gregexpr("\\{", line)))
    closes <- lengths(regmatches(line, gregexpr("\\}", line)))
    brace_depth <- brace_depth + opens - closes
    if (opens > 0) started <- TRUE
    if (started && brace_depth == 0) { end_idx <- i; break }
  }
  eval(parse(text = paste(src[start:end_idx], collapse = "\n")))
  get("inject_careless_extended")
})
cat(sprintf("Loaded inject_careless_extended function\n"))

# ============================================================
# Setup
# ============================================================
NF_LEVELS <- c(8, 16, 30)
IPF_LEVELS <- c(6, 10)
REPS <- 4
N_FIXED <- 500       # bigger sample -> more resp per (pattern, level)
PCT_CARELESS <- 0.50 # 50% careless -> many samples per cell
ITER <- 100
CORPROP <- 0.03
SEED_BASE <- 2026
CARELESS_LEVELS <- seq(0.1, 1.0, 0.1)

# 5 types, 10 levels => each cell gets ~5 respondents per pattern per rep
# Over REPS reps -> 20 resp per (nF,ipf,pattern,level) => enough for smoothing
conditions <- expand.grid(nF = NF_LEVELS, ipf = IPF_LEVELS, stringsAsFactors = FALSE)
conditions$total_items <- conditions$nF * conditions$ipf
TOTAL_CONDITIONS <- nrow(conditions) * REPS

cat(sprintf("Running %d conditions x %d reps = %d RR calls\n",
            nrow(conditions), REPS, TOTAL_CONDITIONS))

# ============================================================
# Main loop: save per-respondent output
# ============================================================
rows <- list(); idx <- 0; t0 <- Sys.time()

for (ci in seq_len(nrow(conditions))) {
  nF <- conditions$nF[ci]; ipf <- conditions$ipf[ci]; tot <- conditions$total_items[ci]
  for (rep_id in seq_len(REPS)) {
    idx <- idx + 1
    set.seed(SEED_BASE + (rep_id - 1) * 31 + ci)

    clean <- simulated_good_responses(nConstructs = nF, nItems = rep(ipf, nF), n = N_FIXED)
    inj <- inject_careless_extended(clean, pct_careless = PCT_CARELESS,
                                    careless_levels = CARELESS_LEVELS)
    rr <- ReReReRe(inj$data_corrupted, corProp = CORPROP, iterations = ITER,
                   align_signs = TRUE, mode = "auto")
    rr_vp <- ReReReRe(inj$data_corrupted, corProp = CORPROP, iterations = ITER,
                      align_signs = TRUE, mode = "auto", variance_penalty = TRUE)

    rows[[length(rows) + 1]] <- data.frame(
      nF = nF, ipf = ipf, total_items = tot, rep = rep_id,
      respondent = seq_len(nrow(inj$data_corrupted)),
      pattern = inj$labels$pattern,
      corruption = round(inj$labels$careless_pct, 2),
      z_standard = rr$z_score,
      z_varpen = rr_vp$z_score
    )

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("[%d/%d] nF=%d ipf=%d rep=%d | %.1f min\n",
                idx, TOTAL_CONDITIONS, nF, ipf, rep_id, elapsed))
  }
}

df <- do.call(rbind, rows)
write.csv(df, "sim_detection_per_resp.csv", row.names = FALSE)

# ============================================================
# Bin corruption into labels
# ============================================================
# "clean" = 0, otherwise the corruption percentage as-is
df <- df %>%
  mutate(corruption_bin = ifelse(pattern == "clean", "clean",
                                  sprintf("%d%%", round(corruption * 100))))

# Bin corruption to clean 10% buckets (fixes rounding artifacts when J*pct is non-integer)
df <- df %>%
  mutate(corruption_bin10 = ifelse(pattern == "clean", NA,
                                   round(corruption * 10) / 10))

# Detection rate per pattern x corruption level (at z_threshold = 1.5)
det_std <- df %>%
  filter(pattern != "clean") %>%
  group_by(pattern, corruption_bin10) %>%
  summarise(n = n(),
            detect_std = mean(z_standard <= 1.5, na.rm = TRUE),
            detect_vp = mean(z_varpen <= 1.5, na.rm = TRUE),
            .groups = "drop") %>%
  rename(corruption = corruption_bin10)

det_clean <- df %>%
  filter(pattern == "clean") %>%
  summarise(detect_std = mean(z_standard <= 1.5, na.rm = TRUE),
            detect_vp = mean(z_varpen <= 1.5, na.rm = TRUE))

cat("\n=== CLEAN false-positive rate at z=1.5 ===\n")
print(det_clean)

cat("\n=== DETECTION RATE per pattern x corruption (at z=1.5) ===\n")
print(as.data.frame(det_std), row.names = FALSE)

# ============================================================
# PLOT 1 — Detection rate vs corruption level
# ============================================================
plot_df <- det_std %>%
  pivot_longer(cols = c(detect_std, detect_vp),
               names_to = "variant", values_to = "detect_rate") %>%
  mutate(variant = recode(variant,
                          "detect_std" = "Standard",
                          "detect_vp" = "Variance-Penalized (α=3, β=0.5)"))

pattern_colors <- c(
  random = "#e41a1c",
  longstring = "#377eb8",
  mixed = "#4daf4a",
  pure_straight = "#984ea3",
  acquiescent = "#ff7f00"
)

p1 <- ggplot(plot_df, aes(x = corruption, y = detect_rate,
                          color = pattern, group = pattern)) +
  geom_hline(yintercept = det_clean$detect_std, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_line(size = 1.1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = pattern_colors, name = "Pattern") +
  scale_x_continuous(breaks = seq(0.1, 1.0, 0.1),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  facet_wrap(~ variant, ncol = 2) +
  labs(
    title = "Detection rate vs corruption level",
    subtitle = "z_threshold = 1.5 | dashed line = false-positive rate on clean respondents",
    x = "Corruption level (% of items replaced)",
    y = "Detection rate (sensitivity per level)",
    caption = sprintf("N=%d, %d%% careless, %d conditions × %d reps. pure_straight/acquiescent: always at 100%%.",
                      N_FIXED, PCT_CARELESS*100, nrow(conditions), REPS)
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "#EEEEEE", color = NA),
        strip.text = element_text(face = "bold"))

ggsave("plot_detection_vs_corruption.png", p1,
       width = 12, height = 6, dpi = 150, bg = "white")
cat("Saved: plot_detection_vs_corruption.png\n")

# ============================================================
# PLOT 2 — z-score distribution by corruption level (standard)
# ============================================================
# Only variable-corruption patterns (random, longstring, mixed)
plot2_df <- df %>%
  filter(pattern %in% c("random", "longstring", "mixed") |
           pattern == "clean") %>%
  mutate(corruption_cat = ifelse(pattern == "clean", "clean",
                                  sprintf("%d%%", round(corruption * 100))),
         corruption_cat = factor(corruption_cat,
                                 levels = c("clean", paste0(seq(10, 100, 10), "%"))))

p2 <- ggplot(plot2_df, aes(x = corruption_cat, y = z_standard, fill = pattern)) +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "red", linewidth = 0.4) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.7) +
  scale_fill_manual(values = c(
    clean = "#cccccc",
    random = pattern_colors["random"],
    longstring = pattern_colors["longstring"],
    mixed = pattern_colors["mixed"]
  ), name = "Pattern") +
  facet_wrap(~ pattern, ncol = 2) +
  labs(
    title = "z-score distribution by corruption level (standard)",
    subtitle = "Red dashed line = z_threshold (1.5). Respondents below the line are flagged careless.",
    x = "Corruption level",
    y = "z-score"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "#EEEEEE", color = NA),
        strip.text = element_text(face = "bold"))

ggsave("plot_zscore_by_corruption.png", p2,
       width = 11, height = 8, dpi = 150, bg = "white")
cat("Saved: plot_zscore_by_corruption.png\n")

# ============================================================
# PLOT 3 — Detection rate by corruption level, all patterns overlay
# ============================================================
p3 <- ggplot(det_std, aes(x = corruption, y = detect_std,
                          color = pattern, group = pattern)) +
  geom_hline(yintercept = det_clean$detect_std, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  annotate("text", x = 0.15, y = det_clean$detect_std + 0.03,
           label = sprintf("clean false-pos = %.1f%%", det_clean$detect_std * 100),
           size = 3.5, color = "grey30") +
  geom_line(size = 1.3) +
  geom_point(size = 3) +
  scale_color_manual(values = pattern_colors, name = "Pattern") +
  scale_x_continuous(breaks = seq(0.1, 1.0, 0.1),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Detection rate vs corruption level — standard method",
    subtitle = sprintf("z_threshold = 1.5 | %d resp per cell | pure_straight & acquiescent are always at 100%%",
                       round(N_FIXED * PCT_CARELESS / 5 * REPS / 10)),
    x = "Corruption level",
    y = "Detection rate (sensitivity at z ≤ 1.5)",
    caption = "Grey dashed line: false-positive rate on attentive respondents."
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

ggsave("plot_detection_summary.png", p3,
       width = 10, height = 6, dpi = 150, bg = "white")
cat("Saved: plot_detection_summary.png\n")

cat(sprintf("\nTotal runtime: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("=== DONE ===\n")
