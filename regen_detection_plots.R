# Regenerate plots from existing sim_detection_per_resp.csv
# (no re-simulation, just re-binned plotting)

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(tidyr)
})

df <- read.csv("sim_detection_per_resp.csv", stringsAsFactors = FALSE)

# Clean 10% bins
df <- df %>%
  mutate(corruption_bin10 = ifelse(pattern == "clean", NA,
                                   round(corruption * 10) / 10))

# Check sample sizes per (pattern, bin)
cat("=== Sample sizes per (pattern, corruption bin) ===\n")
df %>% filter(pattern != "clean") %>%
  group_by(pattern, corruption_bin10) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = corruption_bin10, values_from = n, values_fill = 0) %>%
  as.data.frame() %>% print(row.names = FALSE)

# Detection rates
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
cat("\n=== DETECTION RATE per pattern × corruption bin ===\n")
print(as.data.frame(det_std), row.names = FALSE)

pattern_colors <- c(
  random = "#e41a1c", longstring = "#377eb8", mixed = "#4daf4a",
  pure_straight = "#984ea3", acquiescent = "#ff7f00"
)

# PLOT 1 — Main: detection rate vs corruption (std + varpen side by side)
plot_df <- det_std %>%
  pivot_longer(cols = c(detect_std, detect_vp),
               names_to = "variant", values_to = "detect_rate") %>%
  mutate(variant = recode(variant,
                          "detect_std" = "Standard",
                          "detect_vp" = "Variance-Penalized (α=3, β=0.5)"))

p1 <- ggplot(plot_df, aes(x = corruption, y = detect_rate,
                          color = pattern, group = pattern)) +
  geom_hline(yintercept = det_clean$detect_std, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  annotate("text", x = 0.15, y = det_clean$detect_std + 0.04,
           label = sprintf("clean false-pos = %.1f%%", det_clean$detect_std * 100),
           size = 3.2, color = "grey30") +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = pattern_colors, name = "Pattern") +
  scale_x_continuous(breaks = seq(0.1, 1.0, 0.1),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  facet_wrap(~ variant, ncol = 2) +
  labs(
    title = "Detection rate vs corruption level",
    subtitle = "z_threshold = 1.5 — the method's core challenge: partial corruption (20-60%) barely above false-positive floor",
    x = "Corruption level (% of items replaced)",
    y = "Detection rate (sensitivity per level)",
    caption = "pure_straight / acquiescent: always at 100% corruption by construction."
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "#EEEEEE", color = NA),
        strip.text = element_text(face = "bold"))

ggsave("plot_detection_vs_corruption.png", p1,
       width = 13, height = 6.5, dpi = 150, bg = "white")
cat("Saved: plot_detection_vs_corruption.png\n")

# PLOT 2 — summary single-panel
p3 <- ggplot(det_std, aes(x = corruption, y = detect_std,
                          color = pattern, group = pattern)) +
  geom_hline(yintercept = det_clean$detect_std, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  annotate("text", x = 0.15, y = det_clean$detect_std + 0.03,
           label = sprintf("clean false-pos = %.1f%%", det_clean$detect_std * 100),
           size = 3.5, color = "grey30") +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3) +
  scale_color_manual(values = pattern_colors, name = "Pattern") +
  scale_x_continuous(breaks = seq(0.1, 1.0, 0.1),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Detection rate vs corruption level — standard method",
    subtitle = "z_threshold = 1.5 | pure_straight & acquiescent always at 100%",
    x = "Corruption level",
    y = "Detection rate (sensitivity)",
    caption = "Grey dashed line = false-positive rate on attentive respondents."
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

ggsave("plot_detection_summary.png", p3,
       width = 10, height = 6, dpi = 150, bg = "white")
cat("Saved: plot_detection_summary.png\n")

# PLOT 2b — detection by corruption, split by questionnaire length
det_by_size <- df %>%
  filter(pattern != "clean") %>%
  mutate(size_cat = case_when(
    total_items < 60 ~ "<60 items",
    total_items < 120 ~ "60-119 items",
    total_items < 200 ~ "120-199 items",
    TRUE ~ ">=200 items"
  )) %>%
  mutate(size_cat = factor(size_cat,
                           levels = c("<60 items","60-119 items",
                                      "120-199 items",">=200 items"))) %>%
  group_by(size_cat, pattern, corruption_bin10) %>%
  summarise(detect = mean(z_standard <= 1.5, na.rm = TRUE), n = n(),
            .groups = "drop") %>%
  rename(corruption = corruption_bin10)

p4 <- ggplot(det_by_size, aes(x = corruption, y = detect,
                              color = pattern, group = pattern)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = pattern_colors, name = "Pattern") +
  scale_x_continuous(breaks = seq(0.1, 1.0, 0.2),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  facet_wrap(~ size_cat, nrow = 1) +
  labs(
    title = "Detection rate by questionnaire length and corruption",
    subtitle = "Longer questionnaires detect partial carelessness at lower corruption levels",
    x = "Corruption level",
    y = "Detection rate (z ≤ 1.5)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "#EEEEEE", color = NA),
        strip.text = element_text(face = "bold"))

ggsave("plot_detection_by_size.png", p4,
       width = 13, height = 5, dpi = 150, bg = "white")
cat("Saved: plot_detection_by_size.png\n")

cat("=== DONE ===\n")
