# Regenerate detection plots using auto_z threshold per condition
# (instead of fixed z=1.5)

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(tidyr)
})
source("ReReReRe.R")  # provides get_calibrated_z()

df <- read.csv("sim_detection_per_resp.csv", stringsAsFactors = FALSE)

# Compute auto_z per unique (nF, ipf, total_items) — same threshold for all reps in a cell
cond_thresholds <- df %>%
  distinct(nF, ipf, total_items) %>%
  rowwise() %>%
  mutate(auto_z = round(get_calibrated_z(total_items), 3)) %>%
  ungroup()

cat("=== Auto-z threshold per condition ===\n")
print(as.data.frame(cond_thresholds), row.names = FALSE)

# Merge thresholds back into df
df <- df %>%
  left_join(cond_thresholds, by = c("nF","ipf","total_items")) %>%
  mutate(flagged_autoz = z_standard <= auto_z,
         flagged_autoz_vp = z_varpen <= auto_z,
         corruption_bin10 = ifelse(pattern == "clean", NA,
                                   round(corruption * 10) / 10))

# Detection rate per pattern x corruption bin (auto-z)
det_az <- df %>%
  filter(pattern != "clean") %>%
  group_by(pattern, corruption_bin10) %>%
  summarise(n = n(),
            detect_std = mean(flagged_autoz, na.rm = TRUE),
            detect_vp  = mean(flagged_autoz_vp, na.rm = TRUE),
            .groups = "drop") %>%
  rename(corruption = corruption_bin10)

det_clean_az <- df %>%
  filter(pattern == "clean") %>%
  summarise(detect_std = mean(flagged_autoz, na.rm = TRUE),
            detect_vp  = mean(flagged_autoz_vp, na.rm = TRUE))

cat("\n=== CLEAN false-positive rate at auto-z ===\n")
print(det_clean_az)
cat("\n=== DETECTION RATE per pattern × corruption bin (auto-z) ===\n")
print(as.data.frame(det_az), row.names = FALSE)

pattern_colors <- c(
  random = "#e41a1c", longstring = "#377eb8", mixed = "#4daf4a",
  pure_straight = "#984ea3", acquiescent = "#ff7f00"
)

# ================================================
# PLOT 1 — detection rate vs corruption (std vs varpen, auto-z)
# ================================================
plot_df <- det_az %>%
  pivot_longer(cols = c(detect_std, detect_vp),
               names_to = "variant", values_to = "detect_rate") %>%
  mutate(variant = recode(variant,
                          "detect_std" = "Standard (auto-z)",
                          "detect_vp"  = "Variance-Penalized + auto-z"))

p1 <- ggplot(plot_df, aes(x = corruption, y = detect_rate,
                          color = pattern, group = pattern)) +
  geom_hline(yintercept = det_clean_az$detect_std, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  annotate("text", x = 0.15, y = det_clean_az$detect_std + 0.04,
           label = sprintf("clean false-pos = %.1f%%",
                           det_clean_az$detect_std * 100),
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
    title = "Detection rate vs corruption level — AUTO-Z (2D calibrated)",
    subtitle = "z_threshold depends on total_items (0.53 at 48 items -> 2.18 at 300 items)",
    x = "Corruption level (% of items replaced)",
    y = "Detection rate (sensitivity per level)",
    caption = "pure_straight / acquiescent: always at 100% corruption by construction."
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "#EEEEEE", color = NA),
        strip.text = element_text(face = "bold"))

ggsave("plot_detection_vs_corruption_autoz.png", p1,
       width = 13, height = 6.5, dpi = 150, bg = "white")
cat("Saved: plot_detection_vs_corruption_autoz.png\n")

# ================================================
# PLOT 2 — summary single panel (auto-z)
# ================================================
p2 <- ggplot(det_az, aes(x = corruption, y = detect_std,
                         color = pattern, group = pattern)) +
  geom_hline(yintercept = det_clean_az$detect_std, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  annotate("text", x = 0.15, y = det_clean_az$detect_std + 0.03,
           label = sprintf("clean false-pos = %.1f%%",
                           det_clean_az$detect_std * 100),
           size = 3.5, color = "grey30") +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3) +
  scale_color_manual(values = pattern_colors, name = "Pattern") +
  scale_x_continuous(breaks = seq(0.1, 1.0, 0.1),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Detection rate vs corruption level — AUTO-Z",
    subtitle = "Threshold adapted per questionnaire length",
    x = "Corruption level",
    y = "Detection rate (sensitivity)",
    caption = "Grey dashed line = false-positive rate on attentive respondents."
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

ggsave("plot_detection_summary_autoz.png", p2,
       width = 10, height = 6, dpi = 150, bg = "white")
cat("Saved: plot_detection_summary_autoz.png\n")

# ================================================
# PLOT 3 — detection by questionnaire length (auto-z)
# ================================================
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
  summarise(detect = mean(flagged_autoz, na.rm = TRUE), n = n(),
            .groups = "drop") %>%
  rename(corruption = corruption_bin10)

# Clean FPR per size
clean_fpr_size <- df %>%
  filter(pattern == "clean") %>%
  mutate(size_cat = case_when(
    total_items < 60 ~ "<60 items",
    total_items < 120 ~ "60-119 items",
    total_items < 200 ~ "120-199 items",
    TRUE ~ ">=200 items"
  )) %>%
  group_by(size_cat) %>%
  summarise(fpr = mean(flagged_autoz, na.rm = TRUE), .groups = "drop") %>%
  mutate(size_cat = factor(size_cat,
                           levels = c("<60 items","60-119 items",
                                      "120-199 items",">=200 items")))

p3 <- ggplot(det_by_size, aes(x = corruption, y = detect,
                              color = pattern, group = pattern)) +
  geom_hline(data = clean_fpr_size,
             aes(yintercept = fpr), linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_text(data = clean_fpr_size,
            aes(x = 0.35, y = fpr + 0.05,
                label = sprintf("FPR=%.0f%%", fpr*100)),
            inherit.aes = FALSE, size = 3, color = "grey30") +
  scale_color_manual(values = pattern_colors, name = "Pattern") +
  scale_x_continuous(breaks = seq(0.1, 1.0, 0.2),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  facet_wrap(~ size_cat, nrow = 1) +
  labs(
    title = "Detection rate by questionnaire length — AUTO-Z",
    subtitle = "Auto-z adapts the threshold per questionnaire size, equalizing the false-positive rate across sizes",
    x = "Corruption level",
    y = "Detection rate (z <= auto_z)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "#EEEEEE", color = NA),
        strip.text = element_text(face = "bold"))

ggsave("plot_detection_by_size_autoz.png", p3,
       width = 13, height = 5, dpi = 150, bg = "white")
cat("Saved: plot_detection_by_size_autoz.png\n")

# Print clean FPR per size for reference
cat("\n=== FPR per size category (auto-z) ===\n")
print(as.data.frame(clean_fpr_size), row.names = FALSE)

cat("=== DONE ===\n")
