# Plot_Iterative_Comparison.R
# Confronto iterative vs metodi precedenti
# Due plot: MCC vs total_items e MCC vs nF

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

# Lettura risultati
res <- read.csv("sim_variants_results.csv", stringsAsFactors = FALSE)

# Tieni solo i metodi principali da confrontare
methods_to_plot <- c("std", "efa_d", "iterative", "pf_mean")
method_labels <- c(
  "std" = "Standard (coupled/weighted)",
  "efa_d" = "EFA-D (per-factor)",
  "iterative" = "Iterative EFA",
  "pf_mean" = "Per-factor mean(z)"
)

method_colors <- c(
  "std" = "#377eb8",       # blu
  "efa_d" = "#4daf4a",     # verde
  "iterative" = "#e41a1c", # rosso
  "pf_mean" = "#984ea3"    # viola
)

# Oracle MCC per condizione×rep×metodo (una riga per rep, non per z_threshold)
oracle <- res %>%
  filter(method %in% methods_to_plot) %>%
  group_by(nF, ipf, total_items, rep, method) %>%
  summarise(mcc = max(mcc, na.rm = TRUE), .groups = "drop")

# ============================================================
# PLOT 1 — MCC vs total_items
# ============================================================

summ_items <- oracle %>%
  group_by(total_items, method) %>%
  summarise(
    mean_mcc = mean(mcc, na.rm = TRUE),
    sd_mcc = sd(mcc, na.rm = TRUE),
    n = n(),
    se = sd_mcc / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(method_label = method_labels[method])

p1 <- ggplot(summ_items, aes(x = total_items, y = mean_mcc,
                              color = method, group = method)) +
  geom_ribbon(aes(ymin = mean_mcc - se, ymax = mean_mcc + se, fill = method),
              alpha = 0.15, color = NA) +
  geom_line(size = 1.1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = method_colors, labels = method_labels,
                     name = "Metodo") +
  scale_fill_manual(values = method_colors, guide = "none") +
  scale_x_continuous(breaks = c(12, 24, 40, 48, 72, 80, 90, 120, 180, 200, 300)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(
    title = "Oracle MCC vs Total Items",
    subtitle = "Iterative EFA domina tra 60 e 200 item; Standard vince a >=300",
    x = "Total items (nF x ipf)",
    y = "Oracle MCC (max su z_threshold)",
    caption = "N=300, 15% careless, corruzione 10-100%, GT >50%. Ombra = SE tra reps."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("plot_iterative_vs_total_items.png", p1,
       width = 10, height = 6, dpi = 150, bg = "white")
cat("Saved: plot_iterative_vs_total_items.png\n")

# ============================================================
# PLOT 2 — MCC vs nF (facet per ipf)
# ============================================================

summ_nF <- oracle %>%
  group_by(nF, ipf, method) %>%
  summarise(
    mean_mcc = mean(mcc, na.rm = TRUE),
    sd_mcc = sd(mcc, na.rm = TRUE),
    n = n(),
    se = sd_mcc / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    method_label = method_labels[method],
    ipf_label = paste0("ipf = ", ipf, " (items per factor)")
  )

p2 <- ggplot(summ_nF, aes(x = nF, y = mean_mcc,
                           color = method, group = method)) +
  geom_ribbon(aes(ymin = mean_mcc - se, ymax = mean_mcc + se, fill = method),
              alpha = 0.15, color = NA) +
  geom_line(size = 1.1) +
  geom_point(size = 2.5) +
  facet_wrap(~ ipf_label, ncol = 3) +
  scale_color_manual(values = method_colors, labels = method_labels,
                     name = "Metodo") +
  scale_fill_manual(values = method_colors, guide = "none") +
  scale_x_continuous(breaks = c(4, 8, 12, 20, 30)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(
    title = "Oracle MCC vs Number of Factors (per items-per-factor level)",
    subtitle = "Iterative domina a ipf=6 e 10 (questionari medio-grandi)",
    x = "Number of Factors (nF)",
    y = "Oracle MCC",
    caption = "N=300, 15% careless, corruzione 10-100%, GT >50%. Ombra = SE tra reps."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "#EEEEEE", color = NA),
    strip.text = element_text(face = "bold")
  )

ggsave("plot_iterative_vs_nF.png", p2,
       width = 12, height = 6, dpi = 150, bg = "white")
cat("Saved: plot_iterative_vs_nF.png\n")

# ============================================================
# PLOT 3 — Differenza iterative - std (guadagno visivo)
# ============================================================

diff_df <- oracle %>%
  filter(method %in% c("std", "iterative")) %>%
  pivot_wider(names_from = method, values_from = mcc) %>%
  mutate(delta = iterative - std) %>%
  group_by(total_items, nF, ipf) %>%
  summarise(
    mean_delta = mean(delta, na.rm = TRUE),
    se_delta = sd(delta, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p3 <- ggplot(diff_df, aes(x = total_items, y = mean_delta)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = mean_delta - se_delta, ymax = mean_delta + se_delta),
              fill = "#e41a1c", alpha = 0.2) +
  geom_line(color = "#e41a1c", size = 1.1) +
  geom_point(aes(size = factor(ipf), color = factor(ipf))) +
  scale_color_manual(values = c("3" = "#fdd49e", "6" = "#fc8d59", "10" = "#b30000"),
                     name = "items/factor") +
  scale_size_manual(values = c("3" = 2.5, "6" = 3.5, "10" = 4.5), guide = "none") +
  scale_x_continuous(breaks = c(12, 24, 40, 48, 72, 80, 90, 120, 180, 200, 300)) +
  labs(
    title = "Guadagno Iterative EFA rispetto a Standard",
    subtitle = "Positivo = iterative meglio. Massimo guadagno tra 60-200 items.",
    x = "Total items",
    y = "Delta MCC (iterative - std)",
    caption = "Ombra = SE. Guadagno medio overall: +0.015 MCC."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("plot_iterative_gain.png", p3,
       width = 10, height = 6, dpi = 150, bg = "white")
cat("Saved: plot_iterative_gain.png\n")

# Stampa tabella riassuntiva
cat("\n=== SUMMARY: Mean MCC by method and items bin ===\n")
oracle %>%
  mutate(bin = cut(total_items,
                   breaks = c(0, 30, 60, 100, 200, 1000),
                   labels = c("<30", "30-60", "60-100", "100-200", ">200"))) %>%
  group_by(bin, method) %>%
  summarise(mcc = round(mean(mcc, na.rm = TRUE), 3), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = mcc) %>%
  print()

cat("\nDone.\n")
