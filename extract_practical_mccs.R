# Extract practical (fixed z=1.5) MCCs per cell from sim_variants_results.csv
# Also compute auto-z MCC using 2D lookup (z = 0.505 + 0.0042*total_items clamped)

suppressPackageStartupMessages(library(dplyr))

df <- read.csv("sim_variants_results.csv", stringsAsFactors = FALSE)

# Restrict to 4 main methods we care about
keep <- c("std", "efa_d", "iterative")
df <- df %>% filter(method %in% keep)

# Fixed z=1.5 MCC (closest z in the grid)
z_grid <- unique(df$z_threshold)
target_z <- z_grid[which.min(abs(z_grid - 1.5))]
cat(sprintf("Using z_threshold closest to 1.5: %.1f\n", target_z))

df_z15 <- df %>% filter(z_threshold == target_z) %>%
  group_by(nF, ipf, total_items, method) %>%
  summarise(mcc_z15 = mean(mcc, na.rm = TRUE), .groups = "drop")

# Oracle per cell per method
df_oracle <- df %>%
  group_by(nF, ipf, total_items, rep, method) %>%
  summarise(oracle = max(mcc, na.rm = TRUE), .groups = "drop") %>%
  group_by(nF, ipf, total_items, method) %>%
  summarise(mcc_oracle = mean(oracle, na.rm = TRUE), .groups = "drop")

# Auto-z: for each cell, pick z closest to 0.505 + 0.0042*total_items clamped [0.3, 3.5]
cells <- df %>% distinct(nF, ipf, total_items)
cells$z_auto <- pmin(3.5, pmax(0.3, 0.505 + 0.0042 * cells$total_items))
# Snap to nearest in z_grid
cells$z_auto_snap <- z_grid[apply(outer(cells$z_auto, z_grid, function(a,b) abs(a-b)), 1, which.min)]

df_auto <- df %>%
  inner_join(cells %>% select(nF, ipf, total_items, z_auto_snap),
             by = c("nF", "ipf", "total_items")) %>%
  filter(abs(z_threshold - z_auto_snap) < 0.01) %>%
  group_by(nF, ipf, total_items, method) %>%
  summarise(mcc_auto = mean(mcc, na.rm = TRUE), .groups = "drop")

# Merge
out <- df_oracle %>%
  left_join(df_z15, by = c("nF","ipf","total_items","method")) %>%
  left_join(df_auto, by = c("nF","ipf","total_items","method")) %>%
  arrange(method, total_items)

write.csv(out, "practical_mccs_per_cell.csv", row.names = FALSE)

# Print summary per cell for std method (for Table 2 / peak conditions)
cat("\n=== STD method per-cell summary ===\n")
print(as.data.frame(out %>% filter(method == "std") %>%
  mutate(across(starts_with("mcc"), ~ round(., 3)))), row.names = FALSE)

# Overall summary by bin
cat("\n=== By items bin, mean across methods ===\n")
summ <- out %>%
  mutate(bin = cut(total_items,
                   breaks = c(0, 30, 60, 100, 200, 1000),
                   labels = c("<30", "30-60", "60-100", "100-200", ">200"))) %>%
  group_by(bin, method) %>%
  summarise(oracle = round(mean(mcc_oracle), 3),
            fixed_z15 = round(mean(mcc_z15), 3),
            auto_z = round(mean(mcc_auto), 3),
            .groups = "drop")
print(as.data.frame(summ), row.names = FALSE)

# Global overall means across all cells
cat("\n=== Global overall means (across all 15 cells) ===\n")
glob <- out %>%
  group_by(method) %>%
  summarise(oracle = round(mean(mcc_oracle), 3),
            fixed_z15 = round(mean(mcc_z15), 3),
            auto_z = round(mean(mcc_auto), 3),
            .groups = "drop")
print(as.data.frame(glob), row.names = FALSE)
