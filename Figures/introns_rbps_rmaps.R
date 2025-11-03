library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggrepel)  # For better label placement

# Load data
pladb_int_up <- read.delim("rmaps_out/pladb_int_rmaps/pVal.up.vs.bg.RNAmap.txt", 
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pladb_int_down <- read.delim("rmaps_out/pladb_int_rmaps/pVal.dn.vs.bg.RNAmap.txt", 
                             header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Filter function
filter_by_any_R <- function(df, source_name) {
  r_cols <- grep("^R\\d+$", names(df), value = TRUE)
  df %>%
    filter(if_any(all_of(r_cols), ~ . <= 0.05)) %>%
    mutate(source = source_name)
}

pladb_int_up <- filter_by_any_R(pladb_int_up, "Retained")
pladb_int_down <- filter_by_any_R(pladb_int_down, "Skipped")
pladb_int <- bind_rows(pladb_int_up, pladb_int_down)

# Pivot to long format
pladb_int_long <- pladb_int %>%
  pivot_longer(
    cols = matches("^R\\d+$"),
    names_to = "position",
    values_to = "pvalue"
  ) %>%
  mutate(
    pvalue = as.numeric(pvalue),
    neglog10p = if_else(is.na(pvalue) | pvalue <= 0, NA_real_, -log10(pvalue))
  )

# Define genomic regions
region_map <- tibble(
  position = paste0("R", 1:5),
  region_short = c("Upstream\nExon", "5'SS", "Intron", "3'SS", "Downstream\nExon"),
  region_label = c("Upstream Exon\n(first 50bp)", 
                   "5' Splice Site\n(last 50bp of exon)", 
                   "Intron\n(250bp)", 
                   "3' Splice Site\n(first 50bp of exon)", 
                   "Downstream Exon\n(last 50bp)"),
  length_bp = c(50, 50, 250, 50, 50)
) %>%
  mutate(
    start = cumsum(lag(length_bp, default = 0)),
    mid = start + length_bp / 2,
    end = start + length_bp
  )

pladb_int_long <- pladb_int_long %>%
  left_join(region_map %>% select(position, mid, region_short, region_label), by = "position") %>%
  mutate(
    position = factor(position, levels = paste0("R", 1:5)),
    region_label = factor(region_label, levels = c(
      "Upstream Exon\n(first 50bp)", 
      "5' Splice Site\n(last 50bp of exon)", 
      "Intron\n(250bp)", 
      "3' Splice Site\n(first 50bp of exon)", 
      "Downstream Exon\n(last 50bp)"
    ))
  )

# ============================================================================
# GENOMIC TRACKS PLOT (Top 15 per source)
# ============================================================================

# Get top 15 RBPs per source (Retained and Skipped separately)
top_rbps_per_source <- pladb_int_long %>%
  filter(!is.na(neglog10p), pvalue <= 0.05) %>%
  group_by(RBP, source) %>%
  summarize(max_neglog = max(neglog10p, na.rm = TRUE), .groups = "drop") %>%
  group_by(source) %>%
  slice_max(order_by = max_neglog, n = 15) %>%
  ungroup()

# Filter data to only these top RBPs
track_data <- pladb_int_long %>%
  semi_join(top_rbps_per_source, by = c("RBP", "source")) %>%
  filter(!is.na(neglog10p))

# IMPROVED LABELING: Select top 5 RBPs per source AND per region for labeling
label_data <- pladb_int_long %>%
  filter(!is.na(neglog10p), pvalue <= 0.05) %>%
  semi_join(top_rbps_per_source, by = c("RBP", "source")) %>%
  group_by(source, position) %>%
  slice_max(order_by = neglog10p, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    RBP_clean = sub("\\..*", "", RBP)  # Remove everything after dot
  )

p_tracks <- ggplot() +
  # Individual RBP points
  geom_point(
    data = track_data,
    aes(x = mid, y = neglog10p, color = source),
    size = 1.5, alpha = 0.4, 
    position = position_jitter(width = 8, height = 0, seed = 42)
  ) +
  # Smoothed trend line per source
  geom_smooth(
    data = track_data,
    aes(x = mid, y = neglog10p, color = source, fill = source),
    method = "loess", span = 0.3, alpha = 0.15, linewidth = 1.5, se = TRUE
  ) +
  # Add RBP labels using ggrepel for non-overlapping placement
  geom_text_repel(
    data = label_data,
    aes(x = mid, y = neglog10p, label = RBP_clean, color = source),
    size = 2.2, fontface = "bold", alpha = 0.9,
    max.overlaps = Inf,
    min.segment.length = 0.1,
    segment.size = 0.3,
    segment.alpha = 0.6,
    box.padding = 0.3,
    point.padding = 0.2,
    force = 2,
    seed = 42
  ) +
  # Significance threshold
  geom_hline(
    yintercept = -log10(0.05), 
    linetype = "dashed", 
    color = "grey40", 
    linewidth = 0.4
  ) +
  # Region separators
  geom_vline(
    xintercept = c(50, 100, 350, 400), 
    linetype = "dotted", 
    color = "grey60", 
    linewidth = 0.25
  ) +
  # Scales and labels
  scale_x_continuous(
    breaks = region_map$mid,
    labels = region_map$region_short,
    expand = c(0.02, 0)
  ) +
  scale_color_manual(
    values = c("Retained" = "#D32F2F", "Skipped" = "#1976D2"),
    name = "Splicing Event"
  ) +
  scale_fill_manual(
    values = c("Retained" = "#D32F2F", "Skipped" = "#1976D2"),
    name = "Splicing Event"
  ) +
  labs(
    y = expression(bold(-log[10](italic(p)*"-value"))),
    x = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(face = "bold", size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 9, margin = margin(r = 6)),
    axis.line = element_line(color = "grey30", linewidth = 0.4),
    axis.ticks = element_line(color = "grey30", linewidth = 0.4),
    axis.ticks.length = unit(0.1, "cm"),
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5, margin = margin(b = 8)),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )

# ============================================================================
# BOXPLOT PANEL (All RBPs) with Statistical Tests
# ============================================================================

# Use ALL data (not just top 15)
boxplot_data <- pladb_int_long %>%
  filter(!is.na(neglog10p))

# Calculate Wilcoxon test for each position
stat_tests <- boxplot_data %>%
  group_by(region_label) %>%
  summarize(
    p_value = wilcox.test(neglog10p[source == "Retained"], 
                          neglog10p[source == "Skipped"])$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    label = ifelse(significance != "ns", 
                   paste0("p = ", format.pval(p_value, digits = 2), " ", significance),
                   "ns")
  )

# Get y position for significance labels (top of each facet)
y_positions <- boxplot_data %>%
  group_by(region_label) %>%
  summarize(y_pos = max(neglog10p, na.rm = TRUE) * 1.05, .groups = "drop")

stat_tests <- stat_tests %>%
  left_join(y_positions, by = "region_label")

p_boxplots <- ggplot(boxplot_data, aes(x = source, y = neglog10p, fill = source)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, linewidth = 0.5) +
  geom_jitter(aes(color = source), width = 0.15, alpha = 0.4, size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
             color = "grey40", linewidth = 0.4) +
  # Add significance annotations
  geom_text(
    data = stat_tests,
    aes(x = 1.5, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 2.5, fontface = "bold", color = "black"
  ) +
  facet_wrap(~ region_label, ncol = 5, scales = "free_x") +
  labs(
    y = expression(bold(-log[10](italic(p)*"-value"))),
    x = NULL
  ) +
  scale_fill_manual(
    values = c("Retained" = "#D32F2F", "Skipped" = "#1976D2"),
    name = "Splicing Event"
  ) +
  scale_color_manual(
    values = c("Retained" = "#D32F2F", "Skipped" = "#1976D2"),
    name = "Splicing Event"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    strip.text = element_text(face = "bold", size = 7),
    strip.background = element_rect(fill = "grey90", color = "grey50", linewidth = 0.4),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.4),
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 9, margin = margin(r = 6)),
    axis.line.y = element_line(color = "grey30", linewidth = 0.4),
    axis.ticks.y = element_line(color = "grey30", linewidth = 0.4),
    axis.ticks.length.y = unit(0.1, "cm"),
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5, margin = margin(b = 8)),
    plot.margin = margin(5, 5, 5, 5)
  )

# ============================================================================
# COMBINE PLOTS
# ============================================================================

p_combined <- p_tracks / p_boxplots + 
  plot_layout(heights = c(1, 1.2), guides = "collect") +
  plot_annotation(
    theme = theme(
      legend.position = "top",
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8),
      legend.margin = margin(0, 0, 5, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    )
  )

print(p_combined)

# Save optimized for half A4 (portrait orientation)
# Half A4 portrait: 8.27 inches width, ~5.8 inches height
ggsave("rbps_box_tracks_ssa.pdf", p_combined, 
       width = 8.27, height = 5.8, device = cairo_pdf, dpi = 300)

