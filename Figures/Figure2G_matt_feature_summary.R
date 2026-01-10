library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(showtext)
library(forcats)
library(ggtext)

# ==============================================================================
# SIMPLIFIED DOT MATRIX PLOT - FEATURE EFFECT DIRECTION
# ==============================================================================
#
# PURPOSE: Clean, publication-ready dot matrix showing:
#   - Direction of effect (increase/decrease) for significant features
#   - Organized by: Treatment (columns) × Feature (rows)
#   - Only essential features
#   - Dot size = significance (p-value)
#   - Dot color = direction (red=higher, blue=lower)
#
# ==============================================================================

# ---- Setup fonts and theme ----
showtext::showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)
base_family <- "sans"

# Drug treatment colors
drug_colors <- c(
  "Spliceostatin A" = "#706993",
  "Pladienolide B" = "#70A0AF",
  "Tubercidin" = "#A0C1B9"
)

# Colors for event direction (Included/Skipped, Retained/Excised)
event_colors <- c(
  "Included" = "#D4A574",
  "Skipped" = "#B85450",
  "Retained" = "#D4A574",
  "Excised" = "#B85450"
)

# Colors for effect direction (harmonious with existing palette)
direction_colors <- c(
  "Increase" = "#C85450",    # Beautiful warm red (higher than background)
  "Decrease" = "#5B8FA3"     # Beautiful cool blue (lower than background)
)

# ---- Read data ----
message("Reading data files...")
tub_ex <- read_tsv("Supplementary3_matt_out/tub_ex/tub_exons_with_efeatures.tab", show_col_types = FALSE)
tub_int <- read_tsv("Supplementary3_matt_out/tub_int/tub_introns_with_ifeatures.tab", show_col_types = FALSE)
pladb_ex <- read_tsv("Supplementary3_matt_out/pladb_ex/pladb_exons_with_efeatures.tab", show_col_types = FALSE)
pladb_int <- read_tsv("Supplementary3_matt_out/pladb_int/pladb_introns_with_ifeatures.tab", show_col_types = FALSE)
ssa_ex <- read_tsv("Supplementary3_matt_out/ssa_ex/ssa_exons_with_efeatures.tab", show_col_types = FALSE)
ssa_int <- read_tsv("Supplementary3_matt_out/ssa_int/ssa_introns_with_ifeatures.tab", show_col_types = FALSE)

# Read existing statistical results
stats <- read_csv("feature_comparison_statistics.csv", show_col_types = FALSE)

# ---- Function to calculate effect direction ----
calculate_effect_directions <- function(df, treatment_name, event_type) {
  message("Calculating effect directions for ", treatment_name, " ", event_type, "...")
  
  exclude_cols <- c("START", "END", "SCAFFOLD", "STRAND", "GENEID_ENSEMBL",
                    "DATASET", "EXON_ID", "GENE_BIOTYPE", "EXON_FOUND_IN_GTF",
                    "EXON_FOUND_IN_THESE_TRS", "EXON_COOCCURS_WITH_THESE_EXONS",
                    grep("^SEQ_", colnames(df), value = TRUE))
  
  numeric_cols <- colnames(df)[sapply(df, is.numeric)]
  feature_cols <- setdiff(numeric_cols, exclude_cols)
  
  effect_data <- df %>%
    filter(DATASET %in% c("up", "down", "ndiff")) %>%
    select(DATASET, all_of(feature_cols)) %>%
    pivot_longer(cols = -DATASET, names_to = "feature", values_to = "value") %>%
    filter(!is.na(value)) %>%
    group_by(feature, DATASET) %>%
    summarise(
      median_value = median(value, na.rm = TRUE),
      mean_value = mean(value, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    pivot_wider(
      names_from = DATASET,
      values_from = c(median_value, mean_value)
    ) %>%
    mutate(
      treatment = treatment_name,
      event_type = event_type,
      direction_up_vs_ndiff = median_value_up - median_value_ndiff,
      direction_down_vs_ndiff = median_value_down - median_value_ndiff,
      log2fc_up_vs_ndiff = log2((median_value_up + 0.01) / (median_value_ndiff + 0.01)),
      log2fc_down_vs_ndiff = log2((median_value_down + 0.01) / (median_value_ndiff + 0.01))
    )
  
  return(effect_data)
}

# ---- Calculate effect directions for all datasets ----
message("\n=== Calculating effect directions ===")
all_effects <- bind_rows(
  calculate_effect_directions(tub_ex, "Tubercidin", "exons"),
  calculate_effect_directions(tub_int, "Tubercidin", "introns"),
  calculate_effect_directions(pladb_ex, "Pladienolide B", "exons"),
  calculate_effect_directions(pladb_int, "Pladienolide B", "introns"),
  calculate_effect_directions(ssa_ex, "Spliceostatin A", "exons"),
  calculate_effect_directions(ssa_int, "Spliceostatin A", "introns")
)

write_csv(all_effects, "feature_effect_directions.csv")
message("Saved effect directions to feature_effect_directions.csv")

# ---- Merge with statistical significance ----
stats_filtered <- stats %>%
  filter(group1 %in% c("up", "down"), group2 == "ndiff") %>%
  select(feature, treatment, event_type, comparison = group1, p.value, p.signif)

effects_long <- all_effects %>%
  select(feature, treatment, event_type,
         median_value_up, median_value_down, median_value_ndiff,
         direction_up_vs_ndiff, direction_down_vs_ndiff,
         log2fc_up_vs_ndiff, log2fc_down_vs_ndiff) %>%
  pivot_longer(
    cols = c(direction_up_vs_ndiff, direction_down_vs_ndiff),
    names_to = "comparison_type",
    values_to = "direction"
  ) %>%
  mutate(
    comparison = case_when(
      comparison_type == "direction_up_vs_ndiff" ~ "up",
      comparison_type == "direction_down_vs_ndiff" ~ "down"
    )
  ) %>%
  select(-comparison_type)

log2fc_long <- all_effects %>%
  select(feature, treatment, event_type, log2fc_up_vs_ndiff, log2fc_down_vs_ndiff) %>%
  pivot_longer(
    cols = c(log2fc_up_vs_ndiff, log2fc_down_vs_ndiff),
    names_to = "comparison_type",
    values_to = "log2fc"
  ) %>%
  mutate(
    comparison = case_when(
      comparison_type == "log2fc_up_vs_ndiff" ~ "up",
      comparison_type == "log2fc_down_vs_ndiff" ~ "down"
    )
  ) %>%
  select(feature, treatment, event_type, comparison, log2fc)

effects_long <- effects_long %>%
  left_join(log2fc_long, by = c("feature", "treatment", "event_type", "comparison"))

combined_data <- effects_long %>%
  left_join(stats_filtered, by = c("feature", "treatment", "event_type", "comparison")) %>%
  mutate(
    direction_label = case_when(
      is.na(direction) ~ "No data",
      direction > 0 ~ "Increase",
      direction < 0 ~ "Decrease",
      TRUE ~ "No change"
    ),
    sig_level = case_when(
      is.na(p.value) ~ "ns",
      p.value <= 0.001 ~ "***",
      p.value <= 0.01 ~ "**",
      p.value <= 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    is_significant = p.value <= 0.05
  )

write_csv(combined_data, "feature_consensus_with_direction.csv")
message("Saved combined data to feature_consensus_with_direction.csv")

# ---- Filter for significant features only ----
significant_data <- combined_data %>%
  filter(is_significant == TRUE) %>%
  mutate(
    event_label = case_when(
      event_type == "exons" & comparison == "up" ~ "Included",
      event_type == "exons" & comparison == "down" ~ "Skipped",
      event_type == "introns" & comparison == "up" ~ "Retained",
      event_type == "introns" & comparison == "down" ~ "Excised"
    ),
    feature_display = gsub("_", " ", feature),
    feature_display = gsub("MEDIANLENGTH", "Length", feature_display),
    feature_display = gsub("MAXENTSCR HSAMODEL", "MaxEnt", feature_display),
    feature_display = gsub("UPINTRON", "Up Intron", feature_display),
    feature_display = gsub("DOINTRON", "Down Intron", feature_display)
  )

message("\nTotal significant comparisons: ", nrow(significant_data))

# ---- Publication-ready theme ----
theme_matrix_pub <- function(base_size = 14, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = rel(1.0), hjust = 1, color = "grey20"),
      axis.title = element_blank(),
      axis.ticks = element_line(color = "grey50", linewidth = 0.3),
      axis.ticks.x = element_blank(),
      axis.ticks.length = unit(0.12, "cm"),
      strip.text = element_text(face = "bold", size = rel(1.1)),
      strip.background = element_rect(fill = "grey95", color = "grey70", linewidth = 0.3),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = rel(1.35), hjust = 0.5,
                                margin = margin(b = 8, t = 2)),
      plot.subtitle = element_text(size = rel(1.0), hjust = 0.5,
                                   color = "grey40", margin = margin(b = 10)),
      plot.margin = margin(25, 20, 10, 12)
    )
}

# Prepare data for plotting
plot_data <- significant_data %>%
  mutate(
    treatment = factor(treatment, levels = c("Tubercidin", "Pladienolide B", "Spliceostatin A")),
    event_label = factor(event_label, levels = c("Included", "Skipped", "Retained", "Excised")),
    # Dot size based on significance level (discrete sizes)
    dot_size = case_when(
      sig_level == "***" ~ 6,    # Biggest
      sig_level == "**" ~ 4.5,   # Bigger
      sig_level == "*" ~ 3,      # Small
      TRUE ~ 0                   # No dot (non-significant)
    ),
    # For legend display
    sig_category = factor(sig_level, levels = c("*", "**", "***")),
    # Color based on direction
    direction_color = case_when(
      direction_label == "Increase" ~ "Increase",
      direction_label == "Decrease" ~ "Decrease",
      TRUE ~ "No change"
    )
  )

# ---- Selected features for plotting ----
exon_features_keep <- c(
  "UPINTRON_MEDIANLENGTH",
  "DOINTRON_MEDIANLENGTH",
  "MEDIAN_EXON_NUMBER"
)

intron_features_keep <- c(
  "INTRON_LENGTH",
  "INTRON_GCC",
  "MAXENTSCR_HSAMODEL_5SS"
)

# ---- PLOT 1A: EXONS ----
plot_data_exons <- plot_data %>%
  filter(event_type == "exons", feature %in% exon_features_keep)

exon_display_map <- c(
  "UPINTRON_MEDIANLENGTH" = "Upstream Intron Length",
  "DOINTRON_MEDIANLENGTH" = "Downstream Intron Length",
  "MEDIAN_EXON_NUMBER" = "Exon Number"
)

if (nrow(plot_data_exons) > 0) {
  present_feats <- intersect(exon_features_keep, unique(plot_data_exons$feature))
  plot_data_exons <- plot_data_exons %>%
    mutate(
      feature_display = exon_display_map[feature],
      feature_display = factor(feature_display, levels = rev(exon_display_map[present_feats]))
    )
  
  n_features <- length(unique(plot_data_exons$feature_display))
  
  # Event type annotation bars
  event_bars <- data.frame(
    xmin = c(0.5, 3.5),
    xmax = c(3.5, 6.5),
    ymin = n_features + 0.85,
    ymax = n_features + 1.05,
    event = c("Included", "Skipped"),
    stringsAsFactors = FALSE
  )
  
  # Drug annotation bars
  drug_bars <- data.frame(
    xmin = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5),
    xmax = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5),
    ymin = n_features + 0.55,
    ymax = n_features + 0.75,
    drug = c("Tubercidin", "Pladienolide B", "Spliceostatin A",
             "Tubercidin", "Pladienolide B", "Spliceostatin A"),
    stringsAsFactors = FALSE
  )
  
  p1_exons <- ggplot(plot_data_exons,
                     aes(x = interaction(treatment, event_label, sep = "\n"),
                         y = feature_display)) +
    geom_rect(
      data = drug_bars,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = drug),
      inherit.aes = FALSE,
      color = "white",
      linewidth = 0.5
    ) +
    geom_rect(
      data = event_bars,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = event),
      inherit.aes = FALSE,
      color = "white",
      linewidth = 0.5
    ) +
    geom_text(
      data = event_bars,
      aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = event),
      inherit.aes = FALSE,
      size = 3.2,
      fontface = "bold",
      color = "white",
      family = base_family
    ) +
    geom_point(
      aes(color = direction_color, size = dot_size),
      alpha = 0.85
    ) +
    scale_color_manual(
      name = "Direction",
      values = direction_colors,
      breaks = c("Increase", "Decrease"),
      labels = c("Higher", "Lower")
    ) +
    scale_fill_manual(
      name = NULL,
      values = c(drug_colors, event_colors),
      breaks = c("Included", "Skipped"),
      labels = c("Included", "Skipped")
    ) +
    scale_size_identity() +
    scale_x_discrete(expand = expansion(add = c(0.08, 0.08))) +
    scale_y_discrete(expand = expansion(add = c(0.1, 1.5))) +
    labs(
      title = "Exon Splicing Features",
      subtitle = NULL,
      x = NULL,
      y = NULL
    ) +
    theme_matrix_pub(base_family = base_family) +
    guides(
      fill = guide_legend(
        override.aes = list(shape = 22, size = 3),
        nrow = 1,
        title = NULL,
        order = 1
      ),
      color = guide_legend(
        override.aes = list(size = 3.5),
        nrow = 1,
        order = 2,
        title = "Direction"
      )
    ) +
    # Add manual size legend
    geom_point(
      data = data.frame(
        x = rep(Inf, 3),
        y = rep(Inf, 3),
        size = c(3, 4.5, 6),
        label = c("*", "**", "***")
      ),
      aes(x = x, y = y, size = size),
      inherit.aes = FALSE,
      alpha = 0
    ) +
    scale_size_identity(
      name = "p-value",
      breaks = c(3, 4.5, 6),
      labels = c("* (0.05)", "** (0.01)", "*** (0.001)"),
      guide = guide_legend(
        override.aes = list(alpha = 0.85, color = "grey40"),
        nrow = 1,
        order = 3
      )
    )
  
  # Calculate height accounting for: features + annotation bars (1.05 above top) + margins + legend space
  plot_height <- max(4.5, n_features * 0.3 + 3.5)
  
  ggsave(
    "Figures/dot_matrix_exons.pdf",
    plot = p1_exons,
    width = 5.0,
    height = plot_height,
    device = cairo_pdf,
    bg = "white"
  )
  
  ggsave(
    "Figures/dot_matrix_exons.png",
    plot = p1_exons,
    width = 5.0,
    height = plot_height,
    dpi = 300,
    bg = "white"
  )
  
  message("Saved exon dot matrix to Figures/dot_matrix_exons.pdf/.png")
}

# ---- PLOT 1B: INTRONS ----
plot_data_introns <- plot_data %>%
  filter(event_type == "introns", feature %in% intron_features_keep)

intron_display_map <- c(
  "INTRON_LENGTH" = "Intron Length",
  "INTRON_GCC" = "Intron GC Content",
  "MAXENTSCR_HSAMODEL_5SS" = "MaxEnt 5'SS Score"
)

if (nrow(plot_data_introns) > 0) {
  present_feats_i <- intersect(intron_features_keep, unique(plot_data_introns$feature))
  plot_data_introns <- plot_data_introns %>%
    mutate(
      feature_display = intron_display_map[feature],
      feature_display = factor(feature_display, levels = rev(intron_display_map[present_feats_i]))
    )
  
  n_features_i <- length(unique(plot_data_introns$feature_display))
  
  event_bars_i <- data.frame(
    xmin = c(0.5, 3.5),
    xmax = c(3.5, 6.5),
    ymin = n_features_i + 0.85,
    ymax = n_features_i + 1.05,
    event = c("Retained", "Excised"),
    stringsAsFactors = FALSE
  )
  
  drug_bars_i <- data.frame(
    xmin = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5),
    xmax = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5),
    ymin = n_features_i + 0.55,
    ymax = n_features_i + 0.75,
    drug = c("Tubercidin", "Pladienolide B", "Spliceostatin A",
             "Tubercidin", "Pladienolide B", "Spliceostatin A"),
    stringsAsFactors = FALSE
  )
  
  p1_introns <- ggplot(plot_data_introns,
                       aes(x = interaction(treatment, event_label, sep = "\n"),
                           y = feature_display)) +
    geom_rect(
      data = drug_bars_i,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = drug),
      inherit.aes = FALSE,
      color = "white",
      linewidth = 0.5
    ) +
    geom_rect(
      data = event_bars_i,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = event),
      inherit.aes = FALSE,
      color = "white",
      linewidth = 0.5
    ) +
    geom_text(
      data = event_bars_i,
      aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = event),
      inherit.aes = FALSE,
      size = 3.2,
      fontface = "bold",
      color = "white",
      family = base_family
    ) +
    geom_point(
      aes(color = direction_color, size = dot_size),
      alpha = 0.85
    ) +
    scale_color_manual(
      name = "Direction",
      values = direction_colors,
      breaks = c("Increase", "Decrease"),
      labels = c("Higher", "Lower")
    ) +
    scale_fill_manual(
      name = NULL,
      values = c(drug_colors, event_colors),
      breaks = c("Retained", "Excised"),
      labels = c("Retained", "Excised")
    ) +
    scale_size_identity() +
    scale_x_discrete(expand = expansion(add = c(0.08, 0.08))) +
    scale_y_discrete(expand = expansion(add = c(0.1, 1.5))) +
    labs(
      title = "Intron Splicing Features",
      subtitle = NULL,
      x = NULL,
      y = NULL
    ) +
    theme_matrix_pub(base_family = base_family) +
    guides(
      fill = guide_legend(
        override.aes = list(shape = 22, size = 3),
        nrow = 1,
        title = NULL,
        order = 1
      ),
      color = guide_legend(
        override.aes = list(size = 3.5),
        nrow = 1,
        order = 2,
        title = "Direction"
      )
    ) +
    # Add manual size legend
    geom_point(
      data = data.frame(
        x = rep(Inf, 3),
        y = rep(Inf, 3),
        size = c(3, 4.5, 6),
        label = c("*", "**", "***")
      ),
      aes(x = x, y = y, size = size),
      inherit.aes = FALSE,
      alpha = 0
    ) +
    scale_size_identity(
      name = "p-value",
      breaks = c(3, 4.5, 6),
      labels = c("* (0.05)", "** (0.01)", "*** (0.001)"),
      guide = guide_legend(
        override.aes = list(alpha = 0.85, color = "grey40"),
        nrow = 1,
        order = 3
      )
    )
  
  # Calculate height accounting for: features + annotation bars (1.05 above top) + margins + legend space
  plot_height_i <- max(4.5, n_features_i * 0.3 + 3.5)
  
  ggsave(
    "Figures/dot_matrix_introns.pdf",
    plot = p1_introns,
    width = 5.0,
    height = plot_height_i,
    device = cairo_pdf,
    bg = "white"
  )
  
  ggsave(
    "Figures/dot_matrix_introns.png",
    plot = p1_introns,
    width = 5.0,
    height = plot_height_i,
    dpi = 300,
    bg = "white"
  )
  
  message("Saved intron dot matrix to Figures/dot_matrix_introns.pdf/.png")
}
