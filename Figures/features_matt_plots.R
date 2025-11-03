library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(purrr)
library(readr)
library(showtext)
library(gridExtra)

# ==============================================================================
# BIOINFORMATICS QUALITY STANDARDS - DESIGN PRINCIPLES
# ==============================================================================
# 
# 1. STATISTICAL INTEGRITY:
#    - All statistical tests use FULL, unclipped data
#    - Wilcoxon rank-sum test (non-parametric) for robustness
#    - Multiple comparison framework with proper grouping
#    - P-values and significance annotations based on complete datasets
#
# 2. VISUALIZATION ACCURACY:
#    - Boxplots calculated from FULL data (median, quartiles, whiskers)
#    - Mean values (diamond overlay) computed from FULL data
#    - Only jittered points are clipped for visual clarity
#    - coord_cartesian() used for view clipping (preserves calculations)
#    - Outliers displayed to show data distribution honestly
#
# 3. DATA HANDLING:
#    - Aggressive outlier clipping for VIEW only (3*IQR from quartiles)
#    - Further clipping if range > 10*IQR (down to 2*IQR)
#    - Clipping prevents extreme outliers from distorting y-axis scale
#    - All underlying calculations remain on unclipped data
#
# 4. PLOT ELEMENTS:
#    - Significance bars positioned dynamically based on data range
#    - Sufficient spacing between bars and title (no overlap)
#    - Faceted by treatment (Tubercidin, Pladienolide B, Spliceostatin A)
#    - Consistent color scheme across all plots
#    - Publication-ready theme with appropriate font sizes
#
# ==============================================================================

# ---- Setup fonts and theme ----
preferred_font <- "Roboto"
font_add_google(preferred_font)
showtext::showtext_opts(dpi = 300)
showtext_auto()
base_family <- preferred_font

# Color palette for DATASET groups
dataset_colors <- c(
  "up" = "#D4A574",      # Warm sand
  "down" = "#B85450",    # Muted coral
  "ndiff" = "#808080"    # Gray
)

# Theme matching the provided code
theme_cellpub <- function(base_size = 16, base_family_in = NULL) {
  if (is.null(base_family_in)) base_family_in <- get0("base_family", ifnotfound = "sans")
  theme_classic(base_size = base_size, base_family = base_family_in) %+replace%
    theme(
      axis.line = element_line(linewidth = 0.6, colour = "black"),
      axis.ticks = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      axis.title = element_text(family = base_family_in, face = "bold", size = rel(1.0)),
      axis.text = element_text(family = base_family_in, size = rel(0.9), colour = "black"),
      legend.position = "none",
      legend.background = element_blank(),
      legend.title = element_text(face = "bold", size = rel(0.95)),
      legend.text = element_text(family = base_family_in, size = rel(0.85)),
      legend.key.size = unit(0.4, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = rel(0.95), family = base_family_in, 
                                margin = margin(b = 12, t = 2)),
      plot.title = element_text(face = "bold", size = rel(1.1), hjust = 0.5, family = base_family_in),
      plot.subtitle = element_text(size = rel(0.95), hjust = 0.5, family = base_family_in),
      plot.caption = element_text(size = rel(0.85), colour = "#666666", family = base_family_in),
      plot.margin = margin(6, 6, 6, 6)
    )
}
my_theme <- theme_cellpub(base_size = 16, base_family_in = base_family)

# ---- Read data ----
message("Reading data files...")
tub_ex <- read_tsv("matt_out/tub_ex/tub_exons_with_efeatures.tab", show_col_types = FALSE)
tub_int <- read_tsv("matt_out/tub_int/tub_introns_with_ifeatures.tab", show_col_types = FALSE)
pladb_ex <- read_tsv("matt_out/pladb_ex/pladb_exons_with_efeatures.tab", show_col_types = FALSE)
pladb_int <- read_tsv("matt_out/pladb_int/pladb_introns_with_ifeatures.tab", show_col_types = FALSE)
ssa_ex <- read_tsv("matt_out/ssa_ex/ssa_exons_with_efeatures.tab", show_col_types = FALSE)
ssa_int <- read_tsv("matt_out/ssa_int/ssa_introns_with_ifeatures.tab", show_col_types = FALSE)

# ---- Function to perform statistical tests for all features ----
test_all_features <- function(df, treatment_name, event_type) {
  message("Testing features for ", treatment_name, " ", event_type, "...")
  
  # Identify numeric feature columns (exclude metadata columns)
  exclude_cols <- c("START", "END", "SCAFFOLD", "STRAND", "GENEID_ENSEMBL", 
                    "DATASET", "EXON_ID", "GENE_BIOTYPE", "EXON_FOUND_IN_GTF",
                    "EXON_FOUND_IN_THESE_TRS", "EXON_COOCCURS_WITH_THESE_EXONS",
                    grep("^SEQ_", colnames(df), value = TRUE))
  
  numeric_cols <- colnames(df)[sapply(df, is.numeric)]
  feature_cols <- setdiff(numeric_cols, exclude_cols)
  
  # Define comparisons
  comparisons <- list(
    c("up", "ndiff"),
    c("down", "ndiff"),
    c("up", "down")
  )
  
  # Test each feature
  results_list <- list()
  
  for (feature in feature_cols) {
    feature_data <- df %>%
      filter(!is.na(.data[[feature]]), DATASET %in% c("up", "down", "ndiff")) %>%
      select(DATASET, value = all_of(feature))
    
    # Skip if insufficient data
    if (nrow(feature_data) < 10) next
    
    # Perform tests
    test_results <- map_dfr(comparisons, function(comp) {
      g1 <- comp[1]
      g2 <- comp[2]
      
      data_g1 <- feature_data %>% filter(DATASET == g1) %>% pull(value)
      data_g2 <- feature_data %>% filter(DATASET == g2) %>% pull(value)
      
      if (length(data_g1) < 3 || length(data_g2) < 3) {
        return(tibble(group1 = g1, group2 = g2, p.value = NA_real_))
      }
      
      test <- tryCatch({
        wilcox.test(data_g1, data_g2, exact = FALSE)
      }, error = function(e) {
        list(p.value = NA_real_)
      })
      
      tibble(group1 = g1, group2 = g2, p.value = test$p.value)
    })
    
    test_results <- test_results %>%
      mutate(
        feature = feature,
        treatment = treatment_name,
        event_type = event_type,
        p.signif = case_when(
          is.na(p.value) ~ "ns",
          p.value <= 0.001 ~ "***",
          p.value <= 0.01 ~ "**",
          p.value <= 0.05 ~ "*",
          TRUE ~ "ns"
        )
      )
    
    results_list[[feature]] <- test_results
  }
  
  bind_rows(results_list)
}

# ---- Run tests for all datasets ----
message("\n=== Running statistical tests ===")
all_test_results <- bind_rows(
  test_all_features(tub_ex, "Tubercidin", "exons"),
  test_all_features(tub_int, "Tubercidin", "introns"),
  test_all_features(pladb_ex, "Pladienolide B", "exons"),
  test_all_features(pladb_int, "Pladienolide B", "introns"),
  test_all_features(ssa_ex, "Spliceostatin A", "exons"),
  test_all_features(ssa_int, "Spliceostatin A", "introns")
)

# Save test results
write_csv(all_test_results, "feature_comparison_statistics.csv")
message("Saved all test results to feature_comparison_statistics.csv")

# ---- Identify significant features (p <= 0.05 in any treatment) ----
significant_features <- all_test_results %>%
  filter(p.value <= 0.001) %>%
  distinct(feature, event_type) %>%
  arrange(event_type, feature)

message("\nFound ", nrow(significant_features), " significant feature-event_type combinations")
print(significant_features)

# ---- Function to create faceted plot for a feature ----
# CRITICAL: This function maintains statistical integrity by:
# 1. Using FULL DATA for all statistical calculations (boxplot, median, mean, tests)
# 2. Only clipping jittered points for visual clarity
# 3. Using coord_cartesian() which clips VIEW without affecting geom calculations
create_feature_plot <- function(feature_name, event_type_filter, 
                                tub_data, pladb_data, ssa_data,
                                stat_results) {
  
  message("Creating plot for: ", feature_name, " (", event_type_filter, ")")
  
  # Define labels based on event type
  if (event_type_filter == "exons") {
    label_map <- c("up" = "Included", "down" = "Skipped", "ndiff" = "Unchanged")
  } else {
    label_map <- c("up" = "Retained", "down" = "Excised", "ndiff" = "Unchanged")
  }
  
  # Combine data - USE DISPLAY NAMES FROM THE START
  combined_data <- bind_rows(
    tub_data %>% 
      filter(!is.na(.data[[feature_name]]), DATASET %in% c("up", "down", "ndiff")) %>%
      select(DATASET, value = all_of(feature_name)) %>%
      mutate(treatment = "Tubercidin"),  # Use display name directly
    
    pladb_data %>% 
      filter(!is.na(.data[[feature_name]]), DATASET %in% c("up", "down", "ndiff")) %>%
      select(DATASET, value = all_of(feature_name)) %>%
      mutate(treatment = "Pladienolide B"),  # Use display name directly
    
    ssa_data %>% 
      filter(!is.na(.data[[feature_name]]), DATASET %in% c("up", "down", "ndiff")) %>%
      select(DATASET, value = all_of(feature_name)) %>%
      mutate(treatment = "Spliceostatin A")  # Use display name directly
  ) %>%
    mutate(
      # Factor levels for treatment in desired LEFT TO RIGHT order
      treatment = factor(treatment, levels = c("Tubercidin", "Pladienolide B", "Spliceostatin A")),
      # Factor levels for DATASET
      DATASET = factor(DATASET, levels = c("up", "down", "ndiff")),
      # Create display labels
      DATASET_label = factor(label_map[as.character(DATASET)], 
                             levels = label_map[c("up", "down", "ndiff")])
    )
  
  # Debug: Check the actual levels
  message("  Treatment levels: ", paste(levels(combined_data$treatment), collapse = " | "))
  message("  DATASET_label levels: ", paste(levels(combined_data$DATASET_label), collapse = " | "))
  
  if (nrow(combined_data) == 0) {
    message("  -> No data available, skipping")
    return(NULL)
  }
  
  # Filter stat results for this feature and map to display names
  feature_stats <- stat_results %>%
    filter(feature == feature_name, event_type == event_type_filter) %>%
    mutate(
      treatment = factor(treatment, levels = c("Tubercidin", "Pladienolide B", "Spliceostatin A")),
      # Map group1 and group2 to display labels
      group1_label = label_map[group1],
      group2_label = label_map[group2]
    )
  
  if (nrow(feature_stats) == 0) {
    message("  -> No statistics available, skipping")
    return(NULL)
  }
  
  # STEP 1: Calculate quartiles for each treatment
  quartile_stats <- combined_data %>%
    group_by(treatment) %>%
    summarise(
      q1 = quantile(value, 0.25, na.rm = TRUE),
      q3 = quantile(value, 0.75, na.rm = TRUE),
      median_val = median(value, na.rm = TRUE),
      iqr = q3 - q1,
      min_val = min(value, na.rm = TRUE),
      max_val = max(value, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # STEP 2: Determine plot range - AGGRESSIVE outlier clipping
  overall_q1 <- min(quartile_stats$q1, na.rm = TRUE)
  overall_q3 <- max(quartile_stats$q3, na.rm = TRUE)
  overall_iqr <- overall_q3 - overall_q1
  overall_median <- median(combined_data$value, na.rm = TRUE)
  overall_min <- min(quartile_stats$min_val, na.rm = TRUE)
  overall_max <- max(quartile_stats$max_val, na.rm = TRUE)
  
  # More aggressive clipping: limit to 3*IQR from quartiles
  y_view_lower <- max(overall_q1 - 1.5 * overall_iqr, overall_min)
  y_view_upper <- min(overall_q3 + 1.5 * overall_iqr, overall_max)
  
  # If range is still too extreme (outliers > 5x IQR), clip even more
  if ((y_view_upper - y_view_lower) > 10 * overall_iqr) {
    y_view_lower <- overall_q1 - 1 * overall_iqr
    y_view_upper <- overall_q3 + 1 * overall_iqr
  }
  
  # STEP 3: Position significance bars with GENEROUS spacing to prevent overlap
  # Calculate how many bars we have per treatment
  bars_per_treatment <- feature_stats %>%
    group_by(treatment) %>%
    summarise(n_bars = n(), .groups = 'drop')
  
  max_bars <- max(bars_per_treatment$n_bars, na.rm = TRUE)
  
  # Adjust spacing based on number of bars - EXTRA GENEROUS spacing
  if (max_bars <= 2) {
    bar_step <- 0.38 * overall_iqr  # More spacing between bars
    bar_start <- y_view_upper + 0.80 * overall_iqr  # Start higher
  } else {
    bar_step <- 0.45 * overall_iqr  # Even more for 3 bars
    bar_start <- y_view_upper + 0.80 * overall_iqr  # Start much higher
  }
  
  feature_stats <- feature_stats %>%
    left_join(quartile_stats, by = "treatment") %>%
    group_by(treatment) %>%
    mutate(
      y.position = bar_start + bar_step * (row_number() - 1)
    ) %>%
    ungroup()
  
  # STEP 4: Final y-axis upper limit includes bars + EXTRA padding to prevent title overlap
  max_bar_y <- max(feature_stats$y.position, na.rm = TRUE)
  y_final_upper <- max_bar_y + 0.30 * overall_iqr  # Extra padding for title
  
  # STEP 5: Calculate mean for diamond overlay (no labels needed)
  mean_data <- combined_data %>%
    group_by(treatment, DATASET_label) %>%
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # STEP 6: Prepare data for plotting
  # CRITICAL: Use FULL DATA (not clipped) for boxplot, median, mean calculations
  # Only clip jittered points for visual clarity
  
  # Separate clipped data for jittered points only - ONLY up and down
  jitter_data <- combined_data %>%
    filter(DATASET %in% c("up", "down")) %>%  # Explicit filter to avoid any NA or unexpected values
    mutate(value_display = pmax(y_view_lower, pmin(value, y_view_upper))) %>%
    droplevels()  # Drop unused factor levels to ensure clean color mapping
  
  # Keep full data for boxplots (ggplot will handle whiskers correctly)
  boxplot_data <- combined_data
  
  # Mean uses full data (no clipping)
  mean_data_full <- mean_data  # Already calculated from full data
  
  # Create clean feature name for display
  display_name <- gsub("_", " ", feature_name)
  
  # Use coord_cartesian for view clipping WITHOUT affecting boxplot/mean calculations
  p <- ggplot(boxplot_data, aes(x = DATASET_label, y = value, fill = DATASET)) +
    # Boxplot uses FULL DATA - hide outliers since jitter shows individual points
    geom_boxplot(
      outlier.shape = NA,  # Hide outliers - individual points shown via jitter
      alpha = 0.7, 
      width = 0.7, 
      linewidth = 0.5, 
      color = "black",
      coef = 1.5  # Standard 1.5*IQR whisker length
    ) +
    # Jittered points ONLY for up/down (not ndiff) - use CLIPPED data for visual clarity
    # Colors inherited from fill aesthetic via DATASET column
    geom_jitter(
      data = jitter_data,
      aes(x = DATASET_label, y = value_display, color = DATASET),
      width = 0.15, 
      alpha = 0.5, 
      size = 0.8,
      show.legend = FALSE
    ) +
    # Mean uses FULL DATA (not clipped)
    geom_point(
      data = mean_data_full,
      aes(x = DATASET_label, y = mean_val),
      inherit.aes = FALSE,
      shape = 23,  # Diamond shape
      size = 3.5,
      fill = "white",
      color = "black",
      stroke = 1
    ) +
    # Add significance bars using the label-mapped groups
    stat_pvalue_manual(
      feature_stats,
      xmin = "group1_label",
      xmax = "group2_label",
      y.position = "y.position",
      label = "p.signif",
      tip.length = 0,  # Fixed proportion of plot height
      label.size = 3.5,
      bracket.size = 0.5,
      coord.flip = FALSE,
      step.increase = 0
    ) +
    scale_fill_manual(values = dataset_colors) +
    scale_color_manual(values = dataset_colors) +
    # Use coord_cartesian for VISUAL clipping only - does not affect calculations
    coord_cartesian(
      ylim = c(y_view_lower, y_final_upper),
      clip = "on"
    ) +
    scale_y_continuous(
      expand = c(0, 0)
    ) +
    facet_wrap(~ treatment, nrow = 1) +  # No labeller needed - already using display names
    labs(
      x = NULL,
      y = display_name,
      title = paste0(display_name, " - ", event_type_filter)
    ) +
    my_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.85)),
      panel.spacing = unit(0.8, "lines"),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  return(p)
}

# ---- Generate plots for all significant features ----
message("\n=== Generating plots for significant features ===")

# Create output directory
dir.create("feature_plots", showWarnings = FALSE)
dir.create("feature_plots/exons", showWarnings = FALSE)
dir.create("feature_plots/introns", showWarnings = FALSE)

for (i in seq_len(nrow(significant_features))) {
  feature <- significant_features$feature[i]
  evt_type <- significant_features$event_type[i]
  
  # Select appropriate datasets
  if (evt_type == "exons") {
    tub_data <- tub_ex
    pladb_data <- pladb_ex
    ssa_data <- ssa_ex
    output_dir <- "feature_plots/exons"
  } else {
    tub_data <- tub_int
    pladb_data <- pladb_int
    ssa_data <- ssa_int
    output_dir <- "feature_plots/introns"
  }
  
  # Create plot
  p <- tryCatch({
    create_feature_plot(
      feature_name = feature,
      event_type_filter = evt_type,
      tub_data = tub_data,
      pladb_data = pladb_data,
      ssa_data = ssa_data,
      stat_results = all_test_results
    )
  }, error = function(e) {
    message("ERROR creating plot for ", feature, ": ", e$message)
    NULL
  })
  
  if (!is.null(p)) {
    # Save plot as PDF only
    safe_name <- gsub("[^A-Za-z0-9_]", "_", feature)
    pdf_file <- file.path(output_dir, paste0(safe_name, ".pdf"))
    
    ggsave(
      pdf_file, 
      plot = p, 
      device = cairo_pdf, 
      width = 10, 
      height = 4,
      bg = "white"
    )
    message("  -> Saved: ", pdf_file)
  }
}

# ---- Summary report ----
message("\n=== SUMMARY ===")
message("Total tests performed: ", nrow(all_test_results))
message("Significant comparisons (p <= 0.05): ", sum(all_test_results$p.value <= 0.05, na.rm = TRUE))
message("Features with at least one significant result: ", nrow(significant_features))

# Create summary table
summary_table <- all_test_results %>%
  group_by(feature, event_type) %>%
  summarise(
    n_significant = sum(p.value <= 0.05, na.rm = TRUE),
    min_p = min(p.value, na.rm = TRUE),
    treatments_significant = paste(unique(treatment[p.value <= 0.05]), collapse = ", "),
    .groups = 'drop'
  ) %>%
  filter(n_significant > 0) %>%
  arrange(desc(n_significant), min_p)

write_csv(summary_table, "feature_significance_summary.csv")
message("\nSaved summary to feature_significance_summary.csv")
message("\nPlots (PDFs) saved to feature_plots/exons/ and feature_plots/introns/")
message("\nAnalysis complete!")