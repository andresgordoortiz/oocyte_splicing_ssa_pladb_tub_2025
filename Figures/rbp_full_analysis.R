  library(dplyr)
  library(ggplot2)
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
  library(viridis)
  library(ggridges)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)

  # ============================================================================
  # DEFINE CONSISTENT THEME AND COLORS
  # ============================================================================

  # Drug colors (order: TUB, PLADB, SSA)
  drug_colors <- c(
    "Tubercidin" = "#A0C1B9",        # Light sage green
    "Pladienolide B" = "#70A0AF",    # Teal/blue-green
    "Spliceostatin A" = "#706993"    # Purple/lavender
  )

  # Harmonious splicing colors (complementary to drug palette)
  splicing_colors <- c(
    "Retained" = "#D4A574",  # Warm sand (complementary to cool drug colors)
    "Excised" = "#B85450"    # Muted coral (harmonious contrast)
  )

  # Region colors (matching the FULL region names from region_map)
  region_colors <- c(
    "Upstream Exon" = "#E8E3F5",        # Very light lavender
    "5' Splice Site" = "#C8C3D8",       # Light purple-gray
    "Intron" = "#A8A3BB",               # Medium purple-gray
    "3' Splice Site" = "#4B4856",       # Deeper purple-gray
    "Downstream Exon" = "#372F3F"       # Dark purple-gray
  )

  # Custom theme function matching volcano plot style
  theme_publication <- function(base_size = 10) {
    theme_classic(base_size = base_size) %+replace%
      theme(
        # Axes
        axis.line = element_line(linewidth = 0.6, colour = "black"),
        axis.ticks = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.title = element_text(face = "bold", size = rel(1.0)),
        axis.text = element_text(size = rel(0.9), colour = "black"),

        # Legend
        legend.position = "right",
        legend.background = element_blank(),
        legend.title = element_text(face = "bold", size = rel(0.95)),
        legend.text = element_text(size = rel(0.85)),
        legend.key.size = unit(0.4, "cm"),
        legend.spacing.y = unit(0.1, "cm"),

        # Panel
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),

        # Facets
        strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5),
        strip.text = element_text(face = "bold", size = rel(1.0), margin = margin(5, 5, 5, 5)),

        # Plot
        plot.title = element_text(face = "bold", size = rel(1.1), hjust = 0),
        plot.subtitle = element_text(size = rel(0.9), hjust = 0, color = "grey30"),
        plot.margin = margin(10, 10, 10, 10)
      )
  }

  # ============================================================================
  # DATA LOADING AND PROCESSING FUNCTION
  # ============================================================================

  # ============================================================================
  # DATA LOADING AND PROCESSING FUNCTION
  # ============================================================================
  
  process_drug_data <- function(drug_name) {
    # Load data
    up_file <- paste0("rmaps_out/", drug_name, "_int_rmaps/pVal.up.vs.bg.RNAmap.txt")
    down_file <- paste0("rmaps_out/", drug_name, "_int_rmaps/pVal.dn.vs.bg.RNAmap.txt")
    
    up_data <- read.delim(up_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    down_data <- read.delim(down_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Filter function
    filter_by_any_R <- function(df, source_name) {
      r_cols <- grep("^R\\d+$", names(df), value = TRUE)
      df %>%
        filter(if_any(all_of(r_cols), ~ . <= 0.05)) %>%
        mutate(source = source_name, drug = drug_name)
    }
    
    up_filtered <- filter_by_any_R(up_data, "Retained")      # FIXED: Introns retained
    down_filtered <- filter_by_any_R(down_data, "Excised")   # FIXED: Introns excised (spliced out)
    
    combined <- bind_rows(up_filtered, down_filtered)
    
    # Pivot to long format
    combined_long <- combined %>%
      pivot_longer(
        cols = matches("^R\\d+$"),
        names_to = "position",
        values_to = "pvalue"
      ) %>%
      mutate(
        pvalue = as.numeric(pvalue),
        neglog10p = if_else(is.na(pvalue) | pvalue <= 0, NA_real_, -log10(pvalue)),
        RBP_clean = sub("\\..*", "", RBP)  # Clean RBP names
      )
    
    return(combined_long)
  }

  # ============================================================================
  # LOAD AND COMBINE ALL DRUG DATA
  # ============================================================================

  all_drugs <- map_dfr(c("pladb", "ssa", "tub"), process_drug_data)

  # Define genomic regions with refined labels
  region_map <- tibble(
    position = paste0("R", 1:5),
    region_short = c("Upstream", "5'SS", "Intron", "3'SS", "Downstream"),
    region_full = c("Upstream Exon", "5' Splice Site", "Intron", "3' Splice Site", "Downstream Exon"),
    length_bp = c(50, 50, 250, 50, 50)
  ) %>%
    mutate(
      start = cumsum(lag(length_bp, default = 0)),
      mid = start + length_bp / 2,
      end = start + length_bp
    )

  all_drugs <- all_drugs %>%
    left_join(region_map, by = "position") %>%
    mutate(
      position = factor(position, levels = paste0("R", 1:5)),
      region_full = factor(region_full, levels = c("Upstream Exon", "5' Splice Site", "Intron", "3' Splice Site", "Downstream Exon")),
      drug = factor(drug, levels = c("tub", "pladb", "ssa")),
      drug_label = factor(drug,
                          levels = c("tub", "pladb", "ssa"),
                          labels = c("Tubercidin", "Pladienolide B", "Spliceostatin A"))
    )

  # ============================================================================
  # PANEL A: GENOMIC TRACK OVERVIEW
  # ============================================================================
  
  cat("Creating Panel A: Genomic Track Overview...\n")
  
  # Top RBPs per drug per source (top 15)
  top_rbps_tracks <- all_drugs %>%
    filter(!is.na(neglog10p), pvalue <= 0.05) %>%
    group_by(RBP_clean, drug, source) %>%
    summarize(max_neglog = max(neglog10p, na.rm = TRUE), .groups = "drop") %>%
    group_by(drug, source) %>%
    slice_max(order_by = max_neglog, n = 15) %>%
    ungroup()
  
  track_data <- all_drugs %>%
    inner_join(top_rbps_tracks, by = c("RBP_clean", "drug", "source")) %>%
    filter(!is.na(neglog10p))
  
  # Label top 5 RBPs per region per drug per source
  label_data <- track_data %>%
    group_by(drug_label, position, source) %>%
    slice_max(order_by = neglog10p, n = 3, with_ties = FALSE) %>%
    ungroup()
  
  # Vertical separators at region boundaries
  vline_positions <- region_map$end[-length(region_map$end)]
  
  p_tracks <- ggplot() +
    geom_point(
      data = track_data,
      aes(x = mid, y = neglog10p, color = source),
      size = 2, alpha = 0.5,
      position = position_jitter(width = 8, height = 0, seed = 42)
    ) +
    geom_smooth(
      data = track_data,
      aes(x = mid, y = neglog10p, color = source, fill = source),
      method = "loess", span = 0.3, alpha = 0.15, linewidth = 1.2, se = TRUE
    ) +
    geom_vline(
      xintercept = vline_positions,
      linetype = "dotted",
      color = "grey60",
      linewidth = 0.25
    ) +
    geom_text_repel(
      data = label_data,
      aes(x = mid, y = neglog10p, label = RBP_clean, color = source),
      size = 3, fontface = "bold", alpha = 0.9,
      max.overlaps = Inf,
      min.segment.length = 0.1,
      segment.size = 0.3,
      segment.alpha = 0.6,
      box.padding = 0.3,
      point.padding = 0.2,
      force = 2,
      seed = 42
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.4
    ) +
    facet_wrap(~ drug_label, ncol = 3, scales = "free_x") +
    scale_x_continuous(
      breaks = region_map$mid,
      labels = region_map$region_short,
      expand = c(0.02, 0)
    ) +
    scale_color_manual(values = splicing_colors, name = "Splicing Event") +
    scale_fill_manual(values = splicing_colors, name = "Splicing Event") +
    labs(
      title = "RNA-Binding Protein Enrichment Across Intron Regions",
      y = expression(bold(-log[10](italic(p)*"-value"))),
      x = NULL
    ) +
    theme_publication(base_size = 11) +
    theme(
      legend.position = "top",
      legend.margin = margin(0, 0, 5, 0),
      panel.spacing = unit(1.5, "lines"),
      axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5, vjust = 1, margin = margin(t = 8)),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = rel(1.1), margin = margin(b = 8))
    )
  
  # ============================================================================
  # PANEL B: BOXPLOT COMPARISON WITH SIGNIFICANCE TESTS
  # ============================================================================
  
  cat("Creating Panel B: Statistical Boxplots...\n")
  
  boxplot_data <- all_drugs %>%
    filter(!is.na(neglog10p)) %>%
    mutate(source = factor(source, levels = c("Retained", "Excised")))  # FIXED
  
  # Wilcoxon tests
  stat_tests <- boxplot_data %>%
    group_by(region_full, drug_label) %>%
    summarize(
      p_value = tryCatch({
        wilcox.test(neglog10p[source == "Retained"],  # FIXED
                    neglog10p[source == "Excised"])$p.value  # FIXED
      }, error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    mutate(
      significance = case_when(
        is.na(p_value) ~ "ns",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      label = ifelse(significance != "ns", significance, "")
    )
  
  # Y positions for significance labels
  y_positions <- boxplot_data %>%
    group_by(region_full, drug_label) %>%
    summarize(y_pos = max(neglog10p, na.rm = TRUE) * 1.15, .groups = "drop")
  
  stat_tests <- stat_tests %>%
    left_join(y_positions, by = c("region_full", "drug_label"))
  
  p_boxplots <- ggplot(boxplot_data, aes(x = drug_label, y = neglog10p, fill = source)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, linewidth = 0.5,
                 width = 0.7, position = position_dodge(width = 0.8)) +
    geom_jitter(aes(color = source), alpha = 0.4, size = 1.2,
                position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey40", linewidth = 0.5) +
    geom_text(
      data = stat_tests %>% filter(label != ""),
      aes(x = drug_label, y = y_pos, label = label),
      inherit.aes = FALSE,
      size = 4, fontface = "bold", color = "black"
    ) +
    facet_wrap(~ region_full, ncol = 5, scales = "fixed") +
    labs(
      x = NULL,
      y = expression(bold(-log[10](italic(p)*"-value")))
    ) +
    scale_fill_manual(values = splicing_colors, name = "Splicing Event") +
    scale_color_manual(values = splicing_colors, name = "Splicing Event") +
    theme_publication(base_size = 10) +
    theme(
      legend.position = "top",
      legend.margin = margin(0, 0, 3, 0),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, margin = margin(t = 5)),
      panel.spacing = unit(0.6, "lines"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = rel(0.8), margin = margin(b = 6))
    )
  
  # ============================================================================
  # PANEL C: COMPARATIVE HEATMAP (FIXED - proper column ordering)
  # ============================================================================
  
  cat("Creating Panel C: Comparative Heatmap...\n")
  
  # Collapse multiple entries by taking lowest p-value (highest -log10p)
  all_drugs_collapsed <- all_drugs %>%
    filter(!is.na(neglog10p)) %>%
    group_by(RBP_clean, drug_label, position, source) %>%
    summarize(
      neglog10p = max(neglog10p, na.rm = TRUE),
      pvalue = min(pvalue, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Get top 30 RBPs overall
  top_rbps_heatmap <- all_drugs_collapsed %>%
    filter(pvalue <= 0.05) %>%
    group_by(RBP_clean) %>%
    summarize(max_neglog = max(neglog10p, na.rm = TRUE), .groups = "drop") %>%
    slice_max(order_by = max_neglog, n = 30) %>%
    pull(RBP_clean)
  
  cat(sprintf("Selected %d top RBPs for heatmap\n", length(top_rbps_heatmap)))
  
  heatmap_data <- all_drugs_collapsed %>%
    filter(RBP_clean %in% top_rbps_heatmap) %>%
    dplyr::select(RBP_clean, drug_label, position, source, neglog10p) %>%
    mutate(neglog10p = replace_na(neglog10p, 0))
  
  # FIXED: Use only R1-R5 (the 5 positions that are actually defined in region_map)
  positions <- paste0("R", 1:5)
  drugs <- c("Tubercidin", "Pladienolide B", "Spliceostatin A")
  sources <- c("Retained", "Excised")
  
  # Create ordered column names for R1-R5 only
  col_order <- c()
  for (pos in positions) {
    for (drug in drugs) {
      for (src in sources) {
        col_order <- c(col_order, paste(drug, src, pos, sep = "_"))
      }
    }
  }
  
  # Pivot with explicit factor ordering to ensure correct column creation
  heatmap_data <- heatmap_data %>%
    filter(position %in% positions) %>%  # Keep only R1-R5
    mutate(
      drug_label = factor(drug_label, levels = drugs),
      position = factor(position, levels = positions),
      source = factor(source, levels = sources)
    )
  
  heatmap_wide <- heatmap_data %>%
    unite("col_name", drug_label, source, position, sep = "_", remove = FALSE) %>%
    dplyr::select(RBP_clean, col_name, neglog10p) %>%
    pivot_wider(
      names_from = col_name,
      values_from = neglog10p,
      values_fill = 0
    )
  
  heatmap_matrix <- heatmap_wide %>%
    column_to_rownames("RBP_clean") %>%
    as.matrix()
  
  # Ensure columns are in the correct order
  existing_cols <- intersect(col_order, colnames(heatmap_matrix))
  heatmap_matrix <- heatmap_matrix[, existing_cols, drop = FALSE]
  
  cat(sprintf("Heatmap matrix: %d RBPs x %d conditions\n", 
              nrow(heatmap_matrix), ncol(heatmap_matrix)))
  
  # Color function
  col_fun <- colorRamp2(
    c(0, 2, 4, 6, 8),
    c("white", "#F3E6D0", "#F0A879","#C7664A", "#702632")
  )
  
  # Column annotations - must match column order EXACTLY
  n_positions <- length(positions)  # 5
  n_drugs <- length(drugs)          # 3
  n_sources <- length(sources)      # 2
  # Total: 5 × 3 × 2 = 30 columns
  
  # Create annotation vectors - all should be length 30
  position_ann <- rep(region_map$region_full[1:5], each = n_drugs * n_sources)
  drug_ann <- rep(rep(drugs, each = n_sources), times = n_positions)
  source_ann <- rep(sources, times = n_drugs * n_positions)
  
  # Verify lengths
  cat(sprintf("Annotation lengths - Position: %d, Drug: %d, Source: %d\n",
              length(position_ann), length(drug_ann), length(source_ann)))
  cat(sprintf("Matrix columns: %d\n", ncol(heatmap_matrix)))
  
  position_ann_factor <- factor(position_ann, levels = region_map$region_full[1:5])
  drug_ann_factor <- factor(drug_ann, levels = drugs)
  
  column_ha <- HeatmapAnnotation(
    Position = position_ann_factor,
    Drug = drug_ann_factor,
    Splicing = source_ann,
    col = list(
      Drug = drug_colors,
      Splicing = splicing_colors,
      Position = region_colors
    ),
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    simple_anno_size = unit(5, "mm"),
    gap = unit(1.5, "mm")
  )
  
  # Row clustering
  row_clust <- hclust(dist(heatmap_matrix, method = "euclidean"), method = "complete")
  
  # Create heatmap
  ht <- Heatmap(
    heatmap_matrix,
    name = "-log10(p)",
    col = col_fun,
    cluster_rows = row_clust,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 9, fontface = "bold"),
    column_split = position_ann_factor,
    column_title = NULL,
    top_annotation = column_ha,
    row_title = "RNA-Binding Proteins",
    row_title_gp = gpar(fontsize = 11, fontface = "bold"),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9),
      grid_width = unit(4, "mm"),
      grid_height = unit(4, "mm")
    ),
    width = unit(20, "cm"),
    height = unit(12, "cm"),
    border = TRUE,
    border_gp = gpar(col = "black", lwd = 1)
  )
  
  # ============================================================================
  # PANEL D: SUMMARY STATISTICS
  # ============================================================================

  # Calculate summary statistics
  summary_stats <- all_drugs %>%
    filter(!is.na(neglog10p), pvalue <= 0.05) %>%
    group_by(drug_label, drug, source, region_full) %>%
    summarize(
      n_significant = n_distinct(RBP_clean),
      mean_neglog = mean(neglog10p),
      median_neglog = median(neglog10p),
      .groups = "drop"
    )

  p_summary <- ggplot(summary_stats,
                      aes(x = region_full, y = n_significant, fill = drug_label)) +
    geom_bar(aes(alpha = source), stat = "identity",
             position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = drug_colors, name = "Drug") +
    scale_alpha_manual(values = c("Retained" = 1, "Excised" = 0.6),
                       name = "Splicing") +
    labs(
      title = "Number of Significant RBPs by Region",
      x = "Genomic Region",
      y = "Count of Significant RBPs"
    ) +
    theme_publication(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "right",
      legend.box = "vertical"
    ) +
    guides(
      fill = guide_legend(order = 1),
      alpha = guide_legend(order = 2)
    )

  # ============================================================================
  # COMBINE ALL PANELS WITH PROPER LAYOUT
  # ============================================================================

  # ComplexHeatmap grobs don't work well with patchwork
  # Solution: Use grid.arrange with proper viewport management

  library(gridExtra)
  library(grid)

  # Convert ggplots to grobs
  g_tracks <- ggplotGrob(p_tracks)
  g_boxplots <- ggplotGrob(p_boxplots)

  # Create the combined layout manually using grid
  final_plot_grob <- function() {
    # Create layout: 2 rows, 3 columns
    # Row 1: Panel A (tracks) - spanning all 3 columns
    # Row 2: Panel B (boxplots), Panel C (heatmap), empty space
    grid.newpage()

    # Add title and subtitle
    grid.text("Comprehensive RBP Binding Analysis Across Drug Treatments",
              x = 0.5, y = 0.98,
              gp = gpar(fontsize = 18, fontface = "bold"))
    grid.text("Differential binding patterns at splice sites for retained vs. excised introns",
              x = 0.5, y = 0.95,
              gp = gpar(fontsize = 13, col = "grey30"))

    # Define viewports - boxplot shorter, heatmap bigger
    pushViewport(viewport(layout = grid.layout(3, 2,
                                                heights = unit(c(0.3, 0.25, 0.45), "npc"),
                                                widths = unit(c(0.6, 0.4), "npc"))))

    # Panel A: Tracks (spanning all 2 columns)
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:2))
    grid.rect(gp = gpar(fill = "white", col = NA))
    grid.draw(g_tracks)
    grid.text("A", x = 0.02, y = 0.98, gp = gpar(fontsize = 16, fontface = "bold"))
    popViewport()

    # Panel B: Boxplots (spanning all 2 columns - shorter now)
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:2))
    grid.rect(gp = gpar(fill = "white", col = NA))
    grid.draw(g_boxplots)
    grid.text("B", x = 0.02, y = 0.98, gp = gpar(fontsize = 16, fontface = "bold"))
    popViewport()

    # Panel C: Heatmap (left bottom position, draw directly to avoid viewport conflicts)
    pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
    grid.rect(gp = gpar(fill = "white", col = NA))
    # Draw heatmap directly in this viewport
    pushViewport(viewport(width = 0.95, height = 0.95))
    draw(ht, newpage = FALSE, merge_legend = TRUE)
    upViewport()
    grid.text("C", x = 0.02, y = 0.98, gp = gpar(fontsize = 16, fontface = "bold"))
    popViewport()

    # Panel D position (bottom right) left empty (no drawing)

    # Add caption at bottom
    popViewport()  # Back to root
    grid.text("Significance threshold: p < 0.05 (dashed lines)",
              x = 0.98, y = 0.01,
              just = "right",
              gp = gpar(fontsize = 10, col = "grey50"))
  }  
  # ============================================================================
  # SAVE FINAL FIGURE
  # ============================================================================

  # Suppress LOESS warnings for cleaner output
  suppressWarnings({

    # Publication-ready PDF (cairo_pdf handles fonts better)
    cairo_pdf("rbp_comprehensive_analysis_introns.pdf", width = 18, height = 16,
              family = "sans", pointsize = 10)
    final_plot_grob()
    dev.off()

    # High-quality PNG for presentations
    png("rbp_comprehensive_analysis_introns.png",
        width = 18,
        height = 12,
        units = "in",
        res = 300,
        type = "cairo")
    final_plot_grob()
    dev.off()

    # Individual panel exports for flexibility
    cairo_pdf("panel_A_tracks_introns.pdf", width = 18, height = 5, family = "sans")
    print(p_tracks)
    dev.off()

    # Panel B: Boxplots (was C, now B)
    cairo_pdf("panel_B_boxplots_introns.pdf", width = 6, height = 5, family = "sans")
    print(p_boxplots)
    dev.off()

    # Panel C: Heatmap (was B, now C) - draw directly to avoid viewport conflicts
    cairo_pdf("panel_C_heatmap_introns.pdf", width = 10, height = 10, family = "sans")
    draw(ht, merge_legend = TRUE)
    dev.off()

  })

  cat("\n=== Figure Generation Complete ===\n")
  cat(paste("Matrix dimensions:", nrow(heatmap_matrix), "RBPs x",
            ncol(heatmap_matrix), "conditions\n"))
  cat("\nDrug order (left to right): Tubercidin, Pladienolide B, Spliceostatin A\n")
  cat("Region order (left to right): Upstream Exon, 5'SS, Intron, 3'SS, Downstream Exon\n")
  cat("\nMain figure files created:\n")
  cat("  ✓ rbp_comprehensive_analysis.pdf (18x12 in, publication-ready)\n")
  cat("  ✓ rbp_comprehensive_analysis.png (high-res preview)\n")
  cat("  ✓ rbp_comprehensive_analysis.tiff (journal submission)\n")
  cat("\nIndividual panel files:\n")
  cat("  ✓ panel_A_tracks.pdf (RBP enrichment tracks)\n")
  cat("  ✓ panel_B_boxplots.pdf (statistical comparisons with significance)\n")
  cat("  ✓ panel_C_heatmap.pdf (top 15 RBPs across conditions)\n")

  # ============================================================================
  # VENN DIAGRAM DATA: Common RBP_clean names between 3 drugs
  # ============================================================================

  cat("\n\n=== VENN DIAGRAM DATA FOR RBP_clean NAMES ===\n\n")

  # ============================================================================
  # 1. VENN DIAGRAM: Significant RBPs in ANY POSITION
  # ============================================================================

  cat(strrep("=", 70), "\n")
  cat("1. VENN DIAGRAM FOR RBPs SIGNIFICANT IN ANY POSITION (p <= 0.05)\n")
  cat(strrep("=", 70), "\n\n")

  # Extract unique RBP_clean names per drug (with p-value <= 0.05, any position)
  rbp_tub_any <- all_drugs %>%
    filter(drug == "tub", pvalue <= 0.05) %>%
    pull(RBP_clean) %>%
    unique() %>%
    sort()

  rbp_pladb_any <- all_drugs %>%
    filter(drug == "pladb", pvalue <= 0.05) %>%
    pull(RBP_clean) %>%
    unique() %>%
    sort()

  rbp_ssa_any <- all_drugs %>%
    filter(drug == "ssa", pvalue <= 0.05) %>%
    pull(RBP_clean) %>%
    unique() %>%
    sort()

  # Calculate intersections for ANY POSITION
  # Only in Tubercidin
  only_tub_any <- setdiff(rbp_tub_any, union(rbp_pladb_any, rbp_ssa_any))

  # Only in Pladienolide B
  only_pladb_any <- setdiff(rbp_pladb_any, union(rbp_tub_any, rbp_ssa_any))

  # Only in Spliceostatin A
  only_ssa_any <- setdiff(rbp_ssa_any, union(rbp_tub_any, rbp_pladb_any))

  # Tubercidin & Pladienolide B (not SSA)
  tub_pladb_any <- setdiff(intersect(rbp_tub_any, rbp_pladb_any), rbp_ssa_any)

  # Tubercidin & Spliceostatin A (not PLADB)
  tub_ssa_any <- setdiff(intersect(rbp_tub_any, rbp_ssa_any), rbp_pladb_any)

  # Pladienolide B & Spliceostatin A (not TUB)
  pladb_ssa_any <- setdiff(intersect(rbp_pladb_any, rbp_ssa_any), rbp_tub_any)

  # All three drugs
  all_three_any <- intersect(intersect(rbp_tub_any, rbp_pladb_any), rbp_ssa_any)

  # Print summary counts for ANY POSITION
  cat("SUMMARY COUNTS (for Venn diagram):\n")
  cat("===================================\n")
  cat(sprintf("Total unique RBPs in Tubercidin: %d\n", length(rbp_tub_any)))
  cat(sprintf("Total unique RBPs in Pladienolide B: %d\n", length(rbp_pladb_any)))
  cat(sprintf("Total unique RBPs in Spliceostatin A: %d\n", length(rbp_ssa_any)))
  cat("\n")
  cat(sprintf("Only Tubercidin: %d\n", length(only_tub_any)))
  cat(sprintf("Only Pladienolide B: %d\n", length(only_pladb_any)))
  cat(sprintf("Only Spliceostatin A: %d\n", length(only_ssa_any)))
  cat(sprintf("Tubercidin & Pladienolide B (not SSA): %d\n", length(tub_pladb_any)))
  cat(sprintf("Tubercidin & Spliceostatin A (not PLADB): %d\n", length(tub_ssa_any)))
  cat(sprintf("Pladienolide B & Spliceostatin A (not TUB): %d\n", length(pladb_ssa_any)))
  cat(sprintf("All three drugs: %d\n", length(all_three_any)))
  cat("\n")

  # Print detailed RBP lists for ANY POSITION
  cat("\n\nDETAILED RBP LISTS:\n")
  cat("===================\n\n")

  cat("Only Tubercidin (", length(only_tub_any), " RBPs):\n", sep="")
  cat(paste(only_tub_any, collapse=", "), "\n\n")

  cat("Only Pladienolide B (", length(only_pladb_any), " RBPs):\n", sep="")
  cat(paste(only_pladb_any, collapse=", "), "\n\n")

  cat("Only Spliceostatin A (", length(only_ssa_any), " RBPs):\n", sep="")
  cat(paste(only_ssa_any, collapse=", "), "\n\n")

  cat("Tubercidin & Pladienolide B only (", length(tub_pladb_any), " RBPs):\n", sep="")
  cat(paste(tub_pladb_any, collapse=", "), "\n\n")

  cat("Tubercidin & Spliceostatin A only (", length(tub_ssa_any), " RBPs):\n", sep="")
  cat(paste(tub_ssa_any, collapse=", "), "\n\n")

  cat("Pladienolide B & Spliceostatin A only (", length(pladb_ssa_any), " RBPs):\n", sep="")
  cat(paste(pladb_ssa_any, collapse=", "), "\n\n")

  cat("All three drugs (", length(all_three_any), " RBPs):\n", sep="")
  cat(paste(all_three_any, collapse=", "), "\n\n")

  # ============================================================================
  # 2. VENN DIAGRAM: Significant RBPs in INTRON POSITION ONLY (R3)
  # ============================================================================

  cat("\n\n")
  cat(strrep("=", 70), "\n")
  cat("2. VENN DIAGRAM FOR RBPs SIGNIFICANT IN INTRON POSITION (R3) ONLY\n")
  cat(strrep("=", 70), "\n\n")

  # Extract unique RBP_clean names per drug (with p-value <= 0.05, INTRON position only)
  rbp_tub_intron <- all_drugs %>%
    filter(drug == "tub", pvalue <= 0.05, position == "R3") %>%
    pull(RBP_clean) %>%
    unique() %>%
    sort()

  rbp_pladb_intron <- all_drugs %>%
    filter(drug == "pladb", pvalue <= 0.05, position == "R3") %>%
    pull(RBP_clean) %>%
    unique() %>%
    sort()

  rbp_ssa_intron <- all_drugs %>%
    filter(drug == "ssa", pvalue <= 0.05, position == "R3") %>%
    pull(RBP_clean) %>%
    unique() %>%
    sort()

  # Calculate intersections for INTRON POSITION
  only_tub_intron <- setdiff(rbp_tub_intron, union(rbp_pladb_intron, rbp_ssa_intron))
  only_pladb_intron <- setdiff(rbp_pladb_intron, union(rbp_tub_intron, rbp_ssa_intron))
  only_ssa_intron <- setdiff(rbp_ssa_intron, union(rbp_tub_intron, rbp_pladb_intron))
  tub_pladb_intron <- setdiff(intersect(rbp_tub_intron, rbp_pladb_intron), rbp_ssa_intron)
  tub_ssa_intron <- setdiff(intersect(rbp_tub_intron, rbp_ssa_intron), rbp_pladb_intron)
  pladb_ssa_intron <- setdiff(intersect(rbp_pladb_intron, rbp_ssa_intron), rbp_tub_intron)
  all_three_intron <- intersect(intersect(rbp_tub_intron, rbp_pladb_intron), rbp_ssa_intron)

  # Print summary counts for INTRON POSITION
  cat("SUMMARY COUNTS (for Venn diagram - INTRON only):\n")
  cat("=================================================\n")
  cat(sprintf("Total unique RBPs in Tubercidin (Intron): %d\n", length(rbp_tub_intron)))
  cat(sprintf("Total unique RBPs in Pladienolide B (Intron): %d\n", length(rbp_pladb_intron)))
  cat(sprintf("Total unique RBPs in Spliceostatin A (Intron): %d\n", length(rbp_ssa_intron)))
  cat("\n")
  cat(sprintf("Only Tubercidin: %d\n", length(only_tub_intron)))
  cat(sprintf("Only Pladienolide B: %d\n", length(only_pladb_intron)))
  cat(sprintf("Only Spliceostatin A: %d\n", length(only_ssa_intron)))
  cat(sprintf("Tubercidin & Pladienolide B (not SSA): %d\n", length(tub_pladb_intron)))
  cat(sprintf("Tubercidin & Spliceostatin A (not PLADB): %d\n", length(tub_ssa_intron)))
  cat(sprintf("Pladienolide B & Spliceostatin A (not TUB): %d\n", length(pladb_ssa_intron)))
  cat(sprintf("All three drugs: %d\n", length(all_three_intron)))
  cat("\n")

  # Print detailed RBP lists for INTRON POSITION
  cat("\n\nDETAILED RBP LISTS (INTRON only):\n")
  cat("==================================\n\n")

  cat("Only Tubercidin (", length(only_tub_intron), " RBPs):\n", sep="")
  cat(paste(only_tub_intron, collapse=", "), "\n\n")

  cat("Only Pladienolide B (", length(only_pladb_intron), " RBPs):\n", sep="")
  cat(paste(only_pladb_intron, collapse=", "), "\n\n")

  cat("Only Spliceostatin A (", length(only_ssa_intron), " RBPs):\n", sep="")
  cat(paste(only_ssa_intron, collapse=", "), "\n\n")

  cat("Tubercidin & Pladienolide B only (", length(tub_pladb_intron), " RBPs):\n", sep="")
  cat(paste(tub_pladb_intron, collapse=", "), "\n\n")

  cat("Tubercidin & Spliceostatin A only (", length(tub_ssa_intron), " RBPs):\n", sep="")
  cat(paste(tub_ssa_intron, collapse=", "), "\n\n")

  cat("Pladienolide B & Spliceostatin A only (", length(pladb_ssa_intron), " RBPs):\n", sep="")
  cat(paste(pladb_ssa_intron, collapse=", "), "\n\n")

  cat("All three drugs (", length(all_three_intron), " RBPs):\n", sep="")
  cat(paste(all_three_intron, collapse=", "), "\n\n")

  # Save to a text file for easy reference
  venn_output <- c(
    "=== VENN DIAGRAM DATA FOR RBP_clean NAMES ===",
    "",
    strrep("=", 70),
    "1. VENN DIAGRAM FOR RBPs SIGNIFICANT IN ANY POSITION (p <= 0.05)",
    strrep("=", 70),
    "",
    "SUMMARY COUNTS:",
    "===============",
    sprintf("Total unique RBPs in Tubercidin: %d", length(rbp_tub_any)),
    sprintf("Total unique RBPs in Pladienolide B: %d", length(rbp_pladb_any)),
    sprintf("Total unique RBPs in Spliceostatin A: %d", length(rbp_ssa_any)),
    "",
    sprintf("Only Tubercidin: %d", length(only_tub_any)),
    sprintf("Only Pladienolide B: %d", length(only_pladb_any)),
    sprintf("Only Spliceostatin A: %d", length(only_ssa_any)),
    sprintf("Tubercidin & Pladienolide B (not SSA): %d", length(tub_pladb_any)),
    sprintf("Tubercidin & Spliceostatin A (not PLADB): %d", length(tub_ssa_any)),
    sprintf("Pladienolide B & Spliceostatin A (not TUB): %d", length(pladb_ssa_any)),
    sprintf("All three drugs: %d", length(all_three_any)),
    "",
    "",
    "DETAILED RBP LISTS:",
    "===================",
    "",
    sprintf("Only Tubercidin (%d RBPs):", length(only_tub_any)),
    paste(only_tub_any, collapse=", "),
    "",
    sprintf("Only Pladienolide B (%d RBPs):", length(only_pladb_any)),
    paste(only_pladb_any, collapse=", "),
    "",
    sprintf("Only Spliceostatin A (%d RBPs):", length(only_ssa_any)),
    paste(only_ssa_any, collapse=", "),
    "",
    sprintf("Tubercidin & Pladienolide B only (%d RBPs):", length(tub_pladb_any)),
    paste(tub_pladb_any, collapse=", "),
    "",
    sprintf("Tubercidin & Spliceostatin A only (%d RBPs):", length(tub_ssa_any)),
    paste(tub_ssa_any, collapse=", "),
    "",
    sprintf("Pladienolide B & Spliceostatin A only (%d RBPs):", length(pladb_ssa_any)),
    paste(pladb_ssa_any, collapse=", "),
    "",
    sprintf("All three drugs (%d RBPs):", length(all_three_any)),
    paste(all_three_any, collapse=", "),
    "",
    "",
    "",
    strrep("=", 70),
    "2. VENN DIAGRAM FOR RBPs SIGNIFICANT IN INTRON POSITION (R3) ONLY",
    strrep("=", 70),
    "",
    "SUMMARY COUNTS (INTRON only):",
    "==============================",
    sprintf("Total unique RBPs in Tubercidin (Intron): %d", length(rbp_tub_intron)),
    sprintf("Total unique RBPs in Pladienolide B (Intron): %d", length(rbp_pladb_intron)),
    sprintf("Total unique RBPs in Spliceostatin A (Intron): %d", length(rbp_ssa_intron)),
    "",
    sprintf("Only Tubercidin: %d", length(only_tub_intron)),
    sprintf("Only Pladienolide B: %d", length(only_pladb_intron)),
    sprintf("Only Spliceostatin A: %d", length(only_ssa_intron)),
    sprintf("Tubercidin & Pladienolide B (not SSA): %d", length(tub_pladb_intron)),
    sprintf("Tubercidin & Spliceostatin A (not PLADB): %d", length(tub_ssa_intron)),
    sprintf("Pladienolide B & Spliceostatin A (not TUB): %d", length(pladb_ssa_intron)),
    sprintf("All three drugs: %d", length(all_three_intron)),
    "",
    "",
    "DETAILED RBP LISTS (INTRON only):",
    "==================================",
    "",
    sprintf("Only Tubercidin (%d RBPs):", length(only_tub_intron)),
    paste(only_tub_intron, collapse=", "),
    "",
    sprintf("Only Pladienolide B (%d RBPs):", length(only_pladb_intron)),
    paste(only_pladb_intron, collapse=", "),
    "",
    sprintf("Only Spliceostatin A (%d RBPs):", length(only_ssa_intron)),
    paste(only_ssa_intron, collapse=", "),
    "",
    sprintf("Tubercidin & Pladienolide B only (%d RBPs):", length(tub_pladb_intron)),
    paste(tub_pladb_intron, collapse=", "),
    "",
    sprintf("Tubercidin & Spliceostatin A only (%d RBPs):", length(tub_ssa_intron)),
    paste(tub_ssa_intron, collapse=", "),
    "",
    sprintf("Pladienolide B & Spliceostatin A only (%d RBPs):", length(pladb_ssa_intron)),
    paste(pladb_ssa_intron, collapse=", "),
    "",
    sprintf("All three drugs (%d RBPs):", length(all_three_intron)),
    paste(all_three_intron, collapse=", ")
  )

  writeLines(venn_output, "rbp_venn_diagram_data.txt")
  cat("\n✓ Venn diagram data saved to: rbp_venn_diagram_data.txt\n")
  cat("  (includes BOTH: any position AND intron-specific analysis)\n\n")