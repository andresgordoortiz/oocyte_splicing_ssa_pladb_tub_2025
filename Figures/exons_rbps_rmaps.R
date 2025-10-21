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
library(gridExtra)

# ============================================================================
# DEFINE CONSISTENT THEME AND COLORS
# ============================================================================

# Drug colors (order: TUB, PLADB, SSA)
drug_colors <- c(
  "Tubercidin" = "#A0C1B9",
  "Pladienolide B" = "#70A0AF",
  "Spliceostatin A" = "#706993"
)

# Harmonious splicing colors
splicing_colors <- c(
  "Included" = "#D4A574",
  "Skipped" = "#B85450"
)

# Region colors (8-region exon map)
region_colors <- c(
  "Upstream Exon (last 50bp)" = "#E8E3F5",
  "Upstream Intron (first 50bp)" = "#C8C3D8",
  "Upstream Intron (last 50bp)" = "#A8A3BB",
  "Target Exon (first 50bp)" = "#88839E",
  "Target Exon (last 50bp)" = "#686381",
  "Downstream Intron (first 50bp)" = "#5E5A6A",
  "Downstream Intron (last 50bp)" = "#4B4856",
  "Downstream Exon (first 50bp)" = "#372F3F"
)

theme_publication <- function(base_size = 10) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      axis.line = element_line(linewidth = 0.6, colour = "black"),
      axis.ticks = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      axis.title = element_text(face = "bold", size = rel(1.0)),
      axis.text = element_text(size = rel(0.9), colour = "black"),
      legend.position = "right",
      legend.background = element_blank(),
      legend.title = element_text(face = "bold", size = rel(0.95)),
      legend.text = element_text(size = rel(0.85)),
      legend.key.size = unit(0.4, "cm"),
      legend.spacing.y = unit(0.1, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = rel(1.0), margin = margin(5,5,5,5)),
      plot.title = element_text(face = "bold", size = rel(1.1), hjust = 0),
      plot.subtitle = element_text(size = rel(0.9), hjust = 0, color = "grey30"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# ============================================================================
# DATA LOADING AND PROCESSING FUNCTION
# ============================================================================

process_drug_data_exon <- function(drug_name) {
  # Construct file paths
  up_file <- paste0("rmaps_out/", drug_name, "_ex_rmaps/pVal.up.vs.bg.RNAmap.txt")
  down_file <- paste0("rmaps_out/", drug_name, "_ex_rmaps/pVal.dn.vs.bg.RNAmap.txt")
  
  # Check if files exist
  if (!file.exists(up_file)) {
    stop(paste("File not found:", up_file))
  }
  if (!file.exists(down_file)) {
    stop(paste("File not found:", down_file))
  }
  
  up_data <- read.delim(up_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  down_data <- read.delim(down_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Standardize column names if needed (expecting RBP, R1..R8)
  expected_cols <- c("RBP", paste0("R", 1:8))
  if (!all(expected_cols %in% colnames(up_data))) {
    colnames(up_data)[1:min(length(colnames(up_data)), length(expected_cols))] <- 
      expected_cols[1:min(length(colnames(up_data)), length(expected_cols))]
  }
  if (!all(expected_cols %in% colnames(down_data))) {
    colnames(down_data)[1:min(length(colnames(down_data)), length(expected_cols))] <- 
      expected_cols[1:min(length(colnames(down_data)), length(expected_cols))]
  }
  
  # Filter rows with p <= 0.05 in any position
  filter_by_any_R <- function(df, source_name) {
    r_cols <- grep("^R\\d+$", names(df), value = TRUE)
    df %>%
      filter(if_any(all_of(r_cols), ~ . <= 0.05)) %>%
      mutate(source = source_name, drug = drug_name)
  }
  
  up_filtered <- filter_by_any_R(up_data, "Included")
  down_filtered <- filter_by_any_R(down_data, "Skipped")
  
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
      RBP_clean = sub("\\..*", "", RBP)
    )
  
  return(combined_long)
}

# ============================================================================
# LOAD AND COMBINE ALL DRUG DATA
# ============================================================================

cat("Loading exon data for all drugs...\n")
all_drugs <- map_dfr(c("pladb", "ssa", "tub"), process_drug_data_exon)
cat("Data loaded successfully.\n\n")

# Define 8-region exon-based map
region_map <- tibble(
  position = paste0("R", 1:8),
  region_short = c("Upstream\nExon", "Upstream\nIntron 5'", "Upstream\nIntron 3'",
                   "Target\nExon 5'", "Target\nExon 3'", "Downstream\nIntron 5'",
                   "Downstream\nIntron 3'", "Downstream\nExon"),
  region_full = c("Upstream Exon (last 50bp)",
                  "Upstream Intron (first 50bp)",
                  "Upstream Intron (last 50bp)",
                  "Target Exon (first 50bp)",
                  "Target Exon (last 50bp)",
                  "Downstream Intron (first 50bp)",
                  "Downstream Intron (last 50bp)",
                  "Downstream Exon (first 50bp)"),
  length_bp = c(50, 250, 250, 50, 50, 250, 250, 50)
) %>%
  mutate(
    start = cumsum(lag(length_bp, default = 0)),
    mid = start + length_bp / 2,
    end = start + length_bp
  )

# Join region info and set factor levels
all_drugs <- all_drugs %>%
  left_join(region_map %>% select(position, mid, region_short, region_full), by = "position") %>%
  mutate(
    position = factor(position, levels = paste0("R", 1:8)),
    region_full = factor(region_full, levels = region_map$region_full),
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
    title = "RNA-Binding Protein Enrichment Across Cassette Exon Regions",
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
  mutate(source = factor(source, levels = c("Included", "Skipped")))

# Wilcoxon tests
stat_tests <- boxplot_data %>%
  group_by(region_full, drug_label) %>%
  summarize(
    p_value = tryCatch({
      wilcox.test(neglog10p[source == "Included"], 
                  neglog10p[source == "Skipped"])$p.value
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
  facet_wrap(~ region_full, ncol = 4, scales = "fixed") +
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
    strip.text = element_text(face = "bold", size = rel(1.0), margin = margin(b = 6))
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

# Get top 15 RBPs overall
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

# FIXED: Create columns in the CORRECT order from the start
# Order: Position -> Drug -> Source (matching annotation order)
positions <- paste0("R", 1:8)
drugs <- c("Tubercidin", "Pladienolide B", "Spliceostatin A")
sources <- c("Included", "Skipped")

# Create ordered column names
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
n_positions <- length(positions)
n_drugs <- length(drugs)
n_sources <- length(sources)

position_ann <- rep(region_map$region_full, each = n_drugs * n_sources)
drug_ann <- rep(rep(drugs, each = n_sources), n_positions)
source_ann <- rep(rep(sources, n_drugs), n_positions)

position_ann_factor <- factor(position_ann, levels = region_map$region_full)
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
# COMBINE ALL PANELS
# ============================================================================

cat("Combining panels into final figure...\n")

g_tracks <- ggplotGrob(p_tracks)
g_boxplots <- ggplotGrob(p_boxplots)

final_plot_grob <- function() {
  grid.newpage()
  
  # Title and subtitle
  grid.text("Comprehensive RBP Binding Analysis Across Cassette Exons",
            x = 0.5, y = 0.98,
            gp = gpar(fontsize = 18, fontface = "bold"))
  grid.text("Differential binding patterns at exon/intron regions for included vs. skipped exons",
            x = 0.5, y = 0.95,
            gp = gpar(fontsize = 13, col = "grey30"))
  
  # Layout: 3 rows, 2 columns
  pushViewport(viewport(layout = grid.layout(3, 2,
                                             heights = unit(c(0.3, 0.25, 0.45), "npc"),
                                             widths = unit(c(0.6, 0.4), "npc"))))
  
  # Panel A: Tracks
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:2))
  grid.rect(gp = gpar(fill = "white", col = NA))
  grid.draw(g_tracks)
  grid.text("A", x = 0.02, y = 0.98, gp = gpar(fontsize = 16, fontface = "bold"))
  popViewport()
  
  # Panel B: Boxplots
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:2))
  grid.rect(gp = gpar(fill = "white", col = NA))
  grid.draw(g_boxplots)
  grid.text("B", x = 0.02, y = 0.98, gp = gpar(fontsize = 16, fontface = "bold"))
  popViewport()
  
  # Panel C: Heatmap
  pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
  grid.rect(gp = gpar(fill = "white", col = NA))
  pushViewport(viewport(width = 0.95, height = 0.95))
  draw(ht, newpage = FALSE, merge_legend = TRUE)
  upViewport()
  grid.text("C", x = 0.02, y = 0.98, gp = gpar(fontsize = 16, fontface = "bold"))
  popViewport()
  
  # Caption
  popViewport()
  grid.text("Significance threshold: p < 0.05 (dashed lines)",
            x = 0.98, y = 0.01, just = "right",
            gp = gpar(fontsize = 10, col = "grey50"))
}

# ============================================================================
# SAVE FIGURES
# ============================================================================

cat("\nSaving figures...\n")

suppressWarnings({
  # Main figure
  cairo_pdf("rbp_comprehensive_analysis_exon.pdf", width = 18, height = 16,
            family = "sans", pointsize = 10)
  final_plot_grob()
  dev.off()
  
  png("rbp_comprehensive_analysis_exon.png",
      width = 18, height = 16, units = "in", res = 300, type = "cairo")
  final_plot_grob()
  dev.off()
  
  tiff("rbp_comprehensive_analysis_exon.tiff",
       width = 18, height = 16, units = "in", res = 300, 
       compression = "lzw", type = "cairo")
  final_plot_grob()
  dev.off()
  
  # Individual panels
  cairo_pdf("panel_A_tracks_exon.pdf", width = 18, height = 5, family = "sans")
  print(p_tracks)
  dev.off()
  
  cairo_pdf("panel_B_boxplots_exon.pdf", width = 18, height = 8, family = "sans")
  print(p_boxplots)
  dev.off()
  
  cairo_pdf("panel_C_heatmap_exon.pdf", width = 12, height = 10, family = "sans")
  draw(ht, merge_legend = TRUE)
  dev.off()
})

cat("\n✓ Main figures saved\n")

# ============================================================================
# VENN DIAGRAM DATA
# ============================================================================

cat("\n=== GENERATING VENN DIAGRAM DATA (EXON MODE) ===\n")

# ANY POSITION analysis
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

# Calculate intersections
only_tub_any <- setdiff(rbp_tub_any, union(rbp_pladb_any, rbp_ssa_any))
only_pladb_any <- setdiff(rbp_pladb_any, union(rbp_tub_any, rbp_ssa_any))
only_ssa_any <- setdiff(rbp_ssa_any, union(rbp_tub_any, rbp_pladb_any))
tub_pladb_any <- setdiff(intersect(rbp_tub_any, rbp_pladb_any), rbp_ssa_any)
tub_ssa_any <- setdiff(intersect(rbp_tub_any, rbp_ssa_any), rbp_pladb_any)
pladb_ssa_any <- setdiff(intersect(rbp_pladb_any, rbp_ssa_any), rbp_tub_any)
all_three_any <- intersect(intersect(rbp_tub_any, rbp_pladb_any), rbp_ssa_any)

# TARGET EXON ONLY analysis (R4 & R5)
rbp_tub_target_exon <- all_drugs %>%
  filter(drug == "tub", pvalue <= 0.05, position %in% c("R4", "R5")) %>%
  pull(RBP_clean) %>%
  unique() %>%
  sort()

rbp_pladb_target_exon <- all_drugs %>%
  filter(drug == "pladb", pvalue <= 0.05, position %in% c("R4", "R5")) %>%
  pull(RBP_clean) %>%
  unique() %>%
  sort()

rbp_ssa_target_exon <- all_drugs %>%
  filter(drug == "ssa", pvalue <= 0.05, position %in% c("R4", "R5")) %>%
  pull(RBP_clean) %>%
  unique() %>%
  sort()

only_tub_target_exon <- setdiff(rbp_tub_target_exon, union(rbp_pladb_target_exon, rbp_ssa_target_exon))
only_pladb_target_exon <- setdiff(rbp_pladb_target_exon, union(rbp_tub_target_exon, rbp_ssa_target_exon))
only_ssa_target_exon <- setdiff(rbp_ssa_target_exon, union(rbp_tub_target_exon, rbp_pladb_target_exon))
tub_pladb_target_exon <- setdiff(intersect(rbp_tub_target_exon, rbp_pladb_target_exon), rbp_ssa_target_exon)
tub_ssa_target_exon <- setdiff(intersect(rbp_tub_target_exon, rbp_ssa_target_exon), rbp_pladb_target_exon)
pladb_ssa_target_exon <- setdiff(intersect(rbp_pladb_target_exon, rbp_ssa_target_exon), rbp_tub_target_exon)
all_three_target_exon <- intersect(intersect(rbp_tub_target_exon, rbp_pladb_target_exon), rbp_ssa_target_exon)

# Create output text
venn_output <- c(
  "=================================================================",
  "VENN DIAGRAM DATA FOR RBP ANALYSIS (EXON MODE)",
  "=================================================================",
  "",
  "--- SECTION 1: RBPs SIGNIFICANT IN ANY POSITION ---",
  "",
  sprintf("Total Tubercidin: %d", length(rbp_tub_any)),
  sprintf("Total Pladienolide B: %d", length(rbp_pladb_any)),
  sprintf("Total Spliceostatin A: %d", length(rbp_ssa_any)),
  "",
  sprintf("Only Tubercidin: %d", length(only_tub_any)),
  paste(only_tub_any, collapse = ", "),
  "",
  sprintf("Only Pladienolide B: %d", length(only_pladb_any)),
  paste(only_pladb_any, collapse = ", "),
  "",
  sprintf("Only Spliceostatin A: %d", length(only_ssa_any)),
  paste(only_ssa_any, collapse = ", "),
  "",
  sprintf("Tubercidin & Pladienolide B (not SSA): %d", length(tub_pladb_any)),
  paste(tub_pladb_any, collapse = ", "),
  "",
  sprintf("Tubercidin & Spliceostatin A (not PLADB): %d", length(tub_ssa_any)),
  paste(tub_ssa_any, collapse = ", "),
  "",
  sprintf("Pladienolide B & Spliceostatin A (not TUB): %d", length(pladb_ssa_any)),
  paste(pladb_ssa_any, collapse = ", "),
  "",
  sprintf("All three drugs: %d", length(all_three_any)),
  paste(all_three_any, collapse = ", "),
  "",
  "",
  "--- SECTION 2: RBPs SIGNIFICANT IN TARGET EXON (R4/R5) ONLY ---",
  "",
  sprintf("Total Tubercidin (Target Exon): %d", length(rbp_tub_target_exon)),
  sprintf("Total Pladienolide B (Target Exon): %d", length(rbp_pladb_target_exon)),
  sprintf("Total Spliceostatin A (Target Exon): %d", length(rbp_ssa_target_exon)),
  "",
  sprintf("Only Tubercidin: %d", length(only_tub_target_exon)),
  paste(only_tub_target_exon, collapse = ", "),
  "",
  sprintf("Only Pladienolide B: %d", length(only_pladb_target_exon)),
  paste(only_pladb_target_exon, collapse = ", "),
  "",
  sprintf("Only Spliceostatin A: %d", length(only_ssa_target_exon)),
  paste(only_ssa_target_exon, collapse = ", "),
  "",
  sprintf("Tubercidin & Pladienolide B (not SSA): %d", length(tub_pladb_target_exon)),
  paste(tub_pladb_target_exon, collapse = ", "),
  "",
  sprintf("Tubercidin & Spliceostatin A (not PLADB): %d", length(tub_ssa_target_exon)),
  paste(tub_ssa_target_exon, collapse = ", "),
  "",
  sprintf("Pladienolide B & Spliceostatin A (not TUB): %d", length(pladb_ssa_target_exon)),
  paste(pladb_ssa_target_exon, collapse = ", "),
  "",
  sprintf("All three drugs: %d", length(all_three_target_exon)),
  paste(all_three_target_exon, collapse = ", ")
)

writeLines(venn_output, "rbp_venn_diagram_data_exon.txt")

cat("\n✓ Venn diagram data saved to: rbp_venn_diagram_data_exon.txt\n")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\n=== ANALYSIS SUMMARY ===\n")
cat(sprintf("Total RBPs analyzed: %d\n", length(unique(all_drugs$RBP_clean))))
cat(sprintf("Total significant RBPs (any drug, any position): %d\n", 
            length(unique(all_drugs$RBP_clean[all_drugs$pvalue <= 0.05]))))
cat(sprintf("RBPs in heatmap: %d\n", length(top_rbps_heatmap)))
cat(sprintf("Heatmap dimensions: %d RBPs x %d conditions\n", 
            nrow(heatmap_matrix), ncol(heatmap_matrix)))

cat("\n=== FILES CREATED ===\n")
cat("Main outputs:\n")
cat("  ✓ rbp_comprehensive_analysis_exon.pdf (18x16 in)\n")
cat("  ✓ rbp_comprehensive_analysis_exon.png (high-res)\n")
cat("  ✓ rbp_comprehensive_analysis_exon.tiff (journal submission)\n")
cat("\nIndividual panels:\n")
cat("  ✓ panel_A_tracks_exon.pdf\n")
cat("  ✓ panel_B_boxplots_exon.pdf\n")
cat("  ✓ panel_C_heatmap_exon.pdf\n")
cat("\nData outputs:\n")
cat("  ✓ rbp_venn_diagram_data_exon.txt\n")

cat("\n=== ANALYSIS COMPLETE ===\n")