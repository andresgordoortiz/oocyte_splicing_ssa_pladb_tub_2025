library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(purrr)
library(showtext)
library(readr)
library(ComplexHeatmap)
library(circlize)

# Font setup for Fedora Linux - Liberation Sans is the Arial equivalent
font_add(family = "LiberationSans", 
         regular = "/usr/share/fonts/liberation-sans/LiberationSans-Regular.ttf",
         bold = "/usr/share/fonts/liberation-sans/LiberationSans-Bold.ttf")
showtext::showtext_opts(dpi = 600)
showtext_auto()

# Import betAS output tables
tub_fdr_df <- read_csv("tub_fdr.csv")[,-1]
pladb_fdr_df <- read_csv("pladb_fdr.csv")[,-1]
ssa_fdr_df <- read_csv("ssa_fdr.csv")[,-1]

# Extract significant events
differential_tub <- na.omit(tub_fdr_df[tub_fdr_df$FDR <= 0.05 & abs(tub_fdr_df$deltapsi) >= 0.1,])
differential_pladb <- na.omit(pladb_fdr_df[pladb_fdr_df$FDR <= 0.05 & abs(pladb_fdr_df$deltapsi) >= 0.1,])
differential_ssa <- na.omit(ssa_fdr_df[ssa_fdr_df$FDR <= 0.05 & abs(ssa_fdr_df$deltapsi) >= 0.1,])

# Create list for upset plot
gene_lists <- list(
  "Tubercidin" = unique(differential_tub$GENE),
  "Pladienolide B" = unique(differential_pladb$GENE),
  "Spliceostatin A" = unique(differential_ssa$GENE)
)

# ---- Drug colors ----
ssa_color <- "#706993"    
pladb_color <- "#70A0AF"  
tub_color  <- "#A0C1B9"   

set_colors <- c(
  "Tubercidin" = tub_color,
  "Pladienolide B" = pladb_color,
  "Spliceostatin A" = ssa_color
)

# Create combination matrix
m <- make_comb_mat(gene_lists)

# Get combination sizes for coloring
cs <- comb_size(m)

# Create color function based on which sets are in each combination
comb_colors <- sapply(comb_name(m), function(comb_name) {
  sets_in_comb <- extract_comb(m, comb_name)
  sets_present <- names(sets_in_comb)[sets_in_comb]
  
  if (length(sets_present) == 1) {
    # Single set - use its color
    return(set_colors[sets_present])
  } else if (length(sets_present) == 2) {
    # Two sets - use a darker grey
    return("grey40")
  } else {
    # All three sets - use black
    return("grey20")
  }
})

# Create the upset plot with proper styling
library(ComplexHeatmap)
library(circlize)
library(grid)

# ---------- USER ADJUSTABLE PARAMETERS ----------
# How tall the other barplot is (the one you want to match) and its y max
other_barplot_height_mm <- 80    # mm - height of the other barplot (change if needed)
other_barplot_ylim_max <- 500    # numeric max of the other barplot y axis (you said 500)

# UpSet matrix tuning
mm_per_row <- 7                # mm per row in the matrix body (smaller -> tighter)
# ------------------------------------------------

# compute mm per count from the other barplot so axes are physically comparable
mm_per_count <- other_barplot_height_mm / other_barplot_ylim_max

# set the top annotation target y-range (make it identical numerically)
top_ann_ylim_max <- other_barplot_ylim_max
top_ann_height_mm <- top_ann_ylim_max * mm_per_count # height in mm to match other barplot scale

# matrix / heatmap heights
nrows <- nrow(m)
body_height <- unit(nrows * mm_per_row, "mm")
heatmap_total_height <- unit(top_ann_height_mm, "mm") + body_height

# (optional) comb colors and set colors you already had
# comb_colors <- c(...) # keep your comb_colors as before
# set_colors  <- c(...) # keep your set_colors as before

upset_plot <- UpSet(
  m,
  set_order = c("Tubercidin", "Pladienolide B", "Spliceostatin A"),
  comb_order = order(-comb_size(m)),
  
  # Top bar chart annotation - forcing identical y-axis range and physical height
  top_annotation = upset_top_annotation(
    m,
    add_numbers = TRUE,
    numbers_gp = gpar(fontsize = 12, fontfamily = "sans", col = "black"),
    numbers_rot = 0,
    height = unit(top_ann_height_mm, "mm"),        # physical height matching other barplot
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 9, fontfamily = "sans", fontface = "plain"),
    annotation_name_rot = 90,
    ylim = c(0, top_ann_ylim_max),                 # same numerical range as other barplot
    bar_width = 0.6,
    gp = gpar(fill = comb_colors, col = "#c0c0c0", lwd = 0.5),
    axis_param = list(
      gp = gpar(fontsize = 8, fontfamily = "sans"),
      at = pretty(c(0, top_ann_ylim_max), n = 5),
      side = "left",
      labels_rot = 0
    )
  ),
  
  # Left bar chart annotation (leave as you had it)
  left_annotation = upset_left_annotation(
    m,
    add_numbers = TRUE,
    numbers_gp = gpar(fontsize = 9, fontfamily = "sans", col = "black"),
    width = unit(35, "mm"),
    annotation_name_gp = gpar(fontsize = 9, fontfamily = "sans", fontface = "plain"),
    annotation_name_rot = 0,
    xlim = c(0, max(set_size(m)) * 1.15),
    bar_width = 0.4,
    gp = gpar(fill = set_colors, col = "black", lwd = 0.1),
    axis_param = list(
      gp = gpar(fontsize = 8, fontfamily = "sans"),
      at = pretty(c(0, max(set_size(m))), n = 4),
      side = "bottom",
      labels_rot = 0
    )
  ),
  
  right_annotation = NULL,
  
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 9, fontfamily = "sans", col = "black"),
  
  comb_col = "black",
  bg_col = c("grey95", "white"),
  bg_pt_col = "grey80",
  pt_size = unit(2.4, "mm"),
  lwd = 1.2,
  
  row_names_side = "left",
  column_title = NULL,
  
  # IMPORTANT: set only heatmap_height (absolute); ComplexHeatmap computes internals
  heatmap_height = heatmap_total_height
)

# Draw it
draw(upset_plot)



# Double-column version with better proportions
pdf("upset_genes_figure.pdf", width = 180/25.4, height = 130/25.4)
draw(upset_plot)
dev.off()

png("upset_genes_figure.png", width = 180, height = 130, units = "mm", res = 600)
draw(upset_plot)
dev.off()

# Print to display
draw(upset_plot)