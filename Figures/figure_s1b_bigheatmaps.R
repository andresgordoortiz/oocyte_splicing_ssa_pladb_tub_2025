library(betAS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(showtext)

splicing_data <- getDataset(pathTables = "new_data_INCLUSION_LEVELS_FULL-mm10.tab", tool = "vast-tools")

pladb_data <- splicing_data[,1:18]
ssa_tub_data <- splicing_data[, c(1:6, 19:36)]

pladb_events <- filterEvents(getEvents(pladb_data, tool = "vast-tools"), N = 10)
ssa_tub_events <- filterEvents(getEvents(ssa_tub_data, tool = "vast-tools"), N = 10)

ssa_events <- list()
ssa_events$PSI <- ssa_tub_events$PSI[,c(1:9,13:15)]

tub_events <- list()
tub_events$PSI <- ssa_tub_events$PSI[,c(1:12)]

# Import betAS output tables
tub_fdr_df <- read_csv("tub_fdr.csv")[,-1]
pladb_fdr_df <- read_csv("pladb_fdr.csv")[,-1]
ssa_fdr_df <- read_csv("ssa_fdr.csv")[,-1]

# Extract significant events
differential_tub <- na.omit(tub_fdr_df[tub_fdr_df$FDR <= 0.05 & abs(tub_fdr_df$deltapsi) >= 0.1,])
differential_pladb <- na.omit(pladb_fdr_df[pladb_fdr_df$FDR <= 0.05 & abs(pladb_fdr_df$deltapsi) >= 0.1,])
differential_ssa <- na.omit(ssa_fdr_df[ssa_fdr_df$FDR <= 0.05 & abs(ssa_fdr_df$deltapsi) >= 0.1,])

differential_pladb_exons <- filter(differential_pladb, grepl("EX", differential_pladb$EVENT))
differential_tub_exons <- filter(differential_tub, grepl("EX", differential_tub$EVENT))
differential_ssa_exons <- filter(differential_ssa, grepl("EX", differential_ssa$EVENT))

differential_pladb_introns <- filter(differential_pladb, grepl("INT", differential_pladb$EVENT))
differential_tub_introns <- filter(differential_tub, grepl("INT", differential_tub$EVENT))
differential_ssa_introns <- filter(differential_ssa, grepl("INT", differential_ssa$EVENT))

set.seed(44)

# ---- Font setup ----
font_add(family = "Arial", regular = "/usr/share/fonts/liberation/LiberationSans-Regular.ttf")
showtext::showtext_opts(dpi = 600)
showtext_auto()

# ---- Helper: robust row-wise MAD scaling ----
robust_scale_rows <- function(mat) {
  t(apply(mat, 1, function(x) {
    med <- median(x, na.rm = TRUE)
    mad_val <- mad(x, center = med, na.rm = TRUE)
    if (mad_val == 0) mad_val <- 1
    (x - med) / mad_val
  }))
}

# ---- Function to create publication-ready heatmap ----
make_publication_heatmap <- function(events_obj, diff_table, 
                                     keep_cols = 7:12,
                                     treatment_name = "Treatment",
                                     event_type = "Exons",
                                     color_palette = NULL,
                                     n_clusters = 2) {
  
  # Default color palettes
  if (is.null(color_palette)) {
    if (event_type == "Exons") {
      color_palette <- colorRamp2(c(-2, 0, 2), c("#E0E0E0", "#E0E0E0", "#7A2048"))
    } else {
      color_palette <- colorRamp2(c(-2, 0, 2), c("#E0E0E0", "#E0E0E0", "#2C2C2C"))
    }
  }
  
  # Filter PSI data
  psi_df <- events_obj$PSI %>% filter(complete.cases(dplyr::select(., all_of(keep_cols))))
  uniq <- psi_df %>% distinct(EVENT, .keep_all = TRUE)
  rownames(uniq) <- uniq$EVENT
  
  # Get significant events
  evts <- diff_table$EVENT
  if (length(evts) == 0) {
    stop("No significant events found")
  }
  
  # Extract matrix and scale
  mat <- uniq[rownames(uniq) %in% evts, , drop = FALSE] %>% 
    dplyr::select(all_of(keep_cols)) %>% 
    as.matrix()
  
  if (nrow(mat) == 0) {
    stop("No events found in PSI matrix")
  }
  
  scaled <- robust_scale_rows(mat)
  scaled <- pmax(pmin(scaled, 2), -2)  # Cap at -2 to 2
  
  # Create column annotations
  n_ctrl <- 3
  n_treat <- ncol(scaled) - n_ctrl
  col_split <- factor(c(rep("Control", n_ctrl), rep(treatment_name, n_treat)),
                      levels = c("Control", treatment_name))
  
  # Drug-specific colors
  drug_colors <- c(
    "Tubercidin" = "#A0C1B9",
    "Pladienolide B" = "#70A0AF",
    "Spliceostatin A" = "#706993"
  )
  
  # Column annotation colors - create named vector
  col_colors <- c("Control" = "#D3D3D3")
  col_colors[treatment_name] <- drug_colors[treatment_name]
  
  # Calculate Spearman correlation between samples
  cor_matrix <- cor(scaled, method = "spearman", use = "complete.obs")
  
  # Calculate average within-group and between-group correlations
  ctrl_idx <- 1:n_ctrl
  treat_idx <- (n_ctrl + 1):ncol(scaled)
  
  cor_within_ctrl <- mean(cor_matrix[ctrl_idx, ctrl_idx][upper.tri(cor_matrix[ctrl_idx, ctrl_idx])])
  cor_within_treat <- mean(cor_matrix[treat_idx, treat_idx][upper.tri(cor_matrix[treat_idx, treat_idx])])
  cor_between <- mean(cor_matrix[ctrl_idx, treat_idx])
  
  # Create correlation text annotation
  cor_text <- sprintf("Spearman ρ: Control=%.2f, %s=%.2f, Between=%.2f", 
                      cor_within_ctrl, treatment_name, cor_within_treat, cor_between)
  
  cat(sprintf("  %s\n", cor_text))
  
  column_ha <- HeatmapAnnotation(
    Condition = col_split,
    col = list(Condition = col_colors),
    annotation_name_side = "left",
    annotation_legend_param = list(
      Condition = list(title = "Condition", 
                       title_gp = gpar(fontsize = 10, fontface = "bold"),
                       labels_gp = gpar(fontsize = 9))
    )
  )
  
  # Create heatmap
  ht <- Heatmap(
    scaled,
    name = "Z-score",
    col = color_palette,
    
    # Column settings
    cluster_columns = TRUE,
    cluster_column_slices = FALSE,
    column_split = col_split,
    column_title = NULL,
    column_gap = unit(3, "mm"),
    show_column_names = FALSE,
    top_annotation = column_ha,
    
    # Row settings
    cluster_rows = TRUE,
    row_split = n_clusters,
    cluster_row_slices = FALSE,
    row_title = paste0("Cluster %s"),
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_title_rot = 0,
    row_gap = unit(2, "mm"),
    show_row_names = FALSE,
    show_row_dend = TRUE,
    row_dend_width = unit(10, "mm"),
    
    # Legend settings
    heatmap_legend_param = list(
      title = "Z-score",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9),
      legend_height = unit(4, "cm"),
      legend_direction = "vertical",
      title_position = "topcenter"
    ),
    
    # Appearance
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 0.5),
    
    # Size
    width = unit(8, "cm"),
    height = unit(12, "cm")
  )
  
  # Return both heatmap and correlation text
  return(list(heatmap = ht, cor_text = cor_text))
}

# ---- Create heatmaps for each dataset ----

datasets <- list(
  tub = list(
    events = tub_events,
    diff_exons = differential_tub_exons,
    diff_introns = differential_tub_introns,
    label = "Tubercidin",
    keep_cols = 7:12
  ),
  pladb = list(
    events = pladb_events,
    diff_exons = differential_pladb_exons,
    diff_introns = differential_pladb_introns,
    label = "Pladienolide B",
    keep_cols = 7:12
  ),
  ssa = list(
    events = ssa_events,
    diff_exons = differential_ssa_exons,
    diff_introns = differential_ssa_introns,
    label = "Spliceostatin A",
    keep_cols = 7:12
  )
)

# Define color palettes
exon_palette <- colorRamp2(c(-2, 0, 2), c("#E0E0E0", "#E0E0E0", "#7A2048"))
intron_palette <- colorRamp2(c(-2, 0, 2), c("#E0E0E0", "#E0E0E0", "#2C2C2C"))

# ---- Generate all heatmaps ----
for (nm in names(datasets)) {
  ds <- datasets[[nm]]
  
  cat(sprintf("\n=== Creating heatmaps for %s ===\n", ds$label))
  
  # EXONS heatmap
  if (nrow(ds$diff_exons) > 0) {
    tryCatch({
      pdf(sprintf("%s_exons_heatmap.pdf", nm), width = 8, height = 10)
      
      result <- make_publication_heatmap(
        events_obj = ds$events,
        diff_table = ds$diff_exons,
        keep_cols = ds$keep_cols,
        treatment_name = ds$label,
        event_type = "Exons",
        color_palette = exon_palette,
        n_clusters = 2
      )
      
      ht_exons <- result$heatmap
      cor_text <- result$cor_text
      
      draw(ht_exons, 
           padding = unit(c(5, 2, 2, 2), "mm"))  # Increased top padding
      
      # Add correlation text at bottom
      grid.text(cor_text, x = 0.5, y = 0.02, 
                gp = gpar(fontsize = 9, fontface = "italic", col = "gray30"))
      
      dev.off()
      cat(sprintf("✓ Saved %s_exons_heatmap.pdf (%d events)\n", nm, nrow(ds$diff_exons)))
    }, error = function(e) {
      cat(sprintf("✗ Error creating exons heatmap for %s: %s\n", nm, e$message))
    })
  } else {
    cat(sprintf("✗ No significant exons for %s\n", nm))
  }
  
  # INTRONS heatmap
  if (nrow(ds$diff_introns) > 0) {
    tryCatch({
      pdf(sprintf("%s_introns_heatmap.pdf", nm), width = 8, height = 10)
      
      result <- make_publication_heatmap(
        events_obj = ds$events,
        diff_table = ds$diff_introns,
        keep_cols = ds$keep_cols,
        treatment_name = ds$label,
        event_type = "Introns",
        color_palette = intron_palette,
        n_clusters = 2
      )
      
      ht_introns <- result$heatmap
      cor_text <- result$cor_text
      
      draw(ht_introns,
           padding = unit(c(5, 2, 2, 2), "mm"))  # Increased top padding
      
      # Add correlation text at bottom
      grid.text(cor_text, x = 0.5, y = 0.02, 
                gp = gpar(fontsize = 9, fontface = "italic", col = "gray30"))
      
      dev.off()
      cat(sprintf("✓ Saved %s_introns_heatmap.pdf (%d events)\n", nm, nrow(ds$diff_introns)))
    }, error = function(e) {
      cat(sprintf("✗ Error creating introns heatmap for %s: %s\n", nm, e$message))
    })
  } else {
    cat(sprintf("✗ No significant introns for %s\n", nm))
  }
}

cat("\n=== All heatmaps completed ===\n")