##### RUN AFTER PROTEIN_STABILITY_DELTAPSI_VOLCANO SCRIPT ####
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)

# Get genes for enrichment
genes_for_enrichment <- unique(differential_summary$GENE[grepl("Stability",
                                                               differential_summary$highlight)])

# Perform GO enrichment for all three ontologies
ego_bp <- enrichGO(gene = genes_for_enrichment, 
                   OrgDb = org.Mm.eg.db, 
                   keyType = "SYMBOL",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05)

ego_cc <- enrichGO(gene = genes_for_enrichment, 
                   OrgDb = org.Mm.eg.db, 
                   keyType = "SYMBOL",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05)

ego_mf <- enrichGO(gene = genes_for_enrichment, 
                   OrgDb = org.Mm.eg.db, 
                   keyType = "SYMBOL",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05)

# Enhanced dot plot function
create_enhanced_go_dotplot <- function(ego_result, ontology_name = "BP", 
                                       top_n = 20, output_file = NULL) {
  
  # Prepare data
  ego_df <- as.data.frame(ego_result) %>%
    mutate(gene_count = as.numeric(sapply(strsplit(as.character(GeneRatio), "/"), 
                                          function(x) x[1]))) %>%
    filter(gene_count > 0 & p.adjust <= 0.05) %>%
    arrange(desc(Count)) %>%
    slice_head(n = top_n) %>%
    mutate(
      GeneRatio_numeric = sapply(strsplit(as.character(GeneRatio), "/"), 
                                 function(x) as.numeric(x[1])/as.numeric(x[2])),
      Description = factor(Description, levels = rev(Description)),
      p = p.adjust
    )
  
  if (nrow(ego_df) == 0) {
    message("No significant enrichment terms found for ", ontology_name)
    return(NULL)
  }
  
  # Set ontology full name
  ont_full <- switch(ontology_name,
                     "BP" = "Biological Process",
                     "CC" = "Cellular Component",
                     "MF" = "Molecular Function",
                     ontology_name)
  
  # Create plot with horizontal connecting lines
  p <- ggplot(ego_df, aes(x = GeneRatio_numeric, y = Description)) +
    # Horizontal connecting lines from y-axis to points
    geom_segment(aes(x = 0, xend = GeneRatio_numeric, 
                     y = Description, yend = Description),
                 color = "grey75", linewidth = 0.8) +
    # Points styled like common GO comparison
    geom_point(aes(size = Count, color = p), alpha = 0.85) +
    # Color gradient - lower p-values (more significant) are more red
    scale_color_gradient(low = "#de2d26", high = "#fee5d9",
                         name = "P adj",
                         guide = guide_colorbar(barwidth = 1.2, 
                                                barheight = 8)) +
    # Size scale
    scale_size_continuous(name = "Genes", range = c(4, 14),
                          guide = guide_legend(override.aes = list(color = "black"))) +
    # Expand axes
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.15))) +
    scale_y_discrete(expand = expansion(add = 0.6)) +
    # Labels
    labs(title = paste0("GO ", ont_full),
         x = "Gene Ratio",
         y = NULL) +
    # Theme matching common GO figure
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0,
                                margin = margin(b = 20)),
      axis.title.x = element_text(face = "bold", size = 12, color = "grey20",
                                  margin = margin(t = 12)),
      axis.text.y = element_text(size = 10.5, color = "grey20", hjust = 1,
                                 lineheight = 1.3),
      axis.text.x = element_text(size = 10, color = "grey20"),
      legend.title = element_text(face = "bold", size = 10, 
                                  margin = margin(b = 5)),
      legend.text = element_text(size = 9),
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(color = "grey40", linewidth = 0.5),
      axis.ticks.x = element_line(color = "grey40"),
      plot.margin = margin(20, 20, 20, 20)
    )
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(paste0(output_file, ".pdf"), p, width = 9, height = 7)
    ggsave(paste0(output_file, ".png"), p, width = 9, height = 7, dpi = 300)
    message("Saved: ", output_file, ".pdf and .png")
  }
  
  return(p)
}

# Create plots for all three ontologies
# Biological Process
if (!is.null(ego_bp) && nrow(as.data.frame(ego_bp)) > 0) {
  p_bp <- create_enhanced_go_dotplot(
    ego_result = ego_bp,
    ontology_name = "BP",
    top_n = 20,
    output_file = "GO_BP_enhanced_styled"
  )
  if (!is.null(p_bp)) print(p_bp)
}

# Cellular Component
if (!is.null(ego_cc) && nrow(as.data.frame(ego_cc)) > 0) {
  p_cc <- create_enhanced_go_dotplot(
    ego_result = ego_cc,
    ontology_name = "CC",
    top_n = 20,
    output_file = "GO_CC_enhanced_styled"
  )
  if (!is.null(p_cc)) print(p_cc)
}

# Molecular Function
if (!is.null(ego_mf) && nrow(as.data.frame(ego_mf)) > 0) {
  p_mf <- create_enhanced_go_dotplot(
    ego_result = ego_mf,
    ontology_name = "MF",
    top_n = 20,
    output_file = "GO_MF_enhanced_styled"
  )
  if (!is.null(p_mf)) print(p_mf)
}