# load libraries
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(cowplot)
library(ggrepel)
library(scales)

# Read event info and differential analysis results
tub_fdr_df <- read_csv("tub_fdr.csv")[,-1]
pladb_fdr_df <- read_csv("pladb_fdr.csv")[,-1]
ssa_fdr_df <- read_csv("ssa_fdr.csv")[,-1]

# Filter for significant events
differential_tub <- na.omit(tub_fdr_df[tub_fdr_df$FDR <= 0.05 & abs(tub_fdr_df$deltapsi) >= 0.1,])
differential_pladb <- na.omit(pladb_fdr_df[pladb_fdr_df$FDR <= 0.05 & abs(pladb_fdr_df$deltapsi) >= 0.1,])
differential_ssa <- na.omit(ssa_fdr_df[ssa_fdr_df$FDR <= 0.05 & abs(ssa_fdr_df$deltapsi) >= 0.1,])

datasets <- list(tub = differential_tub,
                 pladb = differential_pladb,
                 ssa = differential_ssa)

map_genes_to_entrez <- function(genes) {
  genes <- unique(na.omit(trimws(as.character(genes))))
  if (length(genes) == 0) return(character(0))
  
  prop_ensembl <- mean(grepl("^ENSMUSG", genes, ignore.case = TRUE))
  if (prop_ensembl > 0.5) {
    mapped <- tryCatch(
      bitr(genes, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb = org.Mm.eg.db),
      error = function(e) data.frame()
    )
  } else {
    mapped <- tryCatch(
      bitr(genes, fromType = "SYMBOL", toType = c("ENTREZID","ENSEMBL"), OrgDb = org.Mm.eg.db),
      error = function(e) data.frame()
    )
    if (nrow(mapped) == 0) {
      mapped <- tryCatch(
        bitr(genes, fromType = "ALIAS", toType = c("ENTREZID","SYMBOL"), OrgDb = org.Mm.eg.db),
        error = function(e) data.frame()
      )
    }
  }
  
  if (nrow(mapped) == 0) {
    warning("No Entrez mappings found; returning original GENE values.")
    return(genes)
  }
  unique(as.character(mapped$ENTREZID))
}

output_dir <- "GO_enrichment_results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Enhanced dot plot function matching the short script style
create_enhanced_go_dotplot <- function(ego_result, ontology_name = "BP", 
                                       treatment_name = "", top_n = 20, 
                                       output_file = NULL) {
  ego_df <- as.data.frame(ego_result) %>%
    mutate(gene_count = as.numeric(sapply(strsplit(as.character(GeneRatio), "/"), 
                                          function(x) x[1]))) %>%
    filter(gene_count > 0 & p.adjust <= 0.05) %>%
    arrange(desc(Count)) %>%
    slice_head(n = top_n) %>%
    mutate(
      GeneRatio_numeric = sapply(strsplit(as.character(GeneRatio), "/"), 
                                 function(x) as.numeric(x[1]) / as.numeric(x[2])),
      Description = factor(Description, levels = rev(Description)),
      p_adj = p.adjust,
      # NEW: color mapped to -log10(p.adjust) so smaller p -> larger numeric -> clearer gradient
      p_col = -log10(p_adj)
    )
  
  if (nrow(ego_df) == 0) {
    message("No significant enrichment terms found for ", treatment_name, " - ", ontology_name)
    return(NULL)
  }
  
  ont_full <- switch(ontology_name,
                     "BP" = "Biological Process",
                     "CC" = "Cellular Component",
                     "MF" = "Molecular Function",
                     ontology_name)
  
  # choose color direction so higher -log10(p) (more significant) = deeper red
  # define breaks for colorbar in terms of original p-values you want labeled
  p_breaks_vals <- c(1e-6, 1e-5,1e-4 ,1e-3, 1e-2, 1e-1, 0.05)
  # convert to -log10 scale for the colorbar positions
  p_breaks_log <- -log10(p_breaks_vals)
  # create labels as scientific notation (or any format you prefer)
  p_breaks_labels <- scales::scientific(p_breaks_vals, digits = 2)
  
  library(scales) # for scientific()
  
  p <- ggplot(ego_df, aes(x = GeneRatio_numeric, y = Description)) +
    geom_segment(aes(x = 0, xend = GeneRatio_numeric, 
                     y = Description, yend = Description),
                 color = "grey75", linewidth = 0.8) +
    # map color to p_col (-log10(p.adjust)), size to Count
    geom_point(aes(size = Count, color = p_col), alpha = 0.85) +
    # color scale: low = light (less significant), high = deep red (more significant)
    scale_color_gradient(
      low = "#fee5d9", high = "#de2d26",
      name = expression(-log[10](adj~P)),
      # set breaks and labels on the -log10 scale but label them with the original p-values:
      breaks = p_breaks_log,
      labels = p_breaks_labels,
      guide = guide_colorbar(barwidth = 1.2, barheight = 8, ticks = TRUE)
    ) +
    scale_size_continuous(name = "Genes", range = c(4, 12),
                          guide = guide_legend(override.aes = list(color = "black"))) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.15))) +
    scale_y_discrete(expand = expansion(add = 0.6)) +
    labs(title = paste0("GO ", ont_full, " — ", treatment_name),
         x = "Gene Ratio",
         y = NULL) +
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
  
  if (!is.null(output_file)) {
    ggsave(paste0(output_file, ".pdf"), p, width = 9, height = 7)
    ggsave(paste0(output_file, ".png"), p, width = 9, height = 7, dpi = 300)
    message("Saved: ", output_file, ".pdf and .png")
  }
  
  return(p)
}


# -------------------- Main loop --------------------
for (nm in names(datasets)) {
  message("---- Processing: ", nm, " ----")
  df <- datasets[[nm]]
  genes_orig <- df$GENE
  
  # Diagnostic: show sample genes
  message("Sample genes from input: ", paste(head(unique(genes_orig), 5), collapse = ", "))
  
  genes_entrez <- map_genes_to_entrez(genes_orig)
  use_entrez <- length(genes_entrez) > 0 && all(grepl("^[0-9]+$", genes_entrez))
  
  message("Unique input genes: ", length(unique(na.omit(genes_orig))))
  message("Mapped Entrez IDs: ", length(genes_entrez))
  message("Mapping success rate: ", 
          round(100 * length(genes_entrez) / length(unique(na.omit(genes_orig))), 1), "%")
  
  # If mapping failed badly, try alternative approaches
  if (length(genes_entrez) < 5) {
    message("WARNING: Very few genes mapped. Trying alternative mapping strategies...")
    
    # Try cleaning gene names
    genes_clean <- unique(na.omit(trimws(toupper(as.character(genes_orig)))))
    genes_clean <- genes_clean[genes_clean != "" & genes_clean != "NA"]
    
    message("Cleaned genes count: ", length(genes_clean))
    message("Sample cleaned genes: ", paste(head(genes_clean, 5), collapse = ", "))
    
    # Try mapping cleaned genes
    genes_entrez <- map_genes_to_entrez(genes_clean)
    message("After cleaning - Mapped Entrez IDs: ", length(genes_entrez))
    
    if (length(genes_entrez) > 0) {
      use_entrez <- all(grepl("^[0-9]+$", genes_entrez))
      genes_orig <- genes_clean  # Use cleaned genes
    }
  }
  
  # Skip if still no genes mapped
  if (length(genes_entrez) == 0) {
    message("ERROR: No genes could be mapped for ", nm, ". Skipping enrichment analysis.")
    next
  }
  
  run_enrichGO_safe <- function(gvec, ont) {
    tryCatch({
      result <- NULL
      if (use_entrez) {
        result <- enrichGO(gene = gvec, OrgDb = org.Mm.eg.db, ont = ont, keyType = "ENTREZID",
                           pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
      } else {
        result <- enrichGO(gene = unique(na.omit(as.character(genes_orig))), OrgDb = org.Mm.eg.db, ont = ont,
                           keyType = "SYMBOL", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
      }
      
      # Check if result is valid
      if (!is.null(result) && nrow(as.data.frame(result)) > 0) {
        result_df <- as.data.frame(result)
        # Check gene counts
        gene_counts <- as.numeric(sapply(strsplit(as.character(result_df$GeneRatio), "/"), function(x) x[1]))
        valid_count <- sum(gene_counts > 0)
        message("  ", ont, ": ", nrow(result_df), " terms found, ", valid_count, " with genes > 0")
        
        if (valid_count == 0) {
          message("  WARNING: All terms have 0 genes! This indicates a mapping problem.")
          return(NULL)
        }
      }
      
      return(result)
    }, error = function(e) {
      message("enrichGO error (", nm, ", ", ont, "): ", e$message)
      NULL
    })
  }
  
  ego_bp <- run_enrichGO_safe(genes_entrez, "BP")
  ego_mf <- run_enrichGO_safe(genes_entrez, "MF")
  ego_cc <- run_enrichGO_safe(genes_entrez, "CC")
  
  # Save CSVs
  save_if <- function(x, suffix) {
    if (!is.null(x) && nrow(as.data.frame(x)) > 0) {
      f <- file.path(output_dir, paste0(nm, "_", suffix, ".csv"))
      write.csv(as.data.frame(x), f, row.names = FALSE)
      message("Saved: ", f)
    }
  }
  save_if(ego_bp, "GO_BP")
  save_if(ego_mf, "GO_MF")
  save_if(ego_cc, "GO_CC")
  
  # -------------------- Create Enhanced Dot Plots for All Ontologies --------------------
  
  # Biological Process
  if (!is.null(ego_bp) && nrow(as.data.frame(ego_bp)) > 0) {
    p_bp <- create_enhanced_go_dotplot(
      ego_result = ego_bp,
      ontology_name = "BP",
      treatment_name = nm,
      top_n = 20,
      output_file = file.path(output_dir, paste0(nm, "_BP_enhanced"))
    )
  } else {
    message("No significant BP enrichment for ", nm)
  }
  
  # Molecular Function
  if (!is.null(ego_mf) && nrow(as.data.frame(ego_mf)) > 0) {
    p_mf <- create_enhanced_go_dotplot(
      ego_result = ego_mf,
      ontology_name = "MF",
      treatment_name = nm,
      top_n = 20,
      output_file = file.path(output_dir, paste0(nm, "_MF_enhanced"))
    )
  } else {
    message("No significant MF enrichment for ", nm)
  }
  
  # Cellular Component
  if (!is.null(ego_cc) && nrow(as.data.frame(ego_cc)) > 0) {
    p_cc <- create_enhanced_go_dotplot(
      ego_result = ego_cc,
      ontology_name = "CC",
      treatment_name = nm,
      top_n = 20,
      output_file = file.path(output_dir, paste0(nm, "_CC_enhanced"))
    )
  } else {
    message("No significant CC enrichment for ", nm)
  }
  
  # Summary
  summary_file <- file.path(output_dir, paste0(nm, "_summary.txt"))
  sink(summary_file)
  cat("Dataset:", nm, "\n")
  cat("Input genes (unique):", length(unique(na.omit(genes_orig))), "\n")
  cat("Mapped Entrez (count):", length(genes_entrez), "\n")
  cat("GO BP terms:", ifelse(is.null(ego_bp), 0, nrow(as.data.frame(ego_bp))), "\n")
  cat("GO MF terms:", ifelse(is.null(ego_mf), 0, nrow(as.data.frame(ego_mf))), "\n")
  cat("GO CC terms:", ifelse(is.null(ego_cc), 0, nrow(as.data.frame(ego_cc))), "\n")
  sink()
  message("Wrote summary: ", summary_file)
}

message("\nCompleted! All outputs in: ", normalizePath(output_dir))
message("Enhanced dot plots generated:")
message("  - BP plots: tub_BP_enhanced.png, pladb_BP_enhanced.png, ssa_BP_enhanced.png")
message("  - MF plots: tub_MF_enhanced.png, pladb_MF_enhanced.png, ssa_MF_enhanced.png")
message("  - CC plots: tub_CC_enhanced.png, pladb_CC_enhanced.png, ssa_CC_enhanced.png")