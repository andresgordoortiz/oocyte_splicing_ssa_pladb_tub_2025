# load libraries
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(pathview)
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

# Enhanced theme for publication-quality plots
pub_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 9),
    panel.grid = element_blank()
  )

output_dir <- "KEGG_enrichment_results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Store KEGG results for comparison
all_kegg_results <- list()

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
  
  # KEGG enrichment requires Entrez IDs
  if (!use_entrez || length(genes_entrez) == 0) {
    message("ERROR: KEGG requires valid Entrez IDs. Skipping ", nm)
    next
  }
  
  # Run KEGG enrichment
  kegg_result <- tryCatch({
    enrichKEGG(gene = genes_entrez,
               organism = 'mmu',  # mouse
               keyType = 'ncbi-geneid',
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05)
  }, error = function(e) {
    message("enrichKEGG error (", nm, "): ", e$message)
    NULL
  })
  
  # Check if result is valid
  if (!is.null(kegg_result) && nrow(as.data.frame(kegg_result)) > 0) {
    kegg_df <- as.data.frame(kegg_result)
    # Check gene counts
    gene_counts <- as.numeric(sapply(strsplit(as.character(kegg_df$GeneRatio), "/"), function(x) x[1]))
    valid_count <- sum(gene_counts > 0)
    message("  KEGG: ", nrow(kegg_df), " pathways found, ", valid_count, " with genes > 0")
    
    if (valid_count > 0) {
      # Filter out invalid entries with 0 genes
      kegg_df <- kegg_df %>%
        mutate(gene_count = as.numeric(sapply(strsplit(as.character(GeneRatio), "/"), function(x) x[1]))) %>%
        filter(gene_count > 0) %>%
        dplyr::select(-gene_count)
      
      all_kegg_results[[nm]] <- kegg_df
      
      # Save results
      write.csv(kegg_df, file.path(output_dir, paste0(nm, "_KEGG.csv")), row.names = FALSE)
      message("Saved: ", file.path(output_dir, paste0(nm, "_KEGG.csv")))
    }
  } else {
    message("No significant KEGG pathways for ", nm)
  }
  
  # -------------------- Enhanced Plots --------------------
  if (!is.null(kegg_result) && nrow(as.data.frame(kegg_result)) > 0) {
    kegg_df <- as.data.frame(kegg_result) %>%
      mutate(gene_count = as.numeric(sapply(strsplit(as.character(GeneRatio), "/"), function(x) x[1]))) %>%
      filter(gene_count > 0) %>%
      arrange(desc(Count)) %>%
      slice_head(n = 20) %>%
      mutate(
        GeneRatio_numeric = sapply(strsplit(as.character(GeneRatio), "/"), 
                                   function(x) as.numeric(x[1])/as.numeric(x[2])),
        Description = factor(Description, levels = rev(Description))
      )
    
    if (nrow(kegg_df) > 0) {
      # Print p-value range for verification
      message("P-value range for ", nm, ": ", 
              formatC(min(kegg_df$p.adjust), format = "e", digits = 2), " to ",
              formatC(max(kegg_df$p.adjust), format = "e", digits = 2))
      
      # Enhanced dot plot
      p_enhanced_dot <- ggplot(kegg_df, aes(x = GeneRatio_numeric, y = Description)) +
        geom_segment(aes(x = 0, xend = GeneRatio_numeric, y = Description, yend = Description),
                     color = "gray70", linewidth = 0.8) +
        geom_point(aes(size = Count, color = p.adjust), alpha = 0.8) +
        scale_color_gradient(low = "#d73027", high = "#4575b4", 
                             name = "Adjusted\np-value",
                             trans = "log10",
                             labels = function(x) formatC(x, format = "e", digits = 1),
                             guide = guide_colorbar(reverse = TRUE)) +
        scale_size_continuous(name = "Gene\nCount", range = c(3, 10)) +
        scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
        labs(title = paste0("KEGG Pathway Enrichment — ", nm),
             x = "Gene Ratio",
             y = NULL) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold", size = 14, hjust = 0),
          axis.title.x = element_text(face = "bold", size = 12),
          axis.text.y = element_text(size = 10, color = "black", face = "plain"),
          axis.text.x = element_text(size = 10, color = "black"),
          legend.title = element_text(face = "bold", size = 10),
          legend.text = element_text(size = 9),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.line.x = element_line(color = "black", linewidth = 0.5),
          axis.ticks.x = element_line(color = "black")
        )
      
      ggsave(file.path(output_dir, paste0(nm, "_KEGG_enhanced_dot.pdf")), 
             p_enhanced_dot, width = 10, height = 7)
      ggsave(file.path(output_dir, paste0(nm, "_KEGG_enhanced_dot.png")), 
             p_enhanced_dot, width = 10, height = 7, dpi = 300)
      
      # Standard bar plot
      p_bar <- barplot(kegg_result, showCategory = 15, title = paste0("KEGG Pathways — ", nm)) + pub_theme
      ggsave(file.path(output_dir, paste0(nm, "_KEGG_bar.pdf")), p_bar, width = 7, height = 5.5)
      ggsave(file.path(output_dir, paste0(nm, "_KEGG_bar.png")), p_bar, width = 7, height = 5.5, dpi = 300)
      
      # Dot plot
      p_dot <- dotplot(kegg_result, showCategory = 15, title = paste0("KEGG Pathways — ", nm)) + pub_theme
      ggsave(file.path(output_dir, paste0(nm, "_KEGG_dot.pdf")), p_dot, width = 7, height = 5.5)
      ggsave(file.path(output_dir, paste0(nm, "_KEGG_dot.png")), p_dot, width = 7, height = 5.5, dpi = 300)
      
      # Gene-concept network
      top_k <- min(8, nrow(as.data.frame(kegg_result)))
      try({
        p_cnet <- cnetplot(kegg_result, showCategory = top_k, foldChange = NULL, 
                           circular = FALSE, colorCategory = FALSE) +
          ggtitle(paste0("Gene-Pathway Network — ", nm)) + pub_theme
        ggsave(file.path(output_dir, paste0(nm, "_KEGG_cnet.pdf")), p_cnet, width = 9, height = 7)
        ggsave(file.path(output_dir, paste0(nm, "_KEGG_cnet.png")), p_cnet, width = 9, height = 7, dpi = 300)
      }, silent = TRUE)
      
      message("Created plots for ", nm)
    }
  }
  
  # Summary
  summary_file <- file.path(output_dir, paste0(nm, "_summary.txt"))
  sink(summary_file)
  cat("Dataset:", nm, "\n")
  cat("Input genes (unique):", length(unique(na.omit(genes_orig))), "\n")
  cat("Mapped Entrez (count):", length(genes_entrez), "\n")
  cat("KEGG pathways:", ifelse(is.null(kegg_result), 0, nrow(as.data.frame(kegg_result))), "\n")
  sink()
  message("Wrote summary: ", summary_file)
}

# -------------------- Common KEGG Pathways Analysis --------------------
message("\n---- Analyzing common KEGG pathways across treatments ----")

if (length(all_kegg_results) >= 2) {
  # Find common pathways
  all_pathways <- lapply(all_kegg_results, function(x) x$Description)
  common_pathways <- Reduce(intersect, all_pathways)
  
  message("Common KEGG pathways found: ", length(common_pathways))
  
  if (length(common_pathways) > 0) {
    # Prepare comparison data
    comparison_data <- do.call(rbind, lapply(names(all_kegg_results), function(nm) {
      df <- all_kegg_results[[nm]]
      df_common <- df[df$Description %in% common_pathways, ]
      if (nrow(df_common) > 0) {
        data.frame(
          Treatment = nm,
          Description = df_common$Description,
          GeneRatio = sapply(strsplit(df_common$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
          Count = df_common$Count,
          p.adjust = df_common$p.adjust,
          stringsAsFactors = FALSE
        )
      }
    }))
    
    # Save common pathways
    write.csv(comparison_data, 
              file.path(output_dir, "common_KEGG_pathways.csv"), 
              row.names = FALSE)
    
    # Select top pathways by average gene ratio or significance
    top_common <- comparison_data %>%
      group_by(Description) %>%
      summarise(avg_ratio = mean(GeneRatio), 
                min_pval = min(p.adjust),
                .groups = "drop") %>%
      arrange(desc(avg_ratio)) %>%
      slice_head(n = min(20, length(common_pathways))) %>%
      pull(Description)
    
    plot_data <- comparison_data %>%
      filter(Description %in% top_common) %>%
      mutate(
        Description = factor(Description, levels = rev(top_common)),
        Treatment = factor(Treatment, levels = c("tub", "pladb", "ssa"))
      )
    
    # Create comparison plot - clean and elegant
    p_comparison <- ggplot(plot_data, aes(x = Treatment, y = Description)) +
      geom_point(aes(size = Count, color = -log10(p.adjust)), alpha = 0.85) +
      scale_color_gradient(low = "#fee5d9", high = "#de2d26",
                           name = "-log10(adj p)") +
      scale_size_continuous(name = "Genes", range = c(4, 14)) +
      scale_x_discrete(expand = expansion(add = 1.2)) +
      scale_y_discrete(expand = expansion(add = 1.5)) +
      labs(title = "Common KEGG Pathways Across Treatments",
           x = NULL,
           y = NULL) +
      theme_void(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5, 
                                  margin = margin(b = 25)),
        axis.text.y = element_text(size = 10.5, color = "grey20", hjust = 1,
                                   margin = margin(r = 12), lineheight = 1.3),
        axis.text.x = element_text(size = 12, color = "grey20", face = "bold",
                                   margin = margin(t = 12)),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 10, margin = margin(b = 5)),
        legend.text = element_text(size = 9),
        legend.spacing.y = unit(0.2, "cm"),
        plot.margin = margin(25, 25, 25, 25)
      )
    
    ggsave(file.path(output_dir, "common_KEGG_comparison.pdf"), 
           p_comparison, width = 8, height = 10)
    ggsave(file.path(output_dir, "common_KEGG_comparison.png"), 
           p_comparison, width = 8, height = 10, dpi = 300)
    
    # Alternative: heatmap-style plot
    p_heatmap <- ggplot(plot_data, aes(x = Treatment, y = Description, fill = GeneRatio)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = Count), color = "white", fontface = "bold", size = 3.5) +
      scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
                           midpoint = median(plot_data$GeneRatio),
                           name = "Gene\nRatio") +
      labs(title = "Common KEGG Pathways - Heatmap View",
           x = NULL, y = NULL) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 11, color = "black", face = "bold", angle = 0),
        legend.title = element_text(face = "bold", size = 10),
        panel.grid = element_blank()
      )
    
    ggsave(file.path(output_dir, "common_KEGG_heatmap.pdf"), 
           p_heatmap, width = 7, height = 8)
    ggsave(file.path(output_dir, "common_KEGG_heatmap.png"), 
           p_heatmap, width = 7, height = 8, dpi = 300)
    
    message("Created comparison plots for common KEGG pathways")
  } else {
    message("No common KEGG pathways found across all treatments")
  }
} else {
  message("Need at least 2 treatments with KEGG results for comparison")
}

message("\nCompleted! All outputs in: ", normalizePath(output_dir))
message("Check for:")
message("  - Enhanced dot plots: *_KEGG_enhanced_dot.png")
message("  - Common pathways comparison: common_KEGG_comparison.png")
message("  - Common pathways heatmap: common_KEGG_heatmap.png")