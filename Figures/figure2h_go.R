# load libraries
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(enrichR)
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

output_dir <- "GO_enrichment_results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

enrichr_dbs <- c("GO_Biological_Process_2021","GO_Molecular_Function_2021","GO_Cellular_Component_2021")

# Store GO results for comparison
all_go_results <- list()

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
  
  # Store BP results for comparison
  if (!is.null(ego_bp) && nrow(as.data.frame(ego_bp)) > 0) {
    ego_bp_df <- as.data.frame(ego_bp)
    # Filter out invalid entries with 0 genes
    ego_bp_df <- ego_bp_df %>%
      mutate(gene_count = as.numeric(sapply(strsplit(as.character(GeneRatio), "/"), function(x) x[1]))) %>%
      filter(gene_count > 0) %>%
      dplyr::select(-gene_count)
    
    if (nrow(ego_bp_df) > 0) {
      all_go_results[[nm]] <- ego_bp_df
    }
  }
  
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
  
  # Enrichr (optional)
  if (length(unique(na.omit(genes_orig))) >= 5) {
    try({
      enr <- enrichr(unique(na.omit(as.character(genes_orig))), enrichr_dbs)
      for (db in names(enr)) {
        if (!is.null(enr[[db]]) && nrow(enr[[db]]) > 0) {
          fdb <- file.path(output_dir, paste0(nm, "_Enrichr_", db, ".csv"))
          write.csv(enr[[db]], fdb, row.names = FALSE)
          message("Saved Enrichr: ", fdb)
        }
      }
    }, silent = TRUE)
  }
  
  # -------------------- Enhanced Plots --------------------
  plots <- list()
  
  if (!is.null(ego_bp) && nrow(as.data.frame(ego_bp)) > 0) {
    # Enhanced dot plot with connecting lines and improved aesthetics
    ego_df <- as.data.frame(ego_bp) %>%
      mutate(gene_count = as.numeric(sapply(strsplit(as.character(GeneRatio), "/"), function(x) x[1]))) %>%
      filter(gene_count > 0) %>%  # Remove entries with 0 genes
      arrange(desc(Count)) %>%
      slice_head(n = 20) %>%
      mutate(
        GeneRatio_numeric = sapply(strsplit(as.character(GeneRatio), "/"), 
                                   function(x) as.numeric(x[1])/as.numeric(x[2])),
        Description = factor(Description, levels = rev(Description))
      )
    
    if (nrow(ego_df) == 0) {
      message("No valid BP enrichment for ", nm, " after filtering")
    } else {
      # Print p-value range for verification
      message("P-value range for ", nm, ": ", 
              formatC(min(ego_df$p.adjust), format = "e", digits = 2), " to ",
              formatC(max(ego_df$p.adjust), format = "e", digits = 2))
      
      p_enhanced_dot <- ggplot(ego_df, aes(x = GeneRatio_numeric, y = Description)) +
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
        labs(title = paste0("GO Biological Process — ", nm),
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
      
      plots$bp_enhanced_dot <- p_enhanced_dot
      ggsave(file.path(output_dir, paste0(nm, "_BP_enhanced_dot.pdf")), 
             p_enhanced_dot, width = 9, height = 7)
      ggsave(file.path(output_dir, paste0(nm, "_BP_enhanced_dot.png")), 
             p_enhanced_dot, width = 9, height = 7, dpi = 300)
      
      # Standard plots
      p_bar <- barplot(ego_bp, showCategory = 12, title = paste0("GO BP (top) — ", nm)) + pub_theme
      plots$bp_bar <- p_bar
      ggsave(file.path(output_dir, paste0(nm, "_BP_bar.pdf")), p_bar, width = 6.5, height = 5)
      ggsave(file.path(output_dir, paste0(nm, "_BP_bar.png")), p_bar, width = 6.5, height = 5, dpi = 300)
      
      top_k <- min(6, nrow(as.data.frame(ego_bp)))
      try({
        p_cnet <- cnetplot(ego_bp, showCategory = top_k, foldChange = NULL, 
                           circular = FALSE, colorCategory = FALSE) +
          ggtitle(paste0("Gene-Concept Network — ", nm)) + pub_theme
        plots$bp_cnet <- p_cnet
        ggsave(file.path(output_dir, paste0(nm, "_BP_cnet.pdf")), p_cnet, width = 8, height = 6)
        ggsave(file.path(output_dir, paste0(nm, "_BP_cnet.png")), p_cnet, width = 8, height = 6, dpi = 300)
      }, silent = TRUE)
    }
  } else {
    message("No significant BP enrichment for ", nm)
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

# -------------------- Common GO Terms Analysis --------------------
message("\n---- Analyzing common GO terms across treatments ----")

if (length(all_go_results) >= 2) {
  # Find common GO terms
  all_terms <- lapply(all_go_results, function(x) x$Description)
  common_terms <- Reduce(intersect, all_terms)
  
  message("Common GO terms found: ", length(common_terms))
  
  if (length(common_terms) > 0) {
    # Prepare comparison data
    comparison_data <- do.call(rbind, lapply(names(all_go_results), function(nm) {
      df <- all_go_results[[nm]]
      df_common <- df[df$Description %in% common_terms, ]
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
    
    # Save common terms
    write.csv(comparison_data, 
              file.path(output_dir, "common_GO_terms.csv"), 
              row.names = FALSE)
    
    # Select top terms by average gene ratio or significance
    top_common <- comparison_data %>%
      group_by(Description) %>%
      summarise(avg_ratio = mean(GeneRatio), 
                min_pval = min(p.adjust),
                .groups = "drop") %>%
      arrange(desc(avg_ratio)) %>%
      slice_head(n = min(15, length(common_terms))) %>%
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
      labs(title = "Common GO Terms Across Treatments",
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
    
    ggsave(file.path(output_dir, "common_GO_comparison.pdf"), 
           p_comparison, width = 8, height = 9)
    ggsave(file.path(output_dir, "common_GO_comparison.png"), 
           p_comparison, width = 8, height = 9, dpi = 300)
    
    # Alternative: heatmap-style plot
    p_heatmap <- ggplot(plot_data, aes(x = Treatment, y = Description, fill = GeneRatio)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = Count), color = "white", fontface = "bold", size = 3.5) +
      scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
                           midpoint = median(plot_data$GeneRatio),
                           name = "Gene\nRatio") +
      labs(title = "Common GO Terms - Heatmap View",
           x = NULL, y = NULL) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 11, color = "black", face = "bold", angle = 0),
        legend.title = element_text(face = "bold", size = 10),
        panel.grid = element_blank()
      )
    
    ggsave(file.path(output_dir, "common_GO_heatmap.pdf"), 
           p_heatmap, width = 7, height = 7)
    ggsave(file.path(output_dir, "common_GO_heatmap.png"), 
           p_heatmap, width = 7, height = 7, dpi = 300)
    
    message("Created comparison plots for common GO terms")
  } else {
    message("No common GO terms found across all treatments")
  }
} else {
  message("Need at least 2 treatments with GO results for comparison")
}

message("\nCompleted! All outputs in: ", normalizePath(output_dir))
message("Check for:")
message("  - Enhanced dot plots: *_BP_enhanced_dot.png")
message("  - Common terms comparison: common_GO_comparison.png")
message("  - Common terms heatmap: common_GO_heatmap.png")