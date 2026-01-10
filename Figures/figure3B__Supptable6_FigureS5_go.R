# ---- 1. Setup and Libraries ----
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(openxlsx) # Added for single Excel output
library(readxl)
library(scales)
library(ggrepel)

# Define IO paths
excel_file <- "Preprocessing/betAS_out/Supplementary2_betAS_splicing_results.xlsx"
output_dir <- "GO_enrichment_results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Helper Functions ----

# Simplified Gene Mapper
map_to_entrez <- function(genes) {
  genes <- unique(na.omit(trimws(as.character(genes))))
  if (length(genes) == 0) return(NULL)
  
  # Try mapping ENSEMBL first, then SYMBOL
  mapped <- tryCatch(
    bitr(genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db),
    error = function(e) NULL
  )
  
  if (is.null(mapped) || nrow(mapped) == 0) {
    mapped <- tryCatch(
      bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db),
      error = function(e) NULL
    )
  }
  
  if (is.null(mapped)) return(NULL)
  return(unique(mapped$ENTREZID))
}

# Function to create the specific "Enhanced" Dotplot style you requested
draw_enhanced_dotplot <- function(enrich_result, title_str) {
  # Convert to DF and calculate numeric ratio
  df <- as.data.frame(enrich_result) %>%
    mutate(
      GeneRatio_numeric = sapply(strsplit(as.character(GeneRatio), "/"), 
                                 function(x) as.numeric(x[1])/as.numeric(x[2]))
    ) %>%
    arrange(desc(Count)) %>%
    slice_head(n = 20) %>% # Top 20 terms
    mutate(Description = factor(Description, levels = rev(Description)))
  
  if(nrow(df) == 0) return(NULL)
  
  ggplot(df, aes(x = GeneRatio_numeric, y = Description)) +
    geom_segment(aes(x = 0, xend = GeneRatio_numeric, y = Description, yend = Description),
                 color = "gray70", linewidth = 0.8) +
    geom_point(aes(size = Count, color = p.adjust), alpha = 0.8) +
    scale_color_gradient(low = "#d73027", high = "#4575b4", 
                         name = "Adj. P-value",
                         trans = "log10",
                         guide = guide_colorbar(reverse = TRUE)) +
    scale_size_continuous(range = c(3, 10)) +
    labs(title = title_str, x = "Gene Ratio", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 10, color = "black"),
      panel.grid.major.y = element_blank()
    )
}

# ---- 3. Load Data ----
raw_data <- list(
  Tub = read_excel(excel_file, sheet = "Tub_FDR"),
  PlaDB = read_excel(excel_file, sheet = "PlaDB_FDR"),
  SSA = read_excel(excel_file, sheet = "SSA_FDR")
)

# Filter for significant events (FDR <= 0.05 & deltaPSI >= 0.1)
datasets <- lapply(raw_data, function(df) {
  na.omit(df[df$FDR <= 0.05 & abs(df$deltapsi) >= 0.1, ])
})

# ---- 4. Main Enrichment Loop ----
excel_sheets <- list()      # To store data for the single Excel file
comparison_list <- list()   # To store results for the final common plot
ontologies <- c("BP", "MF", "CC")

for (drug_name in names(datasets)) {
  message(paste("Processing:", drug_name))
  
  # Map Genes
  genes_entrez <- map_to_entrez(datasets[[drug_name]]$GENE)
  
  if (length(genes_entrez) < 5) {
    message(paste("  Skipping", drug_name, "- insufficient mapped genes."))
    next
  }
  
  # Run Enrichment for each Ontology
  for (ont in ontologies) {
    message(paste("  Running GO:", ont))
    
    ego <- enrichGO(gene = genes_entrez,
                    OrgDb = org.Mm.eg.db,
                    ont = ont,
                    keyType = "ENTREZID",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    readable = TRUE)
    
    # Process if results exist
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      # 1. Add to list for Excel Export
      sheet_name <- paste0(drug_name, "_", ont)
      excel_sheets[[sheet_name]] <- as.data.frame(ego)
      
      # 2. Generate and Save Enhanced Dotplot
      p <- draw_enhanced_dotplot(ego, paste0(drug_name, " - GO ", ont))
      if (!is.null(p)) {
        ggsave(filename = file.path(output_dir, paste0(sheet_name, "_dotplot.png")),
               plot = p, width = 9, height = 7, dpi = 300)
      }
      
      # 3. Store BP results for the Common Plot (usually BP is best for comparison)
      if (ont == "BP") {
        comparison_list[[drug_name]] <- as.data.frame(ego)
      }
    }
  }
}

# ---- 5. Write Single Excel File ----
if (length(excel_sheets) > 0) {
  write.xlsx(excel_sheets, file = file.path(output_dir, "Supplementary6_GO_Results_Combined.xlsx"))
  message("Saved combined Excel file.")
}

# ---- 6. Final Common GO Terms Plot (BP) ----
message("Generating Common Terms Plot...")

if (length(comparison_list) >= 2) {
  # Find intersection of terms across available datasets
  term_lists <- lapply(comparison_list, function(x) x$Description)
  common_terms <- Reduce(intersect, term_lists)
  
  if (length(common_terms) > 0) {
    # Combine data for common terms
    plot_data <- do.call(rbind, lapply(names(comparison_list), function(nm) {
      df <- comparison_list[[nm]]
      df %>% 
        filter(Description %in% common_terms) %>%
        mutate(Treatment = nm)
    }))
    
    # Select top terms to prevent overcrowding (top 15 by count)
    top_terms <- plot_data %>% 
      group_by(Description) %>% 
      summarise(tot = sum(Count)) %>% 
      arrange(desc(tot)) %>% 
      slice_head(n = 15) %>% 
      pull(Description)
    
    final_df <- plot_data %>% 
      filter(Description %in% top_terms) %>%
      mutate(Treatment = factor(Treatment, levels = c("Tub", "PlaDB", "SSA")))
    
    # Comparison Plot
    p_comp <- ggplot(final_df, aes(x = Treatment, y = Description)) +
      geom_point(aes(size = Count, color = -log10(p.adjust)), alpha = 0.9) +
      scale_color_gradient(low = "#fee5d9", high = "#de2d26", name = "-log10(p.adj)") +
      scale_size_continuous(range = c(4, 12)) +
      labs(title = "Common GO Biological Processes", x = NULL, y = NULL) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(size = 10),
        panel.grid = element_line(color = "gray95")
      )
    
    ggsave(file.path(output_dir, "Common_GO_Comparison.png"), p_comp, width = 8, height = 9, dpi = 300)
    ggsave(file.path(output_dir, "Common_GO_Comparison.pdf"), p_comp, width = 8, height = 9)
    message("Saved Common GO Comparison plot.")
    
  } else {
    message("No strictly common GO terms found across all 3 groups.")
  }
}

message("Done! Check 'GO_enrichment_results' folder.")