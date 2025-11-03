library(readxl)
library(tidyverse)
library(ggrepel)
library(biomaRt)
library(dplyr)
library(stringr)
library(org.Mm.eg.db)
library(betAS)

# ============================================================================
# DATA LOADING (Keep your original code)
# ============================================================================

katarina_protein_stable_oocyte <- read_excel("41556_2024_1442_MOESM4_ESM.xlsx", 
                                             sheet = "TableS1")

names(katarina_protein_stable_oocyte)[names(katarina_protein_stable_oocyte) == "gene_name"] <- "GENE"
names(katarina_protein_stable_oocyte)[names(katarina_protein_stable_oocyte) == "mean percent H"] <- "stability"

katarina_protein_stable_oocyte <- katarina_protein_stable_oocyte %>%
  mutate(stability = as.numeric(stability)) %>%
  filter(!is.na(stability))

translated_genes <- read_excel("translated_genes.xlsx")

# Gene symbol mapping
genes_original <- translated_genes$Gene
genes_clean <- str_split(genes_original, " /// ") %>% sapply("[", 1)
is_probe <- str_detect(genes_clean, "^Mm\\.")
genes_to_map <- genes_clean[!is_probe]

ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

res_symbol <- getBM(attributes=c("external_gene_name","mgi_symbol"),
                    filters="external_gene_name",
                    values=genes_to_map,
                    mart=ensembl) %>%
  mutate(input=as.character(external_gene_name))

res_ensembl <- getBM(attributes=c("ensembl_gene_id","mgi_symbol"),
                     filters="ensembl_gene_id",
                     values=genes_to_map,
                     mart=ensembl) %>%
  mutate(input=as.character(ensembl_gene_id))

res_entrez <- getBM(attributes=c("entrezgene_id","mgi_symbol"),
                    filters="entrezgene_id",
                    values=genes_to_map,
                    mart=ensembl) %>%
  mutate(input=as.character(entrezgene_id))

gene_map <- bind_rows(
  res_symbol %>% dplyr::select(input, mgi_symbol),
  res_ensembl %>% dplyr::select(input, mgi_symbol),
  res_entrez %>% dplyr::select(input, mgi_symbol)
) %>% distinct(input, .keep_all = TRUE)

official_symbols_mapped <- sapply(genes_to_map, function(g) {
  res <- gene_map$mgi_symbol[gene_map$input == g]
  if(length(res) > 0) res[1] else NA_character_
})

official_symbols_full <- rep(NA_character_, length(genes_clean))
official_symbols_full[!is_probe] <- official_symbols_mapped

official_symbols_full <- sapply(seq_along(official_symbols_full), function(i) {
  if(!is.na(official_symbols_full[i])) return(official_symbols_full[i])
  match <- str_extract(genes_clean[i], "\\b[A-Z][a-z0-9]+\\b")
  if(!is.na(match)) return(match)
  return(genes_clean[i])
})

translated_genes <- translated_genes %>%
  mutate(symbol = official_symbols_full)

# Clean duplicates
katarina_clean <- katarina_protein_stable_oocyte %>%
  filter(!is.na(GENE)) %>%
  distinct(GENE, .keep_all = TRUE)

translated_clean <- translated_genes %>%
  filter(!is.na(symbol)) %>%
  distinct(symbol, .keep_all = TRUE)

# Import betAS output tables
tub_fdr_df <- read_csv("tub_fdr.csv")[,-1]
pladb_fdr_df <- read_csv("pladb_fdr.csv")[,-1]
ssa_fdr_df <- read_csv("ssa_fdr.csv")[,-1]

differential_tub <- na.omit(tub_fdr_df[tub_fdr_df$FDR <= 0.05 & abs(tub_fdr_df$deltapsi) >= 0.1,])
differential_pladb <- na.omit(pladb_fdr_df[pladb_fdr_df$FDR <= 0.05 & abs(pladb_fdr_df$deltapsi) >= 0.1,])
differential_ssa <- na.omit(ssa_fdr_df[ssa_fdr_df$FDR <= 0.05 & abs(ssa_fdr_df$deltapsi) >= 0.1,])

# ============================================================================
# FIXED DATA PREPARATION - KEEP ALL EVENTS SEPARATE
# ============================================================================

# Combine differential lists - KEEP ALL ROWS (no summarizing!)
differential_all <- bind_rows(
  differential_tub %>% mutate(group = "Tubercidin"),
  differential_pladb %>% mutate(group = "Pladienolide B"),
  differential_ssa %>% mutate(group = "Spliceostatin A")
) 

differential_all$group<-factor(differential_all$group, levels=c("Tubercidin","Pladienolide B","Spliceostatin A"))

# Count in how many groups each gene appears
gene_counts <- differential_all %>%
  group_by(GENE) %>%
  summarize(
    n_groups = n_distinct(group),
    groups_list = paste(sort(unique(group)), collapse = "+"),
    .groups = "drop"
  )

# Join with stability data - KEEP ONE ROW PER GENE-GROUP COMBINATION
differential_summary <- differential_all %>%
  left_join(katarina_clean %>% dplyr::select(GENE, stability), by = "GENE") %>%
  left_join(gene_counts, by = "GENE") %>%
  filter(!is.na(stability))  # Remove genes without stability data

# Define highlighting categories
differential_summary <- differential_summary %>%
  mutate(
    highlight = case_when(
      n_groups == 3 & stability > 10 ~ "All 3 Groups + High Stability",
      n_groups == 3 ~ "All 3 Groups",
      n_groups == 2 & stability > 10 ~ "2 Groups + High Stability",
      n_groups == 2 ~ "2 Groups",
      stability > 10 ~ "High Stability Only",
      TRUE ~ "Other"
    ),
    highlight = factor(highlight, levels = c(
      "All 3 Groups + High Stability",
      "All 3 Groups",
      "2 Groups + High Stability",
      "2 Groups",
      "High Stability Only",
      "Other"
    ))
  )

# Prepare labels: only highly relevant genes
# Label top genes by stability in each facet
label_genes <- differential_summary %>%
  filter(highlight %in% c("All 3 Groups + High Stability", 
                          "2 Groups + High Stability")) %>%
  group_by(group) %>%
  arrange(desc(stability)) %>%
  slice_head(n = 8) %>%  # Top 8 per facet
  ungroup()

# ============================================================================
# PUBLICATION THEME
# ============================================================================

theme_publication <- function(base_size = 11) {
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
# CREATE FACETED VOLCANO PLOT
# ============================================================================

# Color palette
color_palette <- c(
  "All 3 Groups + High Stability" = "#D32F2F",    # Red
  "All 3 Groups" = "#E57373",                      # Light red
  "2 Groups + High Stability" = "#FF6F00",         # Orange
  "2 Groups" = "#FFB74D",                          # Light orange
  "High Stability Only" = "#1976D2",               # Blue
  "Other" = "grey70"                               # Grey
)

# Size mapping
size_palette <- c(
  "All 3 Groups + High Stability" = 3,
  "All 3 Groups" = 2.5,
  "2 Groups + High Stability" = 2.5,
  "2 Groups" = 2,
  "High Stability Only" = 2,
  "Other" = 1.5
)

# Alpha mapping
alpha_palette <- c(
  "All 3 Groups + High Stability" = 0.9,
  "All 3 Groups" = 0.7,
  "2 Groups + High Stability" = 0.8,
  "2 Groups" = 0.6,
  "High Stability Only" = 0.7,
  "Other" = 0.3
)

p_volcano <- ggplot(differential_summary, aes(x = deltapsi, y = stability)) +
  # Background points first (Other)
  geom_point(
    data = differential_summary %>% filter(highlight == "Other"),
    aes(color = highlight, size = highlight, alpha = highlight)
  ) +
  # Foreground points (interesting ones)
  geom_point(
    data = differential_summary %>% filter(highlight != "Other"),
    aes(color = highlight, size = highlight, alpha = highlight)
  ) +
  # Threshold lines
  geom_hline(yintercept = 10, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey60", linewidth = 0.4) +
  # Labels
  geom_text_repel(
    data = label_genes,
    aes(label = GENE),
    size = 2.5,
    fontface = "bold",
    max.overlaps = 15,
    min.segment.length = 0,
    segment.size = 0.3,
    segment.alpha = 0.6,
    box.padding = 0.3,
    point.padding = 0.2,
    force = 2,
    seed = 42
  ) +
  # Facet by group
  facet_wrap(~group, ncol = 3) +
  # Scales
  scale_color_manual(values = color_palette, name = "Category") +
  scale_size_manual(values = size_palette, guide = "none") +
  scale_alpha_manual(values = alpha_palette, guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  # Labels
  labs(
    x = expression(bold(Delta*"PSI")),
    y = "Protein Stability (%)",
    title = "Differential Splicing vs Protein Stability by Treatment",
    subtitle = paste0("Each point represents a significant splicing event (FDR < 0.05, |ΔPSI| ≥ 0.1)")
  ) +
  # Theme
  theme_publication(base_size = 10) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.9)))

# Check data before plotting
cat("\n=== Data Check ===\n")
cat("Total events to plot:", nrow(differential_summary), "\n")
cat("Events per group:\n")
print(table(differential_summary$group))
cat("\nRange of deltapsi:", range(differential_summary$deltapsi, na.rm = TRUE), "\n")
cat("Range of stability:", range(differential_summary$stability, na.rm = TRUE), "\n")
cat("Number of labels:", nrow(label_genes), "\n\n")

# Verify no near-zero deltapsi values
near_zero <- differential_summary %>% filter(abs(deltapsi) < 0.05)
if(nrow(near_zero) > 0) {
  cat("WARNING: Found", nrow(near_zero), "events with |deltapsi| < 0.05\n")
  print(head(near_zero))
} else {
  cat("✓ All deltapsi values are >= 0.1 as expected\n")
}

# Save plot
cat("\nSaving plot to PDF...\n")
ggsave("volcano_stability_deltapsi_faceted.pdf", plot = p_volcano, 
       width = 14, height = 5, device = cairo_pdf)
cat("Plot saved successfully to volcano_stability_deltapsi_faceted.pdf\n\n")

# Try to display
tryCatch({
  print(p_volcano)
}, error = function(e) {
  cat("Note: Could not display plot in RStudio viewer, but PDF was saved successfully.\n")
  cat("Error message:", e$message, "\n")
})

cat("\n=== Plot saved successfully ===\n")

# Summary statistics
cat("\n=== Summary Statistics ===\n")
cat("\nBy highlight category:\n")
print(table(differential_summary$highlight))
cat("\nBy group:\n")
print(table(differential_summary$group))
cat("\nGenes appearing in multiple groups:\n")
print(table(differential_summary$n_groups))
cat("\nGenes with all 3 groups + high stability:\n")
print(differential_summary %>% 
        filter(highlight == "All 3 Groups + High Stability") %>% 
        distinct(GENE) %>% 
        pull(GENE))