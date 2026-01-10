library(ggplot2)
library(ggrepel)
library(dplyr)
library(scales)
library(showtext)
library(readxl) # Required to read the Excel file

# ---- Path to the Excel file created in the previous step ----
excel_path <- "./Preprocessing/gene_expression/Supplementary1_DESeq2_Results.xlsx"

# ---- Load data from Excel tabs ----
# We use the Shrunken tabs as they are better for visualization
ssa_deseq   <- read_excel(excel_path, sheet = "SSA_Shrunken")
pladb_deseq <- read_excel(excel_path, sheet = "PlaDB_Shrunken")
tub_deseq   <- read_excel(excel_path, sheet = "Tub_Shrunken")

# ---- Drug colors ----
ssa_color   <- "#706993"
pladb_color <- "#70A0AF"
tub_color   <- "#A0C1B9"

# ---- Add metadata and combine ----
ssa_deseq$drug   <- "Spliceostatin A"
pladb_deseq$drug <- "Pladienolide B"
tub_deseq$drug   <- "Tubercidin"

combined_df <- bind_rows(ssa_deseq, pladb_deseq, tub_deseq)

# ---- Thresholds ----
logfc_thresh <- 1.5
padj_thresh  <- 0.05

# ---- Handle padj == 0 for plotting ----
min_nonzero_padj <- suppressWarnings(min(combined_df$padj[combined_df$padj > 0], na.rm = TRUE))
if (!is.finite(min_nonzero_padj) || is.na(min_nonzero_padj)) {
  min_nonzero_padj <- 1e-300
}

combined_df <- combined_df %>%
  mutate(padj_for_plot = case_when(
    is.na(padj)     ~ NA_real_,
    padj == 0       ~ min_nonzero_padj / 10,
    TRUE            ~ as.numeric(padj)
  ))

# ---- Define significance and plot colors ----
combined_df <- combined_df %>%
  mutate(
    significant = ifelse(!is.na(padj) & padj <= padj_thresh & abs(log2FoldChange) >= logfc_thresh, "sig", "ns"),
    point_color = ifelse(significant == "sig", drug, "ns"),
    negLog10Padj = -log10(padj_for_plot)
  )

# ---- Clean data for plotting ----
plot_df <- combined_df %>%
  filter(!is.na(log2FoldChange) & !is.na(negLog10Padj) & is.finite(log2FoldChange) & is.finite(negLog10Padj))

# ---- Determine labels (Top 30% of sig genes) ----
label_df <- plot_df %>%
  filter(significant == "sig") %>%
  group_by(drug) %>%
  slice_max(order_by = abs(log2FoldChange), prop = 0.3, with_ties = FALSE) %>%
  ungroup()

# ---- Axis limits with padding ----
x_val <- max(abs(plot_df$log2FoldChange), na.rm = TRUE)
x_lim <- c(-x_val - 0.5, x_val + 0.5) # Symmetrical X axis is usually better
y_lim <- c(0, max(plot_df$negLog10Padj, na.rm = TRUE) * 1.05)

# ---- Count sig genes for legend ----
sig_counts <- combined_df %>%
  filter(significant == "sig") %>%
  group_by(drug) %>%
  summarise(n = n(), .groups = "drop")

all_drugs <- c("Spliceostatin A", "Pladienolide B", "Tubercidin")
sig_counts <- sig_counts %>% 
  right_join(tibble(drug = all_drugs), by = "drug") %>%
  mutate(n = ifelse(is.na(n), 0L, n))

# ---- Plot Setup ----
color_map <- c("ns" = "grey80",
               "Spliceostatin A" = ssa_color,
               "Pladienolide B" = pladb_color,
               "Tubercidin" = tub_color)

legend_labels <- c("Not significant",
                   sapply(all_drugs, function(d) {
                     count <- sig_counts$n[sig_counts$drug == d]
                     paste0(d, " (", count, " sig)")
                   }))

# Font handling
showtext_auto()
showtext.opts(dpi=600)

# ---- Build Volcano Plot ----
volcano_plot_all <- ggplot(plot_df, aes(x = log2FoldChange, y = negLog10Padj)) +
  geom_point(aes(color = point_color), size = 1.5, alpha = 0.5, stroke = 0) +
  # Optional: Uncomment to add labels back
  # geom_text_repel(data = label_df, aes(label = symbol), 
  #                 size = 2, max.overlaps = 50, segment.size = 0.2) +
  scale_color_manual(values = color_map, labels = legend_labels, name = NULL) +
  coord_cartesian(xlim = x_lim, ylim = y_lim, expand = FALSE) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 6),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7)
  ) +
  labs(x = "log2 Fold Change", y = "-log10(adj p-value)")

# ---- Save Plot ----
ggsave("Figures/Figure2b_volcano_combined_expression.png", 
       volcano_plot_all, 
       width = 4, height = 4, units = "in", dpi = 600)