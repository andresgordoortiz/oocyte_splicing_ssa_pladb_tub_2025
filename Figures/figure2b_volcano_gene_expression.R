library(ggplot2)
library(ggrepel)
library(dplyr)
library(scales)
library(showtext)
# ---- If you previously used showtext, disable/remove those calls.
# ---- We will rely on cairo_pdf + family to produce vector text in the PDF.

# ---- Drug colors ----
ssa_color <- "#706993"
pladb_color <- "#70A0AF"
tub_color  <- "#A0C1B9"

# ---- Add a column to each df to track drug/source ----
ssa_deseq$drug <- "Spliceostatin A"
pladb_deseq$drug <- "Pladienolide B"
tub_deseq$drug  <- "Tubercidin"

# ---- Combine all dfs ----
combined_df <- bind_rows(ssa_deseq, pladb_deseq, tub_deseq)

# ---- thresholds ----
logfc_thresh <- 1.5
padj_thresh  <- 0.05

# ---- Protect against padj == 0: create padj_for_plot ----
min_nonzero_padj <- suppressWarnings(min(combined_df$padj[combined_df$padj > 0], na.rm = TRUE))
if (!is.finite(min_nonzero_padj) || is.na(min_nonzero_padj)) {
  min_nonzero_padj <- 1e-300
}
combined_df <- combined_df %>%
  mutate(padj_for_plot = case_when(
    is.na(padj)           ~ NA_real_,
    padj == 0             ~ min_nonzero_padj / 10,
    TRUE                  ~ as.numeric(padj)
  ))

# ---- Define significance using padj (but use padj_for_plot for plotting) ----
combined_df <- combined_df %>%
  mutate(significant = ifelse(!is.na(padj) & padj <= padj_thresh & abs(log2FoldChange) >= logfc_thresh, "sig", "ns"))

# ---- Create point color column: drug name for sig, "ns" for non-sig ----
combined_df <- combined_df %>%
  mutate(point_color = ifelse(significant == "sig", drug, "ns"),
         negLog10Padj = -log10(padj_for_plot))

# ---- Remove rows that cannot be plotted (NA log2FC or NA padj_for_plot or non-finite negLog10Padj) ----
initial_n <- nrow(combined_df)
plot_df <- combined_df %>%
  filter(!is.na(log2FoldChange) & !is.na(negLog10Padj) & is.finite(log2FoldChange) & is.finite(negLog10Padj))
removed_n <- initial_n - nrow(plot_df)
message("Removed ", removed_n, " rows before plotting (NA or non-finite log2FoldChange/padj).")

# ---- Determine labels (top 30% of significant by abs log2FC per drug) ----
label_df <- plot_df %>%
  filter(significant == "sig") %>%
  slice_max(order_by = abs(log2FoldChange), prop = 0.3, with_ties = FALSE) %>%
  ungroup()

# ---- Common axis limits (finite) and a small padding so points at exact limits aren't visually cut ----
x_lim <- range(plot_df$log2FoldChange, na.rm = TRUE)
x_pad <- diff(x_lim) * 0.03
if (!is.finite(x_pad) || x_pad == 0) x_pad <- 0.5
x_lim <- c(x_lim[1] - x_pad, x_lim[2] + x_pad)

y_lim <- c(0, max(plot_df$negLog10Padj, na.rm = TRUE) * 1.05)
if (!is.finite(y_lim[2]) || y_lim[2] <= 0) {
  y_lim <- c(0, 1)
}

# ---- Count significant genes per drug for legend labels ----
sig_counts <- combined_df %>%
  filter(significant == "sig") %>%
  group_by(drug) %>%
  summarise(n = n(), .groups = "drop")

all_drugs <- unique(combined_df$drug)
# ensure all drugs present
sig_counts <- sig_counts %>% right_join(tibble(drug = all_drugs), by = "drug") %>%
  mutate(n = ifelse(is.na(n), 0L, n))

# ---- Build the volcano plot ----
color_map <- c("ns" = "grey80",
               "Spliceostatin A" = ssa_color,
               "Pladienolide B" = pladb_color,
               "Tubercidin" = tub_color)

# prepare legend labels: first "ns", then each drug with count
legend_labels <- c("Not significant",
                   sapply(all_drugs, function(d) {
                     count <- sig_counts$n[sig_counts$drug == d]
                     paste0(d, " (", count, " sig)")
                   }))

font_add(family = "Arial", regular = "/usr/share/fonts/liberation/LiberationSans-Regular.ttf")
showtext::showtext_opts(dpi = 600)
showtext_auto()

volcano_plot_all <- ggplot(plot_df, aes(x = log2FoldChange, y = negLog10Padj)) +
  geom_point(aes(color = point_color), size = 2.5, alpha = 0.6, stroke = 0, show.legend = TRUE) +
  # label only selected genes with repel
  #geom_text_repel(data = label_df, aes(label = symbol),
   #               size = 2.5, box.padding = 0.35, point.padding = 0.3,
    #              segment.color = "grey50", segment.size = 0.4,
     #             max.overlaps = 200) +
  scale_color_manual(values = color_map, labels = legend_labels, name = NULL) +
  coord_cartesian(xlim = x_lim, ylim = y_lim, expand = FALSE, clip = "off") +
  theme_classic(base_size = 14, base_family = "Liberation Sans") +
  theme(
    axis.line = element_line(linewidth = 0.9, colour = "#222222"),
    axis.ticks = element_line(linewidth = 0.9, colour = "#222222"),
    axis.text = element_text(size = rel(1.0), colour = "#111111"),
    axis.title = element_text(face = "plain", size = rel(1.0)),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 3),
    legend.key = element_rect(fill = NA, colour = NA)
    ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)),
         size = "none", alpha = "none") +
  labs(
       x = "log2 Fold Change",
       y = "-log10(adjusted p-value)")

ggsave("volcano_combined_expression.png", volcano_plot_all, device = png, width = 3, height = 3, units = "in", dpi = 600)
