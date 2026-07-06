# =============================================================================
# Reviewer response — GO convergence at pathway vs gene level
# Canonical drug colors used across the manuscript (figure2E, figure2B/C, etc.)
#   Tubercidin      #A0C1B9   sage green
#   Pladienolide B  #70A0AF   teal
#   Spliceostatin A #706993   purple
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ---- canonical colors (mirroring Figures/figure2E_upset_genes.R etc.) -------
tub_color  <- "#A0C1B9"
pladb_color <- "#70A0AF"
ssa_color  <- "#706993"
shared2_color <- "#3D3D3D"   # dark grey for "shared by 2 drugs"
shared3_color <- "#000000"   # black    for "shared by all 3 drugs"

set_colors <- c("Tubercidin" = tub_color,
                "Pladienolide B" = pladb_color,
                "Spliceostatin A" = ssa_color)

# ---- load enrichment --------------------------------------------------------
xl_path <- "Supplementary6_GO_Results_Combined.xlsx"
tub   <- read_excel(xl_path, sheet = "Tub_BP")
pladb <- read_excel(xl_path, sheet = "PlaDB_BP")
ssa   <- read_excel(xl_path, sheet = "SSA_BP")

# ---- keep only terms significant in ALL 3 drugs (the convergent terms) -----
common_terms <- Reduce(intersect, list(tub$Description, pladb$Description, ssa$Description))
message("Terms enriched in all 3 drugs: ", length(common_terms))

# ---- assemble per-term gene sets --------------------------------------------
genes_for <- function(df, term) {
  strsplit(df$geneID[df$Description == term][1], "/", fixed = TRUE)[[1]]
}

rows <- lapply(common_terms, function(term) {
  g_tub   <- genes_for(tub,   term)
  g_pladb <- genes_for(pladb, term)
  g_ssa   <- genes_for(ssa,   term)
  union   <- union(union(g_tub, g_pladb), g_ssa)
  shared3 <- intersect(intersect(g_tub, g_pladb), g_ssa)
  only_t  <- setdiff(setdiff(g_tub, g_pladb), g_ssa)
  only_p  <- setdiff(setdiff(g_pladb, g_tub), g_ssa)
  only_s  <- setdiff(setdiff(g_ssa, g_tub), g_pladb)
  shared2 <- length(union) - length(shared3) - length(only_t) - length(only_p) - length(only_s)
  jaccard_tp <- length(intersect(g_tub, g_pladb)) / length(union(g_tub, g_pladb))
  jaccard_ts <- length(intersect(g_tub, g_ssa))   / length(union(g_tub, g_ssa))
  jaccard_ps <- length(intersect(g_pladb, g_ssa)) / length(union(g_pladb, g_ssa))
  data.frame(
    Description = term,
    union_size  = length(union),
    n_tub       = length(g_tub),
    n_pladb     = length(g_pladb),
    n_ssa       = length(g_ssa),
    only_tub    = length(only_t),
    only_pladb  = length(only_p),
    only_ssa    = length(only_s),
    shared2     = shared2,
    shared3     = length(shared3),
    jaccard_tp  = jaccard_tp,
    jaccard_ts  = jaccard_ts,
    jaccard_ps  = jaccard_ps
  )
})
stats <- do.call(rbind, rows) |> arrange(desc(union_size))

# ---- pick the top 12 convergent terms ---------------------------------------
top <- head(stats, 12)
top$Description <- factor(top$Description, levels = rev(top$Description))

# ---- reshape for stacked composition ----------------------------------------
comp <- top |>
  select(Description, only_tub, only_pladb, only_ssa, shared2, shared3) |>
  rename("Tubercidin-unique"      = only_tub,
         "Pladienolide B-unique"  = only_pladb,
         "Spliceostatin A-unique" = only_ssa,
         "Shared by 2 drugs"      = shared2,
         "Shared by all 3"        = shared3) |>
  pivot_longer(-Description, names_to = "Segment", values_to = "n") |>
  group_by(Description) |>
  mutate(pct = n / sum(n) * 100) |>
  ungroup()

segment_levels <- c("Tubercidin-unique", "Pladienolide B-unique",
                    "Spliceostatin A-unique", "Shared by 2 drugs",
                    "Shared by all 3")
comp$Segment <- factor(comp$Segment, levels = segment_levels)

segment_colors <- c("Tubercidin-unique"      = tub_color,
                    "Pladienolide B-unique"  = pladb_color,
                    "Spliceostatin A-unique" = ssa_color,
                    "Shared by 2 drugs"      = shared2_color,
                    "Shared by all 3"        = shared3_color)

# label suffix showing gene count
comp$label <- ifelse(comp$n > 0, paste0(comp$n), "")

# ---- PLOT 1: 100% stacked horizontal bar -----------------------------------
p1 <- ggplot(comp, aes(x = Description, y = pct, fill = Segment)) +
  geom_bar(stat = "identity", width = 0.72, color = "white", linewidth = 0.4) +
  coord_flip() +
  scale_fill_manual(values = segment_colors, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04)),
                     labels = function(x) paste0(x, "%")) +
  geom_text(data = subset(comp, pct >= 7),
            aes(label = label),
            position = position_stack(vjust = 0.5),
            color = "white", size = 3.2, fontface = "bold") +
  labs(
    title    = "Convergent GO terms: distinct genes per drug within the same pathway",
    subtitle = "Top 12 Biological Processes enriched (FDR < 0.05) in ALL 3 treatments",
    x = NULL, y = "Share of contributing genes"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 11.5, hjust = 0,
                                  margin = margin(b = 4)),
    plot.subtitle = element_text(size = 9, color = "grey30", hjust = 0,
                                  margin = margin(b = 10)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    axis.text.y  = element_text(size = 9.5, color = "black"),
    axis.text.x  = element_text(size = 8.5, color = "grey30"),
    legend.position = "bottom",
    legend.key.size = unit(9, "pt"),
    legend.text     = element_text(size = 9),
    legend.margin   = margin(t = 4),
    plot.margin     = margin(8, 8, 8, 8)
  ) +
  guides(fill = guide_legend(nrow = 2, reverse = FALSE))

# ---- PLOT 2: Jaccard heatmap of gene sets within each term ------------------
jac <- top |>
  select(Description, jaccard_tp, jaccard_ts, jaccard_ps) |>
  rename("Tub vs PlaDB" = jaccard_tp,
         "Tub vs SSA"   = jaccard_ts,
         "PlaDB vs SSA" = jaccard_ps) |>
  pivot_longer(-Description, names_to = "Pair", values_to = "Jaccard")

jac$Description <- factor(jac$Description, levels = levels(top$Description))
jac$Pair <- factor(jac$Pair, levels = c("Tub vs PlaDB", "Tub vs SSA", "PlaDB vs SSA"))

# sequential ramp anchored at SSA purple (the repo's deepest categorical)
p2 <- ggplot(jac, aes(x = Pair, y = Description, fill = Jaccard)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", Jaccard)),
            color = ifelse(jac$Jaccard > 0.18, "white", "black"),
            size = 3) +
  scale_fill_gradient(low = "#F2EFF6", high = ssa_color,
                      name = "Jaccard\n(gene overlap)",
                      limits = c(0, max(jac$Jaccard) * 1.05),
                      breaks = pretty_breaks(4)) +
  labs(
    title    = "Gene-level overlap within each shared pathway is limited",
    subtitle = "Jaccard index of contributing gene sets per drug pair (low = distinct genes per drug)",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 11.5, hjust = 0,
                                  margin = margin(b = 4)),
    plot.subtitle = element_text(size = 9, color = "grey30", hjust = 0,
                                  margin = margin(b = 8)),
    panel.grid    = element_blank(),
    axis.text.x   = element_text(size = 9.5, color = "black", face = "bold"),
    axis.text.y   = element_text(size = 9.5, color = "black"),
    legend.position = "right",
    legend.key.height = unit(0.8, "cm"),
    legend.title      = element_text(size = 9),
    legend.text       = element_text(size = 8.5),
    plot.margin       = margin(8, 8, 8, 8)
  )

# ---- combine & save ---------------------------------------------------------
final <- p1 / p2 + plot_layout(heights = c(1.05, 1))

out_dir <- "Figures/go_convergence"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(out_dir, "Fig_convergence_pathway_vs_gene.png"),
       final, width = 11, height = 9.5, dpi = 300)
ggsave(file.path(out_dir, "Fig_convergence_pathway_vs_gene.pdf"),
       final, width = 11, height = 9.5)

# also save individual panels for flexibility
ggsave(file.path(out_dir, "Fig_convergence_A_stacked.png"),
       p1, width = 10.5, height = 5.5, dpi = 300)
ggsave(file.path(out_dir, "Fig_convergence_B_jaccard.png"),
       p2, width = 8, height = 5.5, dpi = 300)

message("Saved to ", out_dir)
