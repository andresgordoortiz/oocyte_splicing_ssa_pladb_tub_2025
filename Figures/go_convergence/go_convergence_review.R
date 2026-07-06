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

# ---- load enrichment --------------------------------------------------------
go_xl  <- "Supplementary6_GO_Results_Combined.xlsx"
spl_xl <- "Supplementary2_betAS_splicing_results.xlsx"

tub_go   <- read_excel(go_xl,  sheet = "Tub_BP")
pladb_go <- read_excel(go_xl,  sheet = "PlaDB_BP")
ssa_go   <- read_excel(go_xl,  sheet = "SSA_BP")

# splicing events (one row per event, with deltapsi and FDR)
tub_spl   <- read_excel(spl_xl, sheet = "Tub_FDR")
pladb_spl <- read_excel(spl_xl, sheet = "PlaDB_FDR")
ssa_spl   <- read_excel(spl_xl, sheet = "SSA_FDR")

# ---- pick the convergent GO terms (shared across all 3 drugs) -------------
common_terms <- Reduce(intersect, list(tub_go$Description, pladb_go$Description, ssa_go$Description))
message("Terms enriched in all 3 drugs: ", length(common_terms))

# ---- assemble per-term gene sets --------------------------------------------
genes_for <- function(df, term) {
  strsplit(df$geneID[df$Description == term][1], "/", fixed = TRUE)[[1]]
}

per_term_stats <- lapply(common_terms, function(term) {
  g_tub   <- genes_for(tub_go,   term)
  g_pladb <- genes_for(pladb_go, term)
  g_ssa   <- genes_for(ssa_go,   term)
  union   <- union(union(g_tub, g_pladb), g_ssa)
  shared3 <- intersect(intersect(g_tub, g_pladb), g_ssa)
  only_t  <- setdiff(setdiff(g_tub, g_pladb), g_ssa)
  only_p  <- setdiff(setdiff(g_pladb, g_tub), g_ssa)
  only_s  <- setdiff(setdiff(g_ssa, g_tub), g_pladb)
  shared2 <- length(union) - length(shared3) - length(only_t) - length(only_p) - length(only_s)
  data.frame(
    Description = term, union_size = length(union),
    n_tub = length(g_tub), n_pladb = length(g_pladb), n_ssa = length(g_ssa),
    only_tub = length(only_t), only_pladb = length(only_p), only_ssa = length(only_s),
    shared2 = shared2, shared3 = length(shared3)
  )
})
stats <- do.call(rbind, per_term_stats) |> arrange(desc(union_size))

# ---- keep top 12 convergent terms (original, untouched) -------------------
TOP_N_TERMS <- 12
top <- head(stats, TOP_N_TERMS)
top$Description <- factor(top$Description, levels = rev(top$Description))

# ---- compute composition stats for the stacked bar (PLOT 1) ----------------
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
                    "Spliceostatin A-unique",
                    "Shared by 2 drugs",
                    "Shared by all 3")
comp$Segment <- factor(comp$Segment, levels = segment_levels)
segment_colors <- c("Tubercidin-unique"      = tub_color,
                    "Pladienolide B-unique"  = pladb_color,
                    "Spliceostatin A-unique" = ssa_color,
                    "Shared by 2 drugs"      = shared2_color,
                    "Shared by all 3"        = shared3_color)
comp$label <- ifelse(comp$n > 0, as.character(comp$n), "")

# ---- PLOT 1 (unchanged design) — pathway-level convergence -----------------------
p1 <- ggplot(comp, aes(x = Description, y = pct, fill = Segment)) +
  geom_bar(stat = "identity", width = 0.72, color = "white", linewidth = 0.4) +
  coord_flip() +
  scale_fill_manual(values = segment_colors, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04)),
                     labels = function(x) paste0(x, "%")) +
  scale_x_discrete(expand = expansion(add = 0)) +
  geom_text(data = subset(comp, pct >= 7),
            aes(label = label),
            position = position_stack(vjust = 0.5),
            color = "white", size = 3.6, fontface = "bold") +
  labs(
    title    = "Convergent GO terms: distinct genes per drug within the same pathway",
    subtitle = "Top 12 Biological Processes enriched (FDR < 0.05) in ALL 3 treatments",
    x = NULL, y = "Share of contributing genes"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 12.5, hjust = 0,
                                  margin = margin(b = 4)),
    plot.subtitle = element_text(size = 9.5, color = "grey30", hjust = 0,
                                  margin = margin(b = 10)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    axis.text.y  = element_text(size = 10, color = "black", hjust = 1, margin = margin(r = 0)),
    axis.text.x  = element_text(size = 9, color = "grey30"),
    legend.position = "bottom",
    legend.key.size = unit(9, "pt"),
    legend.text     = element_text(size = 9.5),
    legend.margin   = margin(t = 4),
    plot.margin     = margin(8, 8, 8, 0)
  ) +
  guides(fill = guide_legend(nrow = 2, reverse = FALSE))

# ============================================================================
# PLOT 2 — simple gene x drug presence matrix, organised in 3 thematic blocks
#   Block 1: Cell cycle & division (merged cell cycle / nuclear / microtubule terms)
#   Block 2: DNA repair & recombination (merged DNA-repair terms)
#   Block 3: Small GTPase signalling (merged small-GTPase terms)
#   For each block, show the union of contributing genes from all merged
#   terms, then display which drugs each gene is significant in. Genes are
#   grouped by overlap category (shared by all 3, shared by 2, then unique).
#   This is the simple answer to "which genes overlap" — no Venn, no UpSet.
# ============================================================================

# define the 3 thematic blocks
blocks <- list(
  "Cell cycle & division" = c(
    "regulation of cell cycle phase transition",
    "nuclear division",
    "negative regulation of cell cycle",
    "negative regulation of cell cycle process",
    "negative regulation of cell cycle phase transition",
    "regulation of microtubule-based process"
  ),
  "DNA repair & recombination" = c(
    "double-strand break repair",
    "DNA recombination",
    "regulation of DNA repair"
  ),
  "Small GTPase signalling" = c(
    "small GTPase-mediated signal transduction",
    "regulation of small GTPase mediated signal transduction"
  )
)

# cap the number of genes per block so the matrix is readable
N_SHARED3 <- 8    # top shared-by-all-3 per block
N_SHARED2 <- 5    # top shared-by-2 per block
N_UNIQUE  <- 0    # skip drug-unique (clutters; the question is about overlap)

# build per-block gene tables
block_rows <- list()
for (block_name in names(blocks)) {
  block_terms <- blocks[[block_name]]
  g_tub   <- unique(unlist(lapply(block_terms, function(t) genes_for(tub_go,   t))))
  g_pladb <- unique(unlist(lapply(block_terms, function(t) genes_for(pladb_go, t))))
  g_ssa   <- unique(unlist(lapply(block_terms, function(t) genes_for(ssa_go,   t))))
  union_g <- union(union(g_tub, g_pladb), g_ssa)
  if (length(union_g) == 0) next

  in_t <- union_g %in% g_tub
  in_p <- union_g %in% g_pladb
  in_s <- union_g %in% g_ssa
  n    <- in_t + in_p + in_s

  gene_tbl <- data.frame(
    Block   = block_name,
    GENE    = union_g,
    in_tub   = in_t,
    in_pladb = in_p,
    in_ssa   = in_s,
    stringsAsFactors = FALSE
  ) |>
    mutate(Cat = case_when(
      n == 3 ~ "Shared by all 3",
      n == 2 & in_t & in_p ~ "Shared by Tub & PlaDB",
      n == 2 & in_t & in_s ~ "Shared by Tub & SSA",
      n == 2 & in_p & in_s ~ "Shared by PlaDB & SSA",
      in_t ~ "Tubercidin only",
      in_p ~ "Pladienolide B only",
      TRUE  ~ "Spliceostatin A only"
    ))

  # limit per category
  keep_cats <- c("Shared by all 3"      = N_SHARED3,
                 "Shared by Tub & PlaDB" = N_SHARED2,
                 "Shared by Tub & SSA"   = N_SHARED2,
                 "Shared by PlaDB & SSA" = N_SHARED2,
                 "Tubercidin only"        = N_UNIQUE,
                 "Pladienolide B only"    = N_UNIQUE,
                 "Spliceostatin A only"   = N_UNIQUE)
  gene_tbl <- gene_tbl |>
    filter(Cat %in% names(keep_cats)[keep_cats > 0]) |>
    group_by(Cat) |>
    group_modify(~ head(.x, keep_cats[as.character(.y$Cat)[1]])) |>
    ungroup()

  # long form
  long <- bind_rows(
    gene_tbl |> filter(in_tub)   |> mutate(Drug = "Tubercidin"),
    gene_tbl |> filter(in_pladb) |> mutate(Drug = "Pladienolide B"),
    gene_tbl |> filter(in_ssa)   |> mutate(Drug = "Spliceostatin A")
  ) |>
    select(Block, GENE, Drug, Cat)

  block_rows[[length(block_rows) + 1]] <- long
}

p2_long <- do.call(rbind, block_rows) |>
  mutate(
    Drug = factor(Drug, levels = c("Tubercidin", "Pladienolide B", "Spliceostatin A")),
    Block = factor(Block, levels = rev(names(blocks))),
    # ordering: shared-by-2 first (any pair), then shared-by-3
    Cat = factor(as.character(Cat),
                 levels = c("Shared by Tub & PlaDB",
                            "Shared by Tub & SSA",
                            "Shared by PlaDB & SSA",
                            "Shared by all 3",
                            "Tubercidin only",
                            "Pladienolide B only",
                            "Spliceostatin A only"))
  ) |>
  arrange(Block, Cat, GENE)

# one row per (block, gene)
p2_summary <- p2_long |>
  group_by(Block, GENE, Cat) |>
  summarise(
    in_tub   = any(Drug == "Tubercidin"),
    in_pladb = any(Drug == "Pladienolide B"),
    in_ssa   = any(Drug == "Spliceostatin A"),
    .groups  = "drop"
  ) |>
  arrange(Block, Cat, GENE)

# row label = gene name only
p2_summary$RowLabel <- factor(p2_summary$GENE, levels = unique(p2_summary$GENE))

# full matrix
mat_full <- p2_summary |>
  pivot_longer(c(in_tub, in_pladb, in_ssa),
               names_to = "DrugCol", values_to = "InDrug") |>
  mutate(
    Drug = recode(DrugCol,
                  "in_tub"   = "Tubercidin",
                  "in_pladb" = "Pladienolide B",
                  "in_ssa"   = "Spliceostatin A"),
    Drug = factor(Drug, levels = c("Tubercidin", "Pladienolide B", "Spliceostatin A")),
    RowLabel = factor(GENE, levels = levels(p2_summary$RowLabel)),
    Block = factor(Block, levels = rev(names(blocks)))
  )

# colors per drug
drug_dot_colors <- c("Tubercidin"      = tub_color,
                     "Pladienolide B"  = pladb_color,
                     "Spliceostatin A" = ssa_color)

p2 <- ggplot(mat_full, aes(x = Drug, y = RowLabel)) +
  # background grey track
  geom_point(shape = 21, fill = NA, color = "grey85", size = 5, stroke = 0.3) +
  # drug-coloured dots
  geom_point(data = subset(mat_full, InDrug),
             aes(fill = Drug), shape = 21, color = "white",
             size = 5, stroke = 1.0) +
  scale_fill_manual(values = drug_dot_colors, guide = "none") +
  scale_x_discrete(position = "top") +
  facet_grid(Block ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  labs(
    title    = "Which genes overlap across drugs, per convergent pathway block",
    subtitle = "Filled dot in the drug's colour = gene is a significant contributor in that drug's enrichment.\nRows within each block are grouped by overlap category (shared-by-3 first, then shared-by-2)."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12.5, hjust = 0, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 8.5, color = "grey30", hjust = 0, margin = margin(b = 8)),
    panel.grid    = element_blank(),
    axis.text.x.top   = element_text(size = 10, color = "black", face = "bold"),
    axis.text.y   = element_text(size = 8, color = "black", hjust = 1),
    strip.text.y.left = element_text(size = 9.5, face = "bold", angle = 0, hjust = 1),
    strip.background.y = element_rect(fill = "grey95", color = NA),
    plot.margin       = margin(8, 8, 8, 8),
    panel.spacing.y   = unit(0.3, "lines")
  )

# ---- combine & save ---------------------------------------------------------
final <- p1 / p2 + plot_layout(heights = c(0.7, 1.6))

out_dir <- "Figures/go_convergence"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(out_dir, "Fig_convergence_pathway_vs_gene.png"),
       final, width = 14, height = 15, dpi = 300)
ggsave(file.path(out_dir, "Fig_convergence_pathway_vs_gene.pdf"),
       final, width = 14, height = 15)

# also save individual panels
ggsave(file.path(out_dir, "Fig_convergence_A_stacked.png"),
       p1, width = 11, height = 5, dpi = 300)
ggsave(file.path(out_dir, "Fig_convergence_B_gene_membership.png"),
       p2, width = 12, height = 13, dpi = 300)

message("Saved to ", out_dir)
