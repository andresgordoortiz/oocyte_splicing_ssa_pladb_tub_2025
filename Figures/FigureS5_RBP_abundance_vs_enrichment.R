# =============================================================================
# Figure S5: RBP abundance (oocyte RNA-seq) vs motif-enrichment importance
#
# Reviewer point: "check if those proteins [RBPs with higher frequency of
# putative binding sites] are expressed in oocytes, at least at the mRNA level".
#
# Inputs
#   * Supplementary4_rMaps2_out/{drug}_{ex,int}_rmaps/pVal.{up,dn}.vs.bg.RNAmap.txt
#       - pVal.up = "Included" (more inclusion in treatment)
#       - pVal.dn = "Skipped" (more skipping in treatment)
#   * Supplementary1_DESeq2_Results.xlsx
#       - mean_control column = expression in DMSO (oocytes)
#
# For each (drug, direction, event_type) we count, per RBP motif, the number
# of metagene regions (R1..R5 for introns, R1..R8 for exons) where the motif
# p-value is <= 0.05. We then join with the DESeq2 mean_control of the
# corresponding RBP gene to answer: are the RBPs whose motifs are enriched
# in differentially spliced regions actually expressed in oocytes?
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(readxl)
  library(stringr)
  library(patchwork)
})

set.seed(43)

# -----------------------------------------------------------------------------
# Alias map: rMAPS2 RBP labels that are not official MGI gene symbols
# (mostly CISBP / older aliases used in the rMAPS2 motif database).
# -----------------------------------------------------------------------------
rbp_alias <- c(
  "9G8"       = "RBM8A",      # 9G8 = RBM8A (also: alias of CIRBP elsewhere)
  "HNRPLL"    = "HNRNPLL",    # HNRPLL is the legacy / shortened form
  "PTB"       = "PTBP1",      # PTB = PTBP1
  "SF2-ASF"   = "SRSF1",      # SF2/ASF = SRSF1
  "SRp20"     = "SRSF3",      # SRp20 = SRSF3
  "SRp40"     = "SRSF5",      # SRp40 = SRSF5
  "SRp55"     = "SRSF6",      # SRp55 = SRSF6
  "Tra2-beta" = "TRA2B",      # Tra2-beta = TRA2B
  "sf3b1"     = "SF3B1",      # case fix
  "HuR"       = "ELAVL1",     # HuR = ELAVL1
  "RBP"       = "RBPC"        # very ambiguous; map to RBP-cysteine (rare) — will mostly be NA in DESeq2
)

# -----------------------------------------------------------------------------
# 1. Parse rMAPS2 outputs
# -----------------------------------------------------------------------------
# rMAPS2 file naming:
#   pVal.up.vs.bg.RNAmap.txt -> "Included" (PSI goes up)
#   pVal.dn.vs.bg.RNAmap.txt -> "Skipped" (PSI goes down)
#
# Column layout:
#   * exon files: R1..R8 (8 metagene regions around the target exon)
#   * intron files: R1..R5 (5 metagene regions around the intron)
parse_rmaps <- function(drug, event_type, direction) {
  # direction in {"up","dn"} (rMAPS2 naming); event_type in {"ex","int"}
  if (!(event_type %in% c("ex", "int"))) {
    stop("event_type must be 'ex' or 'int'")
  }
  if (!(direction %in% c("up", "dn"))) {
    stop("direction must be 'up' or 'dn'")
  }

  f <- file.path(
    "Supplementary4_rMaps2_out",
    paste0(drug, "_", event_type, "_rmaps"),
    paste0("pVal.", direction, ".vs.bg.RNAmap.txt")
  )
  if (!file.exists(f)) stop("Missing rMAPS2 file: ", f)

  df <- read.delim(f, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE)

  # The exon files have descriptive column names; the intron files use
  # R1..R5 directly. Normalize to R1..Rn.
  if (event_type == "ex") {
    r_actual <- c(
      "smallest_p_in_upstreamExon-3prime",
      "smallest_p_in_upstreamExonIntron",
      "smallest_p_in_upstreamIntron",
      "smallest_p_in_targetExon-5prime",
      "smallest_p_in_targetExon-3prime",
      "smallest_p_in_downstreamIntron",
      "smallest_p_in_downstreamExonIntron",
      "smallest_p_in_downstreamExon-5prime"
    )
    rename_vec <- setNames(paste0("R", seq_along(r_actual)), r_actual)
    df <- df %>% rename_with(~ unname(rename_vec[.x]), .cols = all_of(r_actual))
  }

  r_cols <- grep("^R\\d+$", names(df), value = TRUE)
  stopifnot(length(r_cols) > 0)

  df %>%
    transmute(
      RBP_full   = RBP,
      RBP        = sub("\\..*", "", RBP_full),            # gene symbol
      motif      = sub("^[^.]+\\.", "", RBP_full),        # motif part
      drug       = drug,
      event_type = event_type,
      direction  = ifelse(direction == "up", "Included", "Skipped"),
      best_p     = suppressWarnings(apply(.[, r_cols, drop = FALSE], 1, min,
                                          na.rm = TRUE)),
      n_sig_regions = rowSums(.[, r_cols, drop = FALSE] <= 0.05, na.rm = TRUE),
      # Per-bin minimum p across the row's motifs (used to compute the
      # union of significant bins per (RBP, drug, event, direction) later).
      min_p_R1 = .[, "R1"], min_p_R2 = .[, "R2"],
      min_p_R3 = .[, "R3"], min_p_R4 = .[, "R4"],
      min_p_R5 = .[, "R5"],
      min_p_R6 = if ("R6" %in% names(df)) .[, "R6"] else NA_real_,
      min_p_R7 = if ("R7" %in% names(df)) .[, "R7"] else NA_real_,
      min_p_R8 = if ("R8" %in% names(df)) .[, "R8"] else NA_real_
    ) %>%
    mutate(
      best_p = ifelse(is.finite(best_p) & best_p > 0, best_p, NA_real_),
      neglog10p = -log10(best_p),
      # Map non-canonical rMAPS2 labels to official MGI symbols
      RBP_official = ifelse(RBP %in% names(rbp_alias),
                            unname(rbp_alias[RBP]), RBP),
      drug_label = factor(drug,
                          levels = c("pladb", "ssa", "tub"),
                          labels = c("Pladienolide B",
                                     "Spliceostatin A",
                                     "Tubercidin")),
      event_label = factor(event_type,
                           levels = c("int", "ex"),
                           labels = c("Intron", "Exon"))
    )
}

cat("Loading rMAPS2 results ...\n")
rmaps_all <- bind_rows(
  expand_grid(drug = c("pladb", "ssa", "tub"),
              event_type = c("ex", "int"),
              direction = c("up", "dn")) %>%
    purrr::pmap_dfr(function(drug, event_type, direction) {
      parse_rmaps(drug, event_type, direction)
    })
)
cat("  ->", nrow(rmaps_all), "RBP-motif rows across",
    length(unique(rmaps_all$drug)), "drugs,",
    length(unique(rmaps_all$event_type)), "event types,",
    length(unique(rmaps_all$direction)), "directions\n")

# -----------------------------------------------------------------------------
# 2. Load DESeq2 results (Shrunken, for log2FC stability)
# -----------------------------------------------------------------------------
cat("Loading DESeq2 results ...\n")
deseq_files <- list(
  "Pladienolide B"  = c(sheet = "PlaDB_Shrunken", drug = "pladb"),
  "Spliceostatin A" = c(sheet = "SSA_Shrunken",   drug = "ssa"),
  "Tubercidin"      = c(sheet = "Tub_Shrunken",   drug = "tub")
)

deseq_all <- bind_rows(lapply(names(deseq_files), function(label) {
  info <- deseq_files[[label]]
  read_xlsx("Supplementary1_DESeq2_Results.xlsx", sheet = info[["sheet"]]) %>%
    transmute(
      RBP_candidate = toupper(external_gene_name),
      mgi_symbol    = mgi_symbol,
      baseMean      = baseMean,
      log2FC        = log2FoldChange,
      padj          = padj,
      mean_treatment = mean_treatment,
      mean_control  = mean_control,
      drug          = info[["drug"]],
      drug_label    = label
    )
}))
cat("  ->", nrow(deseq_all), "genes in DESeq2 tables\n")

# -----------------------------------------------------------------------------
# 3. Collapse rMAPS2 to one row per (RBP gene, drug, event_type, direction)
# -----------------------------------------------------------------------------
# For each RBP gene we compute:
#   * best_p: the most significant p-value across all motifs and bins
#   * n_sig_bins: the number of DISTINCT metagene bins (R1..R5 or R1..R8)
#     where at least one motif of this RBP has p <= 0.05. This is the
#     "breadth of enrichment" — bounded by 5 (intron) or 8 (exon) per RBP
#     regardless of how many motifs the RBP has in the database.
#   * n_motifs_sig: how many distinct motifs of this RBP are significant.
bin_cols <- grep("^min_p_R\\d+$", names(rmaps_all), value = TRUE)
rmaps_by_gene <- rmaps_all %>%
  group_by(RBP, RBP_official, drug, drug_label, event_type, event_label, direction) %>%
  summarise(
    best_p        = suppressWarnings(min(best_p, na.rm = TRUE)),
    best_p        = ifelse(is.finite(best_p) & best_p > 0, best_p, NA_real_),
    best_neglog10p = ifelse(is.na(best_p), NA_real_, -log10(best_p)),
    n_motifs      = n_distinct(motif),
    n_motifs_sig  = n_distinct(motif[!is.na(best_p) & best_p <= 0.05]),
    top_motif     = motif[which.max(neglog10p)][1],
    .groups = "drop"
  )

# Per-(RBP, drug, event, direction) bin-level statistics.
# For each metagene bin, we take the minimum p-value across all motifs of
# the RBP (the best evidence the RBP has at that position). We then
#   * n_sig_bins: how many of those minima are <= 0.05
#   * raw_score:   sum of -log10(min p) over the significant bins (this
#                  is biased by n_motifs: an RBP with more motifs in the
#                  CISBP database has more chances to hit a significant
#                  p, so the raw score is inflated for "well-studied" RBPs)
#   * enrichment_score: raw_score / n_motifs  -> a per-motif average.
#                       This is the bias-corrected "importance" metric.
n_motifs_map <- rmaps_all %>%
  group_by(RBP) %>%
  summarise(n_motifs = n_distinct(motif), .groups = "drop")

sig_bin_stats <- rmaps_all %>%
  group_by(RBP, drug, event_type, direction) %>%
  summarise(across(all_of(bin_cols), ~ min(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  left_join(n_motifs_map, by = "RBP") %>%
  rowwise() %>%
  mutate(
    pmat = list(c_across(all_of(bin_cols))),
    n_sig_bins = {
      pmin <- pmat
      pmin[is.infinite(pmin)] <- NA_real_
      sum(pmin <= 0.05, na.rm = TRUE)
    },
    raw_score = {
      pmin <- pmat
      pmin[is.infinite(pmin)] <- NA_real_
      contrib <- ifelse(pmin <= 0.05, -log10(pmax(pmin, 1e-300)), 0)
      sum(contrib, na.rm = TRUE)
    },
    # Bias-corrected: per-motif average. An RBP with 1 motif that hits
    # one bin at p=1e-5 scores 5. An RBP with 5 motifs that hits the
    # same bin at the same p (and 4 other motifs at p>0.5) also scores 5,
    # because we ask: "how strong is the BEST motif evidence per bin,
    # averaged over the RBP's database footprint?".
    enrichment_score = raw_score / n_motifs
  ) %>%
  ungroup() %>%
  select(RBP, drug, event_type, direction,
         n_sig_bins, raw_score,
         n_motifs_total = n_motifs, enrichment_score)

rmaps_by_gene <- rmaps_by_gene %>%
  # Drop the old n_motifs before joining to avoid a name collision
  select(-n_motifs) %>%
  left_join(sig_bin_stats,
            by = c("RBP", "drug", "event_type", "direction")) %>%
  mutate(
    RBP_candidate = toupper(RBP_official),
    is_sig        = n_sig_bins > 0
  )

# -----------------------------------------------------------------------------
# 4. Join with DESeq2 mean_control expression
# -----------------------------------------------------------------------------
joined <- rmaps_by_gene %>%
  left_join(
    deseq_all %>% select(RBP_candidate, drug, mean_control, baseMean,
                         log2FC, padj, mgi_symbol),
    by = c("RBP_candidate", "drug")
  ) %>%
  mutate(
    expressed_in_oocytes = !is.na(mean_control) & mean_control > 0,
    log_mean_control    = log10(pmax(mean_control, 1e-3))
  )

cat("RBPs in rMAPS2 that match a DESeq2 gene (any drug):",
    length(unique(joined$RBP[!is.na(joined$mean_control)])), "\n")
cat("RBPs with significant motif enrichment (any drug/dir/event):",
    sum(joined$is_sig), "of", nrow(joined), "rows\n")

# Save summary table for the supplementary material
out_table <- joined %>%
  arrange(drug_label, event_label, direction, desc(enrichment_score), desc(best_neglog10p)) %>%
  select(drug_label, event_label, direction, RBP, RBP_official, mgi_symbol,
         top_motif, enrichment_score, raw_score, n_motifs_total, n_sig_bins, n_motifs_sig,
         best_p, best_neglog10p,
         mean_control, baseMean, log2FC, padj)
dir.create("Figures/FigureS5_RBP_abundance", showWarnings = FALSE,
           recursive = TRUE)
write.csv(out_table,
          "Figures/FigureS5_RBP_abundance/rbp_abundance_vs_enrichment.csv",
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Plot: RBP abundance (mean_control, oocytes) vs motif-enrichment breadth
# -----------------------------------------------------------------------------
# Use only rows with a DESeq2 match (i.e. an annotated RBP gene in oocytes)
plot_df <- joined %>%
  filter(!is.na(mean_control)) %>%
  mutate(
    # Force UPPERCASE display labels
    RBP = toupper(RBP),
    drug_label = factor(as.character(drug_label),
                        levels = c("Pladienolide B", "Spliceostatin A",
                                   "Tubercidin")),
    direction = factor(direction, levels = c("Included", "Skipped")),
    # Make the facet labels explicit: the rMAPS2 metagene is the
    # differential event PLUS its flanking sequences (50 bp flanking exons
    # for introns; 50-250 bp flanking introns for exons). Don't just say
    # "Intron"/"Exon" - the metagene covers more than the event itself.
    event_label = factor(as.character(event_label),
                         levels = c("Intron", "Exon"),
                         labels = c("Differential introns\n(+ flanking exons)",
                                    "Differential exons\n(+ flanking introns)"))
  )

# Pick RBPs to label per panel:
#  - the top N by enrichment_score (combined across Included/Skipped, one
#    label per RBP), AND
#  - any "anchor" RBPs (e.g. SF3B1) that are present in the panel. Anchors
#    keep their direction, so an anchor that shows up in BOTH Included and
#    Skipped is labelled in both directions.
label_top <- function(d, n = 6, anchors = c("SF3B1")) {
  top <- d %>%
    filter(enrichment_score > 0) %>%
    group_by(RBP) %>%
    summarise(
      score_max    = max(enrichment_score, na.rm = TRUE),
      mean_control = mean_control[which.max(enrichment_score)][1],
      direction    = direction[which.max(enrichment_score)][1],
      .groups = "drop"
    ) %>%
    slice_max(order_by = score_max, n = n, with_ties = FALSE) %>%
    ungroup()

  # For anchors: keep ALL (RBP, direction) rows that are present in the
  # panel, so the same anchor can be labelled twice if it appears in both
  # Included and Skipped. This is important for SF3B1, which is the target
  # of PladB/SSA and so lands in both directions.
  anchor_rows <- d %>%
    filter(RBP %in% anchors, enrichment_score > 0) %>%
    transmute(
      RBP,
      direction,
      score_max    = enrichment_score,
      mean_control = mean_control
    ) %>%
    anti_join(top %>% select(RBP, direction), by = c("RBP", "direction"))

  bind_rows(top, anchor_rows)
}
labels_df <- plot_df %>%
  group_by(drug_label, event_label) %>%
  group_modify(~ label_top(.x, n = 6)) %>%
  ungroup() %>%
  select(drug_label, event_label, RBP, direction,
         mean_control, enrichment_score = score_max)

# A "rug" layer just below the y=0 axis to remind the reader that
# the RBPs at y=0 are not-significant (so they don't visually crowd the panel)
rug_df <- plot_df %>% filter(enrichment_score == 0)

# Add an "abundance score" used to fade low-expressed RBPs. The score is
# log10(mean_control + 1) rescaled to [0, 1] within the full plot, so points
# in the bottom expression quartile are translucent and well-expressed RBPs
# are nearly opaque.
plot_df <- plot_df %>%
  mutate(
    log_abund = log10(mean_control + 1),
    abund_score = scales::rescale(log_abund, to = c(0.35, 0.95))
  )
rug_df <- rug_df %>%
  mutate(
    log_abund = log10(mean_control + 1),
    abund_score = scales::rescale(log_abund, to = c(0.35, 0.95))
  )
labels_df <- labels_df %>%
  left_join(
    plot_df %>% select(drug_label, event_label, RBP, direction, abund_score),
    by = c("drug_label", "event_label", "RBP", "direction")
  )

# Categorical point-size: significant RBPs are large and bolded; non-sig are
# small and translucent. This pulls the eye to the enriched RBPs.
# y-axis = continuous enrichment score (sum of -log10(p) over significant bins)
# x-axis = oocyte mRNA abundance
# point size = n_sig_bins (breadth; secondary visual encoding)
# opacity = oocyte abundance (low-abundance RBPs fade out)
p_scatter <- ggplot(plot_df,
                    aes(x = mean_control, y = enrichment_score)) +
  # Non-significant RBPs as a faint rug just below y=0 (faded by abundance)
  geom_rug(
    data = rug_df,
    aes(x = mean_control),
    sides = "b", outside = TRUE,
    color = "grey70",
    alpha = 0.35 * rug_df$abund_score,
    linewidth = 0.2,
    length = unit(0.03, "npc")
  ) +
  # Non-significant RBPs (enrichment_score == 0) as small open circles on
  # the x-axis line so they still appear in the panel, faded by abundance.
  # Small x-jitter prevents overplotting.
  geom_jitter(
    data = plot_df %>% filter(enrichment_score == 0),
    aes(fill = direction, alpha = abund_score),
    shape = 21, color = "grey55", size = 1.2, stroke = 0.15,
    width = 0.04, height = 0
  ) +
  # Significant RBPs: y = enrichment strength, size = breadth, alpha = abundance
  geom_jitter(
    data = plot_df %>% filter(enrichment_score > 0),
    aes(fill = direction, size = n_sig_bins, alpha = abund_score),
    shape = 21, color = "black", stroke = 0.2,
    width = 0.04, height = 0.04
  ) +
  # Top RBPs labels (white halo + colored text)
  geom_label_repel(
    data = labels_df,
    aes(label = RBP, color = direction),
    size = 2.7, fontface = "bold",
    fill = "white", label.size = 0,
    label.padding = unit(0.1, "lines"),
    max.overlaps = 40, min.segment.length = 0,
    segment.size = 0.3, segment.alpha = 0.6,
    box.padding = 0.3, point.padding = 0.15,
    force = 4, seed = 42
  ) +
  facet_grid(event_label ~ drug_label, scales = "free") +
  scale_x_log10(
    labels = function(x) {
      ifelse(x < 1000, sprintf("%d", x),
             format(x, big.mark = ",", scientific = FALSE))
    },
    breaks = c(1, 10, 100, 1000, 10000)
  ) +
  # Fixed y range so panels are comparable, with a small lower margin for
  # the rug ticks and a generous top margin so labels + reference text don't
  # collide
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 6),
    expand = expansion(mult = c(0.05, 0.18))
  ) +
  scale_size_continuous(range = c(1.4, 3.8), guide = "none") +
  scale_alpha_identity() +
  scale_fill_manual(
    values = c("Included" = "#D4A574", "Skipped" = "#B85450"),
    name = "Splicing event"
  ) +
  scale_color_manual(
    values = c("Included" = "#8A5A2B", "Skipped" = "#7A1F1B"),
    guide = "none"
  ) +
  # Allow the rug below the axis
  coord_cartesian(clip = "off") +
  labs(
    x = "Mean normalized counts in control oocytes (DESeq2 mean_control)",
    y = expression(bold("Per-motif enrichment score") ~
                   " (sum of -log"[10]*"(p) over bins) / motifs in database"),
    title = "RBPs enriched in drug-perturbed splicing are expressed in oocytes",
    subtitle = paste(
      "Y axis = per-motif average of -log10(p) over significant bins (corrects",
      "for the number of motifs the RBP has in the CISBP database);",
      "point area = breadth (significant bins);",
      "opacity = oocyte mRNA abundance (low-abundance RBPs fade out);",
      "dashed line = 1 motif reaching p = 0.05 in 1 bin."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    strip.background = element_rect(fill = "grey95", color = "black",
                                    linewidth = 0.4),
    strip.text = element_text(face = "bold", size = 10),
    axis.title = element_text(face = "bold", size = 11),
    axis.text  = element_text(size = 9, color = "black"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    legend.margin = margin(0, 0, 6, 0),
    plot.title = element_text(face = "bold", size = 13, hjust = 0,
                              margin = margin(b = 4)),
    plot.subtitle = element_text(size = 9, hjust = 0, color = "grey25",
                                 margin = margin(b = 10)),
    plot.margin = margin(12, 14, 12, 12)
  )

# -----------------------------------------------------------------------------
# 6. Save
# -----------------------------------------------------------------------------
dir.create("Figures/FigureS5_RBP_abundance", showWarnings = FALSE,
           recursive = TRUE)
# Use standard pdf() to avoid the cairo dependency (X11 libs unavailable)
pdf("Figures/FigureS5_RBP_abundance/FigureS5_RBP_abundance_vs_enrichment.pdf",
    width = 12, height = 6.5)
print(p_scatter)
dev.off()
ggsave("Figures/FigureS5_RBP_abundance/FigureS5_RBP_abundance_vs_enrichment.png",
       p_scatter, width = 12, height = 6.5, dpi = 300)
cat("Saved FigureS5_RBP_abundance_vs_enrichment.pdf/.png\n")
