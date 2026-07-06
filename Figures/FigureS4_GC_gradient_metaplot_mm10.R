# ============================================================================
# Figure S4 — full CONTINUOUS GC gradient metaplot (mm10)
#
# Reviewer ask (verbatim): "I want to see how the GC varies in the gradient
# between the target, intron/exon upstream/downstream, etc, etc."
#
# Each region of the cassette is broken into N equal proportional bins, GC
# is computed per bin from the actual DNA sequences Matt extracted from
# mm10 (stored in SEQ_UPEXON/SEQ_UPINTRON/SEQ_EXON/SEQ_DOINTRON/SEQ_DOEXON
# for exons and SEQ_LONGESTUPEXON/SEQ_INTRON/SEQ_LONGESTDOEXON for
# introns), and the binned GC is plotted as a continuous line.
#
# Two layouts are produced:
#   - "fine"   : 80 bins — used when we pool across Included+Skipped
#                 (or Retained+Excised) so we have ~1500 events per
#                 category; used for the Differential-vs-Unchanged view.
#   - "coarse" : 30 bins — used when we dissect by Retained/Excised/
#                 Unchanged, where each category has only ~80 events
#                 and 80 bins would be ~1 event / bin (pure noise).
#
# Variable-length introns are stretched proportionally along the
# synthesised 5' → 3' axis. The 5'/3' intron ends are binned at finer
# resolution than the body for both layouts (the splicing decision
# happens there).
#
# Output PNG/PDF files (all saved in this Figures/ directory):
#   - FigureS4_GC_gradient_exons.png/.pdf                    (per-category, fine bins on exons)
#   - FigureS4_GC_gradient_introns.png/.pdf                  (per-category, COARSE bins on introns)
#   - FigureS4_GC_gradient_exons_any_vs_unchanged.png/.pdf   (Differential vs Unchanged, fine bins on exons)
#   - FigureS4_GC_gradient_introns_any_vs_unchanged.png/.pdf (Differential vs Unchanged, fine bins on introns)
#   - FigureS4_GC_gradient_combined_any_vs_unchanged.png/.pdf (top-line reply figure: exons top, introns bottom)
#   - FigureS4_GC_gradient_summary_exons.csv                 (per-bin means + 95% CI)
#   - FigureS4_GC_gradient_summary_introns.csv               (fine 80-bin layout)
#   - FigureS4_GC_gradient_summary_introns_coarse.csv        (coarse 30-bin layout)
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

base_family <- "sans"

# ---- bin layouts ----
# Exon layout (5 regions, 68 bins total)
EXON_LAYOUT <- list(
  upExon    = list(end = "right", n_bins =  8, pad_nt = 200),
  upIntron  = list(mode = "full",  n_bins = 20),
  target    = list(mode = "full",  n_bins = 12),
  doIntron  = list(mode = "full",  n_bins = 20),
  downExon  = list(end = "left",  n_bins =  8, pad_nt = 200)
)

# Two intron layouts — picked by sample size (see file header).
make_intron_layout <- function(level = c("fine", "coarse")) {
  if (level == "fine") {
    list(
      upExon        = list(end = "right", n_bins = 16, pad_nt = 400),
      upIntron5SS   = list(end = "left",  n_bins =  8, pad_nt =  80),
      intronBody    = list(mode = "full",  n_bins = 40),
      doIntron3SS   = list(end = "right", n_bins =  8, pad_nt =  80),
      downExon      = list(end = "left",  n_bins = 16, pad_nt = 400)
    )
  } else {
    list(
      upExon        = list(end = "right", n_bins =  6, pad_nt = 150),
      upIntron5SS   = list(end = "left",  n_bins =  3, pad_nt =  80),
      intronBody    = list(mode = "full",  n_bins = 12),
      doIntron3SS   = list(end = "right", n_bins =  3, pad_nt =  80),
      downExon      = list(end = "left",  n_bins =  6, pad_nt = 150)
    )
  }
}

# Region order used by the plot facets
region_order_exon  <- c("upExon", "upIntron", "Exon", "downIntron", "downExon")
region_order_intron <- c("upExon", "Intron 5' flank", "Intron body",
                         "Intron 3' flank", "downExon")

is_empty_str <- function(s) {
  if (length(s) == 0) return(TRUE)
  if (is.na(s[1]))   return(TRUE)
  nchar(s[1]) == 0
}

# ---- per-bin GC ----
# Average GC of a character vector of nucleotides (skipping N).
gc_of_vec <- function(s) {
  if (length(s) == 0) return(NA_real_)
  s <- s[s %in% c("A","T","G","C")]
  if (length(s) == 0) return(NA_real_)
  100 * sum(s %in% c("G","C")) / length(s)
}

# Break a DNA string into exactly n_bins contiguous chunks (as equal in
# size as possible) and return the GC of each chunk. Returns NA for
# chunks with no sequence.
bin_to_gc <- function(seq, n_bins) {
  if (is_empty_str(seq)) return(rep(NA_real_, n_bins))
  chars <- strsplit(toupper(seq[1]), "")[[1]]
  chars <- chars[chars %in% c("A","T","G","C")]
  n <- length(chars)
  if (n == 0) return(rep(NA_real_, n_bins))
  if (n < n_bins) {
    out <- rep(NA_real_, n_bins)
    out[seq_len(n)] <- vapply(seq_len(n), function(i) gc_of_vec(chars[i]), numeric(1))
    return(out)
  }
  sizes <- rep(floor(n / n_bins), n_bins)
  sizes[seq_len(n %% n_bins)] <- sizes[seq_len(n %% n_bins)] + 1
  ends   <- cumsum(sizes)
  starts <- c(1, head(ends, -1) + 1)
  unname(mapply(function(s, e) gc_of_vec(chars[s:e]), starts, ends))
}

# Take an end-slice of a DNA string (the LAST n nt, or the FIRST n nt).
left_end  <- function(seq, n) {
  if (is_empty_str(seq)) return(NA_character_)
  L <- nchar(seq)
  if (L <= n) return(seq)
  substr(seq, L - n + 1, L)
}
right_end <- function(seq, n) {
  if (is_empty_str(seq)) return(NA_character_)
  L <- nchar(seq)
  if (L <= n) return(seq)
  substr(seq, 1, n)
}

# Peal off n_each_end nucleotides from each end of an intron (used to
# recover the body after the 5'/3' flanks have been binned separately).
strip_intron_ends <- function(seq, n_each_end) {
  if (is_empty_str(seq)) return(NA_character_)
  L <- nchar(seq)
  if (L <= 2 * n_each_end) return(NA_character_)
  substr(seq, n_each_end + 1, L - n_each_end)
}

# Per-event gradient builders: each returns the per-bin GC vector for
# that event (concatenated across the regions).
build_event_exon_gradient <- function(seq_upExon, seq_upIntron, seq_exon,
                                      seq_doIntron, seq_doExon) {
  c(bin_to_gc(left_end(seq_upExon,   EXON_LAYOUT$upExon$pad_nt),
              EXON_LAYOUT$upExon$n_bins),
    bin_to_gc(seq_upIntron,  EXON_LAYOUT$upIntron$n_bins),
    bin_to_gc(seq_exon,      EXON_LAYOUT$target$n_bins),
    bin_to_gc(seq_doIntron,  EXON_LAYOUT$doIntron$n_bins),
    bin_to_gc(right_end(seq_doExon, EXON_LAYOUT$downExon$pad_nt),
              EXON_LAYOUT$downExon$n_bins))
}

build_event_intron_gradient <- function(seq_upExon, seq_intron, seq_doExon,
                                        layout) {
  peel <- layout$upIntron5SS$pad_nt
  c(bin_to_gc(left_end(seq_upExon, layout$upExon$pad_nt),
              layout$upExon$n_bins),
    bin_to_gc(left_end(seq_intron, peel),
              layout$upIntron5SS$n_bins),
    bin_to_gc(strip_intron_ends(seq_intron, peel),
              layout$intronBody$n_bins),
    bin_to_gc(right_end(seq_intron, peel),
              layout$doIntron3SS$n_bins),
    bin_to_gc(right_end(seq_doExon, layout$downExon$pad_nt),
              layout$downExon$n_bins))
}

# Region labels for the per-bin index produced by the gradient builders.
x_anchors_exon <- c(rep("upExon",   EXON_LAYOUT$upExon$n_bins),
                    rep("upIntron", EXON_LAYOUT$upIntron$n_bins),
                    rep("Exon",     EXON_LAYOUT$target$n_bins),
                    rep("downIntron",EXON_LAYOUT$doIntron$n_bins),
                    rep("downExon", EXON_LAYOUT$downExon$n_bins))

x_anchors_intron_from_layout <- function(layout) {
  c(rep("upExon",          layout$upExon$n_bins),
    rep("Intron 5' flank", layout$upIntron5SS$n_bins),
    rep("Intron body",     layout$intronBody$n_bins),
    rep("Intron 3' flank", layout$doIntron3SS$n_bins),
    rep("downExon",        layout$downExon$n_bins))
}

# ---- read Matt tables ----
matt_root <- "../SupplementaryData1_matt_out"
read_matt <- function(drug, etype) {
  f <- if (etype == "exons")
    paste0(matt_root, "/", drug, "_ex/", drug, "_exons_with_efeatures.tab")
  else
    paste0(matt_root, "/", drug, "_int/", drug, "_introns_with_ifeatures.tab")
  read_tsv(f, show_col_types = FALSE)
}

exon_tables <- bind_rows(lapply(c("pladb","ssa","tub"), function(d)
  read_matt(d, "exons") %>% mutate(drug = d)))
intron_tables <- bind_rows(lapply(c("pladb","ssa","tub"), function(d)
  read_matt(d, "introns") %>% mutate(drug = d)))

cat("Loaded:", nrow(exon_tables), "exon rows and",
    nrow(intron_tables), "intron rows\n")

# ---- sub-sample unchanged pool to keep files fast and figures readable ----
set.seed(42)
max_unchanged <- 1500
sample_keep <- function(df) {
  out <- list(); k <- 1
  for (d in unique(df$drug)) {
    for (ds in unique(df$DATASET)) {
      sub <- df %>% filter(drug == d, DATASET == ds)
      n_in <- nrow(sub)
      rows_keep <- if (ds == "ndiff" && n_in > max_unchanged)
        sample.int(n_in, max_unchanged) else seq_len(n_in)
      out[[k]] <- tibble(drug = d, DATASET = ds, row_idx = rows_keep)
      k <- k + 1
    }
  }
  bind_rows(out)
}

km_exon  <- sample_keep(exon_tables)
km_intron <- sample_keep(intron_tables)

exon_tables <- exon_tables %>% group_by(drug, DATASET) %>%
  mutate(row_idx = row_number()) %>% ungroup() %>%
  inner_join(km_exon, by = c("drug","DATASET","row_idx")) %>%
  select(-row_idx)
intron_tables <- intron_tables %>% group_by(drug, DATASET) %>%
  mutate(row_idx = row_number()) %>% ungroup() %>%
  inner_join(km_intron, by = c("drug","DATASET","row_idx")) %>%
  select(-row_idx)

cat("After sampling — exons:\n");   print(exon_tables   %>% count(drug, DATASET))
cat("After sampling — introns:\n"); print(intron_tables %>% count(drug, DATASET))

# ---- build per-event GC gradients ----
build_grad <- function(df, build_fn, layout, x_anchors, n_bins) {
  df %>%
    rowwise() %>%
    mutate(gc_vec = list(build_fn(...))) %>%
    ungroup()
}

# Build helpers need explicit args; we wrap them so they fit the
# rowwise() use above.
build_exon_per_event <- function(df) {
  df %>%
    rowwise() %>%
    mutate(gc_vec = list(build_event_exon_gradient(
      SEQ_UPEXON, SEQ_UPINTRON, SEQ_EXON, SEQ_DOINTRON, SEQ_DOEXON
    ))) %>%
    ungroup() %>%
    mutate(bin = list(seq_len(length(x_anchors_exon)))) %>%
    tidyr::unnest(c(gc_vec, bin)) %>%
    rename(GCC = gc_vec) %>%
    mutate(region = x_anchors_exon[bin])
}

build_intron_per_event <- function(df, layout, x_anchors, label) {
  n_bins <- length(x_anchors)
  df %>%
    rowwise() %>%
    mutate(gc_vec = list(build_event_intron_gradient(
      SEQ_LONGESTUPEXON, SEQ_INTRON, SEQ_LONGESTDOEXON, layout = layout
    ))) %>%
    ungroup() %>%
    mutate(bin = list(seq_len(n_bins))) %>%
    tidyr::unnest(c(gc_vec, bin)) %>%
    rename(GCC = gc_vec) %>%
    mutate(region = x_anchors[bin])
}

exon_grad <- build_exon_per_event(exon_tables) %>% filter(!is.na(GCC))

INTRON_LAYOUT_FINE   <- make_intron_layout("fine")
INTRON_LAYOUT_COARSE <- make_intron_layout("coarse")
x_anchors_intron_fine   <- x_anchors_intron_from_layout(INTRON_LAYOUT_FINE)
x_anchors_intron_coarse <- x_anchors_intron_from_layout(INTRON_LAYOUT_COARSE)

intron_grad <- build_intron_per_event(intron_tables, INTRON_LAYOUT_FINE,
                                     x_anchors_intron_fine, "fine") %>%
  filter(!is.na(GCC))

intron_grad_coarse <- build_intron_per_event(intron_tables, INTRON_LAYOUT_COARSE,
                                           x_anchors_intron_coarse, "coarse") %>%
  filter(!is.na(GCC))

cat("Exon gradient points :", nrow(exon_grad), "\n")
cat("Intron gradient (fine):",   nrow(intron_grad),       "\n")
cat("Intron gradient (coarse):", nrow(intron_grad_coarse),"\n")

# ---- aggregate per (drug, bin, region, category) ----
aggregate_grad <- function(df, x_anchors) {
  df %>%
    group_by(drug, bin, region, category) %>%
    summarise(mean_gcc = mean(GCC),
              se_gcc   = sd(GCC)   / sqrt(sum(!is.na(GCC))),
              n_events = sum(!is.na(GCC)),
              .groups  = "drop") %>%
    group_by(drug, region) %>%
    mutate(bin_in_region = row_number()) %>%
    ungroup() %>%
    mutate(region = as.character(region))
}

# Attach friendly category label
add_exon_category <- function(df) {
  df %>% mutate(category = case_when(
    DATASET == "up"   ~ "Included",
    DATASET == "down" ~ "Skipped",
    DATASET == "ndiff"~ "Unchanged",
    TRUE ~ NA_character_
  ))
}
add_intron_category <- function(df) {
  df %>% mutate(category = case_when(
    DATASET == "up"   ~ "Retained",
    DATASET == "down" ~ "Excised",
    DATASET == "ndiff"~ "Unchanged",
    TRUE ~ NA_character_
  ))
}

exon_grad       <- add_exon_category(exon_grad)
intron_grad     <- add_intron_category(intron_grad)
intron_grad_coarse <- add_intron_category(intron_grad_coarse)

exon_agg          <- aggregate_grad(exon_grad,         x_anchors_exon)
intron_agg        <- aggregate_grad(intron_grad,       x_anchors_intron_fine)
intron_agg_coarse <- aggregate_grad(intron_grad_coarse, x_anchors_intron_coarse)

write_csv(exon_agg,          "FigureS4_GC_gradient_summary_exons.csv")
write_csv(intron_agg,        "FigureS4_GC_gradient_summary_introns.csv")
write_csv(intron_agg_coarse, "FigureS4_GC_gradient_summary_introns_coarse.csv")

# ---- theme ----
theme_grad <- function(base_size = 11, base_family_in = base_family) {
  theme_classic(base_size = base_size, base_family = base_family_in) %+replace%
    theme(
      axis.line     = element_line(linewidth = 0.6, colour = "black"),
      axis.ticks    = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      axis.title    = element_text(face = "bold", size = rel(1.0)),
      axis.text.x   = element_text(size = rel(0.7),  colour = "black", lineheight = 0.95),
      axis.text.y   = element_text(size = rel(0.85), colour = "black"),
      legend.position = "bottom",
      legend.background = element_blank(),
      legend.title = element_text(face = "bold", size = rel(0.85)),
      legend.text  = element_text(size = rel(0.8)),
      legend.key.size = unit(0.4, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5),
      strip.text   = element_text(face = "bold", size = rel(0.9)),
      plot.title   = element_text(face = "bold", size = rel(1.0), hjust = 0.5),
      plot.subtitle = element_text(size = rel(0.75), hjust = 0.5, colour = "grey40",
                                  lineheight = 1.05),
      plot.caption = element_text(size = rel(0.75), colour = "#666666"),
      plot.margin  = margin(12, 8, 12, 8)
    )
}

wrap_subtitle <- function(s, width = 100) {
  paste(strwrap(s, width = width), collapse = "\n")
}

save_plot <- function(plot, base, w, h) {
  ggsave(paste0(base, ".png"), plot = plot, width = w, height = h, dpi = 600, bg = "white")
  ggsave(paste0(base, ".pdf"), plot = plot, width = w, height = h, bg = "white")
  message("Saved: ", base, ".png & .pdf")
}

# ---- gradients ----
splicing_exon_colors  <- c("Included"  = "#D4A574", "Skipped"   = "#B85450",
                            "Unchanged" = "#808080")
splicing_intron_colors <- c("Retained"  = "#D4A574", "Excised"   = "#B85450",
                            "Unchanged" = "#808080")
palette_any <- c("Differential" = "#B85450", "Unchanged" = "#808080")

make_gradient_plot <- function(agg, palette, title, subtitle, region_levels,
                                x_breaks, x_labels, label_size = 8,
                                wrap_text = TRUE) {
  agg <- agg %>% mutate(region = factor(region, levels = region_levels))
  n_max <- max(agg$bin)
  if (wrap_text) subtitle <- wrap_subtitle(subtitle)

  ggplot(agg, aes(x = (bin - 0.5) / n_max, y = mean_gcc,
                  color = category, group = category)) +
    geom_line(linewidth = 1.4, alpha = 0.95) +
    geom_point(size = 2.0, alpha = 0.95) +
    geom_errorbar(aes(ymin = mean_gcc - 1.96 * se_gcc,
                      ymax = mean_gcc + 1.96 * se_gcc),
                  width = 0.0035, linewidth = 0.3, alpha = 0.7) +
    scale_color_manual(values = palette) +
    scale_y_continuous(limits = c(35, 65), breaks = seq(35, 65, 5),
                       expand = expansion(mult = c(0.02, 0.02))) +
    facet_wrap(~ drug, nrow = 1) +
    labs(x = NULL, y = "Mean GC (%)", title = title, subtitle = subtitle) +
    theme_grad() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(size = label_size)) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, limits = c(0, 1),
                        expand = expansion(mult = c(0.01, 0.01)))
}

make_any_plot <- function(agg, palette, title, subtitle, region_levels,
                          x_breaks, x_labels, label_size = 8) {
  df_any <- agg %>%
    mutate(category_any = if_else(category %in% c("Included","Skipped",
                                                   "Retained","Excised"),
                                  "Differential", as.character(category))) %>%
    filter(category_any %in% c("Differential","Unchanged")) %>%
    mutate(region = factor(region, levels = unique(region)))
  n_max <- max(df_any$bin)

  ggplot(df_any, aes(x = (bin - 0.5) / n_max, y = mean_gcc,
                     color = category_any, group = category_any)) +
    geom_line(linewidth = 1.4, alpha = 0.95) +
    geom_point(size = 2.0, alpha = 0.95) +
    geom_errorbar(aes(ymin = mean_gcc - 1.96 * se_gcc,
                      ymax = mean_gcc + 1.96 * se_gcc),
                  width = 0.0035, linewidth = 0.3, alpha = 0.7) +
    scale_color_manual(values = palette) +
    scale_y_continuous(limits = c(35, 65), breaks = seq(35, 65, 5),
                       expand = expansion(mult = c(0.02, 0.02))) +
    facet_wrap(~ drug, nrow = 1) +
    labs(x = NULL, y = "Mean GC (%)", title = title,
         subtitle = wrap_subtitle(subtitle)) +
    theme_grad() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(size = label_size)) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, limits = c(0, 1),
                        expand = expansion(mult = c(0.01, 0.01)))
}

# Nice drug labels for faceting
relabel_drugs <- function(df) {
  df %>% mutate(drug = recode(drug,
                               pladb = "Pladienolide B",
                               ssa   = "Spliceostatin A",
                               tub   = "Tubercidin"),
                drug = factor(drug, levels = c("Tubercidin",
                                                 "Pladienolide B",
                                                 "Spliceostatin A")))
}

# ---- make the four panels ----
exon_agg_l          <- relabel_drugs(exon_agg)
intron_agg_coarse_l <- relabel_drugs(intron_agg_coarse)
intron_agg_l        <- relabel_drugs(intron_agg)

p_exons_v2 <- make_gradient_plot(
  exon_agg_l,
  splicing_exon_colors,
  "Cassette exons — continuous GC gradient",
  "Each region is binned into N equal segments; variable-length introns are stretched proportionally. Error bars = 95% CI of the bin mean.",
  region_order_exon,
  x_breaks = c(0.0875, 0.275, 0.5, 0.725, 0.9125),
  x_labels = c("upExon (200 nt, 8 bins)",
               "upIntron (full, 20 bins)",
               "Target exon (full, 12 bins)",
               "downIntron (full, 20 bins)",
               "downExon (200 nt, 8 bins)"))

p_introns_v2 <- make_gradient_plot(
  intron_agg_coarse_l,
  splicing_intron_colors,
  "Retained introns — continuous GC gradient (per category)",
  "Coarser bins (30 total) because Retained/Excised/Unchanged each have only ~80 events. Error bars = 95% CI of the bin mean.",
  region_order_intron,
  x_breaks = c(0.085, 0.225, 0.55, 0.84, 0.975),
  x_labels = c("upExon", "5' flank", "Intron body", "3' flank", "downExon"),
  label_size = 9)

p_exons_any  <- make_any_plot(
  exon_agg_l,
  palette_any,
  "Cassette exons — gradient GC (any change vs unchanged)",
  "Differential = pooled included + skipped events; Unchanged = ndiff background. Error bars = 95% CI.",
  region_order_exon,
  x_breaks = c(0.0875, 0.275, 0.5, 0.725, 0.9125),
  x_labels = c("upExon (200 nt, 8 bins)",
               "upIntron (full, 20 bins)",
               "Target exon (full, 12 bins)",
               "downIntron (full, 20 bins)",
               "downExon (200 nt, 8 bins)"))

p_introns_any <- make_any_plot(
  intron_agg_l,
  palette_any,
  "Retained introns — gradient GC (any change vs unchanged)",
  "Fine bins at the 5'/3' intron ends where the splicing decision happens (5'SS, branch point, PPT); coarser bins through the body. Error bars = 95% CI.",
  region_order_intron,
  x_breaks = c(0.075, 0.225, 0.55, 0.875, 0.97),
  x_labels = c("upExon (400 nt, 16 bins)",
               "5' flank (80 nt, 8 bins)",
               "Intron body (rest, 40 bins)",
               "3' flank (80 nt, 8 bins)",
               "downExon (400 nt, 16 bins)"))

save_plot(p_exons_v2,  "FigureS4_GC_gradient_exons",  18, 5.0)
save_plot(p_introns_v2,"FigureS4_GC_gradient_introns",18, 5.0)
save_plot(p_exons_any,  "FigureS4_GC_gradient_exons_any_vs_unchanged",  18, 5.0)
save_plot(p_introns_any,"FigureS4_GC_gradient_introns_any_vs_unchanged",18, 5.0)

p_combined_any <- (p_exons_any / p_introns_any) +
  plot_annotation(
    title    = "Continuous GC gradient across cassette exons and retained introns",
    subtitle = paste("Each region binned into N equal segments; the 5'/3'",
                     "intron ends are binned at finer resolution where the",
                     "splicing decision happens. Error bars are 95% CI."),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14, family = base_family),
      plot.subtitle = element_text(size = 11, colour = "grey30",
                                   lineheight = 1.05, family = base_family)
    )
  )
save_plot(p_combined_any, "FigureS4_GC_gradient_combined_any_vs_unchanged",
          18, 10)

message("Done.")
