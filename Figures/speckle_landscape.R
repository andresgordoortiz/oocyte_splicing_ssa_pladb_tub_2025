###############################################################################
# "The Oocyte Speckle-Dependency Landscape"
#
# Three-panel figure testing speckle-dependency predictions for oocyte splicing:
#   Panel A (Małszycki Test): Intron GC vs Exon GC scatter with quadrant stats
#   Panel B (Architecture):   ECDF of intron/exon unit size (3 groups)
#   Panel C (GC Profile):     Position-resolved GC along intron–exon–intron
#                             (Małszycki et al. Fig 6F style, line + heatmap)
#
# Data inputs:
#   - EVENT_INFO-mm10.tab     (VAST-TOOLS event annotation with sequences)
#   - ssa_fdr.csv, tub_fdr.csv, pladb_fdr.csv  (betAS FDR results)
#   - BSgenome.Mmusculus.UCSC.mm10              (intronic sequence extraction)
#
# Install BSgenome if needed:
#   BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
###############################################################################

# ── 0. Libraries ─────────────────────────────────────────────────────────────
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(showtext)
library(patchwork)
library(purrr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)
library(GenomicRanges)

set.seed(42)

# ── 1. Font & Theme ──────────────────────────────────────────────────────────
preferred_font <- "Roboto"
font_add_google(preferred_font)
# showtext is OFF by default for interactive display (avoids viewport errors).
# It is enabled only inside save_plot() for high-quality file output.
showtext_auto(enable = FALSE)
base_family <- "sans"   # fallback for interactive preview

# Helper: save with showtext + high DPI, then disable again
save_plot <- function(filename, plot, width, height, dpi = 600, ...) {
  showtext_auto(enable = TRUE)
  showtext::showtext_opts(dpi = dpi)
  on.exit({
    showtext_auto(enable = FALSE)
    showtext::showtext_opts(dpi = 96)
  }, add = TRUE)
  ggsave(filename, plot = plot, width = width, height = height, dpi = dpi, ...)
}

drug_colors <- c(
  "SSA"   = "#706993",   # Spliceostatin A – purple
  "TUB"   = "#A0C1B9",   # Tubercidin – sage green
  "PLADB" = "#70A0AF"    # Pladienolide B – teal
)

drug_labels <- c(
  "SSA"   = "Spliceostatin A",
  "TUB"   = "Tubercidin",
  "PLADB" = "Pladienolide B"
)

direction_shapes <- c("Included" = 16, "Skipped" = 17)

theme_cellpub <- function(base_size = 14) {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line           = element_line(linewidth = 0.6, colour = "black"),
      axis.ticks          = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks.length   = unit(0.15, "cm"),
      axis.title          = element_text(face = "bold", size = rel(1.0)),
      axis.text           = element_text(size = rel(0.9), colour = "black"),
      legend.background   = element_blank(),
      legend.title        = element_text(face = "bold", size = rel(0.9)),
      legend.text         = element_text(size = rel(0.85)),
      legend.key.size     = unit(0.4, "cm"),
      panel.grid.major    = element_blank(),
      panel.grid.minor    = element_blank(),
      strip.background    = element_rect(fill = "grey95", color = "black",
                                         linewidth = 0.5),
      strip.text          = element_text(face = "bold", size = rel(1.0),
                                         margin = margin(4, 4, 4, 4)),
      plot.title          = element_text(face = "bold", size = rel(1.1),
                                         hjust = 0.5),
      plot.subtitle       = element_text(size = rel(0.85), hjust = 0.5,
                                         color = "grey30"),
      plot.margin         = margin(8, 8, 8, 8)
    )
}
my_theme <- theme_cellpub()

# ── 2. Load Data ─────────────────────────────────────────────────────────────
message("Loading data...")
event_info   <- read.delim("EVENT_INFO-mm10.tab", stringsAsFactors = FALSE)
tub_fdr_df   <- read_csv("tub_fdr.csv",   show_col_types = FALSE)[, -1]
pladb_fdr_df <- read_csv("pladb_fdr.csv", show_col_types = FALSE)[, -1]
ssa_fdr_df   <- read_csv("ssa_fdr.csv",   show_col_types = FALSE)[, -1]

genome <- BSgenome.Mmusculus.UCSC.mm10

# ── 3. Significant Exon Events ──────────────────────────────────────────────
filter_sig <- function(df) {
  na.omit(df[df$FDR <= 0.05 & abs(df$deltapsi) >= 0.1, ])
}

sig_exons <- bind_rows(
  filter_sig(ssa_fdr_df)   %>% filter(grepl("EX", EVENT)) %>% mutate(drug = "SSA"),
  filter_sig(tub_fdr_df)   %>% filter(grepl("EX", EVENT)) %>% mutate(drug = "TUB"),
  filter_sig(pladb_fdr_df) %>% filter(grepl("EX", EVENT)) %>% mutate(drug = "PLADB")
) %>%
  mutate(direction = if_else(deltapsi < 0, "Skipped", "Included"))

message(nrow(sig_exons), " significant exon–drug pairs")



# ── 4. Coordinate Parsing (vectorised) ──────────────────────────────────────
# Handles VAST-TOOLS formats: "chr1:100-200", "chrX:100+150-200", etc.
parse_all_coords <- function(coord_strs) {
  n      <- length(coord_strs)
  chr    <- rep(NA_character_, n)
  starts <- rep(NA_real_, n)
  ends   <- rep(NA_real_, n)

  ok <- !is.na(coord_strs) & coord_strs != ""
  if (any(ok)) {
    chr[ok]   <- sub(":.*", "", coord_strs[ok])
    pos_strs  <- sub(".*:", "", coord_strs[ok])
    nums_list <- regmatches(pos_strs, gregexpr("[0-9]+", pos_strs))
    starts[ok] <- vapply(nums_list, function(x) {
      v <- as.numeric(x); if (length(v) == 0) NA_real_ else min(v)
    }, numeric(1))
    ends[ok] <- vapply(nums_list, function(x) {
      v <- as.numeric(x); if (length(v) == 0) NA_real_ else max(v)
    }, numeric(1))
  }
  tibble(chr = chr, start = starts, end = ends)
}

# ── 5. Batch Intron Feature Extraction ───────────────────────────────────────
# For each exon event, extracts up to `flank_size` nt of intronic sequence on
# each side of the alternative exon, computes GC, and returns intron lengths.
compute_intron_features <- function(df, genome, flank_size = 100) {

  message("  Parsing coordinates...")
  c1 <- parse_all_coords(df$CO_C1)
  a  <- parse_all_coords(df$CO_A)
  c2 <- parse_all_coords(df$CO_C2)
  n  <- nrow(df)

  # Which constitutive exon sits left (lower genomic coord) of A?
  a_mid      <- (a$start + a$end) / 2
  c1_mid     <- (c1$start + c1$end) / 2
  c1_is_left <- c1_mid < a_mid

  left_c_end    <- if_else(c1_is_left, c1$end, c2$end)
  right_c_start <- if_else(c1_is_left, c2$start, c1$start)

  # Full intron boundaries
  li_start <- left_c_end + 1         # left intron start
  li_end   <- a$start - 1            # left intron end
  ri_start <- a$end + 1              # right intron start
  ri_end   <- right_c_start - 1      # right intron end

  li_len <- li_end - li_start + 1
  ri_len <- ri_end - ri_start + 1

  # Flanking intronic region (up to flank_size nt nearest the exon)
  lf_start <- pmax(li_start, a$start - flank_size)
  lf_end   <- a$start - 1
  rf_start <- a$end + 1
  rf_end   <- pmin(ri_end, a$end + flank_size)

  # Validity mask
  valid <- !is.na(a$chr) & !is.na(c1$chr) & !is.na(c2$chr) &
           a$chr == c1$chr & a$chr == c2$chr &
           a$chr %in% seqnames(genome) &
           !is.na(li_len) & !is.na(ri_len) &
           li_len >= 1 & ri_len >= 1 &
           lf_start <= lf_end & rf_start <= rf_end

  # Clamp to chromosome boundaries
  chr_sl    <- seqlengths(genome)
  valid_idx <- which(valid)
  lf_start[valid_idx] <- pmax(1L, lf_start[valid_idx])
  rf_end[valid_idx]   <- pmin(chr_sl[a$chr[valid_idx]], rf_end[valid_idx])
  # Re-check widths after clamping
  valid[valid_idx] <- lf_start[valid_idx] <= lf_end[valid_idx] &
                      rf_start[valid_idx] <= rf_end[valid_idx]
  valid_idx <- which(valid)

  message("  Extracting ", length(valid_idx),
          " flanking sequences from BSgenome.Mmusculus.UCSC.mm10...")

  # Build GRanges and batch-extract sequences
  left_gr  <- GRanges(a$chr[valid_idx],
                       IRanges(lf_start[valid_idx], lf_end[valid_idx]))
  right_gr <- GRanges(a$chr[valid_idx],
                       IRanges(rf_start[valid_idx], rf_end[valid_idx]))

  left_dna  <- getSeq(genome, left_gr)
  right_dna <- getSeq(genome, right_gr)

  # GC content via Biostrings (much faster than character loops)
  left_gc_v  <- as.numeric(letterFrequency(left_dna,  "GC", as.prob = TRUE)) * 100
  right_gc_v <- as.numeric(letterFrequency(right_dna, "GC", as.prob = TRUE)) * 100

  # Map back to full-length vectors
  left_intron_len  <- rep(NA_real_, n)
  right_intron_len <- rep(NA_real_, n)
  left_flank_gc    <- rep(NA_real_, n)
  right_flank_gc   <- rep(NA_real_, n)
  left_flank_seq   <- rep(NA_character_, n)
  right_flank_seq  <- rep(NA_character_, n)

  left_intron_len[valid_idx]  <- li_len[valid_idx]
  right_intron_len[valid_idx] <- ri_len[valid_idx]
  left_flank_gc[valid_idx]    <- left_gc_v
  right_flank_gc[valid_idx]   <- right_gc_v
  left_flank_seq[valid_idx]   <- as.character(left_dna)
  right_flank_seq[valid_idx]  <- as.character(right_dna)

  tibble(
    left_intron_len  = left_intron_len,
    right_intron_len = right_intron_len,
    left_flank_gc    = left_flank_gc,
    right_flank_gc   = right_flank_gc,
    left_flank_seq   = left_flank_seq,
    right_flank_seq  = right_flank_seq,
    avg_intron_gc    = (left_flank_gc + right_flank_gc) / 2,
    avg_intron_len   = (left_intron_len + right_intron_len) / 2
  )
}

# ── GC helper for character sequences (Seq_A column) ────────────────────────
gc_content_vec <- function(seqs) {
  vapply(seqs, function(s) {
    if (is.na(s) || nchar(s) == 0) return(NA_real_)
    ch <- strsplit(toupper(s), "")[[1]]
    100 * sum(ch %in% c("G", "C")) / length(ch)
  }, numeric(1), USE.NAMES = FALSE)
}

# ── GC Profile: per-position GC along intron–exon–intron architecture ────────
# Bins flanking intron sequences (flank_bins windows) and the exon (exon_bins
# windows, length-normalised).  Returns n × total_bins matrix of %GC.
compute_gc_profile <- function(left_seq, exon_seq, right_seq,
                                flank_bins = 20L, exon_bins = 10L) {
  n <- length(left_seq)
  total_bins <- flank_bins + exon_bins + flank_bins
  mat <- matrix(NA_real_, nrow = n, ncol = total_bins)

  bin_gc <- function(chars, n_bins) {
    edges <- round(seq(0, length(chars), length.out = n_bins + 1))
    vapply(seq_len(n_bins), function(b) {
      idx <- seq.int(edges[b] + 1L, edges[b + 1L])
      if (length(idx) == 0) return(NA_real_)
      100 * sum(chars[idx] %in% c("G", "C")) / length(idx)
    }, numeric(1))
  }

  for (i in seq_len(n)) {
    ls <- left_seq[i]
    if (!is.na(ls) && nchar(ls) >= 5)
      mat[i, seq_len(flank_bins)] <-
        bin_gc(strsplit(toupper(ls), "")[[1]], flank_bins)

    es <- exon_seq[i]
    if (!is.na(es) && nchar(es) >= 5)
      mat[i, flank_bins + seq_len(exon_bins)] <-
        bin_gc(strsplit(toupper(es), "")[[1]], exon_bins)

    rs <- right_seq[i]
    if (!is.na(rs) && nchar(rs) >= 5)
      mat[i, flank_bins + exon_bins + seq_len(flank_bins)] <-
        bin_gc(strsplit(toupper(rs), "")[[1]], flank_bins)
  }

  mat
}

# ── 6. Process Significant Events ───────────────────────────────────────────
message("\nProcessing significant events...")
sig_merged <- sig_exons %>%
  left_join(
    event_info %>% select(EVENT, CO_C1, CO_A, CO_C2, Seq_A, LE_o),
    by = "EVENT"
  ) %>%
  filter(!is.na(CO_C1) & !is.na(CO_A) & !is.na(CO_C2))

sig_intron <- compute_intron_features(sig_merged, genome)
sig_full   <- bind_cols(sig_merged, sig_intron) %>%
  mutate(exon_gc = gc_content_vec(Seq_A))

# ── 7. Process Background (non-significant) Events ──────────────────────────
# Background = all tested exon events that are NOT significant in any drug
message("\nProcessing non-significant background events...")
all_sig_ids <- unique(sig_exons$EVENT)

bg_pool <- bind_rows(
  ssa_fdr_df   %>% filter(grepl("EX", EVENT)) %>% select(EVENT, GENE),
  tub_fdr_df   %>% filter(grepl("EX", EVENT)) %>% select(EVENT, GENE),
  pladb_fdr_df %>% filter(grepl("EX", EVENT)) %>% select(EVENT, GENE)
) %>%
  distinct(EVENT, .keep_all = TRUE) %>%
  filter(!EVENT %in% all_sig_ids)
message(nrow(bg_pool), " non-significant exon events in background pool")

bg_n      <- min(3000, nrow(bg_pool))
bg_sample <- bg_pool %>%
  sample_n(bg_n) %>%
  left_join(
    event_info %>% select(EVENT, CO_C1, CO_A, CO_C2, Seq_A, LE_o),
    by = "EVENT"
  ) %>%
  filter(!is.na(CO_C1) & !is.na(CO_A) & !is.na(CO_C2))

bg_intron <- compute_intron_features(bg_sample, genome)
bg_full   <- bind_cols(bg_sample, bg_intron) %>%
  mutate(
    exon_gc   = gc_content_vec(Seq_A),
    drug      = "Background",
    direction = "Unchanged"
  )

# ── Filtered plot-ready data ────────────────────────────────────────────────
sig_plot <- sig_full %>% filter(!is.na(avg_intron_gc) & !is.na(exon_gc))
bg_plot  <- bg_full  %>% filter(!is.na(avg_intron_gc) & !is.na(exon_gc))
message("\n", nrow(sig_plot), " drug-affected events with complete data")
message(nrow(bg_plot),  " non-significant background events with complete data")

###############################################################################
# HELPERS
###############################################################################

# ── GC profile: line + heatmap helper ────────────────────────────────────────
profile_for_group <- function(df, group_name, max_heat = 500,
                               flank_bins = 20L, exon_bins = 10L) {
  total_bins <- flank_bins + exon_bins + flank_bins
  valid <- df %>%
    filter(!is.na(left_flank_seq) & !is.na(Seq_A) & !is.na(right_flank_seq) &
           nchar(left_flank_seq) >= 5 & nchar(Seq_A) >= 5 &
           nchar(right_flank_seq) >= 5)
  if (nrow(valid) == 0)
    return(list(mean = tibble(), heat = tibble(), n = 0L))

  mat <- compute_gc_profile(valid$left_flank_seq, valid$Seq_A,
                             valid$right_flank_seq, flank_bins, exon_bins)

  mean_df <- tibble(bin = seq_len(total_bins),
                     gc  = colMeans(mat, na.rm = TRUE),
                     group = group_name)

  if (nrow(mat) > max_heat) {
    idx <- sort(sample.int(nrow(mat), max_heat))
    mat <- mat[idx, , drop = FALSE]
  }
  mat <- mat[order(rowMeans(mat, na.rm = TRUE)), , drop = FALSE]

  heat_df <- expand.grid(rank = seq_len(nrow(mat)),
                          bin  = seq_len(total_bins)) %>%
    as_tibble() %>%
    mutate(gc = as.vector(mat), group = group_name)

  list(mean = mean_df, heat = heat_df, n = nrow(valid))
}

# Assign quadrants: Q1=leveled (high/high), Q2=low-exon/high-intron,
# Q3=low/low, Q4=high-exon/low-intron
assign_quadrant <- function(exon_gc, intron_gc, med_e, med_i) {
  case_when(
    exon_gc >= med_e & intron_gc >= med_i ~ "Q1: Leveled",
    exon_gc <  med_e & intron_gc >= med_i ~ "Q2: High intron",
    exon_gc <  med_e & intron_gc <  med_i ~ "Q3: Low both",
    exon_gc >= med_e & intron_gc <  med_i ~ "Q4: High exon",
    TRUE ~ NA_character_
  )
}

###############################################################################
# MAIN ANALYSIS — All significant exons combined
###############################################################################

out_dir <- "Figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Export full data table
sig_export <- sig_full %>%
  select(EVENT, GENE, drug, direction, deltapsi, FDR,
         exon_gc, avg_intron_gc,
         left_intron_len, right_intron_len, avg_intron_len,
         left_flank_gc, right_flank_gc)
write_csv(sig_export, file.path(out_dir, "speckle_landscape_data.csv"))

message("\n", strrep("=", 72))
message("  All significant exons (n=", nrow(sig_plot), ")")
message(strrep("=", 72))

# Medians from background
med_exon_gc   <- median(bg_plot$exon_gc,       na.rm = TRUE)
med_intron_gc <- median(bg_plot$avg_intron_gc,  na.rm = TRUE)

sig_plot <- sig_plot %>%
  mutate(quadrant = assign_quadrant(exon_gc, avg_intron_gc,
                                     med_exon_gc, med_intron_gc))
bg_q <- bg_plot %>%
  mutate(quadrant = assign_quadrant(exon_gc, avg_intron_gc,
                                     med_exon_gc, med_intron_gc))

# ── PANEL A — Intron GC vs Exon GC ──────────────────────────────────────
message("\n── Panel A: Intron GC vs Exon GC ──")

quad_counts_sig <- sig_plot %>%
  count(quadrant, name = "n_sig") %>%
  mutate(pct_sig = round(100 * n_sig / sum(n_sig), 1))
quad_counts_bg <- bg_q %>%
  count(quadrant, name = "n_bg") %>%
  mutate(pct_bg = round(100 * n_bg / sum(n_bg), 1))
quad_stats <- left_join(quad_counts_sig, quad_counts_bg, by = "quadrant")
message("  Quadrant distribution:")
for (i in seq_len(nrow(quad_stats)))
  message("    ", quad_stats$quadrant[i],
          ": Affected=", quad_stats$pct_sig[i], "%",
          ", Background=", quad_stats$pct_bg[i], "%")

xr <- range(c(sig_plot$exon_gc, bg_q$exon_gc), na.rm = TRUE)
yr <- range(c(sig_plot$avg_intron_gc, bg_q$avg_intron_gc), na.rm = TRUE)

quad_label_df <- tibble(
  quadrant = c("Q1: Leveled", "Q2: High intron", "Q3: Low both", "Q4: High exon"),
  x = c(med_exon_gc + (xr[2] - med_exon_gc) * 0.55,
        med_exon_gc - (med_exon_gc - xr[1]) * 0.55,
        med_exon_gc - (med_exon_gc - xr[1]) * 0.55,
        med_exon_gc + (xr[2] - med_exon_gc) * 0.55),
  y = c(med_intron_gc + (yr[2] - med_intron_gc) * 0.85,
        med_intron_gc + (yr[2] - med_intron_gc) * 0.85,
        med_intron_gc - (med_intron_gc - yr[1]) * 0.85,
        med_intron_gc - (med_intron_gc - yr[1]) * 0.85)
) %>%
  left_join(quad_stats, by = "quadrant") %>%
  mutate(label = paste0(sub(":.*", "", quadrant),
                        "\nAffected: ", pct_sig, "%",
                        "\nBackground: ", pct_bg, "%"))

p_a <- ggplot() +
  geom_point(data = bg_q,
             aes(x = exon_gc, y = avg_intron_gc),
             color = "grey80", size = 0.6, alpha = 0.3) +
  geom_point(data = sig_plot,
             aes(x = exon_gc, y = avg_intron_gc, color = drug),
             size = 2, alpha = 0.7) +
  geom_vline(xintercept = med_exon_gc,   linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = med_intron_gc, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_text(data = quad_label_df,
            aes(x = x, y = y, label = label),
            size = 2.8, color = "grey35", lineheight = 0.9) +
  scale_color_manual(values = drug_colors, labels = drug_labels) +
  labs(x = "Exon GC Content (%)",
       y = "Flanking Intron GC Content (%)\n(\u00b1100 nt)",
       title = "Intron\u2013Exon GC Architecture",
       color = "Drug") +
  my_theme +
  theme(legend.position = "right")

# ── PANEL B — ECDF of unit size per drug ────────────────────────────────
message("\n── Panel B: ECDF of unit size (per drug) ──")

ecdf_pal <- c(
  "Spliceostatin A" = "#706993",
  "Pladienolide B"  = "#70A0AF",
  "Tubercidin"      = "#A0C1B9",
  "Non-significant" = "grey60"
)

make_arch <- function(df, grp) {
  df %>%
    filter(!is.na(left_intron_len) & !is.na(right_intron_len) & !is.na(LE_o)) %>%
    mutate(unit_size = left_intron_len + LE_o + right_intron_len, group = grp)
}

arch_ssa   <- make_arch(sig_plot %>% filter(drug == "SSA"),   "Spliceostatin A")
arch_pladb <- make_arch(sig_plot %>% filter(drug == "PLADB"), "Pladienolide B")
arch_tub   <- make_arch(sig_plot %>% filter(drug == "TUB"),   "Tubercidin")
arch_bg    <- make_arch(bg_plot, "Non-significant")

arch_df <- bind_rows(arch_ssa, arch_pladb, arch_tub, arch_bg) %>%
  mutate(group = factor(group,
    levels = c("Non-significant", "Spliceostatin A", "Pladienolide B", "Tubercidin")))

for (g in levels(arch_df$group)) {
  n_g <- sum(arch_df$group == g)
  med <- median(arch_df$unit_size[arch_df$group == g], na.rm = TRUE)
  message("  ", g, ": n=", n_g,
          "  median unit = ", format(round(med), big.mark = ","), " nt")
}

# KS tests vs background
ks_results <- list()
ks_lines <- c()
for (dname in c("Spliceostatin A", "Pladienolide B", "Tubercidin")) {
  sub_v <- arch_df$unit_size[arch_df$group == dname]
  bg_v  <- arch_df$unit_size[arch_df$group == "Non-significant"]
  ks <- ks.test(sub_v, bg_v)
  ks_results[[dname]] <- ks
  ks_lines <- c(ks_lines,
    paste0(dname, " vs bg: D=", round(ks$statistic, 3),
           " p=", format.pval(ks$p.value, digits = 2)))
  message("  KS ", dname, " vs bg: D = ", round(ks$statistic, 3),
          ", p = ", format.pval(ks$p.value, digits = 3))
}

p_b <- ggplot(arch_df %>% filter(unit_size > 0),
              aes(x = unit_size, color = group)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  scale_x_log10(labels = scales::label_log(base = 10),
                breaks = 10^seq(2, 6)) +
  scale_color_manual(
    values = ecdf_pal,
    labels = function(x) {
      sapply(x, function(lab) {
        n_g <- sum(arch_df$group == lab)
        paste0(lab, " (n=", format(n_g, big.mark = ","), ")")
      })
    }
  ) +
  annotate("text", x = 200, y = 0.82,
           label = paste(ks_lines, collapse = "\n"),
           hjust = 0, size = 2.6, fontface = "bold", color = "grey30") +
  labs(x = "Intron / Exon Unit Size (nt)", y = "Proportion", color = NULL,
       title = "Intron/Exon Unit Size",
       subtitle = "Per-drug vs non-significant background") +
  my_theme +
  theme(legend.position = c(0.72, 0.22),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA))

# ── PANEL C — GC Content Profile per drug ───────────────────────────────
message("\n── Panel C: GC content profile (per drug) ──")

flank_bins <- 20L; exon_bins <- 10L
total_bins <- flank_bins + exon_bins + flank_bins
exon_left  <- flank_bins + 0.5
exon_right <- flank_bins + exon_bins + 0.5

prof_ssa   <- profile_for_group(sig_plot %>% filter(drug == "SSA"),   "Spliceostatin A")
prof_pladb <- profile_for_group(sig_plot %>% filter(drug == "PLADB"), "Pladienolide B")
prof_tub   <- profile_for_group(sig_plot %>% filter(drug == "TUB"),   "Tubercidin")
prof_bg    <- profile_for_group(bg_plot, "Non-significant")

message("  Events with seqs \u2013 SSA: ", prof_ssa$n,
        "  PLADB: ", prof_pladb$n,
        "  TUB: ", prof_tub$n,
        "  Non-significant: ", prof_bg$n)

grp_lvls <- c("Spliceostatin A", "Pladienolide B", "Tubercidin", "Non-significant")
mean_all <- bind_rows(prof_ssa$mean, prof_pladb$mean, prof_tub$mean, prof_bg$mean) %>%
  mutate(group = factor(group, levels = grp_lvls))

profile_pal <- c(
  "Spliceostatin A" = "#706993",
  "Pladienolide B"  = "#70A0AF",
  "Tubercidin"      = "#A0C1B9",
  "Non-significant" = "grey60"
)

p_gc_line <- ggplot(mean_all, aes(x = bin, y = gc, color = group)) +
  annotate("rect", xmin = exon_left, xmax = exon_right,
           ymin = -Inf, ymax = Inf, fill = "grey92") +
  geom_line(linewidth = 1.3) +
  scale_color_manual(
    values = profile_pal,
    labels = function(x) {
      sapply(x, function(lab) {
        n_ev <- switch(lab,
          "Spliceostatin A" = prof_ssa$n,
          "Pladienolide B"  = prof_pladb$n,
          "Tubercidin"      = prof_tub$n,
          "Non-significant" = prof_bg$n)
        paste0(lab, " (n=", format(n_ev, big.mark = ","), ")")
      })
    }
  ) +
  scale_x_continuous(
    breaks = c(1, flank_bins + exon_bins / 2 + 0.5, total_bins),
    labels = c("-100nt", "", "+100nt"),
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(labels = function(x) paste0(round(x), "%")) +
  labs(x = NULL, y = "GC content",
       title = "GC content profile", color = NULL) +
  my_theme +
  theme(legend.position  = "right",
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank(),
        axis.line.x      = element_blank(),
        plot.margin       = margin(8, 8, 2, 8))

heat_all <- bind_rows(prof_ssa$heat, prof_pladb$heat, prof_tub$heat) %>%
  mutate(group = factor(group,
    levels = c("Spliceostatin A", "Pladienolide B", "Tubercidin")))

p_gc_heat <- ggplot(heat_all, aes(x = bin, y = rank, fill = gc)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradientn(
    colours = c("#313695", "#4575B4", "#ABD9E9", "#FFFFBF",
                "#FEE090", "#F46D43", "#A50026"),
    limits = c(10, 80), na.value = "grey50",
    name = "% GC",
    labels = function(x) paste0(x, "%")
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(
    breaks = c(1, flank_bins + exon_bins / 2 + 0.5, total_bins),
    labels = c("-100nt", "", "+100nt"),
    expand = c(0.01, 0)
  ) +
  labs(x = NULL, y = NULL) +
  my_theme +
  theme(axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank(),
        strip.text.y = element_text(angle = 0, size = rel(0.75)),
        plot.margin  = margin(2, 8, 2, 8))

p_exon_bar <- ggplot() +
  annotate("segment", x = exon_left, xend = exon_right,
           y = 0.5, yend = 0.5,
           linewidth = 8, lineend = "butt", color = "black") +
  annotate("text", x = flank_bins / 2 + 0.5, y = 0.5,
           label = "intron", size = 3.2) +
  annotate("text", x = flank_bins + exon_bins / 2 + 0.5, y = 0.5,
           label = "Exon", size = 3.2, fontface = "bold", color = "white") +
  annotate("text", x = flank_bins + exon_bins + flank_bins / 2 + 0.5, y = 0.5,
           label = "intron", size = 3.2) +
  scale_x_continuous(limits = c(0.5, total_bins + 0.5), expand = c(0.01, 0)) +
  ylim(0, 1) +
  theme_void() +
  theme(plot.margin = margin(0, 8, 4, 8))

p_c <- p_gc_line / p_gc_heat / p_exon_bar +
  plot_layout(heights = c(1.2, 3, 0.2))

# ── COMPOSITE ───────────────────────────────────────────────────────────
message("\n── Assembling composite ──")

composite <- (p_a | p_b) /
             wrap_elements(p_c) +
  plot_annotation(
    tag_levels = "A",
    title = "The Oocyte Speckle-Dependency Landscape",
    theme = theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5,
                                family = base_family)
    )
  ) +
  plot_layout(heights = c(1, 1.5))

# ── SAVE ────────────────────────────────────────────────────────────────
save_plot(file.path(out_dir, "speckle_panelA.png"),
         p_a, width = 7, height = 5.5)
save_plot(file.path(out_dir, "speckle_panelA.pdf"),
         p_a, device = cairo_pdf, width = 7, height = 5.5)

save_plot(file.path(out_dir, "speckle_panelB_ecdf.png"),
         p_b, width = 7, height = 5)
save_plot(file.path(out_dir, "speckle_panelB_ecdf.pdf"),
         p_b, device = cairo_pdf, width = 7, height = 5)

save_plot(file.path(out_dir, "speckle_panelC_gc_profile.png"),
         p_c, width = 7, height = 9)
save_plot(file.path(out_dir, "speckle_panelC_gc_profile.pdf"),
         p_c, device = cairo_pdf, width = 7, height = 9)

save_plot(file.path(out_dir, "speckle_landscape_composite.png"),
         composite, width = 14, height = 14)
save_plot(file.path(out_dir, "speckle_landscape_composite.pdf"),
         composite, device = cairo_pdf, width = 14, height = 14)

message("\n\u2713 All outputs saved to ", out_dir, "/")
message("  Total significant events: ", nrow(sig_plot))
message("\nDone!")
