# ============================================================================
# Figure S4 — Length + GC content panels for included/skipped exons
#                and retained/excised introns, all three drugs.
#
# Purpose: address reviewer comments about
#   (1) GC content (GCC) not being reported for exons in Fig S4
#   (2) confusion about the same features being informative for both
#       included and skipped events.
#
# We re-derive Mann-Whitney U p-values directly on the up_vs_down contrast
# (using ONLY the differentially spliced events in the Matt per-event tables)
# rather than the up_vs_ndiff / down_vs_ndiff tests shown in Matt's overview,
# which the reviewer found confusing.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(ggpubr)
  library(grid)
  library(gridExtra)
})

# ---- fonts (fallback to sans if showtext unavailable) -------------------
preferred_font <- "Roboto"
base_family <- if (requireNamespace("showtext", quietly = TRUE)) {
  tryCatch({
    library(showtext)
    font_add_google(preferred_font)
    showtext::showtext_opts(dpi = 600)
    showtext_auto()
    preferred_font
  }, error = function(e) "sans")
} else {
  "sans"
}

# ---- shared palette + theme ---------------------------------------------
splicing_exon_colors <- c(
  "Included" = "#D4A574",   # warm sand
  "Skipped"  = "#B85450"    # muted coral
)
splicing_intron_colors <- c(
  "Retained" = "#D4A574",
  "Excised"  = "#B85450"
)
drug_colors <- c(
  "Tubercidin"        = "#A0C1B9",
  "Pladienolide B"    = "#70A0AF",
  "Spliceostatin A"   = "#706993"
)

theme_cellpub <- function(base_size = 14, base_family_in = base_family) {
  theme_classic(base_size = base_size, base_family = base_family_in) %+replace%
    theme(
      axis.line     = element_line(linewidth = 0.6, colour = "black"),
      axis.ticks    = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      axis.title    = element_text(face = "bold", size = rel(1.0), family = base_family_in),
      axis.text     = element_text(size = rel(0.85), colour = "black", family = base_family_in),
      legend.position = "none",
      legend.background = element_blank(),
      legend.title = element_text(face = "bold", size = rel(0.95)),
      legend.text  = element_text(size = rel(0.85), family = base_family_in),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5),
      strip.text   = element_text(face = "bold", size = rel(0.95), family = base_family_in),
      plot.title   = element_text(face = "bold", size = rel(1.05), hjust = 0.5, family = base_family_in),
      plot.caption = element_text(size = rel(0.8),  colour = "#666666", family = base_family_in),
      plot.margin  = margin(8, 8, 8, 8)
    )
}
my_theme <- theme_cellpub()

# ---- data ----------------------------------------------------------------
matt_root <- "../SupplementaryData1_matt_out"

read_matt <- function(drug_code, event_type) {
  fname <- if (event_type == "exons") {
    paste0(matt_root, "/", drug_code, "_ex/", drug_code, "_exons_with_efeatures.tab")
  } else {
    paste0(matt_root, "/", drug_code, "_int/", drug_code, "_introns_with_ifeatures.tab")
  }
  read_tsv(fname, show_col_types = FALSE) %>%
    mutate(drug = drug_code, event_type = event_type)
}

events_long <- bind_rows(
  read_matt("pladb", "exons"),
  read_matt("ssa",   "exons"),
  read_matt("tub",   "exons"),
  read_matt("pladb", "introns"),
  read_matt("ssa",   "introns"),
  read_matt("tub",   "introns")
) %>%
  filter(DATASET %in% c("up", "down")) %>%
  mutate(
    event_label = case_when(
      event_type == "exons"  & DATASET == "up"   ~ "Included",
      event_type == "exons"  & DATASET == "down" ~ "Skipped",
      event_type == "introns" & DATASET == "up"  ~ "Retained",
      event_type == "introns" & DATASET == "down"~ "Excised",
      TRUE ~ NA_character_
    ),
    drug_label = recode(drug,
                        "pladb" = "Pladienolide B",
                        "ssa"   = "Spliceostatin A",
                        "tub"   = "Tubercidin"),
    drug_label = factor(drug_label,
                        levels = c("Tubercidin", "Pladienolide B", "Spliceostatin A")),
    event_label = factor(event_label,
                         levels = c("Included",  "Skipped",
                                    "Retained",  "Excised"))
  )

cat("Differential events retained (no ndiff):\n")
print(events_long %>% count(event_type, drug_label, event_label))

# ---- helpers -------------------------------------------------------------
gc_plot_exons   <- function(df) ggplot(df, aes(x = event_label, y = GCC, fill = event_label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.55, linewidth = 0.5, color = "black") +
  geom_jitter(width = 0.18, alpha = 0.55, size = 1.6, color = "black") +
  scale_fill_manual(values = splicing_exon_colors) +
  scale_y_continuous(limits = c(0.35, 0.65), expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = NULL, y = "GCC") +
  my_theme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))

length_plot_exons <- function(df) ggplot(df, aes(x = event_label, y = LEN, fill = event_label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.55, linewidth = 0.5, color = "black") +
  geom_jitter(width = 0.18, alpha = 0.55, size = 1.6, color = "black") +
  scale_fill_manual(values = splicing_exon_colors) +
  scale_y_log10() +
  labs(x = NULL, y = "Exon length (nt, log10)") +
  my_theme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))

gc_plot_introns   <- function(df) ggplot(df, aes(x = event_label, y = GCC, fill = event_label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.55, linewidth = 0.5, color = "black") +
  geom_jitter(width = 0.18, alpha = 0.55, size = 1.6, color = "black") +
  scale_fill_manual(values = splicing_intron_colors) +
  scale_y_continuous(limits = c(0.30, 0.65), expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = NULL, y = "GCC") +
  my_theme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))

length_plot_introns <- function(df) ggplot(df, aes(x = event_label, y = LEN, fill = event_label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.55, linewidth = 0.5, color = "black") +
  geom_jitter(width = 0.18, alpha = 0.55, size = 1.6, color = "black") +
  scale_fill_manual(values = splicing_intron_colors) +
  scale_y_log10() +
  labs(x = NULL, y = "Intron length (nt, log10)") +
  my_theme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))

wilcox_p <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 3 || length(y) < 3) return(NA_real_)
  tryCatch(wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_)
}

fmt_p <- function(p) {
  ifelse(is.na(p), "ns",
         ifelse(p <= 0.001, "***",
                ifelse(p <= 0.01,  "**",
                       ifelse(p <= 0.05,  "*", "ns"))))
}

# ---- safe plot saver (handles missing cairo) ----------------------------
save_plot <- function(plot, base, w, h) {
  png_path <- paste0(base, ".png")
  pdf_path <- paste0(base, ".pdf")
  ggsave(png_path, plot = plot, width = w, height = h, dpi = 600, bg = "white")
  if (capabilities("cairo") && requireNamespace("Cairo", quietly = TRUE)) {
    suppressMessages(ggsave(pdf_path, plot = plot, width = w, height = h,
                            device = cairo_pdf, bg = "white"))
  } else {
    ggsave(pdf_path, plot = plot, width = w, height = h, bg = "white")
  }
  message("Saved: ", png_path, " & ", pdf_path)
}

# ---- Figure S4 panel builder --------------------------------------------
build_exon_panel <- function(drug_label) {
  df_ex <- events_long %>%
    filter(event_type == "exons", drug_label == !!drug_label) %>%
    rename(LEN = EXON_LENGTH, GCC = EXON_GCC)

  p_gc   <- gc_plot_exons(df_ex)
  p_len  <- length_plot_exons(df_ex)

  p_gc_b   <- p_gc   + stat_compare_means(method = "wilcox.test", ref.group = NULL,
                                          comparisons = list(c("Included", "Skipped")),
                                          label = "p.format", label.size = 3.2,
                                          tip.length = 0.02)
  p_len_b  <- p_len  + stat_compare_means(method = "wilcox.test", ref.group = NULL,
                                          comparisons = list(c("Included", "Skipped")),
                                          label = "p.format", label.size = 3.2,
                                          tip.length = 0.02)

  p <- (p_gc_b | p_len_b) +
    plot_annotation(title = drug_label,
                    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5,
                                                             family = base_family, size = 14)))
  p
}

build_intron_panel <- function(drug_label) {
  df_int <- events_long %>%
    filter(event_type == "introns", drug_label == !!drug_label) %>%
    rename(LEN = INTRON_LENGTH, GCC = INTRON_GCC)

  p_gc   <- gc_plot_introns(df_int)
  p_len  <- length_plot_introns(df_int)

  p_gc_b   <- p_gc   + stat_compare_means(method = "wilcox.test",
                                          comparisons = list(c("Retained", "Excised")),
                                          label = "p.format", label.size = 3.2,
                                          tip.length = 0.02)
  p_len_b  <- p_len  + stat_compare_means(method = "wilcox.test",
                                          comparisons = list(c("Retained", "Excised")),
                                          label = "p.format", label.size = 3.2,
                                          tip.length = 0.02)

  p <- (p_gc_b | p_len_b) +
    plot_annotation(title = drug_label,
                    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5,
                                                             family = base_family, size = 14)))
  p
}

exon_panels  <- lapply(levels(events_long$drug_label), build_exon_panel)
intron_panels <- lapply(levels(events_long$drug_label), build_intron_panel)

fig_s4g_exons <- wrap_plots(exon_panels, nrow = 1) +
  plot_annotation(title = "Figure S4G — Exons: GC content and length",
                  subtitle = "Mann–Whitney U tests are run only between Included and Skipped events (no unchanged pool).",
                  theme = theme(plot.title = element_text(face = "bold", family = base_family, size = 16),
                                plot.subtitle = element_text(family = base_family, size = 11, colour = "grey40")))

fig_s4h_introns <- wrap_plots(intron_panels, nrow = 1) +
  plot_annotation(title = "Figure S4H — Introns: GC content and length",
                  subtitle = "Mann–Whitney U tests are run only between Retained and Excised introns (no unchanged pool).",
                  theme = theme(plot.title = element_text(face = "bold", family = base_family, size = 16),
                                plot.subtitle = element_text(family = base_family, size = 11, colour = "grey40")))

save_plot(fig_s4g_exons,   "FigureS4G_exons_length_gcc",   13, 4.2)
save_plot(fig_s4h_introns, "FigureS4H_introns_length_gcc", 13, 4.2)
message("Saved FigureS4G/H panels")

# ---- scatter panels: length vs GCC (the "definition" panel) -------------
# (helper for title case, declared early so build_scatter can use it)
str_to_title <- function(x) {
  parts <- strsplit(x, " ")[[1]]
  paste(toupper(substring(parts, 1, 1)), tolower(substring(parts, 2)), sep = "", collapse = " ")
}

build_scatter <- function(drug_label, event_type_nice, length_col, gcc_col, label_a, label_b) {
  df <- events_long %>%
    filter(drug_label == !!drug_label, event_type == !!event_type_nice)
  colnames(df)[colnames(df) == length_col] <- "LENGTH_VAL"
  colnames(df)[colnames(df) == gcc_col]    <- "GCC_VAL"

  col_pal <- if (event_type_nice == "exons") splicing_exon_colors else splicing_intron_colors

  p <- ggplot(df, aes(x = LENGTH_VAL, y = GCC_VAL, color = event_label)) +
    geom_point(alpha = 0.65, size = 2.4) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.8,
                linetype = "dashed", alpha = 0.6) +
    scale_x_log10() +
    scale_color_manual(values = col_pal) +
    labs(x = paste0(str_to_title(sub("_", " ", length_col)), " (nt, log10)"),
         y = str_to_title(sub("_", " ", gcc_col)),
         title = drug_label) +
    my_theme +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, family = base_family, size = 13),
          legend.position = "bottom",
          legend.title = element_blank())
  p
}

# (helper for title case, dplyr/tidyr already load stringr-free version using base R)
str_to_title <- function(x) {
  parts <- strsplit(x, " ")[[1]]
  paste(toupper(substring(parts, 1, 1)), tolower(substring(parts, 2)), sep = "", collapse = " ")
}

exon_scats  <- lapply(levels(events_long$drug_label),
                      function(d) build_scatter(d, "exons",  "EXON_LENGTH",  "EXON_GCC"))
intron_scats <- lapply(levels(events_long$drug_label),
                       function(d) build_scatter(d, "introns", "INTRON_LENGTH", "INTRON_GCC"))

fig_def_exons <- wrap_plots(exon_scats, nrow = 1) +
  plot_annotation(title = "Exon definition landscape — length vs GC",
                  subtitle = paste0("Included (n=", sum(events_long$event_type == "exons" & events_long$event_label == "Included"), ") ",
                                    "vs Skipped (n=", sum(events_long$event_type == "exons" & events_long$event_label == "Skipped"), "). ",
                                    "Dashed lines: linear fit on log10(length)."),
                  theme = theme(plot.title = element_text(face = "bold", family = base_family, size = 15),
                                plot.subtitle = element_text(family = base_family, size = 10, colour = "grey40")))

fig_def_introns <- wrap_plots(intron_scats, nrow = 1) +
  plot_annotation(title = "Intron definition landscape — length vs GC",
                  subtitle = paste0("Retained (n=", sum(events_long$event_type == "introns" & events_long$event_label == "Retained"), ") ",
                                    "vs Excised (n=", sum(events_long$event_type == "introns" & events_long$event_label == "Excised"), "). ",
                                    "Dashed lines: linear fit on log10(length)."),
                  theme = theme(plot.title = element_text(face = "bold", family = base_family, size = 15),
                                plot.subtitle = element_text(family = base_family, size = 10, colour = "grey40")))

# Safe save — skip cairo_pdf if device unavailable
save_plot <- function(plot, base, w, h) {
  png_path <- paste0(base, ".png")
  pdf_path <- paste0(base, ".pdf")
  ggsave(png_path, plot = plot, width = w, height = h, dpi = 600, bg = "white")
  if (capabilities("cairo") && requireNamespace("Cairo", quietly = TRUE)) {
    suppressMessages(ggsave(pdf_path, plot = plot, width = w, height = h,
                            device = cairo_pdf, bg = "white"))
  } else {
    ggsave(pdf_path, plot = plot, width = w, height = h, bg = "white")
  }
  message("Saved: ", png_path, " & ", pdf_path)
}

save_plot(fig_def_exons,   "Fig_definition_exons_length_vs_gcc",   13, 4.5)
save_plot(fig_def_introns, "Fig_definition_introns_length_vs_gcc", 13, 4.5)
message("Saved definition scatter panels")

# ---- SUMMARY TABLE  ------------------------------------------------------
# For every (drug, event_type, feature) we compute the up_vs_down
# Mann-Whitney p-value directly (no ndiff involved).

features_exon <- c("EXON_LENGTH", "EXON_GCC",
                   "UPEXON_GCC", "DOEXON_GCC",
                   "UPINTRON_GCC", "DOINTRON_GCC",
                   "GCC_5SS_20INT10EX", "GCC_3SS_20INT10EX",
                   "UPINTRON_MEDIANLENGTH", "DOINTRON_MEDIANLENGTH")

features_intron <- c("INTRON_LENGTH", "INTRON_GCC",
                     "UPEXON_GCC", "DOEXON_GCC",
                     "INTRON_5SS_20INT10EX_GCC", "INTRON_3SS_20INT10EX_GCC",
                     "UPEXON_MEDIANLENGTH", "DOEXON_MEDIANLENGTH")

compute_p_row <- function(df, feature) {
  if (!feature %in% colnames(df)) return(NA_real_)
  vals_a <- df %>% filter(DATASET == "up")   %>% pull(all_of(feature))
  vals_b <- df %>% filter(DATASET == "down") %>% pull(all_of(feature))
  wilcox_p(vals_a, vals_b)
}

compute_stats_table <- function(event_type_nice, features) {
  rows <- expand_grid(drug_label = levels(events_long$drug_label),
                      feature = features)
  rows$n_up   <- NA_integer_
  rows$n_down <- NA_integer_
  rows$med_up   <- NA_real_
  rows$med_down <- NA_real_
  rows$p_up_vs_down <- NA_real_

  for (i in seq_len(nrow(rows))) {
    d <- rows$drug_label[i]
    f <- rows$feature[i]
    sub <- events_long %>% filter(drug_label == d,
                                  event_type == event_type_nice,
                                  DATASET %in% c("up", "down"))
    rows$n_up[i]   <- sum(sub$DATASET == "up")
    rows$n_down[i] <- sum(sub$DATASET == "down")
    if (f %in% colnames(events_long)) {
      rows$med_up[i]   <- median(sub[[f]][sub$DATASET == "up"],   na.rm = TRUE)
      rows$med_down[i] <- median(sub[[f]][sub$DATASET == "down"], na.rm = TRUE)
      rows$p_up_vs_down[i] <- wilcox_p(sub[[f]][sub$DATASET == "up"],
                                       sub[[f]][sub$DATASET == "down"])
    }
  }

  rows %>%
    mutate(
      event_type = event_type_nice,
      significance = fmt_p(p_up_vs_down),
      median_delta = med_up - med_down
    )
}

summary_exons  <- compute_stats_table("exons",  features_exon)
summary_introns <- compute_stats_table("introns", features_intron)

summary_all <- bind_rows(summary_exons, summary_introns) %>%
  arrange(event_type, drug_label, feature)

write_csv(summary_all, "summary_up_vs_down_pvalues.csv")
message("Saved summary table → summary_up_vs_down_pvalues.csv  (", nrow(summary_all), " rows)")

# ---- small print of the key rows for the manuscript ---------------------
key_rows <- summary_all %>%
  filter(feature %in% c("EXON_LENGTH", "EXON_GCC",
                        "INTRON_LENGTH", "INTRON_GCC")) %>%
  select(event_type, drug_label, feature, n_up, n_down, med_up, med_down, p_up_vs_down, significance)
print(as.data.frame(key_rows))

cat("\nDone.\n")
