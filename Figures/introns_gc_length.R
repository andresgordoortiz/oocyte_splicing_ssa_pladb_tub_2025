library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(showtext)

# 1) Read event info
event_info <- read.delim("EVENT_INFO-mm10.tab")

set.seed(123)
# Import betAS output tables
tub_fdr_df<-read_csv("tub_fdr.csv")[,-1]
pladb_fdr_df<-read_csv("pladb_fdr.csv")[,-1]
ssa_fdr_df<-read_csv("ssa_fdr.csv")[,-1]

differential_tub<-na.omit(tub_fdr_df[tub_fdr_df$FDR <= 0.05 & abs(tub_fdr_df$deltapsi) >= 0.1,])

differential_pladb<-na.omit(pladb_fdr_df[pladb_fdr_df$FDR <= 0.05 & abs(pladb_fdr_df$deltapsi) >= 0.1,])

differential_ssa<-na.omit(ssa_fdr_df[ssa_fdr_df$FDR <= 0.05 & abs(ssa_fdr_df$deltapsi) >= 0.1,])

# Save separate GC and Length intron plots for SSA, Tubercidin, and Pladb — PNG + PDF (high-res)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)
library(showtext)
library(grid)    # unit()
library(gridExtra)

# ---- showtext / font (high DPI) ----
preferred_font <- "Roboto"
font_add_google(preferred_font)
showtext::showtext_opts(dpi = 600)   # <- high DPI for consistent output
showtext_auto()
base_family <- preferred_font

# ---- styling utilities (harmonized with rbp_full_analysis.R) ----
# Color palette matching rbp_full_analysis.R
splicing_colors <- c(
  "Retained" = "#D4A574",  # Warm sand
  "Excised" = "#B85450",   # Muted coral
  "Unchanged" = "#808080"  # Gray
)
theme_cellpub <- function(base_size = 16, base_family_in = NULL) {
  if (is.null(base_family_in)) base_family_in <- get0("base_family", ifnotfound = "sans")
  theme_classic(base_size = base_size, base_family = base_family_in) %+replace%
    theme(
      # Axes - matching rbp_full_analysis.R
      axis.line = element_line(linewidth = 0.6, colour = "black"),
      axis.ticks = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      axis.title = element_text(family = base_family_in, face = "bold", size = rel(1.0)),
      axis.text = element_text(family = base_family_in, size = rel(0.9), colour = "black"),

      # Legend
      legend.position = "none",
      legend.background = element_blank(),
      legend.title = element_text(face = "bold", size = rel(0.95)),
      legend.text = element_text(family = base_family_in, size = rel(0.85)),
      legend.key.size = unit(0.4, "cm"),

      # Panel
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),

      # Facets - matching rbp_full_analysis.R
      strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = rel(1.0), family = base_family_in, margin = margin(5, 5, 5, 5)),

      # Plot
      plot.title = element_text(face = "bold", size = rel(1.1), hjust = 0.5, family = base_family_in),
      plot.subtitle = element_text(size = rel(0.95), hjust = 0.5, family = base_family_in),
      plot.caption = element_text(size = rel(0.85), colour = "#666666", family = base_family_in),
      plot.margin = margin(6, 6, 6, 6)
    )
}
my_theme <- theme_cellpub(base_size = 16, base_family_in = base_family)

# ---- helper: GC function ----
gc_content <- function(seqs) {
  sapply(seqs, function(seq) {
    if (is.na(seq)) return(NA_real_)
    seq <- toupper(as.character(seq))
    chars <- strsplit(seq, "")[[1]]
    gc  <- sum(chars %in% c("G","C"))
    n   <- length(chars)
    if(n==0) return(NA_real_)
    100 * gc / n
  })
}

# ---- main plotting function (saves separate GC + Length PNG + PDF) ----
make_and_save_introns <- function(dataset_name,
                                  differential_df, fdr_df, event_info,
                                  out_prefix = "introns", sample_max = 10000) {
  ds_label <- toupper(dataset_name)
  message("Processing: ", ds_label)
  # differential intronic events
  introns_gc <- differential_df %>%
    filter(grepl("INT", EVENT)) %>%
    left_join(dplyr::select(event_info, EVENT, Seq_A, LE_o), by = "EVENT") %>%
    mutate(gc = gc_content(Seq_A),
           splicing = if_else(deltapsi < 0, "Excised", "Retained"))
  if(nrow(introns_gc) == 0) {
    message("  -> no differential introns for ", ds_label, "; skipping.")
    return(invisible(NULL))
  }
  # unchanged pool (sample)
  unchanged_pool <- fdr_df %>% filter(grepl("INT", EVENT)) %>%
    filter(!EVENT %in% differential_df$EVENT) %>%
    left_join(dplyr::select(event_info, EVENT, Seq_A, LE_o), by = "EVENT") %>%
    filter(!is.na(Seq_A))
  sample_size <- min(sample_max, nrow(unchanged_pool))
  if(sample_size == 0) {
    message("  -> no unchanged introns to sample for ", ds_label, "; skipping.")
    return(invisible(NULL))
  }
  # sampling reproducibly (global set.seed already set)
  unchanged_gc <- unchanged_pool %>%
    sample_n(size = sample_size) %>%
    mutate(gc = gc_content(Seq_A), splicing = "Unchanged")
  all_introns <- bind_rows(introns_gc, unchanged_gc) %>%
    mutate(splicing = factor(splicing, levels = c("Retained","Excised","Unchanged")))
  if(all(is.na(all_introns$gc))) {
    message("  -> no GC values available for ", ds_label, "; skipping.")
    return(invisible(NULL))
  }
  
  # comparisons & palette
  comparisons <- list(
    c("Retained", "Excised"),
    c("Retained", "Unchanged"),
    c("Excised",  "Unchanged")
  )
  
  # ---- GC plot ----
  p_gc <- ggplot(all_introns, aes(x = splicing, y = gc, fill = splicing)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, linewidth = 0.5) +
    # remove jitter for Unchanged only
    geom_jitter(
      data = dplyr::filter(all_introns, splicing != "Unchanged"),
      aes(x = splicing, y = gc, color = splicing),
      width = 0.15, alpha = 0.5, size = 1.2
    ) +
    scale_color_manual(values = splicing_colors) +
    stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif", label.size = 4) +
    scale_fill_manual(values = splicing_colors) +
    labs(x = "Splicing Outcome", y = "GC Content (%)",
         title = paste0(ds_label, ": Introns — GC content")) +
    my_theme
  
  # ---- Length plot (statistics on original data, display clipped) ----
  # Remove extreme outliers by using 95th percentile clipping for display only
  threshold_95 <- quantile(all_introns$LE_o, 0.95, na.rm = TRUE)
  
  # Create plot data - keep original values for statistics
  plot_data <- all_introns %>%
    filter(!is.na(LE_o)) %>%
    mutate(
      LE_o_display = pmin(LE_o, threshold_95),  # Cap display values at 95th percentile
      is_outlier = LE_o > threshold_95
    )
  
  if(nrow(plot_data) == 0) {
    p_len <- NULL
    message("  -> no length data for ", ds_label, "; skipping length plot.")
  } else {
    message("  -> creating length plot for ", ds_label, " with ", nrow(plot_data), " introns")
    
    # Calculate statistics on ORIGINAL data (including outliers)
    stat_results <- purrr::map_dfr(comparisons, function(comp) {
      g1 <- comp[1]; g2 <- comp[2]
      
      # Use original LE_o values for statistical test
      data_g1 <- plot_data$LE_o[plot_data$splicing == g1]
      data_g2 <- plot_data$LE_o[plot_data$splicing == g2]
      
      if(length(data_g1) < 1 || length(data_g2) < 1) {
        return(tibble(group1 = g1, group2 = g2, p.value = NA_real_))
      }
      
      test_result <- tryCatch({
        wilcox.test(data_g1, data_g2, exact = FALSE)
      }, error = function(e) {
        list(p.value = NA_real_)
      })
      
      tibble(group1 = g1, group2 = g2, p.value = test_result$p.value)
    })
    
    # Convert p-values to significance symbols
    stat_results <- stat_results %>%
      mutate(
        p.signif = case_when(
          is.na(p.value) ~ "ns",
          p.value <= 0.001 ~ "***",
          p.value <= 0.01 ~ "**",
          p.value <= 0.05 ~ "*",
          TRUE ~ "ns"
        )
      )  # Only show significant results
    
    # Position bars based on CLIPPED data range for visibility
    if(nrow(stat_results) > 0) {
      max_display <- max(plot_data$LE_o_display, na.rm = TRUE)
      step_size <- max_display * 0.08  # 8% steps
      
      stat_results <- stat_results %>%
        mutate(
          y.position = max_display + step_size * row_number(),
          xmin = group1,
          xmax = group2,
          label = p.signif
        )
    }
    
    # Create plot with clipped data for display
    p_len <- ggplot(plot_data, aes(x = splicing, y = LE_o_display, fill = splicing)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, linewidth = 0.5) +
      # remove jitter for Unchanged only (display)
      geom_jitter(
        data = dplyr::filter(plot_data, splicing != "Unchanged"),
        aes(x = splicing, y = LE_o_display, color = splicing),
        width = 0.15, alpha = 0.5, size = 1.2
      ) +
      scale_color_manual(values = splicing_colors) +
      scale_fill_manual(values = splicing_colors) +
      labs(x = "Splicing Outcome", y = "Intron Length (nt)",
           title = paste0(ds_label, ": Introns — Length")) +
      my_theme
    
    # Add significance bars if any exist (these were computed on ORIGINAL LE_o)
    if(nrow(stat_results) > 0) {
      p_len <- p_len +
        stat_pvalue_manual(stat_results,
                           xmin = "group1", xmax = "group2",
                           y.position = "y.position", label = "label",
                           tip.length = 0.02, label.size = 4)
      
      # Adjust y-axis to accommodate significance bars
      max_bar_pos <- max(stat_results$y.position)
      p_len <- p_len + coord_cartesian(ylim = c(NA, max_bar_pos * 1.05))
    }
    
    # Add indicator for clipped outliers if any exist
    n_outliers <- sum(plot_data$is_outlier)
    if(n_outliers > 0) {
      p_len <- p_len +
        labs(caption = paste0("Note: ", n_outliers, " outliers > ", round(threshold_95), " nt clipped for display; statistics calculated on full data"))
    }
  }

  # ---- Save files: separate GC and Length PNG + PDF ----
  gc_png  <- paste0(out_prefix, "_", dataset_name, "_gc.png")
  gc_pdf  <- paste0(out_prefix, "_", dataset_name, "_gc.pdf")
  len_png <- paste0(out_prefix, "_", dataset_name, "_length.png")
  len_pdf <- paste0(out_prefix, "_", dataset_name, "_length.pdf")

  # Save GC
  ggsave(gc_png, plot = p_gc, width = 6, height = 5, dpi = 600)
  ggsave(gc_pdf, plot = p_gc, device = cairo_pdf, width = 6, height = 5)
  message("  -> saved: ", gc_png, " and ", gc_pdf)

  # Save Length (if exists)
  if(!is.null(p_len)) {
    ggsave(len_png, plot = p_len, width = 6, height = 5, dpi = 600)
    ggsave(len_pdf, plot = p_len, device = cairo_pdf, width = 6, height = 5)
    message("  -> saved: ", len_png, " and ", len_pdf)
  }

  invisible(list(gc = p_gc, length = p_len))
}

# ---- Run for SSA, tub (tubercidin), pladb ----
# Assumes objects: differential_ssa, differential_tub, differential_pladb,
# and ssa_fdr_df, tub_fdr_df, pladb_fdr_df and event_info exist.
datasets <- list(
  ssa   = list(diff = differential_ssa, fdr = ssa_fdr_df),
  tub   = list(diff = differential_tub, fdr = tub_fdr_df),
  pladb = list(diff = differential_pladb, fdr = pladb_fdr_df)
)

results <- list()
for(nm in names(datasets)) {
  ds <- datasets[[nm]]
  results[[nm]] <- tryCatch({
    make_and_save_introns(dataset_name = nm,
                          differential_df = ds$diff,
                          fdr_df = ds$fdr,
                          event_info = event_info,
                          out_prefix = "introns")
  }, error = function(e) {
    message("ERROR processing ", nm, ": ", e$message)
    NULL
  })
}

message("Done. Individual plots saved for datasets: ", paste(names(Filter(Negate(is.null), results)), collapse = ", "))

# ---- Prepare combined data for unified faceted plots ----
message("\nPreparing combined introns data for unified plot...")

# Drug display names
drug_names <- c(
  "ssa" = "Spliceostatin A",
  "tub" = "Tubercidin",
  "pladb" = "Pladienolide B"
)

# Combine all data
all_drugs_data <- bind_rows(
  # SSA
  differential_ssa %>%
    filter(grepl("INT", EVENT)) %>%
    left_join(dplyr::select(event_info, EVENT, Seq_A, LE_o), by = "EVENT") %>%
    mutate(splicing = if_else(deltapsi < 0, "Excised", "Retained"), drug = "ssa"),

  # TUB
  differential_tub %>%
    filter(grepl("INT", EVENT)) %>%
    left_join(dplyr::select(event_info, EVENT, Seq_A, LE_o), by = "EVENT") %>%
    mutate(splicing = if_else(deltapsi < 0, "Excised", "Retained"), drug = "tub"),

  # PLADB
  differential_pladb %>%
    filter(grepl("INT", EVENT)) %>%
    left_join(dplyr::select(event_info, EVENT, Seq_A, LE_o), by = "EVENT") %>%
    mutate(splicing = if_else(deltapsi < 0, "Excised", "Retained"), drug = "pladb")
) %>%
  mutate(gc = gc_content(Seq_A))

# Add unchanged samples for each drug
sample_max <- 10000
for(drug_code in c("ssa", "tub", "pladb")) {
  fdr_df <- get(paste0(drug_code, "_fdr_df"))
  diff_df <- get(paste0("differential_", drug_code))
  
  unchanged_pool <- fdr_df %>%
    filter(grepl("INT", EVENT)) %>%
    filter(!EVENT %in% diff_df$EVENT) %>%
    left_join(dplyr::select(event_info, EVENT, Seq_A, LE_o), by = "EVENT") %>%
    filter(!is.na(Seq_A))
  
  sample_size <- min(sample_max, nrow(unchanged_pool))
  
  if(sample_size > 0) {
    # sample_n reproducible because of global set.seed above
    unchanged_data <- unchanged_pool %>%
      sample_n(size = sample_size) %>%
      mutate(
        splicing = "Unchanged",
        drug = drug_code,
        gc = gc_content(Seq_A)
      )
    all_drugs_data <- bind_rows(all_drugs_data, unchanged_data)
  }
}

# Set factor levels (same as before)
all_drugs_data <- all_drugs_data %>%
  mutate(
    splicing = factor(splicing, levels = c("Retained", "Excised", "Unchanged")),
    drug = factor(drug, levels = c("tub", "pladb", "ssa")),
    drug_label = factor(drug_names[as.character(drug)], levels = drug_names[c("tub", "pladb", "ssa")])
  )

comparisons <- list(
  c("Retained", "Excised"),
  c("Retained", "Unchanged"),
  c("Excised",  "Unchanged")
)

# ---- Unified GC plot: remove jitter for Unchanged only ----
p_gc_unified <- ggplot(all_drugs_data, aes(x = splicing, y = gc, fill = splicing)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.9, linewidth = 0.5, color = "black") +
  geom_jitter(
    data = dplyr::filter(all_drugs_data, splicing != "Unchanged"),
    aes(x = splicing, y = gc, color = splicing),
    width = 0.15, size = 0.8, alpha = 0.5
  ) +
  scale_color_manual(values = splicing_colors) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = comparisons,
    label = "p.signif",
    label.size = 3.5,
    y.position = 60
  ) +
  scale_fill_manual(values = splicing_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +  # Extra space for significance bars
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  facet_wrap(~ drug_label, nrow = 1) +
  labs(x = "Splicing Outcome", y = "GC Content (%)") +
  my_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.85)),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = rel(0.95),
      margin = margin(b = 12, t = 2)
    ),
    panel.spacing = unit(1.2, "lines")
  )

# ---- Unified Length plot: compute p-values on ORIGINAL LE_o per drug and draw significance on clipped display ----
threshold_95 <- quantile(all_drugs_data$LE_o, 0.95, na.rm = TRUE)
plot_data_len <- all_drugs_data %>%
  filter(!is.na(LE_o)) %>%
  mutate(
    LE_o_display = pmin(LE_o, threshold_95),
    is_outlier = LE_o > threshold_95
  )

# compute per-drug p-values using original LE_o
comparisons <- list(
  c("Retained", "Excised"),
  c("Retained", "Unchanged"),
  c("Excised",  "Unchanged")
)

stat_results_len <- plot_data_len %>%
  group_by(drug_label) %>%
  do({
    df <- .
    purrr::map_dfr(comparisons, function(comp) {
      g1 <- comp[1]; g2 <- comp[2]
      data_g1 <- df$LE_o[df$splicing == g1]
      data_g2 <- df$LE_o[df$splicing == g2]
      pval <- NA_real_
      if(length(data_g1) > 0 && length(data_g2) > 0) {
        pval <- tryCatch(wilcox.test(data_g1, data_g2, exact = FALSE)$p.value, error = function(e) NA_real_)
      }
      tibble(group1 = g1, group2 = g2, p.value = pval)
    }) %>%
      mutate(drug_label = unique(df$drug_label))
  }) %>%
  ungroup() %>%
  mutate(
    p.signif = case_when(
      is.na(p.value) ~ "ns",
      p.value <= 0.001 ~ "***",
      p.value <= 0.01 ~ "**",
      p.value <= 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    label = p.signif
  )

# compute y.positions per drug for display (based on clipped display values)
y_positions <- plot_data_len %>%
  group_by(drug_label) %>%
  summarize(max_display = max(LE_o_display, na.rm = TRUE)) %>%
  arrange(drug_label)

# attach y.positions (one row per comparison within each drug)
stat_results_len <- stat_results_len %>%
  left_join(y_positions, by = "drug_label") %>%
  group_by(drug_label) %>%
  mutate(
    step_size = ifelse(is.finite(max_display) & max_display > 0, max_display * 0.08, 1),
    y.position = max_display + step_size * row_number()
  ) %>%
  ungroup() %>%
  select(group1, group2, drug_label, p.value, label, y.position)

# Build length plot (clipped display), jitter only for non-unchanged
p_len_unified <- ggplot(plot_data_len, aes(x = splicing, y = LE_o_display, fill = splicing)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.9, linewidth = 0.5, color = "black") +
  geom_jitter(
    data = dplyr::filter(plot_data_len, splicing != "Unchanged"),
    aes(x = splicing, y = LE_o_display, color = splicing),
    width = 0.15, size = 0.8, alpha = 0.5
  ) +
  scale_color_manual(values = splicing_colors) +
  scale_fill_manual(values = splicing_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +  # Extra space at top for significance
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  facet_wrap(~ drug_label, nrow = 1) +
  labs(x = "Splicing Outcome", y = "Intron Length (nt)") +
  my_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.85)),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = rel(0.95)),
    panel.spacing = unit(0.5, "lines")
  )

# add manual p-values (computed on original LE_o)
if(nrow(stat_results_len) > 0) {
  p_len_unified <- p_len_unified +
    stat_pvalue_manual(stat_results_len,
                       xmin = "group1", xmax = "group2",
                       y.position = "y.position", label = "label",
                       tip.length = 0.02, label.size = 3.5)
}

# Save unified plots as before
ggsave("introns_unified_gc.png", plot = p_gc_unified, width = 10, height = 3.5, dpi = 600)
ggsave("introns_unified_gc.pdf", plot = p_gc_unified, device = cairo_pdf, width = 10, height = 3.5)
message("  -> saved: introns_unified_gc.png and introns_unified_gc.pdf")

ggsave("introns_unified_length.png", plot = p_len_unified, width = 10, height = 3.5, dpi = 600)
ggsave("introns_unified_length.pdf", plot = p_len_unified, device = cairo_pdf, width = 10, height = 3.5)
message("  -> saved: introns_unified_length.png and introns_unified_length.pdf")