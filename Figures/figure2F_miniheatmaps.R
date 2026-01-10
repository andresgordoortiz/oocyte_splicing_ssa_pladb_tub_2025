library(betAS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggforce)
library(patchwork)
library(grid)

splicing_data <- getDataset(pathTables = "Preprocessing/ssa_pladb_tub_oocyte_INCLUSION_LEVELS_FULL-mm10.tab"
                            , tool = "vast-tools")

pladb_data<-splicing_data[,1:18]
ssa_tub_data<-splicing_data[, c(1:6, 19:36)]

pladb_events <- filterEvents(getEvents(pladb_data, tool = "vast-tools"), N = 10)
ssa_tub_events <- filterEvents(getEvents(ssa_tub_data, tool = "vast-tools"), N = 10)

ssa_events<-list()
ssa_events$PSI<-ssa_tub_events$PSI[,c(1:9,13:15)]

tub_events<-list()
tub_events$PSI<-ssa_tub_events$PSI[,c(1:12)]

exons_pladb <- filterEvents(pladb_events, types = c("C1", "C2", "C3", "S", "MIC"), N = 10)
introns_pladb <- filterEvents(pladb_events, types = c("IR"), N = 10)
alt_pladb <- filterEvents(pladb_events, types = c("Alt5", "Alt3"), N = 10)

exons_ssa_tub <- filterEvents(ssa_tub_events, types = c("C1", "C2", "C3", "S", "MIC"), N = 10)
introns_ssa_tub <- filterEvents(ssa_tub_events, types = c("IR"), N = 10)
alt_ssa_tub <- filterEvents(ssa_tub_events, types = c("Alt5", "Alt3"), N = 10)

# Import betAS output tables
excel_file <- "Preprocessing/betAS_out/Supplementary2_betAS_splicing_results.xlsx"

# Read each sheet (using the names we set in the previous step)
tub_fdr_df   <- read_excel(excel_file, sheet = "Tub_FDR")
pladb_fdr_df <- read_excel(excel_file, sheet = "PlaDB_FDR")
ssa_fdr_df   <- read_excel(excel_file, sheet = "SSA_FDR")

# Extract significant events
differential_tub <- na.omit(tub_fdr_df[tub_fdr_df$FDR <= 0.05 & abs(tub_fdr_df$deltapsi) >= 0.1,])
differential_pladb <- na.omit(pladb_fdr_df[pladb_fdr_df$FDR <= 0.05 & abs(pladb_fdr_df$deltapsi) >= 0.1,])
differential_ssa <- na.omit(ssa_fdr_df[ssa_fdr_df$FDR <= 0.05 & abs(ssa_fdr_df$deltapsi) >= 0.1,])

differential_pladb_exons<-filter(differential_pladb,grepl("EX",differential_pladb$EVENT))
differential_tub_exons<-filter(differential_tub,grepl("EX",differential_tub$EVENT))
differential_ssa_exons<-filter(differential_ssa,grepl("EX",differential_ssa$EVENT))

differential_pladb_introns<-filter(differential_pladb,grepl("INT",differential_pladb$EVENT))
differential_tub_introns<-filter(differential_tub,grepl("INT",differential_tub$EVENT))
differential_ssa_introns<-filter(differential_ssa,grepl("INT",differential_ssa$EVENT))

set.seed(44)

# ---- Updated color palettes ----
# Exons: white -> deep rose (#7A2048)
exon_palette <- c("#E0E0E0", "#E0E0E0", "#7A2048")
# Introns: white -> charcoal black (#2C2C2C)
intron_palette <- c("#E0E0E0", "#E0E0E0", "#2C2C2C")

# ---- Helper: robust row-wise MAD scaling ----
robust_scale_rows <- function(mat) {
  t(apply(mat, 1, function(x) {
    med <- median(x, na.rm = TRUE)
    mad_val <- mad(x, center = med, na.rm = TRUE)
    if (mad_val == 0) mad_val <- 1
    (x - med) / mad_val
  }))
}

# ---- Core: build collapsed summary with Z-scores ----
prepare_collapsed <- function(events_obj, diff_exons, diff_introns,
                              keep_cols = 7:12, pct_threshold = 0.1, n_clusters = 2) {
  if (is.null(events_obj) || is.null(events_obj$PSI)) stop("events_obj with $PSI required")
  psi_df <- events_obj$PSI %>% filter(complete.cases(dplyr::select(., all_of(keep_cols))))
  uniq <- psi_df %>% distinct(EVENT, .keep_all = TRUE)
  rownames(uniq) <- uniq$EVENT
  
  inner <- function(diff_table) {
    if (is.null(diff_table)) return(NULL)
    evts <- diff_table$EVENT[abs(diff_table$deltapsi) >= pct_threshold]
    if (length(evts) == 0) return(NULL)
    mat <- uniq[rownames(uniq) %in% evts, , drop = FALSE] %>% dplyr::select(all_of(keep_cols)) %>% as.matrix()
    if (nrow(mat) == 0) return(NULL)
    scaled <- robust_scale_rows(mat)
    # Scale to -1 to 1 range instead of -2 to 2
    scaled <- pmax(pmin(scaled, 1), -1)
    if (nrow(scaled) > 1) {
      dist_matrix <- as.dist(1 - cor(t(scaled), method = "spearman"))
      hc <- hclust(dist_matrix, method = "ward.D2")
      clusters <- cutree(hc, k = min(n_clusters, nrow(scaled)))
    } else clusters <- rep(1, nrow(scaled))
    collapsed <- data.frame(
      Control = rowMeans(scaled[, 1:3, drop = FALSE], na.rm = TRUE),
      Treat   = rowMeans(scaled[, 4:6, drop = FALSE], na.rm = TRUE),
      Cluster = clusters,
      EVENT   = rownames(scaled),
      stringsAsFactors = FALSE
    )
    collapsed %>%
      group_by(Cluster) %>%
      summarise(
        Control = mean(Control, na.rm = TRUE),
        Treat = mean(Treat, na.rm = TRUE),
        n_events = n(), .groups = "drop"
      ) %>%
      arrange(Cluster)
  }
  
  list(exons = inner(diff_exons), introns = inner(diff_introns))
}

# ---- Plot builder with actual PSI values ----
make_tile_plot <- function(mini_heatmap_data, treat_label, palette_colors, legend_title, show_legend = TRUE) {
  if (is.null(mini_heatmap_data) || nrow(mini_heatmap_data) == 0) {
    return(ggplot() + annotate("text", x = 1, y = 1, label = paste("No events\nfor", treat_label),
                               size = 5, fontface = "bold") + theme_void())
  }
  ggplot_data <- mini_heatmap_data %>%
    pivot_longer(cols = c(Control, Treat), names_to = "Condition", values_to = "PSI") %>%
    mutate(
      Condition = ifelse(Condition == "Control", "Control", treat_label),
      Pattern = factor(paste0("Pattern ", Cluster, "\n(n=", n_events, ")"),
                       levels = paste0("Pattern ", mini_heatmap_data$Cluster, "\n(n=", mini_heatmap_data$n_events, ")")),
      x_num = as.numeric(factor(Condition)),
      y_num = as.numeric(Pattern),
      xmin = x_num - 0.45,
      xmax = x_num + 0.45,
      ymin = y_num - 0.45,
      ymax = y_num + 0.45
    )
  
  p <- ggplot(ggplot_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = PSI)) +
    geom_rect(color = NA, radius = unit(3, "mm")) +
    geom_text(aes(x = x_num, y = y_num, label = sprintf("%.2f", PSI)),
              color = ifelse(abs(ggplot_data$PSI) > 0.75, "white", "black"),
              size = 5, fontface = "bold") +
    scale_fill_gradientn(colors = palette_colors,
                         limits = c(-1, 1), oob = scales::squish,
                         name = legend_title,
                         guide = guide_colorbar(frame.colour = "black", frame.linewidth = 0.6,
                                                ticks.colour = "black", 
                                                barheight = unit(2.5, "cm"),
                                                barwidth = unit(0.4, "cm"))) +
    scale_x_continuous(breaks = unique(ggplot_data$x_num), 
                       labels = levels(factor(ggplot_data$Condition)),
                       expand = c(0.05, 0.05)) +
    scale_y_continuous(breaks = unique(ggplot_data$y_num),
                       labels = levels(ggplot_data$Pattern),
                       expand = c(0.05, 0.05)) +
    coord_fixed(ratio = 1) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 8, face = "bold", vjust = 0.5),
      axis.text.y = element_text(size = 8, face = "bold"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 8),
      legend.position = if (show_legend) "right" else "none",
      panel.grid = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(6, 6, 6, 6)
    )
  
  p
}

showtext::showtext_opts(dpi = 600)
showtext_auto()

# ---- Datasets: order = TUB, PladB, SSA ----
datasets <- list(
  tub  = list(events = if (exists("tub_events")) tub_events else NULL,
              diff_exons = if (exists("differential_tub_exons")) differential_tub_exons else NULL,
              diff_introns = if (exists("differential_tub_introns")) differential_tub_introns else NULL,
              label = "Tubercidin"),
  pladb = list(events = if (exists("pladb_events")) pladb_events else NULL,
               diff_exons = if (exists("differential_pladb_exons")) differential_pladb_exons else NULL,
               diff_introns = if (exists("differential_pladb_introns")) differential_pladb_introns else NULL,
               label = "Pladienolide B"),
  ssa  = list(events = if (exists("ssa_events")) ssa_events else NULL,
              diff_exons = if (exists("differential_ssa_exons")) differential_ssa_exons else NULL,
              diff_introns = if (exists("differential_ssa_introns")) differential_ssa_introns else NULL,
              label = "Spliceostatin A")
)

# Prepare and build plots
exon_plots <- list()
intron_plots <- list()
summaries <- list()

for (nm in names(datasets)) {
  ds <- datasets[[nm]]
  if (is.null(ds$events) || is.null(ds$diff_exons) || is.null(ds$diff_introns)) {
    message(sprintf("Skipping '%s' — required objects not found.", nm))
    exon_plots[[nm]] <- ggplot() + theme_void() + ggtitle(ds$label)
    intron_plots[[nm]] <- ggplot() + theme_void()
    next
  }
  s <- prepare_collapsed(ds$events, ds$diff_exons, ds$diff_introns, n_clusters = 2)
  summaries[[nm]] <- s
  
  # Only show legend on the bottom plot of each column
  show_exon_legend <- (nm == "ssa")
  show_intron_legend <- (nm == "ssa")
  
  p_ex <- make_tile_plot(s$exons, ds$label, exon_palette, "Exons Z-Score", show_legend = show_exon_legend)
  p_in <- make_tile_plot(s$introns, ds$label, intron_palette, "Introns Z-Score", show_legend = show_intron_legend)
  
  # Add title only for the top plot of each column
  p_ex <- p_ex + ggtitle(ds$label) + theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 4)))
  
  exon_plots[[nm]] <- p_ex
  intron_plots[[nm]] <- p_in
}

# Assemble final 3-row x 2-column figure (vertical layout)
left_column <- exon_plots$tub / exon_plots$pladb / exon_plots$ssa
right_column <- intron_plots$tub / intron_plots$pladb / intron_plots$ssa

final_fig <- (left_column | right_column) +
  plot_annotation(
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"))
  )

# Print final figure
print(final_fig)

# ---- Optional: save to file ----
ggsave("Figures/Figure2F_splicing_exon_intron_TUB_PLADB_SSA_vertical_PSI.pdf", final_fig, width = 8, height = 10)