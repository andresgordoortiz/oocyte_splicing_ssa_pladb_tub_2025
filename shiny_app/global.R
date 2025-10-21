# ============================================================================
# Global Configuration for Drug Splicing Effects Explorer
# ============================================================================

# Load required libraries
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(DT)
library(dplyr)
library(ggplot2)
library(ggpubr)  # For stat_compare_means
library(plotly)
library(readr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(scales)
library(viridis)
library(gridExtra)
library(patchwork)

# ============================================================================
# THEME AND COLOR DEFINITIONS
# ============================================================================

# Drug colors (consistent with publication)
DRUG_COLORS <- c(
  "Tubercidin" = "#A0C1B9",
  "Pladienolide B" = "#70A0AF",
  "Spliceostatin A" = "#706993"
)

# Splicing direction colors
SPLICING_COLORS <- c(
  "Included" = "#D4A574",
  "Skipped" = "#B85450",
  "Retained" = "#D4A574",
  "Removed" = "#B85450",
  "Unchanged" = "#808080"
)

# Event type colors
EVENT_COLORS <- c(
  "Exon" = "#7A2048",
  "Intron" = "#2C2C2C",
  "Alt5" = "#4B4856",
  "Alt3" = "#706993"
)

# Professional publication theme
theme_publication <- function(base_size = 12) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      # Axes
      axis.line = element_line(linewidth = 0.6, colour = "black"),
      axis.ticks = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      axis.title = element_text(face = "bold", size = rel(1.0)),
      axis.text = element_text(size = rel(0.9), colour = "black"),
      
      # Legend
      legend.position = "right",
      legend.background = element_blank(),
      legend.title = element_text(face = "bold", size = rel(0.95)),
      legend.text = element_text(size = rel(0.85)),
      legend.key.size = unit(0.4, "cm"),
      
      # Panel
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Facets
      strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = rel(1.0)),
      
      # Plot
      plot.title = element_text(face = "bold", size = rel(1.1), hjust = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Classify event types
classify_event <- function(event_id) {
  ifelse(grepl("EX", event_id), "Exon",
         ifelse(grepl("INT", event_id), "Intron",
                ifelse(grepl("ALTD", event_id), "Alt5",
                       ifelse(grepl("ALTA", event_id), "Alt3", "Other")
                )
         )
  )
}

# Format p-values
format_pval <- function(p) {
  ifelse(p < 0.001, "< 0.001",
         ifelse(p < 0.01, sprintf("%.3f", p),
                sprintf("%.2f", p)))
}

# Calculate summary statistics
calc_summary_stats <- function(data, drug_name, fdr_thresh = 0.05, deltapsi_thresh = 0.1) {
  sig_events <- data %>%
    filter(FDR <= fdr_thresh, abs(deltapsi) >= deltapsi_thresh)
  
  list(
    drug = drug_name,
    total_events = nrow(data),
    significant_events = nrow(sig_events),
    percent_significant = round(100 * nrow(sig_events) / nrow(data), 2),
    upregulated = sum(sig_events$deltapsi > 0),
    downregulated = sum(sig_events$deltapsi < 0),
    median_deltapsi = median(sig_events$deltapsi, na.rm = TRUE),
    max_deltapsi = max(abs(sig_events$deltapsi), na.rm = TRUE)
  )
}

# ============================================================================
# DATA LOADING
# ============================================================================

# Function to load data safely
load_data_safely <- function(file_path, default = NULL) {
  tryCatch({
    if (file.exists(file_path)) {
      if (grepl("\\.csv$", file_path)) {
        data <- read_csv(file_path, show_col_types = FALSE)
        if (ncol(data) > 0 && names(data)[1] == "...1") {
          data <- data[, -1]
        }
        return(data)
      } else if (grepl("\\.tab$", file_path)) {
        return(read_delim(file_path, delim = "\t", show_col_types = FALSE))
      }
    }
    return(default)
  }, error = function(e) {
    warning(paste("Error loading", file_path, ":", e$message))
    return(default)
  })
}

# Load main datasets
load_main_datasets <- function() {
  list(
    tub = load_data_safely("../tub_fdr.csv"),
    pladb = load_data_safely("../pladb_fdr.csv"),
    ssa = load_data_safely("../ssa_fdr.csv")
  )
}

# Load GO enrichment data
load_go_data <- function() {
  drugs <- c("tub", "pladb", "ssa")
  ontologies <- c("BP", "MF", "CC")
  
  go_list <- list()
  for (drug in drugs) {
    go_list[[drug]] <- list()
    for (ont in ontologies) {
      file_path <- paste0("../GO_enrichment_results/", drug, "_GO_", ont, ".csv")
      go_list[[drug]][[ont]] <- load_data_safely(file_path)
    }
  }
  return(go_list)
}

# Load event info
load_event_info <- function() {
  load_data_safely("../EVENT_INFO-mm10.tab")
}

# ============================================================================
# INITIAL DATA LOADING
# ============================================================================

message("Loading application data...")
datasets <- load_main_datasets()
go_data <- load_go_data()
event_info <- load_event_info()

message("Data loading complete!")
