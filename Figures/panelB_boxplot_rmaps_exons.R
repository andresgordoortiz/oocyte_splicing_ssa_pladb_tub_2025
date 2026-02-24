library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggrepel)

# Colors
splicing_colors <- c(
  "Included" = "#D4A574",
  "Skipped" = "#B85450"
)

theme_publication <- function(base_size = 10) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      axis.line = element_line(linewidth = 0.6, colour = "black"),
      axis.ticks = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      axis.title = element_text(face = "bold", size = rel(1.0)),
      axis.text = element_text(size = rel(0.9), colour = "black"),
      legend.position = "right",
      legend.background = element_blank(),
      legend.title = element_text(face = "bold", size = rel(0.95)),
      legend.text = element_text(size = rel(0.85)),
      legend.key.size = unit(0.4, "cm"),
      legend.spacing.y = unit(0.1, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = rel(1.0), margin = margin(5,5,5,5)),
      plot.title = element_text(face = "bold", size = rel(1.1), hjust = 0),
      plot.subtitle = element_text(size = rel(0.9), hjust = 0, color = "grey30"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# ============================================================================
# DATA LOADING
# ============================================================================

process_drug_data_exon <- function(drug_name) {
  up_file <- paste0("Supplementary4_rMaps2_out/", drug_name, "_ex_rmaps/pVal.up.vs.bg.RNAmap.txt")
  down_file <- paste0("Supplementary4_rMaps2_out/", drug_name, "_ex_rmaps/pVal.dn.vs.bg.RNAmap.txt")
  
  up_data <- read.delim(up_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  down_data <- read.delim(down_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Rename descriptive columns to R1-R8
  r_cols_actual <- c(
    "smallest_p_in_upstreamExon.3prime",
    "smallest_p_in_upstreamExonIntron",
    "smallest_p_in_upstreamIntron",
    "smallest_p_in_targetExon.5prime",
    "smallest_p_in_targetExon.3prime",
    "smallest_p_in_downstreamIntron",
    "smallest_p_in_downstreamExonIntron",
    "smallest_p_in_downstreamExon.5prime"
  )
  r_cols_new <- paste0("R", 1:8)
  
  rename_cols <- function(df) {
    rename_vec <- setNames(r_cols_actual, r_cols_new)
    rename(df, all_of(rename_vec))
  }
  
  up_data <- rename_cols(up_data)
  down_data <- rename_cols(down_data)
  
  filter_by_any_R <- function(df, source_name) {
    r_cols <- grep("^R\\d+$", names(df), value = TRUE)
    df %>%
      filter(if_any(all_of(r_cols), ~ . <= 0.05)) %>%
      mutate(source = source_name, drug = drug_name)
  }
  
  combined <- bind_rows(
    filter_by_any_R(up_data, "Included"),
    filter_by_any_R(down_data, "Skipped")
  )
  
  combined %>%
    pivot_longer(
      cols = matches("^R\\d+$"),
      names_to = "position",
      values_to = "pvalue"
    ) %>%
    mutate(
      pvalue = as.numeric(pvalue),
      neglog10p = if_else(is.na(pvalue) | pvalue <= 0, NA_real_, -log10(pvalue)),
      RBP_clean = sub("\\..*", "", RBP)
    )
}

cat("Loading exon data for all drugs...\n")
all_drugs <- map_dfr(c("pladb", "ssa", "tub"), process_drug_data_exon)
cat("Data loaded successfully.\n\n")

# ============================================================================
# REGION MAP - shortened labels for strip text
# ============================================================================

region_map <- tibble(
  position = paste0("R", 1:8),
  region_full = c(
    "Upstream\nExon 3'",
    "Upstream\nIntron 5'",
    "Upstream\nIntron 3'",
    "Target\nExon 5'",
    "Target\nExon 3'",
    "Downstream\nIntron 5'",
    "Downstream\nIntron 3'",
    "Downstream\nExon 5'"
  )
)

all_drugs <- all_drugs %>%
  left_join(region_map, by = "position") %>%
  mutate(
    position = factor(position, levels = paste0("R", 1:8)),
    region_full = factor(region_full, levels = region_map$region_full),
    drug_label = factor(drug,
                        levels = c("tub", "pladb", "ssa"),
                        labels = c("Tubercidin", "Pladienolide B", "Spliceostatin A"))
  )

# ============================================================================
# BOXPLOT
# ============================================================================

cat("Creating boxplot...\n")

boxplot_data <- all_drugs %>%
  filter(!is.na(neglog10p)) %>%
  mutate(source = factor(source, levels = c("Included", "Skipped")))

stat_tests <- boxplot_data %>%
  group_by(region_full, drug_label) %>%
  summarize(
    p_value = tryCatch({
      wilcox.test(neglog10p[source == "Included"],
                  neglog10p[source == "Skipped"])$p.value
    }, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    significance = case_when(
      is.na(p_value) ~ "ns",
      p_value < 0.0001 ~ "****",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    label = ifelse(significance != "ns", significance, "")
  )

y_positions <- boxplot_data %>%
  group_by(region_full, drug_label) %>%
  summarize(y_pos = max(neglog10p, na.rm = TRUE) * 1.15, .groups = "drop")

stat_tests <- stat_tests %>%
  left_join(y_positions, by = c("region_full", "drug_label"))

p_boxplots <- ggplot(boxplot_data, aes(x = drug_label, y = neglog10p, fill = source)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, linewidth = 0.5,
               width = 0.7, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = source), alpha = 0.4, size = 1.2,
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  geom_text(
    data = stat_tests %>% filter(label != ""),
    aes(x = drug_label, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 3, fontface = "bold", color = "black"
  ) +
  facet_wrap(~ region_full, ncol = 8, scales = "fixed") +
  labs(
    x = NULL,
    y = expression(bold(-log[10](italic(p)*"-value")))
  ) +
  scale_fill_manual(values = splicing_colors, name = "Splicing Event") +
  scale_color_manual(values = splicing_colors, name = "Splicing Event") +
  theme_publication(base_size = 9) +
  theme(
    legend.position = "top",
    legend.margin = margin(0, 0, 3, 0),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                               size = 7, margin = margin(t = 3)),
    axis.text.y = element_text(size = 7),
    panel.spacing = unit(0.4, "lines"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 7, 
                              margin = margin(b = 4), lineheight = 0.9)
  )

# ============================================================================
# SAVE
# ============================================================================

pdf("panel_B_boxplots_exon.pdf", width = 9.6, height = 5)
print(p_boxplots)
dev.off()

cat("✓ Saved: panel_B_boxplots_exon.pdf\n")