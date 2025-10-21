library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)
library(showtext)

# Import betAS output tables
tub_fdr_df <- read_csv("tub_fdr.csv")[,-1]
pladb_fdr_df <- read_csv("pladb_fdr.csv")[,-1]
ssa_fdr_df <- read_csv("ssa_fdr.csv")[,-1]

# Base colors
ssa_color   <- "#706993"
pladb_color <- "#70A0AF"
tub_color   <- "#A0C1B9"
base_colors <- tibble(dataset = c("tub", "pladb", "ssa"),
                      base_hex = c(tub_color, pladb_color, ssa_color))

# Filter differentials
differential_tub   <- na.omit(tub_fdr_df[tub_fdr_df$FDR <= 0.05 & abs(tub_fdr_df$deltapsi) >= 0.1,])
differential_pladb <- na.omit(pladb_fdr_df[pladb_fdr_df$FDR <= 0.05 & abs(pladb_fdr_df$deltapsi) >= 0.1,])
differential_ssa   <- na.omit(ssa_fdr_df[ssa_fdr_df$FDR <= 0.05 & abs(ssa_fdr_df$deltapsi) >= 0.1,])

tub_fdr_df$dataset   <- "tub"
pladb_fdr_df$dataset <- "pladb"
ssa_fdr_df$dataset   <- "ssa"

tub_fdr_df$drug   <- "Tubercidin"
pladb_fdr_df$drug <- "Pladienolide B"
ssa_fdr_df$drug   <- "Spliceostatin A"

# combine
combined_df <- bind_rows(tub_fdr_df, pladb_fdr_df, ssa_fdr_df)

# ---- thresholds for significance logic (kept as you used in differential_ dfs) ----
deltapsi_thresh <- 0.1
fdr_thresh      <- 0.05

# ---- Create FDR_for_plot handling zeros ----
# find a small non-zero FDR to replace zero values for plotting
min_nonzero_fdr <- suppressWarnings(min(combined_df$FDR[!is.na(combined_df$FDR) & combined_df$FDR > 0], na.rm = TRUE))
if (!is.finite(min_nonzero_fdr) || is.na(min_nonzero_fdr)) {
  min_nonzero_fdr <- 1e-300
}

combined_df <- combined_df %>%
  mutate(FDR_for_plot = case_when(
    is.na(FDR) ~ NA_real_,
    FDR == 0   ~ min_nonzero_fdr / 10,
    TRUE       ~ as.numeric(FDR)
  ))

# ---- Define significance by presence of EVENT in the differential_* dataframes ----
# we assume each *_fdr_df and differential_* share an 'EVENT' column
combined_df <- combined_df %>%
  mutate(
    significant = case_when(
      dataset == "tub"  & EVENT %in% differential_tub$EVENT  ~ "sig",
      dataset == "pladb" & EVENT %in% differential_pladb$EVENT ~ "sig",
      dataset == "ssa"  & EVENT %in% differential_ssa$EVENT  ~ "sig",
      TRUE ~ "ns"
    ),
    negLog10FDR = -log10(FDR_for_plot)
  )

# ---- Remove rows that cannot be plotted ----
initial_n <- nrow(combined_df)
plot_df <- combined_df %>%
  filter(!is.na(deltapsi) & !is.na(negLog10FDR) & is.finite(deltapsi) & is.finite(negLog10FDR))
removed_n <- initial_n - nrow(plot_df)
message("Removed ", removed_n, " rows before plotting (NA or non-finite deltapsi/FDR).")

# ---- Jitter the top points that came from FDR == 0 so they don't overlap exactly ----
# Only jitter rows where original FDR was exactly zero in the raw data (not the FDR_for_plot).
set.seed(42)
plot_df <- plot_df %>%
  mutate(
    is_fdr_zero = !is.na(FDR) & FDR == 0,
    # small jitter on y only for FDR==0; adjust magnitude if needed
    negLog10FDR_jitter = negLog10FDR + ifelse(is_fdr_zero, runif(n(), -0.08, 0.08), 0)
  )

# ---- Transform y-axis to compress space between 3 and max ----
# Map values: 0-3 stay as is, 3-max get compressed into 3-3.5 range
y_raw_max <- max(plot_df$negLog10FDR, na.rm = TRUE)
break_point <- 3
compressed_max <- 3.5  # Upper limit after compression

transform_y <- function(y) {
  ifelse(y <= break_point, 
         y, 
         break_point + (y - break_point) / (y_raw_max - break_point) * (compressed_max - break_point))
}

inverse_transform_y <- function(y) {
  ifelse(y <= break_point,
         y,
         break_point + (y - break_point) / (compressed_max - break_point) * (y_raw_max - break_point))
}

plot_df <- plot_df %>%
  mutate(negLog10FDR_transformed = transform_y(negLog10FDR_jitter))

label_df_transformed <- plot_df %>%
  filter(significant == "sig") %>%
  slice_max(order_by = abs(deltapsi), prop = 0.005, with_ties = FALSE) %>%
  ungroup()

# Use jittered value for plotting; but keep computation of axis limits on the non-jittered values
# Determine axis limits (cap at least 3 on y)
x_lim <- range(plot_df$deltapsi, na.rm = TRUE)
x_pad <- diff(x_lim) * 0.03
if (!is.finite(x_pad) || x_pad == 0) x_pad <- 0.02
x_lim <- c(x_lim[1] - x_pad, x_lim[2] + x_pad)

y_lim <- c(0, compressed_max)

# ---- color mapping, legend labels with counts ----
color_map <- c("ns" = "grey80",
               "Spliceostatin A" = ssa_color,
               "Pladienolide B" = pladb_color,
               "Tubercidin" = tub_color)

sig_counts <- plot_df %>%
  filter(significant == "sig") %>%
  group_by(drug) %>%
  summarise(n = n(), .groups = "drop")

all_drugs <- c("Spliceostatin A", "Pladienolide B", "Tubercidin")
# ensure all present
sig_counts <- sig_counts %>% right_join(tibble(drug = all_drugs), by = "drug") %>%
  mutate(n = ifelse(is.na(n), 0L, n))

legend_labels <- c("Not significant",
                   sapply(all_drugs, function(d) {
                     count <- sig_counts$n[sig_counts$drug == d]
                     paste0(d, " (", count, " sig)")
                   }))

# Font setup
font_add(family = "Arial", regular = "/usr/share/fonts/liberation/LiberationSans-Regular.ttf")
showtext::showtext_opts(dpi = 600)
showtext_auto()

# ---- Build the volcano plot ----
volcano_plot_all <- ggplot(plot_df, aes(x = deltapsi, y = negLog10FDR_transformed)) +
  # points: color by point_color (drug for sig, "ns" for non-sig)
  geom_point(aes(color = ifelse(significant == "sig", drug, "ns")), size = 2.5, alpha = 0.3, stroke = 0, show.legend = TRUE) +
  # Add break indicator lines
  annotate("segment", x = x_lim[1], xend = x_lim[1] + 0.01, 
           y = break_point, yend = break_point + 0.02, 
           color = "black", linewidth = 0.5) +
  annotate("segment", x = x_lim[1], xend = x_lim[1] + 0.01, 
           y = break_point, yend = break_point - 0.02, 
           color = "black", linewidth = 0.5) +
  scale_color_manual(values = color_map, labels = legend_labels, name = NULL) +
  coord_cartesian(xlim = x_lim, ylim = y_lim, expand = FALSE, clip = "off") +
  scale_y_continuous(
    breaks = c(0, 1, 2, 3, compressed_max),
    labels = c("0", "1", "2", "3", paste0("≥", round(y_raw_max, 1)))
  ) +
  theme_classic(base_size = 14, base_family = "Liberation Sans") +
  theme(
    axis.line = element_line(linewidth = 0.9, colour = "#222222"),
    axis.ticks = element_line(linewidth = 0.9, colour = "#222222"),
    axis.text = element_text(size = rel(1.0), colour = "#111111"),
    axis.title = element_text(face = "plain", size = rel(1.0)),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 3),
    legend.key = element_rect(fill = NA, colour = NA)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)),
         size = "none", alpha = "none") +
  labs(
    x = "Spliced-in rate",
    y = "-log10(adjusted p-value)")


ggsave("volcano_combined_deltapsi.png", volcano_plot_all, device = png, width = 3, height = 3, units = "in", dpi = 600)