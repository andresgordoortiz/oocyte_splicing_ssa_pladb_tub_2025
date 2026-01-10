
# ----- PREPARE DF -------
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(showtext)
library(scales)
library(readr)
library(readxl)


# ---- 1. Import from the combined Excel file ----
excel_file <- "Preprocessing/betAS_out/Supplementary2_betAS_splicing_results.xlsx"

# Read each sheet (using the names we set in the previous step)
tub_fdr_df   <- read_excel(excel_file, sheet = "Tub_FDR")
pladb_fdr_df <- read_excel(excel_file, sheet = "PlaDB_FDR")
ssa_fdr_df   <- read_excel(excel_file, sheet = "SSA_FDR")

# Extract significant (FDR<=0.05 & |deltaPSI| > 10%)
differential_tub <- na.omit(tub_fdr_df[tub_fdr_df$FDR <= 0.05 & abs(tub_fdr_df$deltapsi) >= 0.1,])
differential_pladb <- na.omit(pladb_fdr_df[pladb_fdr_df$FDR <= 0.05 & abs(pladb_fdr_df$deltapsi) >= 0.1,])
differential_ssa <- na.omit(ssa_fdr_df[ssa_fdr_df$FDR <= 0.05 & abs(ssa_fdr_df$deltapsi) >= 0.1,])
showtext::showtext_opts(dpi=600)

# Classification function
classify_event <- function(x) {
  ifelse(grepl("EX", x), "Exon",
         ifelse(grepl("INT", x), "Intron",
                ifelse(grepl("ALTD", x), "Alt5",
                       ifelse(grepl("ALTA", x), "Alt3", NA)
                )
         )
  )
}

# Canonical order
levels_order <- c("Exon", "Intron", "Alt5", "Alt3")

# Build long dataframe
long_df <- bind_rows(
  tibble(Condition = "Tubercidin", Event = classify_event(differential_tub$EVENT)),
  tibble(Condition = "Pladienolide B", Event = classify_event(differential_pladb$EVENT)),
  tibble(Condition = "Spliceostatin A", Event = classify_event(differential_ssa$EVENT))
) %>%
  filter(!is.na(Event)) %>%
  mutate(
    Event = factor(Event, levels = levels_order),
    Condition = factor(Condition, levels = c("Tubercidin", "Pladienolide B", "Spliceostatin A"))
  )

# Summary counts
fig_summary <- long_df %>%
  count(Condition, Event) %>%
  pivot_wider(names_from = Condition, values_from = n, values_fill = 0) %>%
  arrange(match(Event, levels_order))

print(fig_summary)

# ---------- PLOT SETUP ----------

# Font setup
showtext::showtext_opts(dpi = 600)
showtext_auto()

# ---- Drug colors ----
ssa_color <- "#706993"    
pladb_color <- "#70A0AF"  
tub_color  <- "#A0C1B9"   


color_palette <- c(
  "Tubercidin" = tub_color,
  "Pladienolide B" = pladb_color,
  "Spliceostatin A" = ssa_color
)

# Professional publication theme
theme_publication <- function(base_size = 8, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Axes
      axis.line = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks.length = unit(2, "pt"),
      axis.title = element_text(face = "plain", size = base_size),
      axis.text = element_text(size = base_size - 1, colour = "black"),
      axis.text.x = element_text(margin = margin(t = 3)),
      axis.text.y = element_text(margin = margin(r = 3)),
      axis.title.x = element_text(margin = margin(t = 6)),
      axis.title.y = element_text(angle = 90, margin = margin(r = 6)),
      
      # Legend
      legend.position = "top",
      legend.justification = "left",
      legend.direction = "horizontal",
      legend.key.size = unit(8, "pt"),
      legend.key = element_rect(fill = NA, colour = NA),
      legend.background = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 1),
      legend.spacing.x = unit(4, "pt"),
      legend.margin = margin(0, 0, 4, 0),
      
      # Panel
      panel.background = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Plot
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0, margin = margin(b = 6)),
      plot.margin = margin(8, 8, 8, 8)
    )
}

# ---------- CREATE PLOT ----------

p <- ggplot(long_df, aes(x = Event, fill = Condition)) +
  geom_bar(
    position = position_dodge(width = 0.85),
    width = 0.75,
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    position = position_dodge(width = 0.85),
    vjust = -0.4,
    family = "Arial",
    size = 2.2
  ) +
  scale_fill_manual(values = color_palette) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)),
    breaks = pretty_breaks(n = 6)
  ) +
  labs(
    x = "Splicing event type",
    y = "Number of significant events"
  ) +
  coord_cartesian(clip = "off") +
  theme_publication(base_size = 8, base_family = "sans")

print(p)

# ---------- SAVE FIGURE ----------
# For single-column figure (85mm width is standard for many journals)
# For double-column, use width = 180mm (7 inches)

ggsave(
  "Figures/Figure2D_splicing_events_figure.pdf",
  plot = p,
  width = 85,
  height = 65,
  units = "mm",
  device = cairo_pdf
)

ggsave(
  "Figures/Figure2D_splicing_events_figure.png",
  plot = p,
  width = 85,
  height = 65,
  units = "mm",
  dpi = 600
)

# For double-column version
ggsave(
  "Figures/Figure2D_splicing_events_figure_wide.pdf",
  plot = p,
  width = 180,
  height = 100,
  units = "mm",
  device = cairo_pdf
)

