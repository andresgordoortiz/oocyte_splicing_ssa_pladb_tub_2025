# Import betAS output tables

library(readr)
tub_fdr_df   <- read_csv("tub_fdr.csv")[,-1]
pladb_fdr_df <- read_csv("pladb_fdr.csv")[,-1]
ssa_fdr_df   <- read_csv("ssa_fdr.csv")[,-1]

# --- Filter significant events (FDR <= 0.05 & |deltapsi| >= 0.1) ---
differential_tub  <- na.omit(tub_fdr_df[tub_fdr_df$FDR <= 0.05 & abs(tub_fdr_df$deltapsi) >= 0.1, ])
differential_pladb <- na.omit(pladb_fdr_df[pladb_fdr_df$FDR <= 0.05 & abs(pladb_fdr_df$deltapsi) >= 0.1, ])
differential_ssa  <- na.omit(ssa_fdr_df[ssa_fdr_df$FDR <= 0.05 & abs(ssa_fdr_df$deltapsi) >= 0.1, ])

# --- Libraries for plotting / themeing ---
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(showtext)
library(scales)

# --- Event classification function (same logic as your original) ---
classify_event <- function(x) {
  ifelse(grepl("EX", x), "Exon",
         ifelse(grepl("INT", x), "Intron",
                ifelse(grepl("ALTD", x), "Alt5",
                       ifelse(grepl("ALTA", x), "Alt3", NA)
                )
         )
  )
}

# canonical order for event types
levels_order <- c("Exon", "Intron", "Alt5", "Alt3")

# --- Build long dataframe: one row per significant event occurrence ---
long_df <- bind_rows(
  tibble(Condition = "Tubercidin",       Event = classify_event(differential_tub$EVENT),  Inclusion = differential_tub$deltapsi > 0),
  tibble(Condition = "Pladienolide B",   Event = classify_event(differential_pladb$EVENT), Inclusion = differential_pladb$deltapsi > 0),
  tibble(Condition = "Spliceostatin A",  Event = classify_event(differential_ssa$EVENT),  Inclusion = differential_ssa$deltapsi > 0)
) %>%
  filter(!is.na(Event)) %>%
  mutate(
    Event = factor(Event, levels = levels_order),
    Condition = factor(Condition, levels = c("Tubercidin", "Pladienolide B", "Spliceostatin A")),
    Inclusion = ifelse(Inclusion, "Inclusion", "Skipping") %>% factor(levels = c("Skipping", "Inclusion"))
  )

# Optional: quick counts table (wide)
fig_summary <- long_df %>%
  count(Condition, Event) %>%
  pivot_wider(names_from = Condition, values_from = n, values_fill = 0) %>%
  arrange(match(Event, levels_order))

print(fig_summary)  # glance at absolute counts per event type

# create a full grid of Condition × Event × Inclusion
full_grid <- tidyr::expand_grid(
  Condition = unique(long_df$Condition),
  Event = factor(levels_order, levels = levels_order),
  Inclusion = c("Skipping", "Inclusion")
)

summary_df <- long_df %>%
  group_by(Condition, Event, Inclusion) %>%
  summarise(n = n(), .groups = "drop") %>%
  ungroup() %>%
  # left_join to ensure all combinations present
  right_join(full_grid, by = c("Condition", "Event", "Inclusion")) %>%
  mutate(n = replace_na(n, 0)) %>%
  group_by(Condition, Event) %>%
  mutate(
    total = sum(n),
    prop = ifelse(total > 0, n / total, 0),
    pct_label = ifelse(prop > 0, scales::percent(prop, accuracy = 0.1), "")
  ) %>%
  ungroup()


# --- Styling: font + palette + theme (publication-ready) ---
preferred_font <- "Roboto"   # Google Inter for crisp text
font_add_google(preferred_font)
showtext::showtext_opts(dpi = 600)
showtext_auto()

base_family <- preferred_font

# cell-inspired palette function (keeps your original approach)
cell_blue <- "#00A1D7"
make_cell_palette <- function(n, main = cell_blue) {
  anchors <- c(main, "#0073A8", "#9EE8FB", "#4D4D4D", "#E69F00")
  if (n <= length(anchors)) anchors[1:n] else colorRampPalette(anchors)(n)
}

# For Inclusion vs Skipping use a high-contrast pair (Inclusion = blue, Skipping = light grey)
inc_pal <- c("Skipping" = "#DDDDDD", "Inclusion" = cell_blue)

# Publication-like theme (clean, crisp)
theme_cellpub <- function(base_size = 14, base_family = base_family) {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line = element_line(size = 0.8, colour = "#222222"),
      axis.ticks = element_line(size = 0.8, colour = "#222222"),
      axis.ticks.length = unit(3, "pt"),
      axis.title = element_text(face = "plain", size = rel(0.95)),
      axis.text = element_text(size = rel(0.9), colour = "#111111"),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.key.size = unit(10, "pt"),
      legend.background = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = rel(0.9)),
      panel.grid.major.y = element_line(color = alpha("#666666", 0.08), linetype = "dashed", size = 0.4),
      panel.grid.major.x = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = rel(1.0)),
      plot.title = element_text(face = "bold", size = rel(1.05), hjust = 0),
      plot.subtitle = element_text(size = rel(0.95), hjust = 0),
      plot.caption = element_text(size = rel(0.85), colour = "#666666"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

# --- Reorder Event levels optionally by global frequency (most common first) ---
event_order <- long_df %>% count(Event) %>% arrange(desc(n)) %>% pull(Event) %>% as.character()
# Keep canonical order if you prefer; otherwise use event_order
# To keep canonical use: final_event_levels <- levels_order
final_event_levels <- levels_order
summary_df <- summary_df %>% mutate(Event = factor(Event, levels = final_event_levels))

# --- Plot: 100% stacked bars per Event, faceted by Condition ---
p <- ggplot(summary_df, aes(x = Event, y = prop, fill = Inclusion)) +
  geom_col(width = 0.7, colour = "black", linewidth  = 0.25) +
  # percent labels: show only when segment large enough to be legible
  geom_text(
    aes(label = ifelse(prop >= 0.04, pct_label, "")),
    position = position_stack(vjust = 0.5),
    family = base_family,
    size = 3.2
  ) +
  facet_wrap(~ Condition, nrow = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.03))) +
  scale_fill_manual(values = inc_pal) +
  labs(
    x = "Event type",
    y = "Proportion of significant events",
    title = "% Inclusion vs Skipping among significant splicing events",
    subtitle = "For each Condition (facet) and Event type: proportion of events classified as Inclusion (deltapsi > 0) vs Skipping (deltapsi <= 0)",
    caption = "Significant events: FDR <= 0.05 and |deltapsi| >= 0.1"
  ) +
  theme_cellpub(base_size = 14, base_family = base_family) +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5),
    legend.position = "top",
    legend.title = element_blank(),
    strip.text = element_text(size = rel(1.0))
  )

# Print plot into your current device
print(p)

# --- Save high-quality outputs (PDF + PNG) ---
ggsave("figure_pct_inclusion_by_event_and_condition.pdf", plot = p, width = 10, height = 4.0, device = cairo_pdf)
ggsave("figure_pct_inclusion_by_event_and_condition.png", plot = p, width = 10, height = 4.0, dpi = 600)

