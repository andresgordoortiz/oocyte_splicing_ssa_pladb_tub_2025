library(readxl)
library(tidyverse)
library(ggrepel)

# Import betAS output tables
tub_fdr_df<-read_csv("tub_fdr.csv")[,-1]
pladb_fdr_df<-read_csv("pladb_fdr.csv")[,-1]
ssa_fdr_df<-read_csv("ssa_fdr.csv")[,-1]


protein_impact<- read.delim("PROT_IMPACT-mm10-v2.3.tab") %>%
  rename("EVENT"=EventID)


differential_tub<-na.omit(tub_fdr_df[tub_fdr_df$FDR <= 0.05 & abs(tub_fdr_df$deltapsi) >= 0.1,]) %>%
  left_join(protein_impact, by = "EVENT")

differential_pladb<-na.omit(pladb_fdr_df[pladb_fdr_df$FDR <= 0.05 & abs(pladb_fdr_df$deltapsi) >= 0.1,])%>%
  left_join(protein_impact, by = "EVENT")

differential_ssa<-na.omit(ssa_fdr_df[ssa_fdr_df$FDR <= 0.05 & abs(ssa_fdr_df$deltapsi) >= 0.1,])%>%
  left_join(protein_impact, by = "EVENT")


katarina_protein_stable_oocyte <- read_excel("41556_2024_1442_MOESM4_ESM.xlsx", 
                                             sheet = "TableS1", col_types = c("text", 
                                                                              "text", "text", "text", "text", "text", 
                                                                              "text", "numeric", "text", "text", 
                                                                              "text", "numeric", "numeric", "numeric", 
                                                                              "numeric", "numeric", "numeric", 
                                                                              "numeric", "numeric", "numeric", 
                                                                              "numeric", "numeric", "numeric", 
                                                                              "numeric", "numeric", "numeric", 
                                                                              "numeric", "numeric", "numeric", 
                                                                              "numeric", "numeric")) %>%
  rename("GENE"="gene_name", "stability"="mean percent H") 

katarina_protein_stable_oocyte<-filter(katarina_protein_stable_oocyte, 
                                       !is.na(katarina_protein_stable_oocyte$stability))

disruptedtub<-c(differential_tub$GENE[differential_tub$deltapsi>0 & grepl("inclusion",differential_tub$ONTO)],
                differential_tub$GENE[differential_tub$deltapsi<0 & grepl("exclusion",differential_tub$ONTO)])

disruptedpladb<-c(differential_pladb$GENE[differential_pladb$deltapsi>0 & grepl("inclusion",differential_pladb$ONTO)],
                  differential_pladb$GENE[differential_pladb$deltapsi<0 & grepl("exclusion",differential_pladb$ONTO)])

disruptedssa<-c(differential_ssa$GENE[differential_ssa$deltapsi>0 & grepl("inclusion",differential_ssa$ONTO)],
                differential_ssa$GENE[differential_ssa$deltapsi<0 & grepl("exclusion",differential_ssa$ONTO)])





# Publication-ready stability plot with significance asterisks
library(ggplot2)
library(dplyr)
library(tidyr)
library(rstatix)
library(showtext)
library(sysfonts)

# Font & theme setup
preferred_font <- "Roboto"
font_add_google(preferred_font)
showtext::showtext_opts(dpi = 600)
showtext::showtext_auto()
base_family <- preferred_font
cell_blue <- "#00A1D7"

make_cell_palette <- function(n, main = cell_blue) {
  anchors <- c(main, "#0073A8", "#9EE8FB", "#4D4D4D", "#E69F00")
  if (n <= length(anchors)) anchors[1:n] else colorRampPalette(anchors)(n)
}

theme_cellpub <- function(base_size = 18, base_family_in = NULL) {
  if (is.null(base_family_in)) base_family_in <- get0("base_family", ifnotfound = "sans")
  theme_classic(base_size = base_size, base_family = base_family_in) %+replace%
    theme(
      axis.line = element_line(linewidth = 0.9, colour = "#222222"),
      axis.ticks = element_line(linewidth = 0.9, colour = "#222222"),
      axis.ticks.length = unit(3, "pt"),
      axis.title = element_text(family = base_family_in, face = "plain", size = rel(1.0)),
      axis.text = element_text(family = base_family_in, size = rel(0.95), colour = "#111111"),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.key.size = unit(10, "pt"),
      legend.background = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(family = base_family_in, size = rel(0.95)),
      panel.grid.major.y = element_line(color = alpha("#666666", 0.10), linetype = "dashed", linewidth = 0.4),
      panel.grid.major.x = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = rel(1.0), family = base_family_in),
      plot.title = element_text(face = "bold", size = rel(1.05), hjust = 0, family = base_family_in),
      plot.subtitle = element_text(size = rel(0.95), hjust = 0, family = base_family_in),
      plot.caption = element_text(size = rel(0.85), colour = "#666666", family = base_family_in),
      plot.margin = margin(6, 6, 6, 6)
    )
}

# ---------- Data preparation -----------------
df_flags <- katarina_protein_stable_oocyte %>%
  mutate(
    stability = as.numeric(stability),
    disrupted_tub = GENE %in% disruptedtub,
    disrupted_ssa = GENE %in% disruptedssa,
    disrupted_pladb = GENE %in% disruptedpladb
  )

# Create dataset
long_disrupted <- df_flags %>%
  pivot_longer(cols = starts_with("disrupted_"), names_to = "group", values_to = "flag") %>%
  filter(flag == TRUE) %>%
  select(GENE, group, stability)

background <- df_flags %>%
  filter(!disrupted_tub & !disrupted_ssa & !disrupted_pladb) %>%
  select(GENE, stability) %>%
  mutate(group = "background")

df <- bind_rows(long_disrupted, background) %>%
  filter(!is.na(stability) & stability > 0) %>%
  mutate(log_stability = log10(stability))

group_levels <- c("background", "disrupted_tub", "disrupted_pladb", "disrupted_ssa")
df$group <- factor(df$group, levels = group_levels)

# ---------- Statistical tests and significance asterisks -----------------
kw_result <- df %>% kruskal_test(stability ~ group)

# Calculate p-values for each disrupted group vs background
sig_data <- data.frame()
for(grp in c("disrupted_tub", "disrupted_ssa", "disrupted_pladb")) {
  bg_vals <- df %>% filter(group == "background") %>% pull(stability)
  grp_vals <- df %>% filter(group == grp) %>% pull(stability)
  
  if(length(bg_vals) > 0 & length(grp_vals) > 0) {
    p_val <- wilcox.test(bg_vals, grp_vals)$p.value
    # Convert p-value to asterisks
    asterisk <- case_when(
      p_val < 0.001 ~ "***",
      p_val < 0.01 ~ "**", 
      p_val < 0.05 ~ "*",
      TRUE ~ "ns"
    )
    sig_data <- rbind(sig_data, data.frame(group = grp, p_value = p_val, asterisk = asterisk))
  }
}

# Get max values for asterisk positioning
max_vals <- df %>% 
  group_by(group) %>% 
  summarise(max_log = max(log_stability, na.rm = TRUE))

sig_data <- sig_data %>%
  left_join(max_vals, by = "group") %>%
  mutate(y_pos = max_log + 0.1)  # Position asterisks above highest point

# ---------- Create palette and clean labels -----------------
palette <- rev(make_cell_palette(length(group_levels)))
names(palette) <- group_levels

# Clean group labels for display
clean_labels <- c(
  "background" = "Background",
  "disrupted_tub" = "Tubercidin", 
  "disrupted_ssa" = "Spliceostatin A",
  "disrupted_pladb" = "Pladienolide B"
)

# ---------- Create publication-ready plot -----------------
p <- ggplot(df, aes(x = group, y = log_stability, fill = group)) +
  
  # Violin + box plots
  geom_violin(trim = TRUE, width = 0.7, alpha = 0.6, show.legend = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7, show.legend = FALSE) +
  
  # Individual points and labels for disrupted genes only
  geom_jitter(data = df %>% filter(group != "background"),
              aes(color = group, ), 
              width = 0.15, height = 0, size = 2, alpha = 0.8, show.legend = FALSE) +
  geom_text(data = df %>% filter(group != "background"),
            aes(label = GENE),
            position = position_jitter(width = 0.3, height = 0),
            size = 5, family = base_family, show.legend = FALSE, 
            check_overlap = TRUE) +
  
  # Color scheme
  scale_fill_manual(values = palette) +
  scale_color_manual(values = palette) +
  
  # Clean axis labels
  scale_x_discrete(labels = clean_labels) +
  scale_y_continuous(
    name = "Protein Stability (log10)",
    labels = function(x) {
      vals <- 10^x
      ifelse(vals < 0.01 | vals > 1000, 
             sprintf("%.0e", vals),
             sprintf("%.2g", vals))
    },
    expand = expansion(mult = c(0.02, 0.15))  # Extra space at top for asterisks
  ) +
  
  # Add significance asterisks
  geom_text(data = sig_data, 
            aes(x = group, y = y_pos, label = asterisk),
            inherit.aes = FALSE,
            size = 6, 
            family = base_family,
            fontface = "bold") +
  
  # Theme and labels
  theme_cellpub(base_size = 16, base_family_in = base_family) +
  theme(
    axis.text.x = element_text(face = "bold", size = rel(0.9)),
    legend.position = "none"
  ) +
  
  labs(
    x = "",
    title = "Protein Stability Across Disrupted Gene Sets",
    subtitle = paste0("Kruskal-Wallis p = ", 
                      format.pval(kw_result$p, digits = 2))
  )

# Display results
print(p)

# High-resolution save
ggsave("stability_plot_publication.png", p, 
       width = 8, height = 6, dpi = 600, bg = "white")
ggsave("stability_plot_publication.pdf", p, 
       width = 8, height = 6, device = cairo_pdf)

# Print significance results
cat("\nSignificance Results:\n")
print(sig_data)
