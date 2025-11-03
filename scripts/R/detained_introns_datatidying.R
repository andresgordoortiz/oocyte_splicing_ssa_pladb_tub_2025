# https://genesdev.cshlp.org/content/suppl/2014/12/29/29.1.63.DC1/Supplemental_DataS1.txt
# https://genesdev.cshlp.org/content/29/1/63
library(dplyr)
library(tidyr)
library(stringr)
library(fuzzyjoin)

# 1. Read the BED file with detained intron info.
bed <- read.table("detained_introns_mm9.bed", header = FALSE, stringsAsFactors = FALSE)
colnames(bed) <- c("chr", "start", "end", "name", "regulation","strand")

# 2. Read the splicing events file.
events <- read.table("EVENT_INFO-mm9.tab", header = TRUE, sep = "\t")

# 3. Filter events to only include rows where EVENT contains "INT".
events_int <- events %>% 
  filter(grepl("INT", EVENT))

# 4. Extract chromosome, start, and end from the COORD_o column.
# Expected format: "chr9:22253299-22254245"
events_parsed <- events_int %>%
  extract(COORD_o, into = c("event_chr", "event_start", "event_end"),
          regex = "^(chr[^:]+):([0-9]+)-([0-9]+)$", remove = FALSE) %>%
  mutate(
    event_start = as.numeric(event_start),
    event_end   = as.numeric(event_end)
  )

# 5. Fuzzy join the BED data with the parsed events:
# - Chromosome must match exactly.
# - Start and end coordinates are considered a match if they differ by <= 10 bases.
# Note: this takes a lot of time
final_df <- fuzzy_left_join(
  bed, events_parsed,
  by = c("chr" = "event_chr", "start" = "event_start", "end" = "event_end"),
  match_fun = list(`==`, function(x, y) abs(x - y) <= 10, function(x, y) abs(x - y) <= 10)
)

# ---- Exclusion/Inclusion -----

# FMN2
fmndko_data <- getDataset(pathTables = paste0(getwd(),"/notebooks/inclusion_tables/fmndko_INCLUSION_LEVELS_FULL-mm10.tab"), tool = "vast-tools")
fmndko_events <- filterEvents(getEvents(fmndko_data, tool = "vast-tools"), N = 10) # Extract alternative splicing
fmndko_introns <- filterEvents(fmndko_events, types = c("IR"), N = 10)



# Remove any NA rows (adjust as needed)
fmndko_introns_df <- filter(na.omit(fmndko_introns$PSI), EVENT %in% final_df$EVENT)

# Function to compute the PSI difference for each row (KO - Control)
compute_difference <- function(df) {
  avg_KO <- rowMeans(df[, 9:10])
  avg_control <- rowMeans(df[, 7:8])
  avg_KO - avg_control
}

# Calculate the difference and remove any additional NAs
fmndko_introns_df$difference <- compute_difference(fmndko_introns_df)
fmndko_introns_df <- na.omit(fmndko_introns_df)

# Create a group variable based on whether the difference is below or above zero
fmndko_introns_df$group <- ifelse(fmndko_introns_df$difference < 0, "Skipped", "Included")

# Perform a one-sample Wilcoxon signed-rank test (null hypothesis: median = 0)
symmetry_test <- wilcox.test(fmndko_introns_df$difference, mu = 0)

# Compute the percentages (areas under the density curve)
pct_below <- round(mean(fmndko_introns_df$difference < 0) * 100, 1)
pct_above <- round(mean(fmndko_introns_df$difference > 0) * 100, 1)

# ---- Plot -------------------------------------------------------------------

plot_psi_distribution <- ggplot(
  fmndko_introns_df[abs(fmndko_introns_df$difference) > 0, ], 
  aes(x = difference, fill = group)
) +
  geom_histogram(aes(y = ..density..), bins = 200, 
                 position = "identity", alpha = 0.9) +
  labs(
    title = "PSI Difference Distribution of Introns in FMN2",
    subtitle = paste("One-sample Wilcoxon test p-value for Symmetry:", 
                     format(symmetry_test$p.value, digits = 3)),
    x = "ΔPSI (KO - Control)",
    y = "Density",
    fill = "Splicing Direction",
    caption = paste0("Created by AG on ", Sys.Date())
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    legend.position = "right",
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.text.x = element_text(size = 15),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "gray80", size = 0.5),
    panel.grid.minor = element_line(color = "gray95", size = 0.3),
    strip.text = element_text(size = 14, face = "bold"),
    panel.spacing = unit(1.5, "lines")
  ) +
  scale_x_continuous(
    breaks = c(-100, -75, -50, -25, 0, 25, 50, 75, 100)
  ) +
  # Annotate the percentages with polished label boxes
  annotate("label", x = -65, y = 0.1, hjust = 0, vjust = 0,
           label = sprintf("dPSI < 0: %.1f%%", pct_below),
           fill = alpha("white", 0.8), color = "blue", size = 7, fontface = "bold",
           label.size = 0.5) +
  annotate("label", x = 65, y = 0.1, hjust = 1, vjust = 0,
           label = sprintf("dPSI > 0: %.1f%%", pct_above),
           fill = alpha("white", 0.8), color = "red", size = 7, fontface = "bold",
           label.size = 0.5) +
  scale_fill_brewer(palette = "Set1")


print(plot_psi_distribution)




# Spire

spire_data <- getDataset(pathTables = paste0(getwd(),"/notebooks/inclusion_tables/spiredko_INCLUSION_LEVELS_FULL-mm10.tab"), tool = "vast-tools")
spire_events <- filterEvents(getEvents(spire_data, tool = "vast-tools"), N=10) # Extract alternative splicing events
spire_introns <- filterEvents(spire_events, types = c("IR"), N = 10)

# Remove any NA rows (adjust as needed)
spire_introns_df <- filter(na.omit(spire_introns$PSI), EVENT %in% final_df$EVENT)

# Function to compute the PSI difference for each row (KO - Control)
compute_difference <- function(df) {
  avg_KO <- rowMeans(df[, 10:12])
  avg_control <- rowMeans(df[, 7:9])
  avg_KO - avg_control
}

# Calculate the difference and remove any additional NAs
spire_introns_df$difference <- compute_difference(spire_introns_df)
spire_introns_df <- na.omit(spire_introns_df)

# Create a group variable based on whether the difference is below or above zero
spire_introns_df$group <- ifelse(spire_introns_df$difference < 0, "Skipped", "Included")

# Perform a one-sample Wilcoxon signed-rank test (null hypothesis: median = 0)
symmetry_test <- wilcox.test(spire_introns_df$difference, mu = 0)

# Compute the percentages (areas under the density curve)
pct_below <- round(mean(spire_introns_df$difference < 0) * 100, 1)
pct_above <- round(mean(spire_introns_df$difference > 0) * 100, 1)

# ---- Plot -------------------------------------------------------------------

plot_psi_distribution <- ggplot(
  spire_introns_df[abs(spire_introns_df$difference) > 0, ], 
  aes(x = difference, fill = group)
) +
  geom_histogram(aes(y = ..density..), bins = 200, 
                 position = "identity", alpha = 0.9) +
  labs(
    title = "PSI Difference Distribution of Introns in Spire",
    subtitle = paste("One-sample Wilcoxon test p-value for Symmetry:", 
                     format(symmetry_test$p.value, digits = 3)),
    x = "ΔPSI (KO - Control)",
    y = "Density",
    fill = "Splicing Direction",
    caption = paste0("Created by AG on ", Sys.Date())
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    legend.position = "right",
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.text.x = element_text(size = 15),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "gray80", size = 0.5),
    panel.grid.minor = element_line(color = "gray95", size = 0.3),
    strip.text = element_text(size = 14, face = "bold"),
    panel.spacing = unit(1.5, "lines")
  ) +
  scale_x_continuous(
    breaks = c(-100, -75, -50, -25, 0, 25, 50, 75, 100)
  ) +
  # Annotate the percentages with polished label boxes
  annotate("label", x = -65, y = 0.1, hjust = 0, vjust = 0,
           label = sprintf("dPSI < 0: %.1f%%", pct_below),
           fill = alpha("white", 0.8), color = "blue", size = 7, fontface = "bold",
           label.size = 0.5) +
  annotate("label", x = 65, y = 0.1, hjust = 1, vjust = 0,
           label = sprintf("dPSI > 0: %.1f%%", pct_above),
           fill = alpha("white", 0.8), color = "red", size = 7, fontface = "bold",
           label.size = 0.5) +
  scale_fill_brewer(palette = "Set1")


print(plot_psi_distribution)


# Differential Splicing

spire_pdiff_introns <- na.omit(filter(read_csv("notebooks/tables_10n/spire_pdiff_introns.csv")[,-1], abs(deltapsi)>=0.1 & FDR<=0.05))
fmn_pdiff_introns <- na.omit(filter(read_csv("notebooks/tables_10n/fmndko_pdiff_introns.csv")[,-1], abs(deltapsi)>=0.1 & FDR<=0.05))

tao_introns_4d <- na.omit(filter(read_csv("notebooks/tables_10n/tao_pdiff_introns_4d.csv")[,-1], abs(deltapsi)>=0.1 & FDR<=0.05))

