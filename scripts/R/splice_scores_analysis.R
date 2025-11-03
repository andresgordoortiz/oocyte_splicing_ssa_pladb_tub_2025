library(readr)
library(dplyr)
library(tidyverse)
splice_scores <- read.table("notebooks/SPLICE_SITE_SCORES-mm10.tab",sep = "\t")
events_info<-read.table("notebooks/EVENT_INFO-mm10.tab",sep = "\t", header = T)
spire_pdiff_introns <- na.omit(filter(read_csv("notebooks/tables_10n/spire_pdiff_introns.csv")[,-1], abs(deltapsi)>=0.1 & FDR<=0.05))
fmn_pdiff_introns <- na.omit(filter(read_csv("notebooks/tables_10n/fmndko_pdiff_introns.csv")[,-1], abs(deltapsi)>=0.1 & FDR<=0.05))
spire_pdiff_exons <- na.omit(filter(read_csv("notebooks/tables_10n/spire_pdiff_exons.csv")[,-1], abs(deltapsi)>=0.1 & FDR<=0.05))
fmn_pdiff_exons <- na.omit(filter(read_csv("notebooks/tables_10n/fmndko_pdiff_exons.csv")[,-1], abs(deltapsi)>=0.1 & FDR<=0.05))

calculate_GC_percent <- function(df, sequence_column) {
  # Ensure the specified column exists
  if (!sequence_column %in% colnames(df)) {
    stop("The specified column does not exist in the dataframe.")
  }
  
  # Function to calculate GC content for a single sequence
  gc_content <- function(sequence) {
    # Count G and C in the sequence
    gc_count <- sum(tolower(strsplit(sequence, NULL)[[1]]) %in% c('g', 'c'))
    # Calculate the GC percentage
    gc_percentage <- (gc_count / nchar(sequence)) * 100
    return(gc_percentage)
  }
  
  # Apply the GC content function to each row's sequence
  df$GC_percent <- sapply(df[[sequence_column]], gc_content)
  
  return(df)
}


# Function to extract V3 values based on EVENT
extract_V3_for_event <- function(event_value, splice_scores) {
  # Filter rows where V1 matches the event_value
  filtered_rows <- splice_scores %>%
    filter(V1 == event_value)
  
  # Extract the V3 column and return as a list
  return(list(filtered_rows$V3))
}


# Exon List features
fmn_pdiff_introns$sequence<-events_info$Seq_A[match(fmn_pdiff_introns$EVENT, events_info$EVENT)]

fmn_pdiff_introns$splice_scores <- apply(fmn_pdiff_introns, 1, function(row) {
  # Get the event value from the current row
  event_value <- row["EVENT"]
  
  # Call the function to extract the V3 values
  extract_V3_for_event(event_value, splice_scores)
})

# Unnest the list in V3_values into two separate columns
fmn_pdiff_introns <- fmn_pdiff_introns %>%
  unnest(cols = splice_scores) %>%
  separate(splice_scores, into = c("splice_score_3", "splice_score_5"), sep = ",")

# If necessary, clean up the new columns
fmn_pdiff_introns$splice_score_3 <- gsub("[c()]", "", fmn_pdiff_introns$splice_score_3)
fmn_pdiff_introns$splice_score_5 <- gsub("[c()]", "", fmn_pdiff_introns$splice_score_5)
fmn_pdiff_introns$splice_score_3 <- as.numeric(fmn_pdiff_introns$splice_score_3)
fmn_pdiff_introns$splice_score_5 <- as.numeric(fmn_pdiff_introns$splice_score_5)

fmn_pdiff_introns<-calculate_GC_percent(fmn_pdiff_introns,"sequence")

# Exon List features
fmn_pdiff_exons$sequence<-events_info$Seq_A[match(fmn_pdiff_exons$EVENT, events_info$EVENT)]

fmn_pdiff_exons$splice_scores <- apply(fmn_pdiff_exons, 1, function(row) {
  # Get the event value from the current row
  event_value <- row["EVENT"]
  
  # Call the function to extract the V3 values
  extract_V3_for_event(event_value, splice_scores)
})

# Unnest the list in V3_values into two separate columns
fmn_pdiff_exons <- fmn_pdiff_exons %>%
  unnest(cols = splice_scores) %>%
  separate(splice_scores, into = c("splice_score_3", "splice_score_5"), sep = ",")

# If necessary, clean up the new columns
fmn_pdiff_exons$splice_score_3 <- gsub("[c()]", "", fmn_pdiff_exons$splice_score_3)
fmn_pdiff_exons$splice_score_5 <- gsub("[c()]", "", fmn_pdiff_exons$splice_score_5)
fmn_pdiff_exons$splice_score_3 <- as.numeric(fmn_pdiff_exons$splice_score_3)
fmn_pdiff_exons$splice_score_5 <- as.numeric(fmn_pdiff_exons$splice_score_5)

fmn_pdiff_exons<-calculate_GC_percent(fmn_pdiff_exons,"sequence")

# Exon List features
spire_pdiff_introns$sequence<-events_info$Seq_A[match(spire_pdiff_introns$EVENT, events_info$EVENT)]

spire_pdiff_introns$splice_scores <- apply(spire_pdiff_introns, 1, function(row) {
  # Get the event value from the current row
  event_value <- row["EVENT"]
  
  # Call the function to extract the V3 values
  extract_V3_for_event(event_value, splice_scores)
})

# Unnest the list in V3_values into two separate columns
spire_pdiff_introns <- spire_pdiff_introns %>%
  unnest(cols = splice_scores) %>%
  separate(splice_scores, into = c("splice_score_3", "splice_score_5"), sep = ",")

# If necessary, clean up the new columns
spire_pdiff_introns$splice_score_3 <- gsub("[c()]", "", spire_pdiff_introns$splice_score_3)
spire_pdiff_introns$splice_score_5 <- gsub("[c()]", "", spire_pdiff_introns$splice_score_5)
spire_pdiff_introns$splice_score_3 <- as.numeric(spire_pdiff_introns$splice_score_3)
spire_pdiff_introns$splice_score_5 <- as.numeric(spire_pdiff_introns$splice_score_5)

spire_pdiff_introns<-calculate_GC_percent(spire_pdiff_introns,"sequence")

# Exon List features
spire_pdiff_exons$sequence<-events_info$Seq_A[match(spire_pdiff_exons$EVENT, events_info$EVENT)]

spire_pdiff_exons$splice_scores <- apply(spire_pdiff_exons, 1, function(row) {
  # Get the event value from the current row
  event_value <- row["EVENT"]
  
  # Call the function to extract the V3 values
  extract_V3_for_event(event_value, splice_scores)
})

# Unnest the list in V3_values into two separate columns
spire_pdiff_exons <- spire_pdiff_exons %>%
  unnest(cols = splice_scores) %>%
  separate(splice_scores, into = c("splice_score_3", "splice_score_5"), sep = ",")

# If necessary, clean up the new columns
spire_pdiff_exons$splice_score_3 <- gsub("[c()]", "", spire_pdiff_exons$splice_score_3)
spire_pdiff_exons$splice_score_5 <- gsub("[c()]", "", spire_pdiff_exons$splice_score_5)
spire_pdiff_exons$splice_score_3 <- as.numeric(spire_pdiff_exons$splice_score_3)
spire_pdiff_exons$splice_score_5 <- as.numeric(spire_pdiff_exons$splice_score_5)

spire_pdiff_exons<-calculate_GC_percent(spire_pdiff_exons,"sequence")



splice_scores_wide <- splice_scores %>%
  group_by(V1) %>%
  slice(1:2) %>%  # Keep only the first two rows per V1
  mutate(id = row_number()) %>%  # Create an index: 1 for first row, 2 for second row
  ungroup() %>%   # Remove grouping to ensure pivot_wider uses only V1 as key
  pivot_wider(id_cols = V1, names_from = id, values_from = V3, names_prefix = "ss_") %>%
  rename(ss_3 = ss_1, ss_5 = ss_2)

# Define colors
color1 <- "#1f77b4"  # Blue for All Splice Scores (5')
color2 <- "#d62728"  # Red for FMN2 Pdiff Introns (5')
color3 <- "#2ca02c"  # Green for All Splice Scores (3')
color4 <- "#ff7f0e"  # Orange for FMN2 Pdiff Exons (5')
color5 <- "#8c564b"  # Brown for FMN2 Pdiff Exons (3')
color6 <- "#9467bd"  # Purple for FMN2 Pdiff Introns (3')


# 5' Introns Plot: compare All events vs FMN2 Pdiff Introns (5')
plot_5_introns <- ggplot() +
  # Bottom: All Splice Scores (5')
  geom_density(data = filter(splice_scores_wide, grepl("INT",splice_scores_wide$V1)), 
               aes(x = ss_5, y = ..density.., fill = "All Splice Scores (5')"),
               alpha = 0.6, position = "stack", color = NA) +
  # Top: FMN2 Pdiff Introns (5')
  geom_density(data = fmn_pdiff_introns, 
               aes(x = splice_score_5, y = ..density.., fill = "FMN2 Pdiff Introns (5')"),
               alpha = 0.6, position = "stack", color = NA) +
  labs(title = "5' Splice Scores - Introns",
       x = "Score Value",
       y = "Density",
       fill = "Dataset/Region") +
  scale_fill_manual(values = c("All Splice Scores (5')" = color1,
                               "FMN2 Pdiff Introns (5')" = color2)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 16),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# 5' Exons Plot: compare All events vs FMN2 Pdiff Exons (5')
plot_5_exons <- ggplot() +
  # Bottom: All Splice Scores (5')
  geom_density(data = filter(splice_scores_wide, grepl("EX",splice_scores_wide$V1)), 
               aes(x = ss_5, y = ..density.., fill = "All Splice Scores (5')"),
               alpha = 0.6, position = "stack", color = NA) +
  # Top: FMN2 Pdiff Exons (5')
  geom_density(data = fmn_pdiff_exons, 
               aes(x = splice_score_5, y = ..density.., fill = "FMN2 Pdiff Exons (5')"),
               alpha = 0.6, position = "stack", color = NA) +
  labs(title = "5' Splice Scores - Exons",
       x = "Score Value",
       y = "Density",
       fill = "Dataset/Region") +
  scale_fill_manual(values = c("All Splice Scores (5')" = color1,
                               "FMN2 Pdiff Exons (5')" = color4)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 16),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# 3' Introns Plot: compare All events vs FMN2 Pdiff Introns (3')
plot_3_introns <- ggplot() +
  # Bottom: All Splice Scores (3')
  geom_density(data = filter(splice_scores_wide, grepl("INT",splice_scores_wide$V1)), 
               aes(x = ss_3, y = ..density.., fill = "All Splice Scores (3')"),
               alpha = 0.6, position = "stack", color = NA) +
  # Top: FMN2 Pdiff Introns (3')
  geom_density(data = fmn_pdiff_introns, 
               aes(x = splice_score_3, y = ..density.., fill = "FMN2 Pdiff Introns (3')"),
               alpha = 0.6, position = "stack", color = NA) +
  labs(title = "3' Splice Scores - Introns",
       x = "Score Value",
       y = "Density",
       fill = "Dataset/Region") +
  scale_fill_manual(values = c("All Splice Scores (3')" = color1,
                               "FMN2 Pdiff Introns (3')" = color6)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 16),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# 3' Exons Plot: compare All events vs FMN2 Pdiff Exons (3')
plot_3_exons <- ggplot() +
  # Bottom: All Splice Scores (3')
  geom_density(data = filter(splice_scores_wide, grepl("EX",splice_scores_wide$V1)), 
               aes(x = ss_3, y = ..density.., fill = "All Splice Scores (3')"),
               alpha = 0.6, position = "stack", color = NA) +
  # Top: FMN2 Pdiff Exons (3')
  geom_density(data = fmn_pdiff_exons, 
               aes(x = splice_score_3, y = ..density.., fill = "FMN2 Pdiff Exons (3')"),
               alpha = 0.6, position = "stack", color = NA) +
  labs(title = "3' Splice Scores - Exons",
       x = "Score Value",
       y = "Density",
       fill = "Dataset/Region") +
  scale_fill_manual(values = c("All Splice Scores (3')" = color1,
                               "FMN2 Pdiff Exons (3')" = color5)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 16),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# Arrange the four plots in a 2x2 grid
library(gridExtra)
combined_plot_fmn2 <- grid.arrange(plot_5_introns, plot_5_exons, 
                              plot_3_introns, plot_3_exons, ncol = 2)

# Print the combined plot
print(combined_plot_fmn2)


# 5' Introns Plot: compare All events vs SPIRE2 Pdiff Introns (5')
plot_5_introns <- ggplot() +
  # Bottom: All Splice Scores (5')
  geom_density(data = filter(splice_scores_wide, grepl("INT",splice_scores_wide$V1)), 
               aes(x = ss_5, y = ..density.., fill = "All Splice Scores (5')"),
               alpha = 0.6, position = "stack", color = NA) +
  # Top: SPIRE2 Pdiff Introns (5')
  geom_density(data = spire_pdiff_introns, 
               aes(x = splice_score_5, y = ..density.., fill = "SPIRE2 Pdiff Introns (5')"),
               alpha = 0.6, position = "stack", color = NA) +
  labs(title = "5' Splice Scores - Introns",
       x = "Score Value",
       y = "Density",
       fill = "Dataset/Region") +
  scale_fill_manual(values = c("All Splice Scores (5')" = color1,
                               "SPIRE2 Pdiff Introns (5')" = color2)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 16),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# 5' Exons Plot: compare All events vs SPIRE2 Pdiff Exons (5')
plot_5_exons <- ggplot() +
  # Bottom: All Splice Scores (5')
  geom_density(data = filter(splice_scores_wide, grepl("EX",splice_scores_wide$V1)), 
               aes(x = ss_5, y = ..density.., fill = "All Splice Scores (5')"),
               alpha = 0.6, position = "stack", color = NA) +
  # Top: SPIRE2 Pdiff Exons (5')
  geom_density(data = spire_pdiff_exons, 
               aes(x = splice_score_5, y = ..density.., fill = "SPIRE2 Pdiff Exons (5')"),
               alpha = 0.6, position = "stack", color = NA) +
  labs(title = "5' Splice Scores - Exons",
       x = "Score Value",
       y = "Density",
       fill = "Dataset/Region") +
  scale_fill_manual(values = c("All Splice Scores (5')" = color1,
                               "SPIRE2 Pdiff Exons (5')" = color4)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 16),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# 3' Introns Plot: compare All events vs SPIRE2 Pdiff Introns (3')
plot_3_introns <- ggplot() +
  # Bottom: All Splice Scores (3')
  geom_density(data = filter(splice_scores_wide, grepl("INT",splice_scores_wide$V1)), 
               aes(x = ss_3, y = ..density.., fill = "All Splice Scores (3')"),
               alpha = 0.6, position = "stack", color = NA) +
  # Top: SPIRE2 Pdiff Introns (3')
  geom_density(data = spire_pdiff_introns, 
               aes(x = splice_score_3, y = ..density.., fill = "SPIRE2 Pdiff Introns (3')"),
               alpha = 0.6, position = "stack", color = NA) +
  labs(title = "3' Splice Scores - Introns",
       x = "Score Value",
       y = "Density",
       fill = "Dataset/Region") +
  scale_fill_manual(values = c("All Splice Scores (3')" = color1,
                               "SPIRE2 Pdiff Introns (3')" = color6)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 16),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# 3' Exons Plot: compare All events vs SPIRE2 Pdiff Exons (3')
plot_3_exons <- ggplot() +
  # Bottom: All Splice Scores (3')
  geom_density(data = filter(splice_scores_wide, grepl("EX",splice_scores_wide$V1)), 
               aes(x = ss_3, y = ..density.., fill = "All Splice Scores (3')"),
               alpha = 0.6, position = "stack", color = NA) +
  # Top: SPIRE2 Pdiff Exons (3')
  geom_density(data = spire_pdiff_exons, 
               aes(x = splice_score_3, y = ..density.., fill = "SPIRE2 Pdiff Exons (3')"),
               alpha = 0.6, position = "stack", color = NA) +
  labs(title = "3' Splice Scores - Exons",
       x = "Score Value",
       y = "Density",
       fill = "Dataset/Region") +
  scale_fill_manual(values = c("All Splice Scores (3')" = color1,
                               "SPIRE2 Pdiff Exons (3')" = color5)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 16),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# Arrange the four plots in a 2x2 grid
library(gridExtra)
combined_plot_spire <- grid.arrange(plot_5_introns, plot_5_exons, 
                                     plot_3_introns, plot_3_exons, ncol = 2)

# Print the combined plot
print(combined_plot_spire)


# Load required libraries
library(GenomicRanges)
library(Biostrings)

# Convert the sequence column to a DNAStringSet
sequences <- DNAStringSet(fmn_pdiff_introns$sequence)

# Define the pattern
pattern <- "GGAA"


# Function to compute the proportion of pattern in a sequence
compute_proportion <- function(seq, pattern) {
  matches <- matchPattern(pattern, seq)
  # Use nchar(as.character(seq)) to get the length of the sequence
  seq_length <- nchar(as.character(seq))
  proportion <- sum(width(matches)) / seq_length
  return(proportion)
}

# Apply the function to each sequence in the data frame
fmn_pdiff_introns$ggaa_proportion <- sapply(sequences, compute_proportion, pattern = pattern)
sequences_all <- DNAStringSet(sample(filter(na.omit(events_info), grepl("INT", EVENT))$Seq_A, 5000))
ggaa_all <- data.frame(pattern=sapply(sequences_all, compute_proportion, pattern = pattern))
# Print results
library(ggplot2)

ggplot() +
  geom_boxplot(data = fmn_pdiff_introns, aes(x = factor(1), y = ggaa_proportion), fill = "blue", alpha = 0.5) +
  geom_boxplot(data = ggaa_all, aes(x = factor(2), y = pattern), fill = "red", alpha = 0.5) +
  scale_x_discrete(labels = c("fmn_pdiff_exons", "ggaa_all")) +
  labs(x = "Dataset", y = "Value", title = "Comparison of Distributions") +
  theme_minimal()

