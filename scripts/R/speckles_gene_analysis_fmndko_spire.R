#----
# Speckles gene analysis
#----

# Nuclear Speckle Protein Splicing Correlation { .tabset}

Nuclear speckles genes were obtained from the Human Protein Atlas
Differentially spliced events were filtered to obtain those proteins localizing within speckles. I performed pearson correlation between those selected events. In adittion to the FDR threshold, they also need to have a |dPSI|>=0.1.

## Spire

```{r nuclear spckeles data hadnling spire}
# Load required libraries
library(tidyverse)
library(ComplexHeatmap)

# Set seed for reproducibility
set.seed(43)

# Load nuclear speckles protein atlas data
speckles_proteins <- read_tsv("nuclear_speckles_protein_human_atlas.tsv")

# Filter for spire events associated with nuclear speckle proteins
spire_speckle_prots <- spire_events$PSI %>%
  filter((EVENT %in% differential_spire_exons$EVENT[
    tolower(differential_spire_exons$GENE) %in% tolower(speckles_proteins$Gene) & 
      abs(differential_spire_exons$deltapsi) >= 0.1]) |
      (EVENT %in% differential_spire_introns$EVENT[
        tolower(differential_spire_introns$GENE) %in% tolower(speckles_proteins$Gene) & 
          abs(differential_spire_introns$deltapsi) >= 0.1]) |
      (EVENT %in% differential_spire_alt$EVENT[
        tolower(differential_spire_alt$GENE) %in% tolower(speckles_proteins$Gene) & 
          abs(differential_spire_alt$deltapsi) >= 0.1]))



# Add corresponding gene names next to event names
spire_speckle_prots <- spire_speckle_prots %>%
  mutate(EVENT_GENE = paste(EVENT, GENE, sep = " | "))

rownames(spire_speckle_prots) <- spire_speckle_prots$EVENT_GENE

# Transpose selected columns for correlation analysis
spire_speckle_prots_t <- t(spire_speckle_prots[, 7:12])

# Compute correlation matrix
cor_matrix <- cor(spire_speckle_prots_t, method = "pearson")

# Generate heatmap with enhanced visualization
heatmap <- Heatmap(cor_matrix,
                   name = "Correlation",                      # Color legend title
                   clustering_distance_rows = "euclidean",   # Row clustering method
                   clustering_method_rows = "complete",      # Clustering method for rows
                   show_row_dend = TRUE,                      # Show row dendrogram
                   show_column_dend = TRUE,                   # Show column dendrogram
                   column_title = "Gene Correlation Heatmap", # Title for the heatmap
                   row_title = "Genes",                       # Title for rows
                   column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                   row_title_gp = gpar(fontsize = 14, fontface = "bold"),
                   heatmap_legend_param = list(title_gp = gpar(fontsize = 12))) # Legend formatting

# Draw the heatmap
draw(heatmap)

```

```{r table speckles spire}
# Render DataTable with enhancements

datatable(
  spire_speckle_prots,
  options = list(
    pageLength = 10,                      # Rows per page
    autoWidth = TRUE,                     # Adjust column widths automatically
    dom = 'Bfrtip',                       # Add buttons for export
    buttons = c("copy", "csv", "excel", "pdf"),  # Simplified button definitions
    rownames = FALSE,                       # Disable row names
    extensions = "Buttons"                  # Enable export buttons
  ))
```



## FMNDKO

```{r nuclear spckeles data hadnling fmndko}
set.seed(43)

fmndko_speckle_prots <- fmndko_events$PSI %>%
  filter((EVENT %in% differential_fmndko_exons$EVENT[
    tolower(differential_fmndko_exons$GENE) %in% tolower(speckles_proteins$Gene) & 
      abs(differential_fmndko_exons$deltapsi) >= 0.1]) |
      (EVENT %in% differential_fmndko_introns$EVENT[
        tolower(differential_fmndko_introns$GENE) %in% tolower(speckles_proteins$Gene) & 
          abs(differential_fmndko_introns$deltapsi) >= 0.1]) |
      (EVENT %in% differential_fmndko_alt$EVENT[
        tolower(differential_fmndko_alt$GENE) %in% tolower(speckles_proteins$Gene) & 
          abs(differential_fmndko_alt$deltapsi) >= 0.1]))



# Add corresponding gene names next to event names
fmndko_speckle_prots <- fmndko_speckle_prots %>%
  mutate(EVENT_GENE = paste(EVENT, GENE, sep = " | "))

rownames(fmndko_speckle_prots) <- fmndko_speckle_prots$EVENT_GENE

# Transpose selected columns for correlation analysis
fmndko_speckle_prots_t <- t(fmndko_speckle_prots[, 7:10])

# Compute correlation matrix
cor_matrix <- cor(fmndko_speckle_prots_t, method = "pearson")

# Generate heatmap with enhanced visualization
heatmap <- Heatmap(cor_matrix,
                   name = "Correlation",                      # Color legend title
                   clustering_distance_rows = "euclidean",   # Row clustering method
                   clustering_method_rows = "complete",      # Clustering method for rows
                   show_row_dend = TRUE,                      # Show row dendrogram
                   show_column_dend = TRUE,                   # Show column dendrogram
                   column_title = "Gene Correlation Heatmap", # Title for the heatmap
                   row_title = "Genes",                       # Title for rows
                   column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                   row_title_gp = gpar(fontsize = 14, fontface = "bold"),
                   heatmap_legend_param = list(title_gp = gpar(fontsize = 12))) # Legend formatting

# Draw the heatmap
draw(heatmap)


```

```{r table speckles fmndko}
# Render DataTable with enhancements


datatable(
  fmndko_speckle_prots,
  options = list(
    pageLength = 10,                      # Rows per page
    autoWidth = TRUE,                     # Adjust column widths automatically
    dom = 'Bfrtip',                       # Add buttons for export
    buttons = c("copy", "csv", "excel", "pdf"),  # Simplified button definitions
    rownames = FALSE,                       # Disable row names
    extensions = "Buttons"                  # Enable export buttons
  ))
```