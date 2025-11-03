# Nuclear Speckle Protein Splicing Correlation

Nuclear speckles genes were obtained from the Human Protein Atlas.
Differentially spliced events were filtered to obtain those proteins experimentally localizing within speckles.

```{r nuclear spckeles data hadnling spire}
set.seed(43)

library(readr)
library(dplyr)
library(ComplexHeatmap)
library(grid)  # for gpar()

# Read in the speckles proteins list
speckles_proteins <- read_tsv("nuclear_speckles_protein_human_atlas.tsv")

# Filter events based on gene names and PSI change threshold
pladb_speckle_prots <- pladb_events$PSI %>%
  filter(
    (EVENT %in% differential_pladb10mm_exons$EVENT[
      tolower(differential_pladb10mm_exons$GENE) %in% tolower(speckles_proteins$Gene) & 
        abs(differential_pladb10mm_exons$deltapsi) >= 0.1]) |
      
      (EVENT %in% differential_pladb10mm_introns$EVENT[
        tolower(differential_pladb10mm_introns$GENE) %in% tolower(speckles_proteins$Gene) & 
          abs(differential_pladb10mm_introns$deltapsi) >= 0.1]) |
      
      (EVENT %in% differential_pladb1mm_introns$EVENT[
        tolower(differential_pladb1mm_introns$GENE) %in% tolower(speckles_proteins$Gene) & 
          abs(differential_pladb1mm_introns$deltapsi) >= 0.1]) |
      
      (EVENT %in% differential_pladb1mm_exons$EVENT[
        tolower(differential_pladb1mm_exons$GENE) %in% tolower(speckles_proteins$Gene) & 
          abs(differential_pladb1mm_exons$deltapsi) >= 0.1]) |
      
      # Added alt dataset filter:
      (EVENT %in% differential_pladb10mm_alt$EVENT[
        tolower(differential_pladb10mm_alt$GENE) %in% tolower(speckles_proteins$Gene) &
          abs(differential_pladb10mm_alt$deltapsi) >= 0.1])
  ) %>%
  na.omit()

# Update row names to include gene labels
rownames(pladb_speckle_prots) <- paste(pladb_speckle_prots$EVENT, pladb_speckle_prots$GENE, sep = " | ")

# Transpose the filtered matrix for correlation calculation
# (Assuming columns 7:15 are the numeric PSI values)
pladb_speckle_prots_t <- t(pladb_speckle_prots[, 7:15])

# Calculate the correlation matrix using pairwise complete observations
cor_matrix <- cor(pladb_speckle_prots_t, use = "pairwise.complete.obs")

# Create an empty matrix to store p-values
n <- ncol(pladb_speckle_prots_t)
p_matrix <- matrix(NA, nrow = n, ncol = n)
colnames(p_matrix) <- colnames(pladb_speckle_prots_t)
rownames(p_matrix) <- colnames(pladb_speckle_prots_t)

# Compute p-values for each pair using cor.test()
for (i in 1:n) {
  for (j in i:n) {
    x <- pladb_speckle_prots_t[, i]
    y <- pladb_speckle_prots_t[, j]
    
    # Use only complete cases for this pair
    valid <- complete.cases(x, y)
    if (sum(valid) > 2) {  # require at least 3 pairs to compute correlation
      test <- cor.test(x[valid], y[valid])
      p_matrix[i, j] <- test$p.value
      p_matrix[j, i] <- test$p.value
    } else {
      p_matrix[i, j] <- NA
      p_matrix[j, i] <- NA
    }
  }
}

# Set the significance threshold (alpha)
alpha <- 0.05


# Perform hierarchical clustering on rows
row_dend <- hclust(dist(cor_matrix, method = "euclidean"), method = "complete")

# Extract the order from the row dendrogram
row_order <- row_dend$order

# Apply the same order to the columns
heatmap <- Heatmap(
  cor_matrix[row_order, row_order],   # Reorder both rows and columns
  name = "Correlation",
  show_row_dend = TRUE,               # Show row dendrogram
  show_column_dend = FALSE,            # Hide column dendrogram to maintain order
  column_title = "Correlation Heatmap",
  row_title = "Variables",
  column_title_gp = gpar(fontsize = 14),
  row_title_gp = gpar(fontsize = 14),
  column_names_rot = 45, 
  column_names_gp = gpar(fontsize = 5),
  row_names_gp = gpar(fontsize = 5),
  na_col = "white"
)
draw(heatmap)


```

```{r table speckles pladb}
# Render DataTable with enhancements
datatable(
  pladb_speckle_prots              # Enable export buttons
)
```

# RNA Binding Protein Splicing Correlation

This plot displays splicing events that are included or excluded in a similar pattern. **Why is this important?** PladB alters the splicing patterns of numerous genes. However, among these genes, could there be other splicing or transcription regulators that might contribute to the observed oocyte phenotype?
  
  ```{r rna binding data hadnling spire, fig.width=14, fig.height=10}
set.seed(43)


# Filter events based on gene names and PSI change threshold
pladb_rna_prots <- pladb_events$PSI %>%
  filter(
    (EVENT %in% differential_pladb10mm_exons$EVENT[
      tolower(differential_pladb10mm_exons$GENE) %in% tolower(rna_genes) & 
        abs(differential_pladb10mm_exons$deltapsi) >= 0.1]) |
      
      (EVENT %in% differential_pladb10mm_introns$EVENT[
        tolower(differential_pladb10mm_introns$GENE) %in% tolower(rna_genes) & 
          abs(differential_pladb10mm_introns$deltapsi) >= 0.1]) |
      
      (EVENT %in% differential_pladb1mm_introns$EVENT[
        tolower(differential_pladb1mm_introns$GENE) %in% tolower(rna_genes) & 
          abs(differential_pladb1mm_introns$deltapsi) >= 0.1]) |
      
      (EVENT %in% differential_pladb1mm_exons$EVENT[
        tolower(differential_pladb1mm_exons$GENE) %in% tolower(rna_genes) & 
          abs(differential_pladb1mm_exons$deltapsi) >= 0.1]) |
      
      # Added alt dataset filter:
      (EVENT %in% differential_pladb10mm_alt$EVENT[
        tolower(differential_pladb10mm_alt$GENE) %in% tolower(rna_genes) &
          abs(differential_pladb10mm_alt$deltapsi) >= 0.1])
  ) %>%
  na.omit()

# Update row names to include gene labels
rownames(pladb_rna_prots) <- paste(pladb_rna_prots$EVENT, pladb_rna_prots$GENE, sep = " | ")

# Transpose the filtered matrix for correlation calculation
# (Assuming columns 7:15 are the numeric PSI values)
pladb_rna_prots_t <- t(pladb_rna_prots[, 7:15])

# Calculate the correlation matrix using pairwise complete observations
cor_matrix <- cor(pladb_rna_prots_t, use = "pairwise.complete.obs")

# Create an empty matrix to store p-values
n <- ncol(pladb_rna_prots_t)
p_matrix <- matrix(NA, nrow = n, ncol = n)
colnames(p_matrix) <- colnames(pladb_rna_prots_t)
rownames(p_matrix) <- colnames(pladb_rna_prots_t)

# Compute p-values for each pair using cor.test()
for (i in 1:n) {
  for (j in i:n) {
    x <- pladb_rna_prots_t[, i]
    y <- pladb_rna_prots_t[, j]
    
    # Use only complete cases for this pair
    valid <- complete.cases(x, y)
    if (sum(valid) > 2) {  # require at least 3 pairs to compute correlation
      test <- cor.test(x[valid], y[valid])
      p_matrix[i, j] <- test$p.value
      p_matrix[j, i] <- test$p.value
    } else {
      p_matrix[i, j] <- NA
      p_matrix[j, i] <- NA
    }
  }
}

# Set the significance threshold (alpha)
alpha <- 0.05


# Perform hierarchical clustering on rows
row_dend <- hclust(dist(cor_matrix, method = "euclidean"), method = "complete")

# Extract the order from the row dendrogram
row_order <- row_dend$order

# Apply the same order to the columns
heatmap <- Heatmap(
  cor_matrix[row_order, row_order],   # Reorder both rows and columns
  name = "Correlation",
  show_row_dend = TRUE,               # Show row dendrogram
  show_column_dend = FALSE,            # Hide column dendrogram to maintain order
  column_title = "Correlation Heatmap",
  row_title = "Variables",
  column_title_gp = gpar(fontsize = 14),
  row_title_gp = gpar(fontsize = 14),
  column_names_rot = 45, 
  column_names_gp = gpar(fontsize = 6),
  row_names_gp = gpar(fontsize = 6),
  na_col = "white"
)
draw(heatmap)
```