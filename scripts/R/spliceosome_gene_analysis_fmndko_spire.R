#-----
# Spliceosome analysis
#-----


# Spliceosome Ontology Genes { .tabset}


## Spire

```{r, splicesosomo ontology behaviour}

enrichdata<-enrichr(unique(c(differential_spire_exons$GENE, differential_spire_introns$GENE, differential_spire_alt$GENE,differential_fmndko_exons$GENE, differential_fmndko_introns$GENE, differential_fmndko_alt$GENE)), databases = c("GO_Biological_Process_2023","GO_Cellular_Component_2023","GO_Molecular_Function_2023"))

# Extract and format gene names
genes_splicing_1 <- enrichdata[[1]] %>%
  filter(str_detect(Term, "Splic")) %>%  # Select rows with 'splic' in Term
  pull(Genes) %>%                         # Extract Genes column
  str_split(";") %>%                      # Split gene names
  unlist() %>%                            # Flatten list
  str_to_lower() %>%                      # Convert to lowercase
  str_replace("^.", toupper) 
# Capitalize first letter
genes_splicing_2 <- enrichdata[[2]] %>%
  filter(str_detect(Term, "splic") | str_detect(Term, "Splic")) %>%  # Select rows with 'splic' in Term
  pull(Genes) %>%                         # Extract Genes column
  str_split(";") %>%                      # Split gene names
  unlist() %>%                            # Flatten list
  str_to_lower() %>%                      # Convert to lowercase
  str_replace("^.", toupper) 

genes_splicing<-unique(c(genes_splicing_1,genes_splicing_2))

combined_df <- bind_rows(
  exons_spire = differential_spire_exons,
  introns_spire = differential_spire_introns,
  .id = "source"  # This will be the new column indicating the origin
)
splicing_df<-filter(combined_df, GENE %in% genes_splicing)
splicing_df$protein_impact<-protein_impact$ONTO[match(splicing_df$EVENT, protein_impact$EventID)]



# Prepare the summary data frame, ensuring proper sorting and factor ordering
splicing_df_summary <- splicing_df %>%
  # Remove duplicate rows for the same gene and event (if applicable)
  distinct(GENE, EVENT, .keep_all = TRUE) %>%
  group_by(GENE) %>%
  summarise(
    nEvents   = n(),                                        # count unique events per gene
    meanDelta = mean(abs(deltapsi), na.rm = TRUE),          # mean of absolute deltapsi
    sdDelta   = sd(abs(deltapsi), na.rm = TRUE)             # standard deviation (for reference)
  ) %>%
  # First sort by the number of events (desc), then by mean absolute deltapsi (desc)
  arrange(desc(nEvents), desc(meanDelta)) %>%
  # Set the factor levels for GENE according to the sorted order
  mutate(GENE = factor(GENE, levels = unique(GENE)))

# Update the original data frame:
# - Create a column for the absolute value of deltapsi.
# - Ensure that GENE is a factor with levels set by the summary.
# - Create a new column 'impact' that indicates whether protein_impact contains "ORF".
splicing_df <- splicing_df %>%
  mutate(
    absDelta = abs(deltapsi),
    GENE = factor(GENE, levels = levels(splicing_df_summary$GENE)),
    impact = ifelse(grepl("ORF", protein_impact, ignore.case = TRUE), "ORF", "Non-ORF")
  )

# Create the publication-ready plot:
p <- ggplot() +
  # Bar plot showing mean absolute deltapsi per gene, colored by the number of events
  geom_bar(
    data = splicing_df_summary, 
    aes(x = GENE, y = meanDelta, fill = factor(nEvents)),
    stat = "identity", 
    width = 0.7, 
    alpha = 0.8
  ) +
  # Overlay individual data points with jitter, mapping shape by the 'impact' column
  geom_jitter(
    data = splicing_df,
    aes(x = GENE, y = absDelta, shape = impact),
    width = 0.15,    # small horizontal jitter to reduce overplotting
    size = 2.5, 
    alpha = 0.6,     # adjust transparency
    color = "black"  # point border color
  ) +
  # Use a pleasing color palette for the bars (from RColorBrewer)
  scale_fill_brewer(palette = "Set2", name = "Number of Events") +
  # Manually set the shapes: for example, 16 (circle) for Non-ORF and 17 (triangle) for ORF
  scale_shape_manual(values = c("Non-ORF" = 16, "ORF" = 17), name = "Protein Impact") +
  # Publication-quality theme settings
  theme_bw(base_size = 14) +
  labs(
    x = "Gene ID",
    y = "|dPSI|",
    title = "Splicing related Genes in Spire"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# Display the plot
print(p)


```


```{r splicing plots spire, warning=FALSE, message=FALSE, fig.width=16, fig.height=8}
spire_shared<-filter(spire_events$PSI, EVENT %in% splicing_df$EVENT)

# Add GENE names to the EVENT column
spire_shared <- spire_shared %>%
  mutate(EVENT = paste0(EVENT, " (", GENE, ")"))

spire_shared <- spire_shared %>%
  pivot_longer(
    cols = 7:10,  # Only selecting columns 7 to 10
    names_to = "Condition",
    values_to = "Value"
  ) %>%
  mutate(
    Group = case_when(
      Condition %in% names(spire_shared)[7:9] ~ "Control",
      Condition %in% names(spire_shared)[10:12] ~ "spireKO"
    )
  )


spire_shared$Group<-factor(spire_shared$Group, levels = c("Control","spireKO"))
# Plot with updated EVENT names
ggplot(spire_shared, aes(x = Group, y = Value, color = EVENT)) +
  geom_boxplot(size = 1,show.legend = FALSE) +
  labs(
    title = "Events of Splicing Proteins Spire",
    x = "",
    y = "PSI"
  ) +
  theme_minimal() +
  facet_wrap(~EVENT)
```



## FMNDKO

```{r splicing fmndko}

enrichdata<-enrichr(unique(c(differential_spire_exons$GENE, differential_spire_introns$GENE, differential_spire_alt$GENE,differential_fmndko_exons$GENE, differential_fmndko_introns$GENE, differential_fmndko_alt$GENE)), databases = c("GO_Biological_Process_2023","GO_Cellular_Component_2023","GO_Molecular_Function_2023"))

# Extract and format gene names
genes_splicing_1 <- enrichdata[[1]] %>%
  filter(str_detect(Term, "Splic")) %>%  # Select rows with 'splic' in Term
  pull(Genes) %>%                         # Extract Genes column
  str_split(";") %>%                      # Split gene names
  unlist() %>%                            # Flatten list
  str_to_lower() %>%                      # Convert to lowercase
  str_replace("^.", toupper) 
# Capitalize first letter
genes_splicing_2 <- enrichdata[[2]] %>%
  filter(str_detect(Term, "splic") | str_detect(Term, "Splic")) %>%  # Select rows with 'splic' in Term
  pull(Genes) %>%                         # Extract Genes column
  str_split(";") %>%                      # Split gene names
  unlist() %>%                            # Flatten list
  str_to_lower() %>%                      # Convert to lowercase
  str_replace("^.", toupper) 

genes_splicing<-unique(c(genes_splicing_1,genes_splicing_2))

combined_df <- bind_rows(
  exons_fmndko = differential_fmndko_exons,
  introns_fmndko = differential_fmndko_introns,
  .id = "source"  # This will be the new column indicating the origin
)
splicing_df<-filter(combined_df, GENE %in% genes_splicing)
splicing_df$protein_impact<-protein_impact$ONTO[match(splicing_df$EVENT, protein_impact$EventID)]



# Prepare the summary data frame, ensuring proper sorting and factor ordering
splicing_df_summary <- splicing_df %>%
  # Remove duplicate rows for the same gene and event (if applicable)
  distinct(GENE, EVENT, .keep_all = TRUE) %>%
  group_by(GENE) %>%
  summarise(
    nEvents   = n(),                                        # count unique events per gene
    meanDelta = mean(abs(deltapsi), na.rm = TRUE),          # mean of absolute deltapsi
    sdDelta   = sd(abs(deltapsi), na.rm = TRUE)             # standard deviation (for reference)
  ) %>%
  # First sort by the number of events (desc), then by mean absolute deltapsi (desc)
  arrange(desc(nEvents), desc(meanDelta)) %>%
  # Set the factor levels for GENE according to the sorted order
  mutate(GENE = factor(GENE, levels = unique(GENE)))

# Update the original data frame:
# - Create a column for the absolute value of deltapsi.
# - Ensure that GENE is a factor with levels set by the summary.
# - Create a new column 'impact' that indicates whether protein_impact contains "ORF".
splicing_df <- splicing_df %>%
  mutate(
    absDelta = abs(deltapsi),
    GENE = factor(GENE, levels = levels(splicing_df_summary$GENE)),
    impact = ifelse(grepl("ORF", protein_impact, ignore.case = TRUE), "ORF", "Non-ORF")
  )

# Create the publication-ready plot:
p <- ggplot() +
  # Bar plot showing mean absolute deltapsi per gene, colored by the number of events
  geom_bar(
    data = splicing_df_summary, 
    aes(x = GENE, y = meanDelta, fill = factor(nEvents)),
    stat = "identity", 
    width = 0.7, 
    alpha = 0.8
  ) +
  # Overlay individual data points with jitter, mapping shape by the 'impact' column
  geom_jitter(
    data = splicing_df,
    aes(x = GENE, y = absDelta, shape = impact),
    width = 0.15,    # small horizontal jitter to reduce overplotting
    size = 2.5, 
    alpha = 0.6,     # adjust transparency
    color = "black"  # point border color
  ) +
  # Use a pleasing color palette for the bars (from RColorBrewer)
  scale_fill_brewer(palette = "Set2", name = "Number of Events") +
  # Manually set the shapes: for example, 16 (circle) for Non-ORF and 17 (triangle) for ORF
  scale_shape_manual(values = c("Non-ORF" = 16, "ORF" = 17), name = "Protein Impact") +
  # Publication-quality theme settings
  theme_bw(base_size = 14) +
  labs(
    x = "Gene ID",
    y = "|dPSI|",
    title = "Splicing related Genes in FMNDKO"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# Display the plot
print(p)


```

```{r splicing plots, warning=FALSE, message=FALSE, fig.width=16, fig.height=8}
fmndko_shared<-filter(fmndko_events$PSI, EVENT %in% splicing_df$EVENT)

# Add GENE names to the EVENT column
fmndko_shared <- fmndko_shared %>%
  mutate(EVENT = paste0(EVENT, " (", GENE, ")"))

fmndko_shared <- fmndko_shared %>%
  pivot_longer(
    cols = 7:10,  # Only selecting columns 7 to 10
    names_to = "Condition",
    values_to = "Value"
  ) %>%
  mutate(
    Group = case_when(
      Condition %in% names(fmndko_shared)[7:8] ~ "Control",
      Condition %in% names(fmndko_shared)[9:10] ~ "fmndko"
    )
  )


fmndko_shared$Group<-factor(fmndko_shared$Group, levels = c("Control","fmndko"))
# Plot with updated EVENT names
ggplot(fmndko_shared, aes(x = Group, y = Value, color = EVENT)) +
  geom_boxplot(size = 1,show.legend = FALSE) +
  labs(
    title = "Events of Splicing Proteins FMNDKO",
    x = "",
    y = "PSI"
  ) +
  theme_minimal() +
  facet_wrap(~EVENT)
```