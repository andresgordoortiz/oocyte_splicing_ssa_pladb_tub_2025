#----
# Mtascape Analysis
#----

# Metascape Ontology Analysis { .tabset}

## Cell Cycle related

set.seed(44)
metascape_alltocontrol <- read_excel("results/metascape_n10/pladb_all_n10.xlsx")
metascape_allgenes_10mmvscontrol <- read_excel("results/metascape_n10/pladb_10vscontrol_n10.xlsx")
metascape_allgenes_1mmvscontrol <- read_excel("results/metascape_n10/pladb_1vscontrol_n10.xlsx")
metascape_allgenes_10mmvs1mm <- read_excel("results/metascape_n10/pladb_10vs1_n10.xlsx")

cell_cycle_genes<-unique(c(
  metascape_alltocontrol$MyList[which(metascape_alltocontrol$`GO:0010564 regulation of cell cycle proce`=="1.0")],
  metascape_allgenes_10mmvscontrol$MyList[which(metascape_allgenes_10mmvscontrol$`GO:0010564 regulation of cell cycle proce`=="1.0" | metascape_allgenes_10mmvscontrol$`GO:1903047 mitotic cell cycle process`=="1.0" )],
  metascape_allgenes_1mmvscontrol$MyList[which(metascape_allgenes_1mmvscontrol$`GO:1901987 regulation of cell cycle phase`=="1.0")])
)

combined_df <- bind_rows(
  exons_10mm = differential_pladb10mm_exons,
  exons_1mm = differential_pladb1mm_exons,
  introns_10mm = differential_pladb10mm_introns,
  introns_1mm = differential_pladb1mm_introns,
  alt_10mm = differential_pladb10mm_alt,
  alt_1mm = differential_pladb1mm_alt,
  .id = "source"  # This will be the new column indicating the origin
)
cell_cycle_df<-filter(combined_df, GENE %in% cell_cycle_genes)
cell_cycle_df$protein_impact<-protein_impact$ONTO[match(cell_cycle_df$EVENT, protein_impact$EventID)]



# Prepare the summary data frame, ensuring proper sorting and factor ordering
cell_cycle_df_summary <- cell_cycle_df %>%
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
cell_cycle_df <- cell_cycle_df %>%
  mutate(
    absDelta = abs(deltapsi),
    GENE = factor(GENE, levels = levels(cell_cycle_df_summary$GENE)),
    impact = ifelse(grepl("ORF", protein_impact, ignore.case = TRUE), "ORF", "Non-ORF")
  )

# Create the publication-ready plot:
p <- ggplot() +
  # Bar plot showing mean absolute deltapsi per gene, colored by the number of events
  geom_bar(
    data = cell_cycle_df_summary, 
    aes(x = GENE, y = meanDelta, fill = factor(nEvents)),
    stat = "identity", 
    width = 0.7, 
    alpha = 0.8
  ) +
  # Overlay individual data points with jitter, mapping shape by the 'impact' column
  geom_jitter(
    data = cell_cycle_df,
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
    title = "Cell Cycle related Genes"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# Display the plot
print(p)





event_names <- cell_cycle_df$EVENT[cell_cycle_df$impact=="ORF"]
pladb_events_shared<-filter(pladb_events$PSI, EVENT %in% event_names)

# Add GENE names to the EVENT column
pladb_events_shared <- pladb_events_shared %>%
  mutate(EVENT = paste0(EVENT, " (", GENE, ")"))

# Reshape data for plotting
pladb_events_shared <- pladb_events_shared %>%
  pivot_longer(
    cols = 7:15,
    names_to = "Condition",
    values_to = "Value"
  ) %>%
  mutate(
    Group = case_when(
      grepl("038|039|040", Condition) ~ "Control",
      grepl("041|042|043", Condition) ~ "Pladb1mM",
      grepl("044|045|046", Condition) ~ "Pladb10mM"
    )
  )

pladb_events_shared$Group<-factor(pladb_events_shared$Group, levels = c("Control","Pladb1mM","Pladb10mM"))
# Plot with updated EVENT names
ggplot(pladb_events_shared, aes(x = Group, y = Value, color = EVENT)) +
  geom_boxplot(size = 1,show.legend = FALSE) +
  labs(
    title = "ORF Disrupted Events Cell Cycle Related",
    x = "",
    y = "PSI"
  ) +
  theme_minimal() +
  facet_wrap(~EVENT)

## Cytoskeleton related


set.seed(44)



cytoskeleton_genes<-unique(c(
  metascape_alltocontrol$MyList[which(metascape_alltocontrol$`GO:0097749 membrane tubulation`=="1.0")],
  metascape_allgenes_10mmvscontrol$MyList[which(metascape_allgenes_10mmvscontrol$`GO:0032465 regulation of cytokinesis`=="1.0" | metascape_allgenes_10mmvscontrol$`GO:0010638 positive regulation of organel`=="1.0")],
  metascape_allgenes_10mmvs1mm$MyList[which(metascape_allgenes_10mmvs1mm$`GO:0022027 interkinetic nuclear migration`=="1.0" | metascape_allgenes_10mmvs1mm$`GO:0000226 microtubule cytoskeleton organ`=="1.0" | metascape_allgenes_10mmvs1mm$`GO:0010638 positive regulation of organel`=="1.0" | metascape_allgenes_10mmvs1mm$`GO:0046599 regulation of centriole replic`=="1.0")],
  metascape_allgenes_1mmvscontrol$MyList[which(metascape_allgenes_1mmvscontrol$`GO:0051640 organelle localization`=="1.0" | metascape_allgenes_1mmvscontrol$`GO:0051168 nuclear export`=="1.0" | metascape_allgenes_1mmvscontrol$`GO:0030036 actin cytoskeleton organizatio`=="1.0" | metascape_allgenes_1mmvscontrol$`GO:0097749 membrane tubulation`=="1.0")]))
  

combined_df <- bind_rows(
  exons_10mm = differential_pladb10mm_exons,
  exons_1mm = differential_pladb1mm_exons,
  introns_10mm = differential_pladb10mm_introns,
  introns_1mm = differential_pladb1mm_introns,
  alt_10mm = differential_pladb10mm_alt,
  alt_1mm = differential_pladb1mm_alt,
  .id = "source"  # This will be the new column indicating the origin
)

cytoskeleton_df<-filter(combined_df, GENE %in% cytoskeleton_genes)
cytoskeleton_df$protein_impact<-protein_impact$ONTO[match(cytoskeleton_df$EVENT, protein_impact$EventID)]



# Prepare the summary data frame, ensuring proper sorting and factor ordering
cytoskeleton_df_summary <- cytoskeleton_df %>%
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
cytoskeleton_df <- cytoskeleton_df %>%
  mutate(
    absDelta = abs(deltapsi),
    GENE = factor(GENE, levels = levels(cytoskeleton_df_summary$GENE)),
    impact = ifelse(grepl("ORF", protein_impact, ignore.case = TRUE), "ORF", "Non-ORF")
  )

# Create the publication-ready plot:
p <- ggplot() +
  # Bar plot showing mean absolute deltapsi per gene, colored by the number of events
  geom_bar(
    data = cytoskeleton_df_summary, 
    aes(x = GENE, y = meanDelta, fill = factor(nEvents)),
    stat = "identity", 
    width = 0.7, 
    alpha = 0.8
  ) +
  # Overlay individual data points with jitter, mapping shape by the 'impact' column
  geom_jitter(
    data = cytoskeleton_df,
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
    title = "Cytoskeleton related Genes"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# Display the plot
print(p)




event_names <- cytoskeleton_df$EVENT[cytoskeleton_df$impact=="ORF"]
pladb_events_shared<-filter(pladb_events$PSI, EVENT %in% event_names)

# Add GENE names to the EVENT column
pladb_events_shared <- pladb_events_shared %>%
  mutate(EVENT = paste0(EVENT, " (", GENE, ")"))

# Reshape data for plotting
pladb_events_shared <- pladb_events_shared %>%
  pivot_longer(
    cols = 7:15,
    names_to = "Condition",
    values_to = "Value"
  ) %>%
  mutate(
    Group = case_when(
      grepl("038|039|040", Condition) ~ "Control",
      grepl("041|042|043", Condition) ~ "Pladb1mM",
      grepl("044|045|046", Condition) ~ "Pladb10mM"
    )
  )

pladb_events_shared$Group<-factor(pladb_events_shared$Group, levels = c("Control","Pladb1mM","Pladb10mM"))
# Plot with updated EVENT names
ggplot(pladb_events_shared, aes(x = Group, y = Value, color = EVENT)) +
  geom_boxplot(size = 1,show.legend = FALSE) +
  labs(
    title = "ORF Disrupted Events Cytoskeleton Related",
    x = "",
    y = "PSI"
  ) +
  theme_minimal() +
  facet_wrap(~EVENT)



## Chromosome related


chromosome_genes<-unique(c(
  metascape_alltocontrol$MyList[which(metascape_alltocontrol$`GO:0051053 negative regulation of DNA met`=="1.0" | metascape_alltocontrol$`GO:0033044 regulation of chromosome organ`=="1.0" | metascape_alltocontrol$`R-MMU-73894 DNA Repair`=="1.0" | metascape_alltocontrol$`GO:0006974 DNA damage response`=="1.0")],
  metascape_allgenes_10mmvscontrol$MyList[which(metascape_allgenes_10mmvscontrol$`R-MMU-6804115 TP53 regulates transcription o`=="1.0" | metascape_allgenes_10mmvscontrol$`GO:0033044 regulation of chromosome organ`=="1.0" | metascape_allgenes_10mmvscontrol$`GO:2000104 negative regulation of DNA-tem`=="1.0" | metascape_allgenes_10mmvscontrol$`GO:0051053 negative regulation of DNA met`=="1.0")],
  metascape_allgenes_10mmvs1mm$MyList[which(metascape_allgenes_10mmvs1mm$`GO:0034502 protein localization to chromo`=="1.0" | metascape_allgenes_10mmvs1mm$`GO:0006307 DNA alkylation repair`=="1.0" | metascape_allgenes_10mmvs1mm$`GO:0006974 DNA damage response`=="1.0")],
  metascape_allgenes_1mmvscontrol$MyList[which(metascape_allgenes_1mmvscontrol$`GO:0006271 DNA strand elongation involved`=="1.0" | metascape_allgenes_1mmvscontrol$`GO:0006281 DNA repair`=="1.0" )]))

chromosome_df<-filter(combined_df, GENE %in% chromosome_genes)
chromosome_df$protein_impact<-protein_impact$ONTO[match(chromosome_df$EVENT, protein_impact$EventID)]



# Prepare the summary data frame, ensuring proper sorting and factor ordering
chromosome_df_summary <- chromosome_df %>%
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
chromosome_df <- chromosome_df %>%
  mutate(
    absDelta = abs(deltapsi),
    GENE = factor(GENE, levels = levels(chromosome_df_summary$GENE)),
    impact = ifelse(grepl("ORF", protein_impact, ignore.case = TRUE), "ORF", "Non-ORF")
  )

# Create the publication-ready plot:
p <- ggplot() +
  # Bar plot showing mean absolute deltapsi per gene, colored by the number of events
  geom_bar(
    data = chromosome_df_summary, 
    aes(x = GENE, y = meanDelta, fill = factor(nEvents)),
    stat = "identity", 
    width = 0.7, 
    alpha = 0.8
  ) +
  # Overlay individual data points with jitter, mapping shape by the 'impact' column
  geom_jitter(
    data = chromosome_df,
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
    title = "Chromosome related Genes"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# Display the plot
print(p)


event_names <- chromosome_df$EVENT[chromosome_df$impact=="ORF"]
pladb_events_shared<-filter(pladb_events$PSI, EVENT %in% event_names)

# Add GENE names to the EVENT column
pladb_events_shared <- pladb_events_shared %>%
  mutate(EVENT = paste0(EVENT, " (", GENE, ")"))

# Reshape data for plotting
pladb_events_shared <- pladb_events_shared %>%
  pivot_longer(
    cols = 7:15,
    names_to = "Condition",
    values_to = "Value"
  ) %>%
  mutate(
    Group = case_when(
      grepl("038|039|040", Condition) ~ "Control",
      grepl("041|042|043", Condition) ~ "Pladb1mM",
      grepl("044|045|046", Condition) ~ "Pladb10mM"
    )
  )

pladb_events_shared$Group<-factor(pladb_events_shared$Group, levels = c("Control","Pladb1mM","Pladb10mM"))
# Plot with updated EVENT names
ggplot(pladb_events_shared, aes(x = Group, y = Value, color = EVENT)) +
  geom_boxplot(size = 1,show.legend = FALSE) +
  labs(
    title = "ORF Disrupted Events Chromosome Related",
    x = "",
    y = "PSI"
  ) +
  theme_minimal() +
  facet_wrap(~EVENT)



## RNA related 



splicesome_genes<-unique(c(
  metascape_alltocontrol$MyList[which(metascape_alltocontrol$`GO:0006366 transcription by RNA polymeras`=="1.0" | metascape_alltocontrol$`R-MMU-3700989 Transcriptional Regulation by`=="1.0" | metascape_alltocontrol$`GO:0000432 positive regulation of transcr`=="1.0" | metascape_alltocontrol$`WP310 mRNA processing`=="1.0")],
  metascape_allgenes_10mmvscontrol$MyList[which(metascape_allgenes_10mmvscontrol$`GO:0006417 regulation of translation`=="1.0" | metascape_allgenes_10mmvscontrol$`GO:0016071 mRNA metabolic process`=="1.0")],
  metascape_allgenes_10mmvs1mm$MyList[which(metascape_allgenes_10mmvs1mm$`GO:0006417 regulation of translation`=="1.0" | metascape_allgenes_10mmvs1mm$`GO:0016071 mRNA metabolic process`=="1.0" | metascape_allgenes_10mmvs1mm$`GO:1900368 regulation of post-transcripti`=="1.0")],
  metascape_allgenes_1mmvscontrol$MyList[which(metascape_allgenes_1mmvscontrol$`GO:0009451 RNA modification`=="1.0" | metascape_allgenes_1mmvscontrol$`WP310 mRNA processing`=="1.0" )]))

spliceosome_df<-filter(combined_df, GENE %in% splicesome_genes)
spliceosome_df$protein_impact<-protein_impact$ONTO[match(spliceosome_df$EVENT, protein_impact$EventID)]



# Prepare the summary data frame, ensuring proper sorting and factor ordering
spliceosome_df_summary <- spliceosome_df %>%
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
spliceosome_df <- spliceosome_df %>%
  mutate(
    absDelta = abs(deltapsi),
    GENE = factor(GENE, levels = levels(spliceosome_df_summary$GENE)),
    impact = ifelse(grepl("ORF", protein_impact, ignore.case = TRUE), "ORF", "Non-ORF")
  )

# Create the publication-ready plot:
p <- ggplot() +
  # Bar plot showing mean absolute deltapsi per gene, colored by the number of events
  geom_bar(
    data = spliceosome_df_summary, 
    aes(x = GENE, y = meanDelta, fill = factor(nEvents)),
    stat = "identity", 
    width = 0.7, 
    alpha = 0.8
  ) +
  # Overlay individual data points with jitter, mapping shape by the 'impact' column
  geom_jitter(
    data = spliceosome_df,
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
    title = "RNA related Genes"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# Display the plot
print(p)


event_names <- spliceosome_df$EVENT[spliceosome_df$impact=="ORF"]
pladb_events_shared<-filter(pladb_events$PSI, EVENT %in% event_names)

# Add GENE names to the EVENT column
pladb_events_shared <- pladb_events_shared %>%
  mutate(EVENT = paste0(EVENT, " (", GENE, ")"))

# Reshape data for plotting
pladb_events_shared <- pladb_events_shared %>%
  pivot_longer(
    cols = 7:15,
    names_to = "Condition",
    values_to = "Value"
  ) %>%
  mutate(
    Group = case_when(
      grepl("038|039|040", Condition) ~ "Control",
      grepl("041|042|043", Condition) ~ "Pladb1mM",
      grepl("044|045|046", Condition) ~ "Pladb10mM"
    )
  )

pladb_events_shared$Group<-factor(pladb_events_shared$Group, levels = c("Control","Pladb1mM","Pladb10mM"))
# Plot with updated EVENT names
ggplot(pladb_events_shared, aes(x = Group, y = Value, color = EVENT)) +
  geom_boxplot(size = 1,show.legend = FALSE) +
  labs(
    title = "ORF Disrupted Events RNA Related",
    x = "",
    y = "PSI"
  ) +
  theme_minimal() +
  facet_wrap(~EVENT)


