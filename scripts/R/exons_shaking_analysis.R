vastdb_exons<-arrange(vastdb_exons, GENE, EVENT)

spire_constitutive <- spire_exons$PSI %>%
  filter(
    (if_all(7:9, ~ . <= 5) | if_all(7:9, ~ . >= 95)) |  # Columns 7-9 all 0 or all 100
      (if_all(10:12, ~ . <= 5) | if_all(10:12, ~ . >= 95)) # Columns 10-12 all 0 or all 100
  )


vastdb_exons$constitutive<-ifelse(vastdb_exons$EVENT %in% spire_constitutive$EVENT, yes = T, no = F)
vastdb_exons_model <- vastdb_exons %>%
  arrange(GENE, EVENT) %>%
  group_by(GENE) %>%
  summarise(GENE=GENE,pattern = paste(constitutive, collapse = ",")) %>%
  unique()


final_alt_exons_gene_model<-na.omit(spire_exons$PSI)
final_alt_exons_gene_model <- final_alt_exons_gene_model %>%
  left_join(vastdb_exons_model, by = "GENE")


# Compute the mean of columns 10 to 12 for each row and the mean of columns 7 to 9 for each row
mean_10_12 <- rowMeans(final_alt_exons_gene_model[, 10:12], na.rm = TRUE)
mean_7_9  <- rowMeans(final_alt_exons_gene_model[, 7:9], na.rm = TRUE)

# Compute the difference: (mean of columns 10:12) - (mean of columns 7:9)
final_alt_exons_gene_model$mean_diff <- mean_10_12 - mean_7_9

final_alt_exons_gene_model$included_control<-ifelse(final_alt_exons_gene_model$mean_diff < 0, "Skipped", "Included")

final_alt_exons_gene_model$constitutive<-ifelse(final_alt_exons_gene_model$EVENT %in% spire_constitutive$EVENT, yes = T, no = F)


plot_psi_distribution <- ggplot(
  final_alt_exons_gene_model[abs(final_alt_exons_gene_model$mean_diff)>=1,], 
  aes(y=..density..,x = mean_diff, fill = ifelse(mean_diff < 0, "Skipped", "Included"))
) +
  geom_histogram(bins = 100, position = "identity") +
  labs(
    title = "PSI Difference Distribution Grouped by Number of Exons present in a Gene",
    subtitle = "Spire, dPSI = KO - Control",
    x = "ΔPSI (DKO - WT)",
    y = "# of Exons",
    fill = "|ΔPSI| >= 10",
    caption = paste0("Created by AG on ", Sys.Date())
  ) +
  theme_minimal(base_family = font) +
  theme(
    legend.position = "right",
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.text.x = element_text(size = 15),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "gray80", size = 0.5),
    panel.grid.minor = element_line(color = "gray95", size = 0.3),
    strip.text = element_text(size = 14, face = "bold"),
    panel.spacing = unit(1.5, "lines")  # Increase the distance between facets
  ) +
  scale_fill_brewer(palette = "Set1")

plot_psi_distribution

final_alt_exons_gene_model$pattern[final_alt_exons_gene_model$mean_diff<(-10)]