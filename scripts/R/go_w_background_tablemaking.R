# Gene Ontology BG creation
spire_bg<-data.frame(bg=unique(c(
  spire_pdiff_exons[(spire_pdiff_exons$FDR > 0.05 | abs(spire_pdiff_exons$deltapsi) < 0.1), "GENE"],
  spire_pdiff_introns[(spire_pdiff_introns$FDR > 0.05 | abs(spire_pdiff_introns$deltapsi) < 0.1), "GENE"],
  spire_pdiff_alt[(spire_pdiff_alt$FDR > 0.05 | abs(spire_pdiff_alt$deltapsi) < 0.1), "GENE"]
)))
tmp<-spire_bg[!spire_bg$bg %in% c(differential_spire_exons$GENE,differential_spire_introns$GENE,differential_spire_alt$GENE),]
write.csv(tmp, "bg_spire.csv", row.names = FALSE)

fmn2_bg<-data.frame(bg=unique(c(
  fmndko_pdiff_exons[(fmndko_pdiff_exons$FDR > 0.05 | abs(fmndko_pdiff_exons$deltapsi) < 0.1), "GENE"],
  fmndko_pdiff_introns[(fmndko_pdiff_introns$FDR > 0.05 | abs(fmndko_pdiff_introns$deltapsi) < 0.1), "GENE"],
  fmndko_pdiff_alt[(fmndko_pdiff_alt$FDR > 0.05 | abs(fmndko_pdiff_alt$deltapsi) < 0.1), "GENE"]
)))
tmp<-fmn2_bg[!fmn2_bg$bg %in% c(differential_fmndko_exons$GENE,differential_fmndko_introns$GENE,differential_fmndko_alt$GENE),]
write.csv(tmp, "bg_fmndko.csv", row.names = FALSE)

write.csv(unique(c(differential_fmndko_exons$GENE,differential_fmndko_introns$GENE,differential_fmndko_alt$GENE)), "diff_fmndko.csv", row.names = FALSE)
write.csv(unique(c(differential_spire_exons$GENE,differential_spire_introns$GENE,differential_spire_alt$GENE)), "diff_spire.csv", row.names = FALSE)

library(readxl)
library(dplyr)

# List of Excel files
files <- c("notebooks/tables_10n/pladb_pdiff_alt.csv",
           "notebooks/tables_10n/pladb_pdiff_exons.csv",
           "notebooks/tables_10n/pladb_pdiff_introns.csv")


# Initialize vectors to store unique gene names
significant_genes <- c()
nonsignificant_genes <- c()

# Loop through each file
for (file in files) {
  df <- read.csv(file)  # Read CSV file
  
  # Check if required columns exist
  if (!all(c("GENE", "FDR", "deltapsi") %in% colnames(df))) {
    stop(paste("Missing required columns in file:", file))
  }
  
  # Get all genes from the file
  all_genes_file <- unique(df$GENE)
  
  # Filter significant genes (FDR <= 0.05 and DeltaPSI >= 0.1) for the file
  sig_genes_file <- df %>%
    filter(FDR <= 0.05 & deltapsi >= 0.1) %>%
    pull(GENE)
  
  # Define non-significant genes for this file as those that are in the file but not in the file's significant list
  nonsig_genes_file <- setdiff(all_genes_file, sig_genes_file)
  
  # Update the global lists with a union to keep unique gene names
  significant_genes <- union(significant_genes, sig_genes_file)
  nonsignificant_genes <- union(nonsignificant_genes, nonsig_genes_file)
}

# Now remove from the nonsignificant list any gene that ever appeared as significant
nonsignificant_genes <- setdiff(nonsignificant_genes, significant_genes)

# Optionally, save the results to CSV files
write.csv(significant_genes, "significant_genes.csv", row.names = FALSE, quote = FALSE)
write.csv(nonsignificant_genes, "nonsignificant_genes.csv", row.names = FALSE, quote = FALSE)


pladb_bg<-data.frame(bg=unique(c(
  pladb_pdiff_exons[(fmndko_pdiff_exons$FDR > 0.05 | abs(fmndko_pdiff_exons$deltapsi) < 0.1), "GENE"],
  fmndko_pdiff_introns[(fmndko_pdiff_introns$FDR > 0.05 | abs(fmndko_pdiff_introns$deltapsi) < 0.1), "GENE"],
  fmndko_pdiff_alt[(fmndko_pdiff_alt$FDR > 0.05 | abs(fmndko_pdiff_alt$deltapsi) < 0.1), "GENE"]
)))
tmp<-fmn2_bg[!fmn2_bg$bg %in% c(differential_fmndko_exons$GENE,differential_fmndko_introns$GENE,differential_fmndko_alt$GENE),]
write.csv(tmp, "bg_fmndko.csv", row.names = FALSE)

write.csv(unique(c(differential_fmndko_exons$GENE,differential_fmndko_introns$GENE,differential_fmndko_alt$GENE)), "diff_fmndko.csv", row.names = FALSE)
write.csv(unique(c(differential_spire_exons$GENE,differential_spire_introns$GENE,differential_spire_alt$GENE)), "diff_spire.csv", row.names = FALSE)

# GO analysis

library(clusterProfiler)
library(org.Mm.eg.db) 

# Spire

gene_ids <- bitr(unique(c(differential_spire_exons$GENE,differential_spire_introns$GENE,differential_spire_alt$GENE)), fromType = "SYMBOL", 
                 toType = "ENTREZID", OrgDb = org.Mm.eg.db)
background_ids <- bitr(spire_bg$bg, fromType = "SYMBOL", 
                       toType = "ENTREZID", OrgDb = org.Mm.eg.db)

go_bp <- enrichGO(gene = gene_ids$ENTREZID, 
                  universe = background_ids$ENTREZID,
                  OrgDb = org.Mm.eg.db, 
                  ont = "BP", 
                  pAdjustMethod = "BH", 
                  pvalueCutoff = 0.05)

go_mf <- enrichGO(gene = gene_ids$ENTREZID, 
                  universe = background_ids$ENTREZID,
                  OrgDb = org.Mm.eg.db, 
                  ont = "MF", 
                  pAdjustMethod = "BH", 
                  pvalueCutoff = 0.05)

go_cc <- enrichGO(gene = gene_ids$ENTREZID, 
                  universe = background_ids$ENTREZID,
                  OrgDb = org.Mm.eg.db, 
                  ont = "CC", 
                  pAdjustMethod = "BH", 
                  pvalueCutoff = 0.05)

library(ggplot2)
dotplot(go_bp, showCategory = 10, title = "Biological Process")
dotplot(go_mf, showCategory = 10, title = "Molecular Function")
dotplot(go_cc, showCategory = 10, title = "Cellular Component")


# FMN2

gene_ids <- bitr(unique(c(differential_fmndko_exons$GENE,differential_fmndko_introns$GENE,differential_fmndko_alt$GENE)), fromType = "SYMBOL", 
                 toType = "ENTREZID", OrgDb = org.Mm.eg.db)
background_ids <- bitr(fmn2_bg$bg, fromType = "SYMBOL", 
                       toType = "ENTREZID", OrgDb = org.Mm.eg.db)

go_bp <- enrichGO(gene = gene_ids$ENTREZID, 
                  universe = background_ids$ENTREZID,
                  OrgDb = org.Mm.eg.db, 
                  ont = "BP", 
                  pAdjustMethod = "BH", 
                  pvalueCutoff = 0.05)

go_mf <- enrichGO(gene = gene_ids$ENTREZID, 
                  universe = background_ids$ENTREZID,
                  OrgDb = org.Mm.eg.db, 
                  ont = "MF", 
                  pAdjustMethod = "BH", 
                  pvalueCutoff = 0.05)

go_cc <- enrichGO(gene = gene_ids$ENTREZID, 
                  universe = background_ids$ENTREZID,
                  OrgDb = org.Mm.eg.db, 
                  ont = "CC", 
                  pAdjustMethod = "BH", 
                  pvalueCutoff = 0.05)

library(ggplot2)
dotplot(go_bp, showCategory = 10, title = "Biological Process")
dotplot(go_mf, showCategory = 10, title = "Molecular Function")
dotplot(go_cc, showCategory = 10, title = "Cellular Component")
