library(DESeq2)
library(biomaRt)
library(apeglm)  # Required for proper shrinkage

# ===== SETTINGS =====
samples_dir <- "./Preprocessing/gene_expression/star_out"
sample_names <- c("reads_10_ctl_13_03_25_merged", "reads_11_ctl_13_03_25_merged", 
                  "reads_12_ctl_13_03_25_merged", "reads_13_pladb_3mm_13_03_25_merged",
                  "reads_14_pladb_3mm_13_03_25_merged","reads_15_pladb_3mm_13_03_25_merged",
                  "reads_1_ctl_11_03_25_merged","reads_2_ctl_11_03_25_merged",
                  "reads_3_ctl_11_03_25_merged","reads_4_tub_10mm_11_03_25_merged",
                  "reads_5_tub_10mm_11_03_25_merged","reads_6_tub_10mm_11_03_25_merged",
                  "reads_7_ssa_1mm_11_03_25_merged","reads_8_ssa_1mm_11_03_25_merged",
                  "reads_9_ssa_1mm_11_03_25_merged")

condition <- factor(c(rep("control_pladb",3), rep("pladb",3), rep("control",3), 
                      rep("tub",3), rep("ssa",3)))

stranded_column <- 2
output_dir <- "./Preprocessing/gene_expression"
alpha <- 0.05

# ===== READ COUNTS =====
read_star_counts <- function(fp, col = 2) {
  df <- read.delim(fp, header = FALSE, stringsAsFactors = FALSE)
  counts <- as.integer(df[[col]])
  names(counts) <- df[[1]]
  counts
}

counts_list <- lapply(sample_names, function(s) {
  fp <- file.path(samples_dir, s, "ReadsPerGene.out.tab")
  read_star_counts(fp, col = stranded_column)
})
names(counts_list) <- sample_names

common_genes <- Reduce(intersect, lapply(counts_list, names))
count_mat <- sapply(counts_list, function(x) as.integer(x[common_genes]))
rownames(count_mat) <- common_genes
colnames(count_mat) <- sample_names

count_mat <- count_mat[grepl("^ENSMUSG", rownames(count_mat)), , drop = FALSE]
message(sprintf("Initial genes: %d", nrow(count_mat)))

# ===== MINIMAL PRE-FILTER =====
keep_expressed <- rowSums(count_mat >= 1) >= 3
count_mat <- count_mat[keep_expressed, , drop = FALSE]
message(sprintf("After minimal filter: %d genes", nrow(count_mat)))

# ===== DESEQ2 =====
coldata <- data.frame(row.names = sample_names, condition = condition)
dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, 
                              design = ~ condition)
dds <- DESeq(dds)
norm_counts <- counts(dds, normalized = TRUE)

# Save transformed data
vst_mat <- vst(dds, blind = FALSE)
write.csv(as.data.frame(norm_counts), file.path(output_dir, "normalized_counts.csv"))
write.csv(as.data.frame(assay(vst_mat)), file.path(output_dir, "vst_transformed.csv"))

# ===== ANNOTATION =====
mart <- tryCatch(
  useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl"),
  error = function(e) useMart("ensembl", dataset = "mmusculus_gene_ensembl")
)

gene_ids <- sub("\\..*$", "", rownames(dds))
ann_df <- getBM(attributes = c("ensembl_gene_id","mgi_symbol","external_gene_name"),
                filters = "ensembl_gene_id", values = unique(gene_ids), mart = mart)
ann_df$symbol <- ifelse(!is.na(ann_df$mgi_symbol) & ann_df$mgi_symbol != "", 
                        ann_df$mgi_symbol, ann_df$external_gene_name)
ann_df <- ann_df[!duplicated(ann_df$ensembl_gene_id), ]

# ===== CONTRAST FUNCTION =====
run_contrast <- function(dds, norm_counts, ann_df, cond1, cond2, outfile, 
                         basemean_cutoff = 10) {
  message(sprintf("\n=== %s vs %s ===", cond1, cond2))
  
  idx1 <- which(colData(dds)$condition == cond1)
  idx2 <- which(colData(dds)$condition == cond2)
  
  # Get unshrunken results
  res_raw <- results(dds, contrast = c("condition", cond1, cond2), alpha = alpha)
  
  # Apply apeglm shrinkage (best method)
  coef_name <- resultsNames(dds)[grep(paste0("condition_", cond1), resultsNames(dds))]
  res <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)
  
  # Build results table
  res_df <- as.data.frame(res)
  res_df$ensembl_gene_id_version <- rownames(res_df)
  res_df$ensembl_gene_id <- sub("\\..*$", "", rownames(res_df))
  res_df <- merge(res_df, ann_df, by = "ensembl_gene_id", all.x = TRUE)
  
  # Add raw (unshrunken) LFC for comparison
  res_df$log2FC_raw <- res_raw$log2FoldChange[match(res_df$ensembl_gene_id_version, 
                                                    rownames(res_raw))]
  
  # Add mean counts
  res_df$mean_treatment <- rowMeans(norm_counts[res_df$ensembl_gene_id_version, idx1, drop=FALSE])
  res_df$mean_control <- rowMeans(norm_counts[res_df$ensembl_gene_id_version, idx2, drop=FALSE])
  
  # Filter by baseMean (removes noisy low-count genes)
  res_df <- res_df[res_df$baseMean >= basemean_cutoff, ]
  res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
  
  # Summary
  n_sig <- sum(res_df$padj < alpha, na.rm = TRUE)
  n_up <- sum(res_df$padj < alpha & res_df$log2FoldChange > 0, na.rm = TRUE)
  n_down <- sum(res_df$padj < alpha & res_df$log2FoldChange < 0, na.rm = TRUE)
  
  message(sprintf("Significant (padj < %.2f, baseMean >= %d): %d", 
                  alpha, basemean_cutoff, n_sig))
  message(sprintf("  Upregulated: %d | Downregulated: %d", n_up, n_down))
  
  if (n_sig > 0) {
    sig_genes <- res_df[!is.na(res_df$padj) & res_df$padj < alpha, ]
    message(sprintf("  Shrunken log2FC range: [%.2f, %.2f]", 
                    min(sig_genes$log2FoldChange), max(sig_genes$log2FoldChange)))
    
    # Show top genes with both shrunken and raw LFC
    message("\n  Top 5 genes (shrunken LFC | raw LFC):")
    top5 <- head(sig_genes, 5)
    for (i in 1:nrow(top5)) {
      sym <- ifelse(is.na(top5$symbol[i]), top5$ensembl_gene_id[i], top5$symbol[i])
      message(sprintf("    %s | LFC: %.2f (raw: %.2f) | padj: %.2e | counts: %.1f vs %.1f",
                      sym, top5$log2FoldChange[i], top5$log2FC_raw[i], top5$padj[i],
                      top5$mean_treatment[i], top5$mean_control[i]))
    }
  }
  
  write.csv(res_df, file.path(output_dir, outfile), row.names = FALSE)
  message(sprintf("Saved: %s", outfile))
  
  return(res_df)
}

# ===== DIAGNOSTIC: CHECK EXTREME LFC GENES =====
check_extreme_lfc <- function(dds, norm_counts, condition_name, control_name) {
  message(sprintf("\n=== DIAGNOSTIC: Extreme LFC genes (%s vs %s) ===", 
                  condition_name, control_name))
  
  res_raw <- results(dds, contrast = c("condition", condition_name, control_name))
  extreme <- res_raw[!is.na(res_raw$log2FoldChange) & abs(res_raw$log2FoldChange) > 15, ]
  
  if (nrow(extreme) == 0) {
    message("No genes with |log2FC| > 15")
    return(NULL)
  }
  
  message(sprintf("Found %d genes with |log2FC| > 15\n", nrow(extreme)))
  
  idx1 <- which(colData(dds)$condition == condition_name)
  idx2 <- which(colData(dds)$condition == control_name)
  
  for (i in 1:min(10, nrow(extreme))) {
    gene_id <- rownames(extreme)[i]
    counts_treat <- norm_counts[gene_id, idx1]
    counts_ctrl <- norm_counts[gene_id, idx2]
    
    message(sprintf("Gene %d: %s | log2FC = %.2f", i, gene_id, extreme$log2FoldChange[i]))
    message(sprintf("  Treatment counts: %s (mean: %.2f)", 
                    paste(round(counts_treat, 1), collapse=", "), mean(counts_treat)))
    message(sprintf("  Control counts: %s (mean: %.2f)", 
                    paste(round(counts_ctrl, 1), collapse=", "), mean(counts_ctrl)))
    message(sprintf("  baseMean: %.2f | padj: %.2e\n", 
                    extreme$baseMean[i], extreme$padj[i]))
  }
  
  return(extreme)
}

# ===== RUN DIAGNOSTICS FIRST =====
message("\n========== CHECKING FOR EXTREME LFC (BEFORE SHRINKAGE) ==========")
check_extreme_lfc(dds, norm_counts, "pladb", "control_pladb")
check_extreme_lfc(dds, norm_counts, "tub", "control")
check_extreme_lfc(dds, norm_counts, "ssa", "control")

# ===== RUN CONTRASTS WITH SHRINKAGE =====
message("\n========== RUNNING CONTRASTS WITH APEGLM SHRINKAGE ==========")
pladb_deseq <- run_contrast(dds, norm_counts, ann_df, "pladb", "control_pladb", 
                            "pladb_deseq_shrunken.csv")
tub_deseq <- run_contrast(dds, norm_counts, ann_df, "tub", "control", 
                          "tub_deseq_shrunken.csv")
ssa_deseq <- run_contrast(dds, norm_counts, ann_df, "ssa", "control", 
                          "ssa_deseq_shrunken.csv")

message("\n=== COMPLETE ===")
message("Files contain both shrunken and raw log2FC for comparison")