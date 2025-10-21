# ---- USER SETTINGS (edit only these) ---------------------------------------
samples_dir <- "./Preprocessing/gene_expression/star_out"
sample_names <- c("reads_10_ctl_13_03_25_merged", "reads_11_ctl_13_03_25_merged", "reads_12_ctl_13_03_25_merged", 
                  "reads_13_pladb_3mm_13_03_25_merged","reads_14_pladb_3mm_13_03_25_merged","reads_15_pladb_3mm_13_03_25_merged",
                  "reads_1_ctl_11_03_25_merged","reads_2_ctl_11_03_25_merged","reads_3_ctl_11_03_25_merged",
                  "reads_4_tub_10mm_11_03_25_merged","reads_5_tub_10mm_11_03_25_merged","reads_6_tub_10mm_11_03_25_merged",
                  "reads_7_ssa_1mm_11_03_25_merged","reads_8_ssa_1mm_11_03_25_merged","reads_9_ssa_1mm_11_03_25_merged"
)
condition <- factor(rep(c("control_pladb","pladb","control","tub","ssa"),each=3))
stranded_column <- 2   # STAR unstranded counts
output_dir <- "./Preprocessing/gene_expression"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- dependencies -----------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2","biomaRt"), ask = FALSE, update = FALSE)
library(DESeq2)
library(biomaRt)

# ---- helpers ----------------------------------------------------------------
read_star_counts <- function(fp, col = 2) {
  df <- read.delim(fp, header = FALSE, stringsAsFactors = FALSE)
  counts <- as.integer(df[[col]])
  names(counts) <- df[[1]]
  counts
}

# Read counts
counts_list <- lapply(sample_names, function(s) {
  fp <- file.path(samples_dir, s, "ReadsPerGene.out.tab")
  if (!file.exists(fp)) stop("Missing: ", fp)
  read_star_counts(fp, col = stranded_column)
})
names(counts_list) <- sample_names

common_genes <- Reduce(intersect, lapply(counts_list, names))
message("Genes retained after intersecting all samples: ", length(common_genes))
count_mat <- sapply(counts_list, function(x) as.integer(x[common_genes]))
rownames(count_mat) <- common_genes
colnames(count_mat) <- sample_names

# Keep only Ensembl genes
keep_genes <- grepl("^ENSMUSG", rownames(count_mat))
count_mat <- count_mat[keep_genes, , drop = FALSE]
message("Genes retained after keeping ENSMUSG only: ", nrow(count_mat))

# ---- DESeq2 -----------------------------------------------------------------
coldata <- data.frame(row.names = sample_names, condition = condition)
dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = ~ condition)

# Filter low counts: keep genes with >=10 reads total
dds <- dds[rowSums(counts(dds)) >= 10, ]
message("Genes retained after filtering low counts: ", nrow(dds))

dds <- DESeq(dds)

# Normalized counts and VST
norm_counts <- counts(dds, normalized = TRUE)
vst_mat <- vst(dds, blind = FALSE)  # blind=FALSE to account for design
write.csv(as.data.frame(norm_counts), file.path(output_dir, "normalized_counts.csv"))
write.csv(as.data.frame(assay(vst_mat)), file.path(output_dir, "vst_transformed.csv"))

# ---- annotation (biomaRt) --------------------------------------------------
mart <- NULL
try({ mart <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl") }, silent = TRUE)
if (is.null(mart)) mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

attrs <- c("ensembl_gene_id","mgi_symbol","external_gene_name")
gene_ids_versioned <- rownames(dds)
gene_ids <- sub("\\..*$", "", gene_ids_versioned)

ann_df <- getBM(attributes = attrs, filters = "ensembl_gene_id", values = unique(gene_ids), mart = mart)
# Safe symbol assignment
ann_df$symbol <- ifelse(!is.na(ann_df$mgi_symbol) & ann_df$mgi_symbol != "", ann_df$mgi_symbol, ann_df$external_gene_name)
ann_df <- ann_df[!duplicated(ann_df$ensembl_gene_id), c("ensembl_gene_id","symbol")]

# ---- DESeq2 contrasts -------------------------------------------------------
run_contrast <- function(dds, ann_df, cond1, cond2, outfile) {
  res <- results(dds, contrast = c("condition", cond1, cond2))
  res <- lfcShrink(dds, contrast = c("condition", cond1, cond2), type = "ashr")
  res_df <- as.data.frame(res)
  res_df$ensembl_gene_id_version <- rownames(res_df)
  res_df$ensembl_gene_id <- sub("\\..*$", "", rownames(res_df))
  # Merge annotation safely
  res_df <- merge(res_df, ann_df, by = "ensembl_gene_id", all.x = TRUE)
  write.csv(res_df, file.path(output_dir, outfile), row.names = FALSE)
  return(res_df)
}

pladb_deseq <- run_contrast(dds, ann_df, "control_pladb", "pladb", "pladb_deseq.csv")
tub_deseq   <- run_contrast(dds, ann_df, "control", "tub", "tub_deseq.csv")
ssa_deseq   <- run_contrast(dds, ann_df, "control", "ssa", "ssa_deseq.csv")
