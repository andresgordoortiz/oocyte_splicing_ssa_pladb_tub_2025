library(DESeq2)
library(biomaRt)
library(apeglm)
library(openxlsx) # New dependency for Excel output

# ===== SETTINGS =====
# (Keeping your original settings)
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

# Create directory if it doesn't exist
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

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

# ===== PRE-FILTER =====
keep_expressed <- rowSums(count_mat >= 1) >= 3
count_mat <- count_mat[keep_expressed, , drop = FALSE]

# ===== DESEQ2 =====
coldata <- data.frame(row.names = sample_names, condition = condition)
dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, 
                              design = ~ condition)
dds <- DESeq(dds)
norm_counts <- counts(dds, normalized = TRUE)

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

# ===== UPDATED CONTRAST FUNCTION =====
run_contrast_tables <- function(dds, norm_counts, ann_df, cond1, cond2, basemean_cutoff = 10) {
  message(sprintf("\n=== Processing %s vs %s ===", cond1, cond2))
  
  idx1 <- which(colData(dds)$condition == cond1)
  idx2 <- which(colData(dds)$condition == cond2)
  
  # 1. Get RAW results
  res_raw_obj <- results(dds, contrast = c("condition", cond1, cond2), alpha = alpha)
  
  # 2. Get SHRUNKEN results
  coef_name <- resultsNames(dds)[grep(paste0("condition_", cond1), resultsNames(dds))]
  res_shrunk_obj <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)
  
  # Helper to format dataframes
  format_df <- function(res_obj) {
    df <- as.data.frame(res_obj)
    df$ensembl_gene_id_version <- rownames(df)
    df$ensembl_gene_id <- sub("\\..*$", "", rownames(df))
    df <- merge(df, ann_df, by = "ensembl_gene_id", all.x = TRUE)
    df$mean_treatment <- rowMeans(norm_counts[df$ensembl_gene_id_version, idx1, drop=FALSE])
    df$mean_control <- rowMeans(norm_counts[df$ensembl_gene_id_version, idx2, drop=FALSE])
    df <- df[df$baseMean >= basemean_cutoff, ]
    df <- df[order(df$padj, na.last = TRUE), ]
    return(df)
  }
  
  return(list(
    shrunken = format_df(res_shrunk_obj),
    raw = format_df(res_raw_obj)
  ))
}

# ===== GENERATE DATA =====
comparisons <- list(
  list(c1="pladb", c2="control_pladb", name="PlaDB"),
  list(c1="tub",   c2="control",       name="Tub"),
  list(c1="ssa",   c2="control",       name="SSA")
)

# Initialize Excel Workbook
wb <- createWorkbook()

for (comp in comparisons) {
  results_list <- run_contrast_tables(dds, norm_counts, ann_df, comp$c1, comp$c2)
  
  # Create Tab Names (Excel tabs max 31 chars)
  shrunk_name <- paste0(comp$name, "_Shrunken")
  raw_name <- paste0(comp$name, "_Raw")
  
  # Add Shrunken Sheet
  addWorksheet(wb, shrunk_name)
  writeData(wb, shrunk_name, results_list$shrunken)
  
  # Add Raw Sheet
  addWorksheet(wb, raw_name)
  writeData(wb, raw_name, results_list$raw)
}

# ===== SAVE OUTPUT =====
excel_file <- file.path(output_dir, "Supplementary1_DESeq2_Results.xlsx")
saveWorkbook(wb, excel_file, overwrite = TRUE)

message(sprintf("\nAnalysis Complete. Results saved to: %s", excel_file))