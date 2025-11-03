library(dplyr)
library(stringr)
library(readr)

# Helper functions ------------------------------------------------------------
extract_genes_from_gtf <- function(gtf_file, keep_version = TRUE) {
  cat("Extracting gene symbols from GTF:", gtf_file, "\n")
  if (!file.exists(gtf_file)) stop("GTF file not found: ", gtf_file)
  
  gtf_lines <- readLines(gtf_file)
  gtf_lines <- gtf_lines[!grepl("^#", gtf_lines)]
  
  gene_name_to_id <- list()
  
  for (i in seq_along(gtf_lines)) {
    line <- gtf_lines[i]
    if (grepl("\\tgene\\t", line)) {
      gene_id_match <- str_extract(line, 'gene_id[= ]"?([^";\\s]+)"?')
      gene_name_match <- str_extract(line, 'gene_name[= ]"?([^";\\s]+)"?')
      
      if (!is.na(gene_id_match) && !is.na(gene_name_match)) {
        gene_id <- str_replace(gene_id_match, 'gene_id[= ]"?([^";\\s]+)"?', '\\1')
        gene_name <- str_replace(gene_name_match, 'gene_name[= ]"?([^";\\s]+)"?', '\\1')
        
        if (!gene_name %in% names(gene_name_to_id)) {
          gene_name_to_id[[gene_name]] <- if (keep_version) gene_id else sub("\\.\\d+$", "", gene_id)
        }
      }
    }
    if (i %% 50000 == 0) cat("  Processed", i, "lines...\n")
  }
  
  data.frame(
    mgi_symbol = names(gene_name_to_id),
    ensembl_gene_id = unlist(gene_name_to_id, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

build_rmaps_table <- function(joined_df, dataset_label, gene_map) {
  if (nrow(joined_df) == 0) return(tibble())
  
  gene_col <- if ("GENE" %in% names(joined_df)) "GENE" else "GENE.x"
  
  # Parse CO_A: "chr16:18795825-18795944"
  coords <- str_match(joined_df$CO_A, "^(.+):(\\d+)-(\\d+)$")
  
  # Parse REF_CO for strand: "chr16:18797349,18795825-18795944,18795226:-"
  strand <- str_extract(joined_df$REF_CO, "[+-]$")
  
  joined_df %>%
    mutate(
      SCAFFOLD = coords[, 2],
      START = as.integer(coords[, 3]),
      END = as.integer(coords[, 4]),
      STRAND = strand,
      DATASET = dataset_label,
      GENE = .data[[gene_col]]
    ) %>%
    select(GENE, START, END, SCAFFOLD, STRAND, DATASET) %>%
    left_join(gene_map, by = c("GENE" = "mgi_symbol")) %>%
    rename(GENEID_ENSEMBL = ensembl_gene_id) %>%
    select(START, END, SCAFFOLD, STRAND, GENEID_ENSEMBL, DATASET) %>%
    filter(!is.na(SCAFFOLD) & !is.na(START) & !is.na(END) & !is.na(STRAND) & !is.na(GENEID_ENSEMBL))
}

# Load data -------------------------------------------------------------------
cat("Loading data...\n")

if (!exists("tub_fdr_df")) {
  if (file.exists("tub_fdr.csv")) tub_fdr_df <- read_csv("tub_fdr.csv", show_col_types = FALSE)[, -1]
  else stop("tub_fdr_df not found")
}
if (!exists("pladb_fdr_df")) {
  if (file.exists("pladb_fdr.csv")) pladb_fdr_df <- read_csv("pladb_fdr.csv", show_col_types = FALSE)[, -1]
  else stop("pladb_fdr_df not found")
}
if (!exists("ssa_fdr_df")) {
  if (file.exists("ssa_fdr.csv")) ssa_fdr_df <- read_csv("ssa_fdr.csv", show_col_types = FALSE)[, -1]
  else stop("ssa_fdr_df not found")
}
if (!exists("event_info")) {
  if (file.exists("EVENT_INFO-mm10.tab")) event_info <- read.delim("EVENT_INFO-mm10.tab")
  else stop("event_info not found")
}

# Get gene mapping ------------------------------------------------------------
if (file.exists("mm10_gene_mapping.csv")) {
  cat("Loading existing gene mapping from mm10_gene_mapping.csv\n")
  gene_map <- read_csv("mm10_gene_mapping.csv", show_col_types = FALSE)
  cat("Loaded", nrow(gene_map), "gene mappings\n\n")
} else {
  cat("Gene mapping file not found. Extracting from GTF...\n")
  gtf_file <- NULL
  possible_gtf <- c(
    "reference_genome/gencode.vM10.annotation.gtf",
    "gencode.vM10.annotation.gtf",
    "mm10.gtf"
  )
  for (path in possible_gtf) {
    if (file.exists(path)) {
      gtf_file <- path
      break
    }
  }
  if (is.null(gtf_file)) stop("No GTF file found")
  
  gene_map <- extract_genes_from_gtf(gtf_file, keep_version = TRUE)
  write_csv(gene_map, "mm10_gene_mapping.csv")
  cat("Saved gene mapping to mm10_gene_mapping.csv\n")
  cat("Gene mapping:", nrow(gene_map), "genes\n\n")
}

# Process each treatment ------------------------------------------------------
treatments <- list(tub = tub_fdr_df, pladb = pladb_fdr_df, ssa = ssa_fdr_df)

for (trname in names(treatments)) {
  cat("Processing:", trname, "\n")
  df <- treatments[[trname]]
  
  # Separate exon and intron events
  exon_events <- df %>% filter(grepl("EX", EVENT))
  intron_events <- df %>% filter(grepl("INT", EVENT))
  
  cat("  Total - Exons:", nrow(exon_events), "| Introns:", nrow(intron_events), "\n")
  
  # Define groups for EXONS
  exon_up <- exon_events %>% filter(!is.na(deltapsi) & deltapsi >= 0.1 & !is.na(FDR) & FDR <= 0.05)
  exon_down <- exon_events %>% filter(!is.na(deltapsi) & deltapsi <= -0.1 & !is.na(FDR) & FDR <= 0.05)
  exon_nonsig <- exon_events %>% filter(is.na(FDR) | FDR > 0.05 | is.na(deltapsi) | abs(deltapsi) < 0.1)
  
  # Define groups for INTRONS
  intron_up <- intron_events %>% filter(!is.na(deltapsi) & deltapsi >= 0.1 & !is.na(FDR) & FDR <= 0.05)
  intron_down <- intron_events %>% filter(!is.na(deltapsi) & deltapsi <= -0.1 & !is.na(FDR) & FDR <= 0.05)
  intron_nonsig <- intron_events %>% filter(is.na(FDR) | FDR > 0.05 | is.na(deltapsi) | abs(deltapsi) < 0.1)
  
  # Sample background if needed
  max_bg <- 10000
  if (nrow(exon_nonsig) > max_bg) {
    set.seed(42)
    exon_nonsig <- exon_nonsig %>% sample_n(max_bg)
  }
  if (nrow(intron_nonsig) > max_bg) {
    set.seed(42)
    intron_nonsig <- intron_nonsig %>% sample_n(max_bg)
  }
  
  cat("  Exons - up:", nrow(exon_up), "| down:", nrow(exon_down), "| nonsig:", nrow(exon_nonsig), "\n")
  cat("  Introns - up:", nrow(intron_up), "| down:", nrow(intron_down), "| nonsig:", nrow(intron_nonsig), "\n")
  
  # Join with event_info to get coordinates
  exon_up_info <- exon_up %>% inner_join(event_info, by = "EVENT")
  exon_down_info <- exon_down %>% inner_join(event_info, by = "EVENT")
  exon_nonsig_info <- exon_nonsig %>% inner_join(event_info, by = "EVENT")
  
  intron_up_info <- intron_up %>% inner_join(event_info, by = "EVENT")
  intron_down_info <- intron_down %>% inner_join(event_info, by = "EVENT")
  intron_nonsig_info <- intron_nonsig %>% inner_join(event_info, by = "EVENT")
  
  # Build tables using CO_A coordinates and REF_CO strand
  exon_table <- bind_rows(
    build_rmaps_table(exon_up_info, "up", gene_map),
    build_rmaps_table(exon_down_info, "down", gene_map),
    build_rmaps_table(exon_nonsig_info, "ndiff", gene_map)
  )
  
  intron_table <- bind_rows(
    build_rmaps_table(intron_up_info, "up", gene_map),
    build_rmaps_table(intron_down_info, "down", gene_map),
    build_rmaps_table(intron_nonsig_info, "ndiff", gene_map)
  )
  
  # Write output
  if (nrow(exon_table) > 0) {
    write_tsv(exon_table, paste0(trname, "_exons.tab"))
    cat("  Wrote", nrow(exon_table), "exons to", paste0(trname, "_exons.tab\n"))
  }
  
  if (nrow(intron_table) > 0) {
    write_tsv(intron_table, paste0(trname, "_introns.tab"))
    cat("  Wrote", nrow(intron_table), "introns to", paste0(trname, "_introns.tab\n"))
  }
  cat("\n")
}

cat("Done! Files created with DATASET values: up, down, ndiff\n")