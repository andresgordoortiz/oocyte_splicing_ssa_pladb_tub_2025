# Export exons and introns for rMAPS RBP binding analysis (mm10 mouse genome)
# Creates .tab files with required columns for matt rna_maps_cisbp
# 
# KEY DIFFERENCE FROM tableformatt.R:
# - This script adds UPSTRM_EX_BORDER and DOSTRM_EX_BORDER columns needed for RBP analysis
# - Exon and intron events are processed SEPARATELY with different logic
# 
# For EXONS (EX events only - excludes ALTA/ALTD):
#   START, END, SCAFFOLD, STRAND, UPSTRM_EX_BORDER, DOSTRM_EX_BORDER, GENEID_ENSEMBL, GROUP
#   - START/END: alternative exon coordinates (from CO_A)
#   - UPSTRM_EX_BORDER: end coordinate of upstream constitutive exon (from CO_C1)
#   - DOSTRM_EX_BORDER: start coordinate of downstream constitutive exon (from CO_C2)
# 
# For INTRONS (INT events only - excludes ALTA/ALTD):
#   UPSTRM_EX_BORDER, START, END, DOSTRM_EX_BORDER, SCAFFOLD, STRAND, GENEID_ENSEMBL, GROUP
#   - START/END: intron coordinates (from C1 end + 1 to C2 start - 1)
#   - UPSTRM_EX_BORDER: end coordinate of upstream exon (C1)
#   - DOSTRM_EX_BORDER: start coordinate of downstream exon (C2)
# 
# For each treatment (tub, pladb, ssa), creates separate files for exons and introns
# with up (deltapsi >= 0.1), down (deltapsi <= -0.1), and nonsig (background) events

library(dplyr)
library(stringr)
library(readr)

# --- Helper functions -------------------------------------------------------
# Parse coordinate string like "chr16:18797349-18797399"
parse_coord_range <- function(x) {
  if (is.na(x) || x == "") return(list(chr = NA_character_, start = NA_integer_, end = NA_integer_))
  parts <- str_split_fixed(x, ":", 2)
  if (ncol(parts) < 2) return(list(chr = NA_character_, start = NA_integer_, end = NA_integer_))
  chr <- parts[1]
  range <- parts[2]
  if (!str_detect(range, "-")) return(list(chr = chr, start = as.integer(range), end = as.integer(range)))
  se <- str_split_fixed(range, "-", 2)
  start <- suppressWarnings(as.integer(se[1]))
  end   <- suppressWarnings(as.integer(se[2]))
  list(chr = chr, start = start, end = end)
}

# Parse strand from REF_CO
parse_strand <- function(refco) {
  if (is.na(refco) || refco == "") return(NA_character_)
  m <- str_match(refco, ":([+-])$")
  if (!is.na(m[1,2])) return(m[1,2])
  lastchar <- substr(refco, nchar(refco), nchar(refco))
  if (lastchar %in% c("+", "-")) return(lastchar)
  NA_character_}

# Extract gene symbols and IDs directly from GTF file
extract_genes_from_gtf <- function(gtf_file, keep_version = TRUE) {
  cat("Extracting gene symbols and IDs from GTF file:", gtf_file, "\n")
  
  if (!file.exists(gtf_file)) {
    stop("GTF file not found: ", gtf_file)
  }
  
  cat("Reading GTF file...\n")
  gtf_lines <- readLines(gtf_file)
  gtf_lines <- gtf_lines[!grepl("^#", gtf_lines)]
  
  cat("Processing", length(gtf_lines), "GTF lines...\n")
  
  gene_name_to_id <- list()
  gene_id_to_name <- list()
  gene_biotypes <- list()
  
  for (i in seq_along(gtf_lines)) {
    line <- gtf_lines[i]
    
    if (grepl("\\tgene\\t", line)) {
      gene_id_match <- str_extract(line, 'gene_id[= ]"?([^";\\s]+)"?')
      gene_name_match <- str_extract(line, 'gene_name[= ]"?([^";\\s]+)"?')
      gene_biotype_match <- str_extract(line, 'gene_type[= ]"?([^";\\s]+)"?|gene_biotype[= ]"?([^";\\s]+)"?')
      
      has_biotype <- !is.na(gene_biotype_match)
      
      if (!is.na(gene_id_match) && !is.na(gene_name_match)) {
        gene_id <- str_replace(gene_id_match, 'gene_id[= ]"?([^";\\s]+)"?', '\\1')
        gene_name <- str_replace(gene_name_match, 'gene_name[= ]"?([^";\\s]+)"?', '\\1')
        
        gene_id_versioned <- gene_id
        gene_id_unversioned <- sub("\\.\\d+$", "", gene_id)
        
        if (!gene_name %in% names(gene_name_to_id)) {
          if (keep_version) {
            gene_name_to_id[[gene_name]] <- gene_id_versioned
          } else {
            gene_name_to_id[[gene_name]] <- gene_id_unversioned
          }
          gene_biotypes[[gene_name]] <- has_biotype
        }
        
        gene_id_to_name[[gene_id_versioned]] <- gene_name
      }
    }
    
    if (i %% 50000 == 0) {
      cat("  Processed", i, "lines, found", length(gene_name_to_id), "genes so far...\n")
    }
  }
  
  mapping <- data.frame(
    mgi_symbol = names(gene_name_to_id),
    ensembl_gene_id = unlist(gene_name_to_id, use.names = FALSE),
    has_biotype = unlist(gene_biotypes, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  
  cat("Successfully extracted", nrow(mapping), "unique gene symbol -> Ensembl ID mappings\n")
  
  return(mapping[, c("mgi_symbol", "ensembl_gene_id")])
}

# Build rMAPS table for exons with flanking exon borders for RBP analysis
build_exon_rbp_table <- function(joined_df, group_label, gene_map) {
  if (nrow(joined_df) == 0) return(tibble::tibble())
  
  gene_col <- if ("GENE" %in% names(joined_df)) {
    "GENE"
  } else if ("GENE.x" %in% names(joined_df)) {
    cat("    Using GENE.x column (from FDR data)\n")
    "GENE.x"
  } else if ("GENE.y" %in% names(joined_df)) {
    cat("    Using GENE.y column (from event_info)\n")
    "GENE.y"
  } else {
    stop("No GENE column found in joined data")
  }
  
  # Parse coordinates for alternative exon (CO_A) and flanking exons (C1, C2)
  parsed_a <- lapply(joined_df$CO_A, parse_coord_range)
  parsed_c1 <- lapply(joined_df$CO_C1, parse_coord_range)
  parsed_c2 <- lapply(joined_df$CO_C2, parse_coord_range)
  
  out <- joined_df %>%
    mutate(
      START = sapply(parsed_a, function(x) x$start),
      END = sapply(parsed_a, function(x) x$end),
      SCAFFOLD = sapply(parsed_a, function(x) x$chr),
      STRAND = sapply(REF_CO, parse_strand),
      # For upstream/downstream borders, use end of C1 and start of C2
      UPSTRM_EX_BORDER = sapply(parsed_c1, function(x) x$end),
      DOSTRM_EX_BORDER = sapply(parsed_c2, function(x) x$start),
      GROUP = group_label,
      GENE = .data[[gene_col]]
    ) %>%
    select(GENE, START, END, SCAFFOLD, STRAND, UPSTRM_EX_BORDER, DOSTRM_EX_BORDER, GROUP) %>%
    left_join(gene_map, by = c("GENE" = "mgi_symbol")) %>%
    rename(GENEID_ENSEMBL = ensembl_gene_id) %>%
    select(START, END, SCAFFOLD, STRAND, UPSTRM_EX_BORDER, DOSTRM_EX_BORDER, GENEID_ENSEMBL, GROUP)
  
  # Diagnostic: count filtering reasons
  n_before <- nrow(out)
  n_missing_scaffold <- sum(is.na(out$SCAFFOLD))
  n_missing_coords <- sum(is.na(out$START) | is.na(out$END))
  n_missing_strand <- sum(is.na(out$STRAND))
  n_missing_gene <- sum(is.na(out$GENEID_ENSEMBL))
  n_missing_borders <- sum(is.na(out$UPSTRM_EX_BORDER) | is.na(out$DOSTRM_EX_BORDER))
  
  out <- out %>%
    filter(!is.na(SCAFFOLD) & !is.na(START) & !is.na(END) & !is.na(STRAND) & 
             !is.na(GENEID_ENSEMBL) & !is.na(UPSTRM_EX_BORDER) & !is.na(DOSTRM_EX_BORDER))
  
  n_after <- nrow(out)
  if (n_before > n_after) {
    cat("    [EXON] Filtered out", n_before - n_after, "events:",
        "scaffold=", n_missing_scaffold,
        "coords=", n_missing_coords,
        "strand=", n_missing_strand,
        "gene=", n_missing_gene,
        "borders=", n_missing_borders, "\n")
  }
  
  out
}

# Build rMAPS table for introns with flanking exon borders for RBP analysis
build_intron_rbp_table <- function(joined_df, group_label, gene_map) {
  if (nrow(joined_df) == 0) return(tibble::tibble())
  
  gene_col <- if ("GENE" %in% names(joined_df)) {
    "GENE"
  } else if ("GENE.x" %in% names(joined_df)) {
    cat("    Using GENE.x column (from FDR data)\n")
    "GENE.x"
  } else if ("GENE.y" %in% names(joined_df)) {
    cat("    Using GENE.y column (from event_info)\n")
    "GENE.y"
  } else {
    stop("No GENE column found in joined data")
  }
  
  # Parse coordinates for flanking exons (C1, C2) and alternative exon (A)
  # For intron retention (INT) events:
  # CO_A contains the RETAINED INTRON coordinates (the alternative region)
  # CO_C1 and CO_C2 are the flanking exons
  parsed_c1 <- lapply(joined_df$CO_C1, parse_coord_range)
  parsed_c2 <- lapply(joined_df$CO_C2, parse_coord_range)
  parsed_a <- lapply(joined_df$CO_A, parse_coord_range)
  
  out <- joined_df %>%
    mutate(
      # Intron coordinates come directly from CO_A (the retained intron)
      START = sapply(parsed_a, function(x) x$start),
      END = sapply(parsed_a, function(x) x$end),
      SCAFFOLD = sapply(parsed_a, function(x) x$chr),
      STRAND = sapply(REF_CO, parse_strand),
      # Flanking exon borders: use end of C1 and start of C2
      UPSTRM_EX_BORDER = sapply(parsed_c1, function(x) x$end),
      DOSTRM_EX_BORDER = sapply(parsed_c2, function(x) x$start),
      GROUP = group_label,
      GENE = .data[[gene_col]]
    ) %>%
    select(GENE, START, END, SCAFFOLD, STRAND, UPSTRM_EX_BORDER, DOSTRM_EX_BORDER, GROUP) %>%
    left_join(gene_map, by = c("GENE" = "mgi_symbol")) %>%
    rename(GENEID_ENSEMBL = ensembl_gene_id) %>%
    select(UPSTRM_EX_BORDER, START, END, DOSTRM_EX_BORDER, SCAFFOLD, STRAND, GENEID_ENSEMBL, GROUP)
  
  # Diagnostic: count filtering reasons
  n_before <- nrow(out)
  n_missing_scaffold <- sum(is.na(out$SCAFFOLD))
  n_missing_coords <- sum(is.na(out$START) | is.na(out$END))
  n_missing_strand <- sum(is.na(out$STRAND))
  n_missing_gene <- sum(is.na(out$GENEID_ENSEMBL))
  n_missing_borders <- sum(is.na(out$UPSTRM_EX_BORDER) | is.na(out$DOSTRM_EX_BORDER))
  n_invalid_coords <- sum(!is.na(out$START) & !is.na(out$END) & out$START >= out$END)
  
  out <- out %>%
    filter(!is.na(SCAFFOLD) & !is.na(START) & !is.na(END) & !is.na(STRAND) & 
             !is.na(GENEID_ENSEMBL) & !is.na(UPSTRM_EX_BORDER) & !is.na(DOSTRM_EX_BORDER) &
             START < END)
  
  n_after <- nrow(out)
  if (n_before > n_after) {
    cat("    [INTRON] Filtered out", n_before - n_after, "events:",
        "scaffold=", n_missing_scaffold,
        "coords=", n_missing_coords,
        "strand=", n_missing_strand,
        "gene=", n_missing_gene,
        "borders=", n_missing_borders,
        "invalid_coords=", n_invalid_coords, "\n")
  }
  
  out
}

# --- Load data --------------------------------------------------------------
cat("========================================\n")
cat("Loading input data...\n")
cat("========================================\n")

if (!exists("tub_fdr_df")) {
  if (file.exists("tub_fdr.csv")) tub_fdr_df <- read_csv("tub_fdr.csv")[,-1]
  else stop("tub_fdr_df not found and tub_fdr.csv not present")
}
if (!exists("pladb_fdr_df")) {
  if (file.exists("pladb_fdr.csv")) pladb_fdr_df <- read_csv("pladb_fdr.csv")[,-1]
  else stop("pladb_fdr_df not found and pladb_fdr.csv not present")
}
if (!exists("ssa_fdr_df")) {
  if (file.exists("ssa_fdr.csv")) ssa_fdr_df <- read_csv("ssa_fdr.csv")[,-1]
  else stop("ssa_fdr_df not found and ssa_fdr.csv not present")
}
if (!exists("event_info")) {
  if (file.exists("EVENT_INFO-mm10.tab")) event_info <- read.delim("EVENT_INFO-mm10.tab")
  else stop("event_info not found")
}

if (!"EVENT" %in% names(event_info)) stop("event_info must contain column 'EVENT'")
if (!"GENE" %in% names(event_info)) stop("event_info must contain column 'GENE'")

cat("Data loaded successfully.\n")

# --- Get gene symbol to Ensembl ID mapping ----------------------------------
cat("\n========================================\n")
cat("Getting gene symbol to Ensembl ID mapping\n")
cat("========================================\n")

all_genes <- unique(event_info$GENE)
all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]

cat("Found", length(all_genes), "unique gene symbols in event_info\n")

if (file.exists("mm10_gene_mapping.csv")) {
  cat("Loading existing gene mapping from mm10_gene_mapping.csv\n")
  gene_map <- read_csv("mm10_gene_mapping.csv")
  cat("Loaded", nrow(gene_map), "gene mappings\n")
} else {
  gtf_file <- NULL
  possible_gtf_paths <- c(
    "reference_genome/gencode.vM10.annotation.gtf",
    "gencode.vM10.annotation.gtf",
    "reference_genome/gencode.vM10.gtf",
    "mm10.gtf",
    "GRCm38.gtf"
  )
  
  for (path in possible_gtf_paths) {
    if (file.exists(path)) {
      gtf_file <- path
      break
    }
  }
  
  if (is.null(gtf_file)) {
    stop("No GTF file found. Please provide one of:\n  ",
         paste(possible_gtf_paths, collapse = "\n  "))
  }
  
  gene_map <- extract_genes_from_gtf(gtf_file, keep_version = TRUE)
  
  if (nrow(gene_map) > 0) {
    write_csv(gene_map, "mm10_gene_mapping.csv")
    cat("\nSaved gene mapping to mm10_gene_mapping.csv for future use\n")
  }
}

if (nrow(gene_map) == 0) {
  stop("No gene mappings available. Check your GTF file.")
}

genes_in_map <- sum(all_genes %in% gene_map$mgi_symbol)
cat("\nMapping coverage:", genes_in_map, "/", length(all_genes), 
    "(", round(100 * genes_in_map / length(all_genes), 1), "%)\n")

# --- Process each treatment -------------------------------------------------
treatments <- list(
  tub = tub_fdr_df,
  pladb = pladb_fdr_df,
  ssa = ssa_fdr_df
)

for (trname in names(treatments)) {
  cat("\n========================================\n")
  cat("Processing treatment:", trname, "\n")
  cat("========================================\n")
  df <- treatments[[trname]]
  
  # ===== SEPARATE EXON AND INTRON EVENTS FIRST =====
  # Exon events: EX only (not ALTA or ALTD)
  # Intron events: INT only
  exon_events <- df %>% filter(grepl("EX", EVENT) & !grepl("ALTA|ALTD", EVENT))
  intron_events <- df %>% filter(grepl("INT", EVENT) & !grepl("ALTA|ALTD", EVENT))
  
  cat("  Total events - Exons:", nrow(exon_events), "| Introns:", nrow(intron_events), "\n")
  
  # Define groups for EXONS
  exon_up <- exon_events %>% filter(!is.na(deltapsi) & deltapsi >= 0.1 & !is.na(FDR) & FDR <= 0.05)
  exon_down <- exon_events %>% filter(!is.na(deltapsi) & deltapsi <= -0.1 & !is.na(FDR) & FDR <= 0.05)
  exon_nonsig <- exon_events %>% filter(is.na(FDR) | FDR > 0.05 | is.na(deltapsi) | abs(deltapsi) < 0.1)
  
  # Define groups for INTRONS
  intron_up <- intron_events %>% filter(!is.na(deltapsi) & deltapsi >= 0.1 & !is.na(FDR) & FDR <= 0.05)
  intron_down <- intron_events %>% filter(!is.na(deltapsi) & deltapsi <= -0.1 & !is.na(FDR) & FDR <= 0.05)
  intron_nonsig <- intron_events %>% filter(is.na(FDR) | FDR > 0.05 | is.na(deltapsi) | abs(deltapsi) < 0.1)
  
  # Sample background events to max 10,000 for each type
  max_background <- 10000
  if (nrow(exon_nonsig) > max_background) {
    cat("  Sampling", max_background, "exon background events from", nrow(exon_nonsig), "(seed=42)\n")
    set.seed(42)
    exon_nonsig <- exon_nonsig %>% sample_n(max_background)
  }
  if (nrow(intron_nonsig) > max_background) {
    cat("  Sampling", max_background, "intron background events from", nrow(intron_nonsig), "(seed=42)\n")
    set.seed(42)
    intron_nonsig <- intron_nonsig %>% sample_n(max_background)
  }
  
  cat("  Exon events - up:", nrow(exon_up), "| down:", nrow(exon_down), "| nonsig:", nrow(exon_nonsig), "\n")
  cat("  Intron events - up:", nrow(intron_up), "| down:", nrow(intron_down), "| nonsig:", nrow(intron_nonsig), "\n")
  
  # Join with event_info to get coordinates
  cat("\n  Joining with event_info...\n")
  # EXONS
  exon_up_info <- exon_up %>% inner_join(event_info, by = "EVENT")
  exon_down_info <- exon_down %>% inner_join(event_info, by = "EVENT")
  exon_nonsig_info <- exon_nonsig %>% inner_join(event_info, by = "EVENT")
  
  # INTRONS
  intron_up_info <- intron_up %>% inner_join(event_info, by = "EVENT")
  intron_down_info <- intron_down %>% inner_join(event_info, by = "EVENT")
  intron_nonsig_info <- intron_nonsig %>% inner_join(event_info, by = "EVENT")
  
  # Build tables for exons with RBP format
  cat("\n  Building exon RBP tables...\n")
  cat("    [INPUT] up:", nrow(exon_up_info), "| down:", nrow(exon_down_info), "| ndiff:", nrow(exon_nonsig_info), "\n")
  exon_up <- build_exon_rbp_table(exon_up_info, "up", gene_map)
  exon_down <- build_exon_rbp_table(exon_down_info, "down", gene_map)
  exon_nonsig <- build_exon_rbp_table(exon_nonsig_info, "ndiff", gene_map)
  cat("    [OUTPUT] up:", nrow(exon_up), "| down:", nrow(exon_down), "| ndiff:", nrow(exon_nonsig), "\n")
  
  # Build tables for introns with RBP format
  cat("\n  Building intron RBP tables...\n")
  cat("    [INPUT] up:", nrow(intron_up_info), "| down:", nrow(intron_down_info), "| ndiff:", nrow(intron_nonsig_info), "\n")
  intron_up <- build_intron_rbp_table(intron_up_info, "up", gene_map)
  intron_down <- build_intron_rbp_table(intron_down_info, "down", gene_map)
  intron_nonsig <- build_intron_rbp_table(intron_nonsig_info, "ndiff", gene_map)
  cat("    [OUTPUT] up:", nrow(intron_up), "| down:", nrow(intron_down), "| ndiff:", nrow(intron_nonsig), "\n")
  
  # Combine all groups (ndiff last as reference group)
  exon_combined <- bind_rows(exon_up, exon_down, exon_nonsig)
  intron_combined <- bind_rows(intron_up, intron_down, intron_nonsig)
  
  # Write output files
  if (nrow(exon_combined) > 0) {
    exon_file <- paste0(trname, "_exons_rbp.tab")
    write_tsv(exon_combined, exon_file)
    cat("\n  -> Wrote", nrow(exon_combined), "exons to", exon_file, "\n")
    cat("     up:", nrow(exon_up), "| down:", nrow(exon_down), "| ndiff:", nrow(exon_nonsig), "\n")
    cat("     Columns:", paste(names(exon_combined), collapse=", "), "\n")
  } else {
    cat("\n  -> No exon data to write\n")
  }
  
  if (nrow(intron_combined) > 0) {
    intron_file <- paste0(trname, "_introns_rbp.tab")
    write_tsv(intron_combined, intron_file)
    cat("  -> Wrote", nrow(intron_combined), "introns to", intron_file, "\n")
    cat("     up:", nrow(intron_up), "| down:", nrow(intron_down), "| ndiff:", nrow(intron_nonsig), "\n")
    cat("     Columns:", paste(names(intron_combined), collapse=", "), "\n")
  } else {
    cat("  -> No intron data to write\n")
  }
}

cat("\n========================================\n")
cat("All rMAPS RBP files created successfully.\n")
cat("========================================\n")
cat("\nGROUP column contains:\n")
cat("  - 'up': events with deltapsi >= 0.1 and FDR <= 0.05\n")
cat("  - 'down': events with deltapsi <= -0.1 and FDR <= 0.05\n")
cat("  - 'ndiff': non-significant background (reference group, max 10,000)\n")
cat("\nExample matt rna_maps_cisbp command for exons:\n")
cat("matt rna_maps_cisbp tub_exons_rbp.tab UPSTRM_EX_BORDER START END DOSTRM_EX_BORDER \\\n")
cat("  SCAFFOLD STRAND GROUP[up,down,ndiff] 31 35 135 \\\n")
cat("  reference_genome/GRCm38.primary_assembly.genome.fa cisbprna_regexps \\\n")
cat("  -d tub_exons_rbp_maps -p 0.05\n")
cat("\nExample matt rna_maps_cisbp command for introns:\n")
cat("matt rna_maps_cisbp tub_introns_rbp.tab UPSTRM_EX_BORDER START END DOSTRM_EX_BORDER \\\n")
cat("  SCAFFOLD STRAND GROUP[up,down,ndiff] 31 35 135 \\\n")
cat("  reference_genome/GRCm38.primary_assembly.genome.fa cisbprna_regexps \\\n")
cat("  -d tub_introns_rbp_maps -p 0.05\n")
cat("\nNote: Use -p or -fdr arguments to highlight significant enrichment/depletion.\n")
cat("The last group (ndiff) is the reference group for statistical testing.\n")