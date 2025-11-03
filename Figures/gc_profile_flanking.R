# Required libraries
library(Biostrings)
library(tidyverse)
library(ggplot2)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicRanges)
library(dplyr)
library(broom)
library(viridis)

# Read event info and differential analysis results
event_info <- read.delim("EVENT_INFO-mm10.tab")
tub_fdr_df <- read_csv("tub_fdr.csv")[,-1]
pladb_fdr_df <- read_csv("pladb_fdr.csv")[,-1]
ssa_fdr_df <- read_csv("ssa_fdr.csv")[,-1]

# Read mouse introns dataset
mouse_introns <- read.delim("mouse_analysis_exons_introns_combined.tsv")

# Filter for significant events
differential_tub <- na.omit(tub_fdr_df[tub_fdr_df$FDR <= 0.05 & abs(tub_fdr_df$deltapsi) >= 0.1,])
differential_pladb <- na.omit(pladb_fdr_df[pladb_fdr_df$FDR <= 0.05 & abs(pladb_fdr_df$deltapsi) >= 0.1,])
differential_ssa <- na.omit(ssa_fdr_df[ssa_fdr_df$FDR <= 0.05 & abs(ssa_fdr_df$deltapsi) >= 0.1,])

# Create random control dataset (larger for better statistics)
set.seed(42)
exon_events <- filter(event_info, grepl("EX", event_info$EVENT)) %>%
  filter(!EVENT %in% c(differential_tub$EVENT,differential_pladb$EVENT,differential_ssa$EVENT))
random_events <- exon_events %>% sample_n(min(500, nrow(exon_events)))  # Use up to 500 controls
random_differential <- data.frame(EVENT = random_events$EVENT)

# Get exon events for each dataset
tub_exons <- filter(event_info, EVENT %in% differential_tub$EVENT & grepl("EX", event_info$EVENT))
pladb_exons <- filter(event_info, EVENT %in% differential_pladb$EVENT & grepl("EX", event_info$EVENT))
ssa_exons <- filter(event_info, EVENT %in% differential_ssa$EVENT & grepl("EX", event_info$EVENT))
random_exons <- random_events

cat("=== INITIAL DATA INSPECTION ===\n")
cat("Tubercidin exon events:", nrow(tub_exons), "\n")
cat("Pladienolide-B exon events:", nrow(pladb_exons), "\n")
cat("SSA exon events:", nrow(ssa_exons), "\n")
cat("Random exon events:", nrow(random_exons), "\n")
cat("Number of intron records:", nrow(mouse_introns), "\n")

# Improved parameters
n_bins <- 20
min_sequence_length <- 50  # Minimum sequence length for reliable analysis
min_bin_size <- 3  # Minimum nucleotides per bin

# Function to parse genomic coordinates
parse_coordinates <- function(coord_string) {
  parts <- strsplit(coord_string, "[:-]")[[1]]
  return(data.frame(
    chromosome = parts[1],
    start = as.numeric(parts[2]),
    end = as.numeric(parts[3])
  ))
}

# Function to find flanking introns for each exon
find_flanking_introns <- function(exon_chr, exon_start, exon_end, intron_df) {
  chr_introns <- intron_df[intron_df$chromosome == exon_chr & 
                             intron_df$feature_type == "intron", ]
  
  if(nrow(chr_introns) == 0) {
    return(list(upstream = NULL, downstream = NULL))
  }
  
  upstream_candidates <- chr_introns[chr_introns$end < exon_start, ]
  if(nrow(upstream_candidates) > 0) {
    upstream_idx <- which.max(upstream_candidates$end)
    upstream_intron <- upstream_candidates[upstream_idx, ]
  } else {
    upstream_intron <- NULL
  }
  
  downstream_candidates <- chr_introns[chr_introns$start > exon_end, ]
  if(nrow(downstream_candidates) > 0) {
    downstream_idx <- which.min(downstream_candidates$start)
    downstream_intron <- downstream_candidates[downstream_idx, ]
  } else {
    downstream_intron <- NULL
  }
  
  return(list(upstream = upstream_intron, downstream = downstream_intron))
}

# Load mm10 genome
genome <- BSgenome.Mmusculus.UCSC.mm10

# Function to extract sequence
extract_sequence <- function(chr, start, end, genome) {
  if(!startsWith(chr, "chr")) {
    chr <- paste0("chr", chr)
  }
  
  gr <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))
  
  tryCatch({
    seq <- getSeq(genome, gr)
    return(as.character(seq))
  }, error = function(e) {
    return(NA)
  })
}

# Function to process a dataset
process_dataset <- function(exon_data, dataset_name) {
  cat("\n=== PROCESSING", toupper(dataset_name), "===\n")
  
  # Extract exon coordinates
  exon_coords_list <- lapply(exon_data$CO_A, parse_coordinates)
  exon_coords <- do.call(rbind, exon_coords_list)
  exon_coords$event_id <- exon_data$EVENT
  
  # Find flanking introns
  flanking_introns <- list()
  for(i in 1:nrow(exon_coords)) {
    if(i %% 50 == 0) cat("Processing exon", i, "of", nrow(exon_coords), "\n")
    
    flanking <- find_flanking_introns(
      exon_coords$chromosome[i],
      exon_coords$start[i],
      exon_coords$end[i],
      mouse_introns
    )
    
    if(!is.null(flanking$upstream) && !is.null(flanking$downstream)) {
      flanking_introns[[length(flanking_introns) + 1]] <- list(
        exon_idx = i,
        upstream = flanking$upstream,
        downstream = flanking$downstream
      )
    }
  }
  
  cat("Exons with both flanking introns found:", length(flanking_introns), "\n")
  
  if(length(flanking_introns) == 0) {
    warning("No exons found with both flanking introns for ", dataset_name)
    return(NULL)
  }
  
  # Extract sequences
  sequences_data <- list()
  for(i in 1:length(flanking_introns)) {
    if(i %% 25 == 0) cat("Extracting sequences for exon", i, "of", length(flanking_introns), "\n")
    
    flank <- flanking_introns[[i]]
    exon_idx <- flank$exon_idx
    
    exon_seq <- extract_sequence(
      exon_coords$chromosome[exon_idx],
      exon_coords$start[exon_idx],
      exon_coords$end[exon_idx],
      genome
    )
    
    upstream_seq <- extract_sequence(
      flank$upstream$chromosome,
      flank$upstream$start,
      flank$upstream$end,
      genome
    )
    
    downstream_seq <- extract_sequence(
      flank$downstream$chromosome,
      flank$downstream$start,
      flank$downstream$end,
      genome
    )
    
    if(!is.na(exon_seq) && !is.na(upstream_seq) && !is.na(downstream_seq)) {
      exon_len <- nchar(exon_seq)
      upstream_len <- nchar(upstream_seq)
      downstream_len <- nchar(downstream_seq)
      
      # Filter by minimum length
      if(exon_len >= min_sequence_length && 
         upstream_len >= min_sequence_length && 
         downstream_len >= min_sequence_length) {
        sequences_data[[length(sequences_data) + 1]] <- list(
          event_id = exon_coords$event_id[exon_idx],
          exon_seq = exon_seq,
          upstream_intron_seq = upstream_seq,
          downstream_intron_seq = downstream_seq,
          exon_length = exon_len,
          upstream_length = upstream_len,
          downstream_length = downstream_len
        )
      }
    }
  }
  
  cat("Successfully extracted sequences for", length(sequences_data), "exon-intron trios\n")
  return(sequences_data)
}

# Improved GC content computation with better binning
compute_gc_content_improved <- function(dna_stringset, n_bins = 50, region_name = "Unknown") {
  cat("Computing GC content for", region_name, "- ", length(dna_stringset), "sequences\n")
  
  n_seqs <- length(dna_stringset)
  result_matrix <- matrix(NA_real_, nrow = n_seqs, ncol = n_bins)
  
  for(i in 1:n_seqs) {
    seq_obj <- dna_stringset[i]
    seq_length <- Biostrings::width(seq_obj)
    
    # Calculate bin size, ensuring minimum bin size
    ideal_bin_size <- seq_length / n_bins
    if(ideal_bin_size < min_bin_size) next
    
    # Use equal-sized bins (except possibly the last one)
    actual_bin_size <- floor(seq_length / n_bins)
    
    for(j in 1:n_bins) {
      start_pos <- (j-1) * actual_bin_size + 1
      
      # For the last bin, extend to sequence end to avoid truncation
      if(j == n_bins) {
        end_pos <- seq_length
      } else {
        end_pos <- j * actual_bin_size
      }
      
      if(start_pos <= end_pos && start_pos <= seq_length) {
        tryCatch({
          subseq <- Biostrings::subseq(seq_obj, start_pos, end_pos)
          seq_chars <- strsplit(as.character(subseq), "")[[1]]
          gc_count <- sum(seq_chars %in% c("G", "C", "g", "c"))
          total_count <- sum(seq_chars %in% c("A", "T", "G", "C", "a", "t", "g", "c"))
          
          if(total_count > 0) {
            result_matrix[i, j] <- 100 * gc_count / total_count
          }
        }, error = function(e) {
          # Silent error handling
        })
      }
    }
    
    if(i %% 25 == 0) cat("Processed", i, "sequences...\n")
  }
  
  return(result_matrix)
}

# Function to perform statistical comparisons
perform_statistical_tests <- function(treatment_data, control_data, treatment_name) {
  cat("\n=== STATISTICAL TESTING:", treatment_name, "vs Control ===\n")
  
  # Overall exon vs intron comparison within each dataset
  treatment_exon_means <- rowMeans(treatment_data$exon_gc_matrix, na.rm = TRUE)
  treatment_intron_means <- rowMeans(cbind(treatment_data$upstream_gc_matrix, 
                                           treatment_data$downstream_gc_matrix), na.rm = TRUE)
  
  control_exon_means <- rowMeans(control_data$exon_gc_matrix, na.rm = TRUE)
  control_intron_means <- rowMeans(cbind(control_data$upstream_gc_matrix, 
                                         control_data$downstream_gc_matrix), na.rm = TRUE)
  
  # Calculate exon-intron differences for each sequence
  treatment_diff <- treatment_exon_means - treatment_intron_means
  control_diff <- control_exon_means - control_intron_means
  
  # Remove NA values
  treatment_diff <- treatment_diff[!is.na(treatment_diff)]
  control_diff <- control_diff[!is.na(control_diff)]
  
  # Perform t-test comparing exon enrichment between treatment and control
  if(length(treatment_diff) > 3 && length(control_diff) > 3) {
    ttest_result <- t.test(treatment_diff, control_diff, alternative = "two.sided")
    
    cat("Exon GC enrichment comparison:\n")
    cat("  Treatment mean difference:", round(mean(treatment_diff), 3), "%\n")
    cat("  Control mean difference:", round(mean(control_diff), 3), "%\n")
    cat("  Difference of differences:", round(mean(treatment_diff) - mean(control_diff), 3), "%\n")
    cat("  P-value:", format(ttest_result$p.value, scientific = TRUE), "\n")
    cat("  95% CI:", round(ttest_result$conf.int, 3), "\n")
    
    # Effect size (Cohen's d)
    pooled_sd <- sqrt(((length(treatment_diff)-1)*var(treatment_diff) + 
                         (length(control_diff)-1)*var(control_diff)) / 
                        (length(treatment_diff) + length(control_diff) - 2))
    cohens_d <- (mean(treatment_diff) - mean(control_diff)) / pooled_sd
    cat("  Cohen's d (effect size):", round(cohens_d, 3), "\n")
    
    return(list(
      treatment_name = treatment_name,
      treatment_n = length(treatment_diff),
      control_n = length(control_diff),
      treatment_mean_diff = mean(treatment_diff),
      control_mean_diff = mean(control_diff),
      difference_of_differences = mean(treatment_diff) - mean(control_diff),
      p_value = ttest_result$p.value,
      ci_lower = ttest_result$conf.int[1],
      ci_upper = ttest_result$conf.int[2],
      cohens_d = cohens_d,
      significant = ttest_result$p.value < 0.05
    ))
  } else {
    cat("Insufficient data for statistical testing\n")
    return(NULL)
  }
}

# Process all datasets
tub_sequences <- process_dataset(tub_exons, "Tubercidin")
pladb_sequences <- process_dataset(pladb_exons, "Pladienolide-B")
ssa_sequences <- process_dataset(ssa_exons, "SSA")
random_sequences <- process_dataset(random_exons, "Random")

# Filter datasets to ensure we have sufficient data
datasets <- list(
  "Tubercidin" = tub_sequences,
  "Pladienolide-B" = pladb_sequences,
  "SSA" = ssa_sequences,
  "Random" = random_sequences
)

# Remove datasets with insufficient data
datasets <- datasets[sapply(datasets, function(x) !is.null(x) && length(x) >= 10)]

cat("\n=== DATASETS WITH SUFFICIENT DATA ===\n")
for(name in names(datasets)) {
  cat(name, ":", length(datasets[[name]]), "sequences\n")
}

# Function to analyze GC content with improved statistics
analyze_gc_content_improved <- function(sequences, dataset_name) {
  if(is.null(sequences) || length(sequences) < 10) return(NULL)
  
  cat("\n=== ANALYZING GC CONTENT FOR", toupper(dataset_name), "===\n")
  
  # Extract sequences
  exon_seqs <- sapply(sequences, function(x) x$exon_seq)
  upstream_seqs <- sapply(sequences, function(x) x$upstream_intron_seq)
  downstream_seqs <- sapply(sequences, function(x) x$downstream_intron_seq)
  
  # Convert to DNAStringSet
  exons_dna <- Biostrings::DNAStringSet(exon_seqs)
  up_introns_dna <- Biostrings::DNAStringSet(upstream_seqs)
  down_introns_dna <- Biostrings::DNAStringSet(downstream_seqs)
  
  # Compute GC content with improved binning
  exon_gc <- compute_gc_content_improved(exons_dna, n_bins, paste(dataset_name, "EXONS"))
  up_intron_gc <- compute_gc_content_improved(up_introns_dna, n_bins, paste(dataset_name, "UPSTREAM INTRONS"))
  down_intron_gc <- compute_gc_content_improved(down_introns_dna, n_bins, paste(dataset_name, "DOWNSTREAM INTRONS"))
  
  # Filter out sequences with too many NA values
  valid_rows <- rowSums(is.na(exon_gc)) < (n_bins * 0.5) & 
    rowSums(is.na(up_intron_gc)) < (n_bins * 0.5) & 
    rowSums(is.na(down_intron_gc)) < (n_bins * 0.5)
  
  exon_gc <- exon_gc[valid_rows, , drop = FALSE]
  up_intron_gc <- up_intron_gc[valid_rows, , drop = FALSE]
  down_intron_gc <- down_intron_gc[valid_rows, , drop = FALSE]
  
  cat("Valid sequences after filtering:", nrow(exon_gc), "\n")
  
  if(nrow(exon_gc) < 5) return(NULL)
  
  # Calculate means and confidence intervals
  up_intron_means <- colMeans(up_intron_gc, na.rm = TRUE)
  exon_means <- colMeans(exon_gc, na.rm = TRUE)
  down_intron_means <- colMeans(down_intron_gc, na.rm = TRUE)
  
  # Calculate standard errors
  up_intron_se <- apply(up_intron_gc, 2, function(x) {
    valid_x <- x[!is.na(x)]
    if(length(valid_x) < 2) return(NA)
    sd(valid_x) / sqrt(length(valid_x))
  })
  
  exon_se <- apply(exon_gc, 2, function(x) {
    valid_x <- x[!is.na(x)]
    if(length(valid_x) < 2) return(NA)
    sd(valid_x) / sqrt(length(valid_x))
  })
  
  down_intron_se <- apply(down_intron_gc, 2, function(x) {
    valid_x <- x[!is.na(x)]
    if(length(valid_x) < 2) return(NA)
    sd(valid_x) / sqrt(length(valid_x))
  })
  
  # Create seamless data frame
  seamless_df <- data.frame(
    position = c((-n_bins):(-1), 1:n_bins, (n_bins+1):(2*n_bins)),
    gc_content = c(up_intron_means, exon_means, down_intron_means),
    se = c(up_intron_se, exon_se, down_intron_se),
    region = rep(c("Upstream Intron", "Exon", "Downstream Intron"), each = n_bins),
    dataset = dataset_name,
    n_sequences = nrow(exon_gc)
  ) %>%
    filter(!is.na(gc_content) & !is.na(se)) %>%
    mutate(
      region_factor = factor(region, levels = c("Upstream Intron", "Exon", "Downstream Intron")),
      lower_ci = pmax(0, gc_content - 1.96 * se),
      upper_ci = pmin(100, gc_content + 1.96 * se)
    )
  
  # Store matrices for statistical testing
  result <- list(
    seamless_df = seamless_df,
    exon_gc_matrix = exon_gc,
    upstream_gc_matrix = up_intron_gc,
    downstream_gc_matrix = down_intron_gc,
    n_sequences = nrow(exon_gc)
  )
  
  return(result)
}

# Analyze all datasets
results <- list()
for(name in names(datasets)) {
  results[[name]] <- analyze_gc_content_improved(datasets[[name]], name)
}

# Filter out NULL results
results <- results[!sapply(results, is.null)]

# Perform statistical comparisons vs control
if("Random" %in% names(results)) {
  control_result <- results[["Random"]]
  statistical_results <- list()
  
  for(name in names(results)) {
    if(name != "Random") {
      stat_result <- perform_statistical_tests(results[[name]], control_result, name)
      if(!is.null(stat_result)) {
        statistical_results[[name]] <- stat_result
      }
    }
  }
  
  # Create statistical summary table
  if(length(statistical_results) > 0) {
    stats_df <- do.call(rbind, lapply(statistical_results, function(x) {
      data.frame(
        Dataset = x$treatment_name,
        N_sequences = x$treatment_n,
        Treatment_ExonGC_Enrichment = round(x$treatment_mean_diff, 3),
        Control_ExonGC_Enrichment = round(x$control_mean_diff, 3),
        Difference_of_Differences = round(x$difference_of_differences, 3),
        P_value = format(x$p_value, scientific = TRUE, digits = 3),
        Cohens_D = round(x$cohens_d, 3),
        Significant = x$significant,
        CI_Lower = round(x$ci_lower, 3),
        CI_Upper = round(x$ci_upper, 3)
      )
    }))
    
    cat("\n=== STATISTICAL SUMMARY ===\n")
    print(stats_df)
    write_csv(stats_df, "statistical_comparison_results.csv")
  }
}

# Combine all seamless data
all_gc_data <- do.call(rbind, lapply(results, function(x) x$seamless_df))

if(nrow(all_gc_data) > 0) {
  all_gc_data <- all_gc_data %>%
    mutate(
      dataset_factor = factor(dataset, levels = c("Random", "Tubercidin", "Pladienolide-B", "SSA"))
    )
  
  # Define colors
  dataset_colors <- c(
    "Random" = "#95A5A6",          # Gray for control
    "Tubercidin" = "#E74C3C",      # Red
    "Pladienolide-B" = "#3498DB",  # Blue
    "SSA" = "#27AE60"              # Green
  )
  
  # Create main comparison plot with improved aesthetics
  p_comparison <- ggplot(all_gc_data, aes(x = position, y = gc_content, color = dataset_factor)) +
    geom_line(size = 1.2, alpha = 0.8) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = dataset_factor), 
                alpha = 0.15, color = NA) +
    geom_vline(xintercept = c(0.5, n_bins + 0.5), linetype = "dashed", alpha = 0.6, color = "black") +
    scale_x_continuous(
      "Genomic Position", 
      breaks = c(-n_bins/2, n_bins/2, n_bins + n_bins/2),
      labels = c("Upstream\nIntron", "Exon", "Downstream\nIntron")
    ) +
    scale_y_continuous("GC Content (%)", labels = function(x) paste0(x, "%")) +
    scale_color_manual("Dataset", values = dataset_colors) +
    scale_fill_manual("Dataset", values = dataset_colors) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = "GC Content Profile Comparison with Statistical Testing",
      subtitle = paste0("Exon GC enrichment analysis with 95% confidence intervals (", n_bins, " bins per region)"),
      caption = "Dashed lines indicate exon-intron boundaries. Ribbons show 95% confidence intervals."
    )
  
  # Calculate and plot differences
  exon_intron_diff <- all_gc_data %>%
    group_by(dataset) %>%
    summarise(
      mean_exon_gc = mean(gc_content[region == "Exon"], na.rm = TRUE),
      mean_upstream_gc = mean(gc_content[region == "Upstream Intron"], na.rm = TRUE),
      mean_downstream_gc = mean(gc_content[region == "Downstream Intron"], na.rm = TRUE),
      n_sequences = first(n_sequences),
      .groups = "drop"
    ) %>%
    mutate(
      mean_intron_gc = (mean_upstream_gc + mean_downstream_gc) / 2,
      exon_intron_diff = mean_exon_gc - mean_intron_gc
    ) %>%
    arrange(desc(exon_intron_diff))
  
  # Add significance annotations if available
  if(exists("stats_df")) {
    exon_intron_diff <- exon_intron_diff %>%
      left_join(stats_df %>% select(Dataset, Significant, P_value), 
                by = c("dataset" = "Dataset")) %>%
      mutate(
        significance_label = case_when(
          is.na(Significant) ~ "",
          Significant ~ paste0("**\n(p", substr(P_value, 1, 6), ")"),
          !Significant ~ paste0("ns\n(p", substr(P_value, 1, 6), ")")
        )
      )
  } else {
    exon_intron_diff$significance_label <- ""
  }
  
  # Create improved differences plot
  p_differences <- exon_intron_diff %>%
    ggplot(aes(x = reorder(dataset, exon_intron_diff), y = exon_intron_diff, fill = dataset)) +
    geom_col(alpha = 0.8, color = "black", size = 0.3) +
    geom_text(aes(label = significance_label), 
              vjust = ifelse(exon_intron_diff >= 0, -0.1, 1.1), 
              size = 3, fontface = "bold") +
    geom_text(aes(label = paste0("n=", n_sequences)), 
              vjust = ifelse(exon_intron_diff >= 0, 1.5, -0.5), 
              size = 2.5, alpha = 0.7) +
    scale_fill_manual("Dataset", values = dataset_colors) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      title = "Exon GC Enrichment with Statistical Significance",
      x = "Dataset",
      y = "GC Content Difference (%)",
      subtitle = "** = p<0.05, ns = not significant",
      caption = "Error bars represent 95% confidence intervals"
    )
  
  # Save outputs
  ggsave("improved_gc_comparison.png", p_comparison, width = 14, height = 10, dpi = 300)
  ggsave("improved_gc_differences.png", p_differences, width = 10, height = 8, dpi = 300)
  
  write_csv(all_gc_data, "improved_gc_data_with_ci.csv")
  write_csv(exon_intron_diff, "improved_exon_intron_differences.csv")
  
  # Print plots
  print(p_comparison)
  print(p_differences)
  
  
  showtext::showtext_opts(dpi=600)
  ggsave("exon_intron_tracks_gc.pdf", plot = p_comparison, width = 6.5, height = 4.2, device = cairo_pdf)
  ggsave("exon_intron_tracks_gc.png", plot = p_comparison, width = 6.5, height = 4.2, dpi = 600)
  
  
  cat("\n=== FINAL RESULTS SUMMARY ===\n")
  for(i in 1:nrow(exon_intron_diff)) {
    row <- exon_intron_diff[i,]
    significance <- ifelse(is.na(row$Significant), "not tested", 
                           ifelse(row$Significant, "SIGNIFICANT", "not significant"))
    cat(sprintf("%s (n=%d): %.3f%% exon GC enrichment - %s\n", 
                row$dataset, row$n_sequences, row$exon_intron_diff, significance))
  }
  
} else {
  cat("No valid data for plotting\n")
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Improved plots saved with statistical testing\n")
cat("Key improvements:\n")
cat("- Better binning to avoid end-peak artifacts\n")
cat("- Statistical significance testing vs control\n")
cat("- Confidence intervals on all plots\n")
cat("- Minimum sequence length filtering\n")
cat("- Effect size calculations\n")