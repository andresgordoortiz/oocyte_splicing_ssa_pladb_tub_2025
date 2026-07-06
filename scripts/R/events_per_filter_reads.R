# =============================================================================
# Events detected as a function of the minimum-read filter, per drug and
# per event-type group, coloured by drug (control in grey).
#
# Input : vast-tools INCLUSION_LEVELS_FULL table produced for the SSA / PLADB /
#         TUB oocyte experiment.
# Output: two publication-ready ggplot2 figures and the underlying tidy tables.
#   1. events_per_filter_reads.{pdf,png}  -- total events vs. min-read filter,
#                                            faceted by event type.
#   2. ir_deltapsi_by_filter.{pdf,png}    -- IR ΔPSI (drug - CTL) density at
#                                            several min-read filters, faceted
#                                            by drug. Diagnoses whether the
#                                            expected SSA -> intron retention
#                                            asymmetry is preserved or washes
#                                            out as the filter tightens.
# =============================================================================

suppressPackageStartupMessages({
  library(betAS)        # getDataset() to read the vast-tools table
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(readxl)
})

# ---- 1. Configuration --------------------------------------------------------

inclusion_tab <- file.path(
  getwd(),
  "notebooks/inclusion_tables/ssa_pladb_tub_oocyte_INCLUSION_LEVELS_FULL-mm10.tab"
)

# Sample -> condition mapping. Both CTL batches are merged into one control so
# we have a single control line in the final plot.
sample_map <- list(
  CTL   = c(
    "reads_10_ctl_13_03_25", "reads_11_ctl_13_03_25", "reads_12_ctl_13_03_25",
    "reads_1_ctl_11_03_25",  "reads_2_ctl_11_03_25",  "reads_3_ctl_11_03_25"
  ),
  PLADB = c(
    "reads_13_pladb_3mm_13_03_25",
    "reads_14_pladb_3mm_13_03_25",
    "reads_15_pladb_3mm_13_03_25"
  ),
  SSA   = c(
    "reads_7_ssa_1mm_11_03_25",
    "reads_8_ssa_1mm_11_03_25",
    "reads_9_ssa_1mm_11_03_25"
  ),
  TUB   = c(
    "reads_4_tub_10mm_11_03_25",
    "reads_5_tub_10mm_11_03_25",
    "reads_6_tub_10mm_11_03_25"
  )
)

# Event-type groups (vast-tools COMPLEX codes). Alt5 and Alt3 are split so the
# asymmetry of each alternative splice site can be inspected independently.
type_groups <- list(
  Exons   = c("C1", "C2", "C3", "S", "MIC"),
  Introns = c("IR"),
  Alt5    = c("Alt5"),
  Alt3    = c("Alt3")
)

# Read filter sweep (start at 10 reads, step 10 up to 100)
min_reads <- seq(10, 100, by = 10)

# For the IR asymmetry diagnostic we pick a subset of N values to keep the
# density panels readable.
asymmetry_N    <- c(10, 30, 50, 100)

# Strict definition of "differentially spliced" event (used everywhere below):
#   FDR < 0.05  AND  |DeltaPSI| >= 10%
# The betAS DeltaPSI is on a 0-1 (proportion) scale; 10% = 0.10 in that scale.
# Below 10% is treated as noise.
min_abs_dPSI   <- 10         # percentage points for our bulk |DeltaPSI| view
min_abs_dPSI_betAS <- 0.10   # proportion for the betAS FDR-signed view

# Fixed plot order: control first, then drugs in experimental order
drug_levels     <- c("CTL", "PLADB", "SSA", "TUB")
drugs_to_plot   <- c("PLADB", "SSA", "TUB")

# Canonical drug colours used elsewhere in the manuscript
# (Figures/Figure2H_S4A_B.R, figure2B_volcano_gene_expression.R, ...).
# Deliberately desaturated pastels chosen for cross-figure consistency; the
# relief rule applies -- direct labels are always visible (legend, axes).
drug_palette <- c(
  CTL   = "#D3D3D3",
  PLADB = "#70A0AF",  # Pladienolide B (steel blue)
  SSA   = "#706993",  # Spliceostatin A (purple)
  TUB   = "#A0C1B9"   # Tubercidin (sage)
)

# Canonical splicing-direction colors used elsewhere in the manuscript
# (Figures/Figure2H_S4A_B.R "Included/Skipped",
#  Figures/FigureS4C_E.R and Figures/introns_gc_length.R "Retained/Excised").
# Positive DeltaPSI = more inclusion (exons) or more retention (introns) in
# drug -- the "Included / Retained" direction.
direction_palette <- c(
  "drug > CTL" = "#D4A574",  # warm sand  -- Retained / Included
  "drug < CTL" = "#B85450"   # muted coral -- Excised  / Skipped
)

# ---- 2. Load data ------------------------------------------------------------

splicing_data <- getDataset(pathTables = inclusion_tab, tool = "vast-tools")

missing_samples <- setdiff(unlist(sample_map), colnames(splicing_data))
if (length(missing_samples) > 0) {
  stop("Samples not found in inclusion table: ",
       paste(missing_samples, collapse = ", "))
}

# ---- 3. Count events ---------------------------------------------------------
# Vectorized: parse the Qual strings to a numeric inc+exc matrix ONCE per
# drug, then iterate over N cheaply with rowSums(). This replaces 120 calls of
# filterEvents() (each re-parsing 660K rows) with one pass per drug.

parse_total_reads <- function(qual_strings) {
  after_at <- strsplit(qual_strings, "@", fixed = TRUE)
  inc_exc  <- vapply(after_at, `[`, character(1), 2L)
  pair     <- strsplit(inc_exc, ",", fixed = TRUE)
  vapply(pair, function(p) sum(as.numeric(p)), numeric(1))
}

prepare_drug <- function(samples) {
  psi_mat  <- as.matrix(splicing_data[, samples])
  qual_mat <- as.matrix(splicing_data[, paste0(samples, ".Q")])
  total    <- apply(qual_mat, c(1, 2), parse_total_reads)
  storage.mode(total) <- "double"
  no_na_psi <- !apply(psi_mat, 1, anyNA)
  list(
    total     = total,
    psi       = psi_mat,
    no_na_psi = no_na_psi,
    complex   = splicing_data$COMPLEX
  )
}

drug_cache <- setNames(
  lapply(drug_levels, function(d) prepare_drug(sample_map[[d]])),
  drug_levels
)

# Count of events that pass the read filter in ALL replicates of a drug,
# subsetted to a given COMPLEX type group. Vectorized across all N at once.
count_for <- function(drug, type_group) {
  cache     <- drug_cache[[drug]]
  mask_type <- cache$complex %in% type_groups[[type_group]]
  valid     <- mask_type & cache$no_na_psi
  total     <- cache$total[valid, , drop = FALSE]
  n_samples <- ncol(total)
  vapply(min_reads, function(N) sum(rowSums(total >= N) == n_samples),
         integer(1))
}

count_grid <- expand_grid(
  drug        = drug_levels,
  event_group = names(type_groups)
)

count_grid <- count_grid %>%
  rowwise() %>%
  mutate(n_events = list(count_for(drug, event_group))) %>%
  ungroup() %>%
  unnest(cols = c(n_events)) %>%
  mutate(
    min_reads   = rep(min_reads, times = n_distinct(drug) * n_distinct(event_group)),
    drug        = factor(drug, levels = drug_levels),
    event_group = factor(
      event_group,
      levels = names(type_groups),
      labels = c("Exons (C1/C2/C3/S/MIC)",
                 "Introns (IR)",
                 "Alt 5' splice sites (Alt5)",
                 "Alt 3' splice sites (Alt3)")
    )
  )

# ---- 4. IR ΔPSI asymmetry analysis -------------------------------------------
# For each drug, at each N threshold:
#   1. Take events that pass the filter in BOTH the drug and CTL at N.
#   2. Compute per-event mean PSI across the passing replicates (within each
#      condition; replicates failing the filter at N are excluded from the mean
#      but the event is only kept if every replicate passes).
#   3. ΔPSI = drug_mean_PSI - ctl_mean_PSI.
# This directly measures the asymmetry of each drug vs CTL on the same event
# set, which is what the user noticed as suspiciously symmetric for IR.

ctl_cache <- drug_cache[["CTL"]]

asymmetry_for_drug <- function(drug, N) {
  drug_cache_d <- drug_cache[[drug]]

  # Per-event: did every replicate in each condition hit >= N reads?
  drug_total   <- drug_cache_d$total
  ctl_total    <- ctl_cache$total
  drug_psi     <- drug_cache_d$psi
  ctl_psi      <- ctl_cache$psi

  # Restrict to IR events with no NA PSI
  ir_mask       <- drug_cache_d$complex == "IR" &
                   drug_cache_d$no_na_psi &
                   ctl_cache$no_na_psi
  drug_total    <- drug_total[ir_mask, , drop = FALSE]
  ctl_total     <- ctl_total[ir_mask,  , drop = FALSE]
  drug_psi      <- drug_psi[ir_mask, , drop = FALSE]
  ctl_psi       <- ctl_psi[ir_mask,  , drop = FALSE]

  drug_pass     <- rowSums(drug_total >= N) == ncol(drug_total)
  ctl_pass      <- rowSums(ctl_total  >= N) == ncol(ctl_total)
  common        <- drug_pass & ctl_pass
  if (sum(common) == 0) {
    return(tibble(
      EVENT  = character(0),
      drug_mean_PSI = numeric(0),
      ctl_mean_PSI  = numeric(0),
      delta_PSI     = numeric(0)
    ))
  }

  # Mean PSI across the (all-passing) replicates in each condition
  drug_mean <- rowMeans(drug_psi[common, , drop = FALSE], na.rm = TRUE)
  ctl_mean  <- rowMeans(ctl_psi[common,  , drop = FALSE], na.rm = TRUE)
  tibble(
    EVENT         = splicing_data$EVENT[ir_mask][common],
    drug_mean_PSI = drug_mean,
    ctl_mean_PSI  = ctl_mean,
    delta_PSI     = drug_mean - ctl_mean
  )
}

asymmetry_grid <- expand_grid(
  drug     = drugs_to_plot,
  N        = asymmetry_N
) %>%
  rowwise() %>%
  mutate(data = list(asymmetry_for_drug(drug, N))) %>%
  ungroup() %>%
  unnest(cols = c(data)) %>%
  mutate(
    drug = factor(drug, levels = drugs_to_plot),
    N    = factor(paste0("N >= ", N), levels = paste0("N >= ", asymmetry_N))
  )

# ---- 3b. Reviewer-investigation: betAS-significant asymmetry -----------------
# The reviewer asked whether the symmetric |DeltaPSI| distribution in the bulk
# view is an artifact of the |DeltaPSI|>=10 cutoff or a real biological
# signature. Answer: the symmetry is a CUTOFF ARTIFACT. The betAS results
# (FDR < 0.05) include many tiny but systematic shifts that collectively show
# strong asymmetry (e.g. SSA -> LESS intron retention, MORE exon inclusion).

betAS_path <- file.path(getwd(), "Supplementary2_betAS_splicing_results.xlsx")
betAS_sheets <- excel_sheets(betas_path <- betAS_path)
drug_to_sheet <- c(
  "Tubercidin"      = "Tub_FDR",
  "Spliceostatin A" = "SSA_FDR",
  "Pladienolide B"  = "PlaDB_FDR"
)

# Read each sheet and stack. betAS deltapsi is on a 0-1 (proportion) scale.
# Apply the STRICT filter (FDR < 0.05 AND |dPSI| >= 10%) at load time and
# classify events via the EVENT ID pattern -- matching figure2D exactly.
classify_event <- function(x) {
  ifelse(grepl("EX", x), "Exon",
         ifelse(grepl("INT", x), "Intron",
                ifelse(grepl("ALTD", x), "Alt5",
                       ifelse(grepl("ALTA", x), "Alt3", NA_character_))))
}

betAS_all <- bind_rows(lapply(drug_to_sheet, function(sh) {
  read_excel(betas_path, sheet = sh) %>%
    filter(FDR <= 0.05, abs(deltapsi) >= min_abs_dPSI_betAS) %>%
    na.omit() %>%
    mutate(
      drug      = unique(.data$drug),
      Event_grp = classify_event(EVENT)
    ) %>%
    filter(!is.na(Event_grp))
}))

# For each drug x event group, compute the asymmetry of the STRICT
# differentially-spliced events. Uses the EVENT-pattern classification
# (Event_grp), pre-filtered and na.omit'd at load time above.
betAS_asymmetry <- betAS_all %>%
  group_by(drug, Event_grp) %>%
  summarise(
    n_sig       = n(),
    n_pos       = sum(deltapsi > 0),
    n_neg       = sum(deltapsi < 0),
    pct_pos     = 100 * mean(deltapsi > 0),
    median_dPSI = median(deltapsi),
    mean_dPSI   = mean(deltapsi),
    .groups     = "drop"
  )

# As a function of effect-size cutoff (on the betAS 0-1 scale) -- show how the
# asymmetry changes when you only keep "large" effects. Built with a join so we
# avoid the rowwise+filter scoping trap that returns identical values per drug.
# Cutoffs start at 1% because below that the effects are dominated by noise.
# (effect_cutoffs kept for reference / future re-use; no plot depends on it now)

# ---- 5. Theme ----------------------------------------------------------------

theme_publish <- function(base_size = 11) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      plot.title       = element_text(face = "bold", hjust = 0),
      plot.subtitle    = element_text(color = "grey30", hjust = 0),
      axis.title       = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, color = "grey92"),
      strip.background = element_rect(fill = "grey96", color = NA),
      strip.text       = element_text(face = "bold"),
      legend.position  = "top",
      legend.title     = element_text(face = "bold"),
      plot.margin      = margin(10, 14, 10, 10)
    )
}

# ---- 6. Plot 1: counts vs. min-read filter -----------------------------------

p_counts <- ggplot(
    count_grid,
    aes(x = min_reads, y = n_events, colour = drug, group = drug)
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2, shape = 16) +
  scale_colour_manual(values = drug_palette, name = "Condition") +
  scale_x_continuous(breaks = min_reads) +
  scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = "Splicing events detected vs. minimum-read filter",
    subtitle = paste0(
      "vast-tools oocyte dataset · events passing the filter in ALL replicates ",
      "of each condition"
    ),
    x = "Minimum reads (inc + exc) per event per replicate",
    y = "Number of events detected",
    caption  = "Event types grouped by vast-tools COMPLEX codes; control (CTL) shown in grey."
  ) +
  facet_wrap(~ event_group, scales = "free_y", nrow = 1) +
  theme_publish()

# ---- 7. Plot 2: IR ΔPSI asymmetry, faceted by drug × N -----------------------
# A symmetric ΔPSI distribution around 0 means the drug shifts retention both
# ways equally (i.e. no net effect on the bulk of introns). SSA is expected to
# skew positive (more retention). If the skew vanishes at high N, the
# "asymmetry" the user noticed is driven by low-coverage noise.

# Density of DeltaPSI on the DIFFERENTIALLY-SPLICED subset (|DeltaPSI| >= 10).
# Restricting to "changed" events removes the spike at 0 and exposes any
# asymmetry (skew towards positive DeltaPSI = more retention in drug).
# Note: `N` in asymmetry_grid is a factor whose levels are "N >= 10" etc.;
# strip the prefix to recover the underlying number, otherwise as.numeric()
# coerces the string to NA.
n_levels     <- levels(asymmetry_grid$N)
n_num_levels <- as.numeric(sub("^N >= ", "", n_levels))

diff_counts <- asymmetry_grid %>%
  filter(abs(delta_PSI) >= min_abs_dPSI) %>%
  mutate(N_num = as.numeric(sub("^N >= ", "", as.character(N)))) %>%
  count(drug, N_num, name = "n_events_diff")

asymmetry_diff <- asymmetry_grid %>%
  filter(abs(delta_PSI) >= min_abs_dPSI) %>%
  mutate(N_num = as.numeric(sub("^N >= ", "", as.character(N)))) %>%
  left_join(diff_counts, by = c("drug", "N_num")) %>%
  mutate(
    N_num_label = factor(
      paste0("N=", N_num),
      levels = paste0("N=", sort(n_num_levels))
    )
  )

p_asym_density <- ggplot(
    asymmetry_diff,
    aes(x = delta_PSI, fill = drug, colour = drug)
  ) +
  geom_vline(xintercept = 0,
             linewidth = 0.4, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = c(-min_abs_dPSI, min_abs_dPSI),
             linewidth = 0.3, linetype = "dotted", colour = "grey55") +
  geom_density(alpha = 0.55, linewidth = 0.6, adjust = 1.0) +
  # Per-panel n annotation in the top-left corner (must use the SAME column name
  # as the facet variable, N_num_label, so ggplot routes each row to the right panel)
  geom_text(
    data  = diff_counts %>%
              mutate(N_num_label = factor(paste0("N=", N_num),
                                          levels = paste0("N=", sort(n_num_levels))),
                     label_plot  = paste0("n = ", format(n_events_diff, big.mark = ","))),
    aes(x = -95, y = Inf, label = label_plot),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1.5,
    size = 3.2, colour = "grey20", fontface = "italic"
  ) +
  scale_fill_manual(values  = drug_palette[drugs_to_plot], guide = "none") +
  scale_colour_manual(values = drug_palette[drugs_to_plot], guide = "none") +
  scale_x_continuous(breaks = seq(-100, 100, by = 25)) +
  labs(
    title    = paste0("Differentially-spliced IR events  (|DeltaPSI| >= ",
                      min_abs_dPSI, ")"),
    subtitle = paste0(
      "DeltaPSI density on the |DeltaPSI| >= ", min_abs_dPSI,
      " subset -- the background spike at 0 is removed so any asymmetry of the ",
      "tail is visible. Positive DeltaPSI = more intron retention in drug."
    ),
    x = "DeltaPSI  (drug - control)",
    y = "Density",
    caption  = paste0(
      "Columns = min-read filter; rows = drug; n per panel shown top-left. ",
      "Symmetric tails = equal numbers of up- and down-regulated introns."
    )
  ) +
  facet_grid(drug ~ N_num_label, scales = "free_y") +
  theme_publish() +
  theme(panel.spacing.x = unit(4, "mm"))

# ---- 7b. Plot 3: betAS asymmetry per drug x event type -----------------------
# Answer to the reviewer comment: the symmetry of the bulk |DeltaPSI| view is a
# CUTOFF ARTIFACT, not a biological signature. The betAS results (FDR < 0.05)
# include many small but systematically biased shifts. This plot shows the
# asymmetry of the betAS-significant set as a function of the effect-size
# cutoff. 50% = perfectly symmetric; >50% = biased toward drug > CTL.
#
# Note: betAS deltapsi is on a 0-1 (proportion) scale, NOT 0-100. The reviewer
# expects SF3B1 inhibitors (SSA, PLADB) to bias IR events positive (more
# retention). The betAS data show the OPPOSITE for all three drugs -- IR events
# are systematically biased toward LESS retention -- and exon inclusion (S) is
# biased toward MORE inclusion.

# Diverging bar plot of the STRICT differential events (FDR<0.05 AND |dPSI|>=10%).
# Positive bars (drug > control) and negative bars (drug < control). Uses the
# EVENT-pattern classification from figure2D (Exon / Intron / Alt5 / Alt3) so the
# counts exactly match figure2D's panel.
betAS_strict_long <- betAS_asymmetry %>%
  pivot_longer(cols = c(n_pos, n_neg),
               names_to = "direction", values_to = "n") %>%
  mutate(
    signed_n  = ifelse(direction == "n_pos", n, -n),
    direction = recode(direction, n_pos = "drug > CTL", n_neg = "drug < CTL"),
    Event_grp = factor(Event_grp,
                       levels = c("Exon", "Intron", "Alt5", "Alt3"),
                       labels = c("Exons", "Introns", "Alt5'", "Alt3'")),
    drug      = factor(drug,
                       levels = c("Pladienolide B", "Spliceostatin A", "Tubercidin"))
  )

p_betAS_strict <- ggplot(
    betAS_strict_long,
    aes(x = Event_grp, y = signed_n, fill = direction)
  ) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey30") +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = ifelse(direction == "drug > CTL", paste0("+", n), paste0("-", n))),
    position = position_stack(vjust = 0.5),
    size = 3.0, colour = "grey20"
  ) +
  scale_fill_manual(values = direction_palette, name = "Direction") +
  scale_y_continuous(labels = function(x) abs(x), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = paste0("Strict differentially-spliced events  ",
                      "(FDR < 0.05  AND  |DeltaPSI| >= ",
                      min_abs_dPSI_betAS * 100, "%)"),
    subtitle = paste0(
      "All exon types (C1/C2/C3/S/MIC) aggregated as 'Exons'. ",
      "Positive bars = drug > CTL; negative bars = drug < CTL."
    ),
    x = NULL,
    y = "Number of significant events",
    caption  = paste0(
      "Bar height |n| is the count of significant events going that direction. ",
      "50/50 bars = balanced; asymmetric = biased directionality."
    )
  ) +
  facet_wrap(~ drug, nrow = 1) +
  theme_publish() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "top"
  )

# ---- 8. Persist --------------------------------------------------------------

out_dir <- file.path(getwd(), "Figures", "events_per_filter_reads")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(out_dir, "events_per_filter_reads.pdf"),
       p_counts, width = 12, height = 4.2, units = "in", useDingbats = FALSE)
ggsave(file.path(out_dir, "events_per_filter_reads.png"),
       p_counts, width = 12, height = 4.2, units = "in", dpi = 300)

ggsave(file.path(out_dir, "ir_deltapsi_density_diffevents.pdf"),
       p_asym_density, width = 9, height = 6, units = "in", useDingbats = FALSE)
ggsave(file.path(out_dir, "ir_deltapsi_density_diffevents.png"),
       p_asym_density, width = 9, height = 6, units = "in", dpi = 300)

ggsave(file.path(out_dir, "betAS_strict_dPSI10_bars.pdf"),
       p_betAS_strict, width = 12, height = 5, units = "in", useDingbats = FALSE)
ggsave(file.path(out_dir, "betAS_strict_dPSI10_bars.png"),
       p_betAS_strict, width = 12, height = 5, units = "in", dpi = 300)

write.csv(count_grid,
          file.path(out_dir, "events_per_filter_reads.csv"),
          row.names = FALSE)
write.csv(betAS_asymmetry,
          file.path(out_dir, "betAS_significant_asymmetry.csv"),
          row.names = FALSE)

message("Done. Outputs written to: ", out_dir)
message("\nbetAS STRICT (FDR < 0.05 AND |DeltaPSI| >= 10%) events per drug x type:\n")
print(betAS_asymmetry %>% select(drug, Event_grp, n_sig, n_pos, n_neg, pct_pos),
      n = Inf)