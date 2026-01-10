library(tidyverse)
library(scales)
library(ggpubr)
library(rstatix)
library(readxl)


# ---- 1. Import from the combined Excel file ----
excel_file <- "Preprocessing/betAS_out/Supplementary2_betAS_splicing_results.xlsx"

# Read each sheet (using the names we set in the previous step)
tub_fdr_df   <- read_excel(excel_file, sheet = "Tub_FDR")
pladb_fdr_df <- read_excel(excel_file, sheet = "PlaDB_FDR")
ssa_fdr_df   <- read_excel(excel_file, sheet = "SSA_FDR")

# Base colors
ssa_color   <- "#706993"
pladb_color <- "#70A0AF"
tub_color   <- "#A0C1B9"
base_colors <- tibble(dataset = c("tub", "pladb", "ssa"),
                      base_hex = c(tub_color, pladb_color, ssa_color))

# Filter differentials
differential_tub   <- na.omit(tub_fdr_df[tub_fdr_df$FDR <= 0.05 & abs(tub_fdr_df$deltapsi) >= 0.1,])
differential_pladb <- na.omit(pladb_fdr_df[pladb_fdr_df$FDR <= 0.05 & abs(pladb_fdr_df$deltapsi) >= 0.1,])
differential_ssa   <- na.omit(ssa_fdr_df[ssa_fdr_df$FDR <= 0.05 & abs(ssa_fdr_df$deltapsi) >= 0.1,])

# Protein impact lookup
# Define the filename
file_path <- "PROT_IMPACT-mm10-v2.3.tab"

# Check if the file exists
if (!file.exists(file_path)) {
  stop(paste("File '", file_path, "' not found. ",
             "Please download it from the VASTDB website: https://vastdb.crg.eu/wiki/Main_Page", 
             sep = ""))
}

# Protein impact lookup (only runs if file exists)
protein_impact <- read.table(file_path, sep = "\t",
                             header = TRUE, stringsAsFactors = FALSE, 
                             fill = TRUE, quote = "") %>%
  dplyr::rename("EVENT" = EventID)



delta_col <- "deltapsi"
event_col <- "EVENT"
onto_col  <- "ONTO"

simplify_onto <- function(onto_vec) {
  onto_vec[is.na(onto_vec)] <- ""
  out <- rep("Unknown/noncoding", length(onto_vec))
  out[grepl("inclusion", onto_vec, ignore.case = TRUE)] <- "ORF disrupted upon inclusion"
  out[grepl("exclusion", onto_vec, ignore.case = TRUE)]  <- "ORF disrupted upon exclusion"
  out[grepl("alternative|isoform", onto_vec, ignore.case = TRUE)] <- "Alternative isoform"
  out[grepl("\\bUTR\\b|UTR", onto_vec, ignore.case = TRUE)] <- "Regulatory"
  out
}

protein_impact_lookup <- protein_impact %>%
  dplyr::select(EVENT, !!sym(onto_col)) %>%
  dplyr::rename(ONTO = !!sym(onto_col)) %>%
  mutate(Impact = simplify_onto(ONTO)) %>%
  dplyr::select(EVENT, Impact)

# Combine datasets
differential_list <- list(tub = differential_tub, pladb = differential_pladb, ssa = differential_ssa)

combined <- imap_dfr(differential_list, ~ .x %>%
                       left_join(protein_impact_lookup, by = "EVENT") %>%
                       mutate(Impact = coalesce(Impact, "Unknown/noncoding"),
                              dataset = .y,
                              deltapsi = as.numeric(deltapsi)))

impact_levels <- c(
  "ORF disrupted upon inclusion",
  "ORF disrupted upon exclusion",
  "Alternative isoform",
  "Regulatory",
  "Unknown/noncoding"
)

combined <- combined %>%
  mutate(Impact = factor(Impact, levels = impact_levels),
         dataset = factor(dataset, levels = c("tub", "pladb", "ssa")))

# Summary counts & proportions
summary_counts <- combined %>%
  group_by(dataset, Impact) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(dataset) %>%
  mutate(total = sum(n),
         prop = n / total)

library(colorspace)

# Replace only this function: pale -> intense along the same hue
make_harmonic_palette_v2 <- function(base_hex, n) {
  # Validate n and base_hex
  if (is.null(n) || !is.finite(n) || n <= 0) return(character(0))
  n <- as.integer(n)
  base_hex <- if (!is.null(base_hex) && length(base_hex) > 0) trimws(as.character(base_hex)[1]) else NA_character_
  valid_hex <- function(x) is.character(x) && length(x) == 1 && !is.na(x) && grepl("^#[0-9A-Fa-f]{6}$", x)
  if (n == 1L) return(ifelse(valid_hex(base_hex), base_hex, "#808080"))
  n_colored <- n - 1L
  if (n_colored <= 0L) return("#808080")
  if (!valid_hex(base_hex)) {
    vals <- seq(0.9, 0.35, length.out = n_colored)
    cols <- vapply(vals, function(v) grDevices::rgb(v, v, v), character(1))
    return(c(cols, "#808080"))
  }
  
  # Convert hex -> rgb (3 x 1 matrix)
  rgb_mat <- tryCatch(grDevices::col2rgb(base_hex) / 255, error = function(e) NULL)
  if (is.null(rgb_mat) || !is.matrix(rgb_mat) || any(!is.finite(rgb_mat))) {
    vals <- seq(0.9, 0.35, length.out = n_colored)
    cols <- vapply(vals, function(v) grDevices::rgb(v, v, v), character(1))
    return(c(cols, "#808080"))
  }
  
  base_hsv <- tryCatch(grDevices::rgb2hsv(rgb_mat[1,], rgb_mat[2,], rgb_mat[3,]), error = function(e) NULL)
  if (is.null(base_hsv) || !is.matrix(base_hsv) || any(is.na(base_hsv))) {
    vals <- seq(0.9, 0.35, length.out = n_colored)
    cols <- vapply(vals, function(v) grDevices::rgb(v, v, v), character(1))
    return(c(cols, "#808080"))
  }
  
  # Extract HSV components reliably
  base_h <- as.numeric(base_hsv["h", 1])
  base_s <- as.numeric(base_hsv["s", 1])
  base_v <- as.numeric(base_hsv["v", 1])
  
  clamp <- function(x, lo = 0, hi = 1) pmax(lo, pmin(hi, x))
  
  # --- NEW: Pale -> Intense ---
  # Pale start: low saturation, high brightness (light pastel)
  s_start <- clamp(base_s * 0.18, 0.05, 0.45)    # quite desaturated for pale
  # Intense end: somewhat boosted saturation (but not oversaturated)
  s_end   <- clamp(base_s * 1.15 + 0.05, 0.25, 1)
  
  # Brightness: start brighter (paler look), end slightly lower (richer)
  v_start <- clamp(base_v * 1.06, 0.85, 1.00)    # bright / pale
  v_end   <- clamp(max(base_v, 0.60), 0.60, 0.98) # keep end reasonably bright (avoid dark)
  
  # Safety: make sure ranges are ordered sensibly
  if (s_end < s_start) { tmp <- s_start; s_start <- s_end; s_end <- tmp }
  if (v_start < v_end) { tmp <- v_start; v_start <- v_end; v_end <- tmp }
  
  # Build sequences: saturation increases, brightness slightly decreases
  s_vals <- seq(s_start, s_end, length.out = n_colored)
  v_vals <- seq(v_start, v_end, length.out = n_colored)
  
  # Generate colours keeping hue fixed
  colors <- mapply(function(s, v) grDevices::hsv(base_h, s, v),
                   s_vals, v_vals, SIMPLIFY = TRUE, USE.NAMES = FALSE)
  
  # Safety: if any color is accidentally near-black, nudge brightness and recompute
  is_blackish <- function(hex) {
    ok <- tryCatch({ m <- grDevices::col2rgb(hex)/255; max(m) < 0.06 }, error = function(e) TRUE)
    isTRUE(ok)
  }
  if (any(vapply(colors, is_blackish, logical(1)))) {
    v_vals <- pmax(v_vals, 0.55)
    colors <- mapply(function(s, v) grDevices::hsv(base_h, s, v),
                     s_vals, v_vals, SIMPLIFY = TRUE, USE.NAMES = FALSE)
  }
  
  # Append grey for Unknown/noncoding
  c(colors, "#808080")
}



# Grid for plotting
grid <- tibble(x = rep(1:10, times = 10), y = rep(10:1, each = 10), pos = 1:100)

# Map impacts to full dots (100 total)
prop_to_100dots_round <- function(props) {
  raw <- props * 100
  floored <- floor(raw)
  remainder <- 100 - sum(floored)
  if (remainder > 0) {
    fracs <- raw - floored
    ord <- order(fracs, decreasing = TRUE)
    for (i in seq_len(remainder)) floored[ord[i]] <- floored[ord[i]] + 1L
  }
  floored
}

# ensure base_colors are characters
base_colors <- base_colors %>% mutate(dataset = as.character(dataset),
                                      base_hex = as.character(base_hex))

plot_prep <- summary_counts %>%
  group_by(dataset) %>%
  arrange(Impact) %>%  
  summarise(Impacts = list(as.character(Impact)),
            Props   = list(prop),
            .groups = "drop") %>%
  left_join(base_colors, by = "dataset") %>%
  mutate(counts = map(Props, prop_to_100dots_round),
         tones  = map2(base_hex, Impacts, ~ make_harmonic_palette_v2(.x, length(.y))),
         impact_df = pmap(list(Impacts, Props, counts, tones),
                          function(Impacts, Props, counts, tones) {
                            tibble(Impact = Impacts, prop = Props, count = counts, tone = tones)
                          })) %>%
  select(dataset, impact_df) %>%
  unnest(cols = c(impact_df)) %>%
  mutate(tone = as.character(tone))


# Expand to one row per dot
full_rows <- plot_prep %>%
  filter(count > 0) %>%
  uncount(weights = count, .remove = FALSE) %>%
  mutate(type = "full") %>%
  select(dataset, Impact, prop, tone, type)

# Assign positions
full_rows <- full_rows %>%
  group_by(dataset) %>%
  mutate(pos = row_number()) %>%
  ungroup() %>%
  left_join(grid, by = c("pos" = "pos"))

# DEBUG - Check for issues
message("Number of rows in full_rows: ", nrow(full_rows))
message("Number of rows with NA x: ", sum(is.na(full_rows$x)))
message("Number of rows with NA y: ", sum(is.na(full_rows$y)))
message("Number of rows with NA tone: ", sum(is.na(full_rows$tone)))
message("Sample of tones: ", paste(head(unique(full_rows$tone), 10), collapse = ", "))

# Check if tones are valid hex colors
check_hex <- function(x) {
  grepl("^#[0-9A-Fa-f]{6}$", x)
}
invalid_tones <- full_rows$tone[!check_hex(full_rows$tone)]
if(length(invalid_tones) > 0) {
  message("WARNING: Invalid hex colors detected: ", paste(unique(invalid_tones), collapse = ", "))
}

# MAIN FIX: Use explicit fill mapping instead of identity scale
# Create a named vector for the color palette
all_tones <- unique(full_rows$tone)
names(all_tones) <- all_tones

# Alternative approach if tone column has issues
if(any(is.na(full_rows$tone)) || any(!check_hex(full_rows$tone))) {
  message("Fixing invalid tones...")
  full_rows <- full_rows %>%
    mutate(tone = ifelse(is.na(tone) | !check_hex(tone), "#808080", tone))
}

# Create the plot with explicit color mapping
p <- ggplot(full_rows) +
  geom_point(aes(x = x, y = y, fill = tone),
             shape = 21, size = 25, stroke = 0) +
  scale_fill_identity() +  # This should work if tones are valid hex codes
  facet_wrap(~ dataset, nrow = 1, strip.position = "top") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 14),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(8, 10, 8, 10),
        plot.title = element_text(face = "bold", hjust = 0.5))

# Add labels
labels_df <- full_rows %>%
  group_by(dataset, Impact) %>%
  summarise(median_y = median(y), prop = dplyr::first(prop), tone = dplyr::first(tone), .groups = "drop") %>%
  arrange(dataset, Impact) %>%
  mutate(pct_text = paste0(percent(prop, accuracy = 0.1)),
         name_text = Impact,
         # place labels a bit to the right of the grid maximum (grid x runs 1:10)
         label_x = max(grid$x) + 2.2)   # <- move label further right than before

# --- plot: increase panel spacing, enlarge right plot margin, allow drawing outside panels ---
p <- ggplot(full_rows) +
  geom_point(aes(x = x, y = y, fill = tone),
             shape = 21, size = 25, stroke = 0.3) +
  scale_fill_identity() +
  facet_wrap(~ dataset, nrow = 1, strip.position = "top") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 14),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    # more space between facets so labels sit in the gutter and don't overlap next plot
    panel.spacing = unit(6.0, "lines"),
    # larger right margin so long labels are not clipped by the edge of the canvas
    plot.margin = margin(8, 80, 8, 10),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  # draw text outside panels (requires ggplot2 that supports coord clip argument)
  geom_text(data = labels_df, aes(x = label_x, y = median_y + 0.28, label = name_text, color = tone),
            hjust = 0, size = 5, fontface = "bold") +
  geom_text(data = labels_df, aes(x = label_x, y = median_y - 0.28, label = pct_text, color = tone),
            hjust = 0, size = 8) +
  scale_color_identity() +
  # keep aspect ratio, extend limits a bit, and allow drawing outside the panel area
  coord_fixed(clip = "off", expand = FALSE, xlim = c(-1.8, max(grid$x) + 4.5), ylim = c(0.15, 10.95))

# Save
ggsave("Figures/Figure3A_functional_impact_dotplot.pdf", p, device = cairo_pdf, width = 40, height = 6, units = "in", dpi = 600)


# ---- 2. Export Differential Events with Protein Impact to Excel ----

# Ensure you have the library installed
# install.packages("writexl")
library(writexl)

# 1. Prepare the data for export
# We use the 'combined' dataframe which already has the protein impacts joined
export_list <- combined %>%
  # Remove the helper 'dataset' column and relocate 'Impact' next to 'EVENT'
  dplyr::select(EVENT, Impact, everything(), -dataset) %>%
  # Convert factors to characters for cleaner Excel display
  mutate(across(where(is.factor), as.character)) %>%
  # Split into a named list based on the drug/dataset names
  # (Since 'combined$dataset' was a factor, the sheets will be in order: tub, pladb, ssa)
  split(combined$dataset)

# 2. Define the output path
output_excel <- "Supplementary5_Differential_Splicing_Protein_Impact_Results.xlsx"

# 3. Write to Excel
write_xlsx(export_list, path = output_excel)

message("Excel file successfully created at: ", output_excel)