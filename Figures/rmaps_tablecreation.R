# Export significant and non-significant events into tab-separated files
# For each treatment (tub, pladb, ssa) the script will produce files by event type
# and by significance class (up, down, nonsig).
# Expected objects / CSVs:
# - tub_fdr_df, pladb_fdr_df, ssa_fdr_df  (or tub_fdr.csv, pladb_fdr.csv, ssa_fdr.csv)
# - event_info (or event_info.csv) containing columns: EVENT, CO_C1, CO_A, CO_C2, REF_CO

library(dplyr)
library(stringr)
library(readr)

# --- helper functions -------------------------------------------------------
# parse a single coordinate string like "chr16:18797349-18797399" -> list(chr, start, end)
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

# parse strand from REF_CO like "chr16:18797349,18795825-18795944,18795226:-" -> + or -
parse_strand <- function(refco) {
  if (is.na(refco) || refco == "") return(NA_character_)
  # take last character after last ':' which should be + or -
  m <- str_match(refco, ":([+-])$")
  if (!is.na(m[1,2])) return(m[1,2])
  # fallback: last character
  lastchar <- substr(refco, nchar(refco), nchar(refco))
  if (lastchar %in% c("+", "-")) return(lastchar)
  NA_character_
}

# build final table given joined df containing EVENT and event_info columns
build_table <- function(joined_df) {
  if (nrow(joined_df) == 0) return(tibble::tibble())
  out <- joined_df %>%
    rowwise() %>%
    mutate(
      .A = list(parse_coord_range(CO_A)),
      .C1 = list(parse_coord_range(CO_C1)),
      .C2 = list(parse_coord_range(CO_C2)),
      chr = .A$chr,
      exonStart = .A$start,
      exonEnd = .A$end,
      firstExonStart = .C1$start,
      firstExonEnd   = .C1$end,
      secondExonStart = .C2$start,
      secondExonEnd   = .C2$end,
      strand = parse_strand(REF_CO)
    ) %>%
    ungroup() %>%
    select(chr, strand, exonStart, exonEnd, firstExonStart, firstExonEnd, secondExonStart, secondExonEnd)
  out
}

# --- load data (use existing objects if present, otherwise try to read CSVs) ---
if (!exists("tub_fdr_df")) {
  if (file.exists("tub_fdr.csv")) tub_fdr_df <- read_csv("tub_fdr.csv")[,-1]
  else stop("tub_fdr_df not found and tub_fdr.csv not present in working directory")
}
if (!exists("pladb_fdr_df")) {
  if (file.exists("pladb_fdr.csv")) pladb_fdr_df <- read_csv("pladb_fdr.csv")[,-1]
  else stop("pladb_fdr_df not found and pladb_fdr.csv not present in working directory")
}
if (!exists("ssa_fdr_df")) {
  if (file.exists("ssa_fdr.csv")) ssa_fdr_df <- read_csv("ssa_fdr.csv")[,-1]
  else stop("ssa_fdr_df not found and ssa_fdr.csv not present in working directory")
}

if (!exists("event_info")) {
  if (file.exists("EVENT_INFO-mm10.tab")) event_info <- read.delim("EVENT_INFO-mm10.tab")
  else stop("event_info not found and event_info.csv not present in working directory")
}

# Ensure EVENT column exists in event_info
if (!"EVENT" %in% names(event_info)) stop("event_info must contain column 'EVENT'")

# Create differential sets if they don't already exist (user mentioned they have them)
if (!exists("differential_tub")) {
  differential_tub <- tub_fdr_df %>% filter(FDR <= 0.05 & abs(deltapsi) >= 0.1) %>% na.omit()
}
if (!exists("differential_pladb")) {
  differential_pladb <- pladb_fdr_df %>% filter(FDR <= 0.05 & abs(deltapsi) >= 0.1) %>% na.omit()
}
if (!exists("differential_ssa")) {
  differential_ssa <- ssa_fdr_df %>% filter(FDR <= 0.05 & abs(deltapsi) >= 0.1) %>% na.omit()
}

# Non-significant = full set minus differential events
nonsig_tub <- anti_join(tub_fdr_df, differential_tub, by = "EVENT")
nonsig_pladb <- anti_join(pladb_fdr_df, differential_pladb, by = "EVENT")
nonsig_ssa <- anti_join(ssa_fdr_df, differential_ssa, by = "EVENT")

# Pack into a list for iteration
treatments <- list(
  tub = list(full = tub_fdr_df, diff = differential_tub, nonsig = nonsig_tub),
  pladb = list(full = pladb_fdr_df, diff = differential_pladb, nonsig = nonsig_pladb),
  ssa = list(full = ssa_fdr_df, diff = differential_ssa, nonsig = nonsig_ssa)
)

# Definitions for event types (patterns in EVENT column)
event_types <- list(
  EX = "EX",
  INT = "INT",
  ALTD = "ALTD", # A5SS
  ALTA = "ALTA"  # A3SS
)

# significance classes and how to get the corresponding EVENT lists
# diff_up: events in differential with deltapsi > 0
# diff_down: deltapsi < 0
# nonsig: the nonsig df

for (trname in names(treatments)) {
  cat("Processing treatment:", trname, "\n")
  t <- treatments[[trname]]
  
  # compute lists of EVENTs
  diff_df <- t$diff
  up_events <- diff_df %>% filter(deltapsi > 0) %>% pull(EVENT) %>% unique()
  down_events <- diff_df %>% filter(deltapsi < 0) %>% pull(EVENT) %>% unique()
  nonsig_events <- t$nonsig %>% pull(EVENT) %>% unique()
  
  classes <- list(
    up = up_events,
    down = down_events,
    nonsig = nonsig_events
  )
  
  for (cls in names(classes)) {
    evlist <- classes[[cls]]
    if (length(evlist) == 0) {
      cat("  - No events for class:", cls, "in treatment:", trname, "\n")
      next
    }
    
    # subset event_info to these events
    sub_info <- event_info %>% filter(EVENT %in% evlist)
    
    for (etype in names(event_types)) {
      pat <- event_types[[etype]]
      sel <- sub_info %>% filter(str_detect(EVENT, fixed(pat)))
      if (nrow(sel) == 0) {
        cat("    - No", etype, "events for class", cls, "(treatment", trname, ")\n")
        next
      }
      
      out_table <- build_table(sel)
      
      # some rows might have NA chr if CO_A missing; drop those (can't export coordinates)
      out_table <- out_table %>% filter(!is.na(chr) & !is.na(exonStart) & !is.na(exonEnd))
      
      if (nrow(out_table) == 0) {
        cat("    - After parsing coordinates, nothing to write for", etype, cls, trname, "\n")
        next
      }
      
      fname <- paste0(trname, "_", etype, "_", cls, ".txt")
      write_tsv(out_table, fname)
      cat("    -> Wrote", nrow(out_table), "rows to", fname, "\n")
    }
  }
}

cat("All done.\n")
