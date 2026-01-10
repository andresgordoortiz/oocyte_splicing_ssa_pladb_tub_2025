# LOAD LIBRARIES
library(betAS)
library(dplyr)
library(readxl)
library(writexl)

# 1. IMPORT DATA
inclusion_file <- "Preprocessing/ssa_pladb_tub_oocyte_INCLUSION_LEVELS_FULL-mm10.tab"
splicing_data <- getDataset(pathTables = inclusion_file, tool = "vast-tools")

# 2. SUBSET SAMPLES & FILTER EVENTS
# Split data into the two experimental batches, since they were sequenced in different
# days and therefore must be filtered separetly 

pladb_data    <- splicing_data[, 1:18]
ssa_tub_data  <- splicing_data[, c(1:6, 19:36)]

# Filter for events with sufficient coverage (N=10)
pladb_events   <- filterEvents(getEvents(pladb_data, tool = "vast-tools"), N = 10)
ssa_tub_events <- filterEvents(getEvents(ssa_tub_data, tool = "vast-tools"), N = 10)

# 3. DEFINE SAMPLE COLUMNS
# Pladienolide B Comparison
cols_ctrl_pladb <- 7:9
cols_pladb      <- 10:12

# Spliceostatin A & Tubercidin Comparison
cols_ctrl_common <- 7:9
cols_tub         <- 10:12
cols_ssa         <- 13:15

# 4. RUN BETAS FDR CALCULATIONS
set.seed(42) # Ensure reproducibility for simulations

message("Calculating FDR for Tubercidin...")
tub_fdr_df <- prepareTableVolcanoFDR(
  psitable = ssa_tub_events$PSI, qualtable = ssa_tub_events$Qual,
  colsA = cols_ctrl_common, colsB = cols_tub,
  maxDevTable = maxDevSimulationN100, nsim = 1000, seed = TRUE, npoints = 500
)

message("Calculating FDR for Pladienolide B...")
pladb_fdr_df <- prepareTableVolcanoFDR(
  psitable = pladb_events$PSI, qualtable = pladb_events$Qual,
  colsA = cols_ctrl_pladb, colsB = cols_pladb,
  maxDevTable = maxDevSimulationN100, nsim = 1000, seed = TRUE
)

message("Calculating FDR for Spliceostatin A...")
ssa_fdr_df <- prepareTableVolcanoFDR(
  psitable = ssa_tub_events$PSI, qualtable = ssa_tub_events$Qual,
  colsA = cols_ctrl_common, colsB = cols_ssa,
  maxDevTable = maxDevSimulationN100, nsim = 1000, seed = TRUE
)

# 5. COMBINE AND SAVE TO EXCEL

sheets <- list(
  "Tub_FDR"   = as.data.frame(tub_fdr_df),
  "PlaDB_FDR" = as.data.frame(pladb_fdr_df),
  "SSA_FDR"   = as.data.frame(ssa_fdr_df)
)

write_xlsx(sheets, path = "Preprocessing/betAS_out/Supplementary2_betAS_splicing_results.xlsx")

message("Done! Results saved to Preprocessing/betAS_out/Supplementary2_betAS_splicing_results.xlsx")