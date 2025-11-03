#!/usr/bin/env Rscript
# Setup script for betAS local installation
# Run this script to ensure betAS is properly installed from the local submodule

cat("=== betAS Local Setup ===\n")

# Check if renv is available
if (!require("renv", quietly = TRUE)) {
  stop("renv is not installed. Please install renv first: install.packages('renv')")
}

# Install betAS from local submodule
cat("Installing betAS from local submodule...\n")
renv::install("./betAS", prompt = FALSE)

# Test installation
cat("Testing betAS installation...\n")
library(betAS)

# Check if getDataset is available
if (exists("getDataset")) {
  cat("✓ getDataset function is available\n")
  
  # Test with internal data
  tryCatch({
    testData <- getDataset(pathTables = NULL, tool = "vast-tools")
    cat("✓ getDataset test successful - data dimensions:", dim(testData), "\n")
  }, error = function(e) {
    cat("✗ getDataset test failed:", e$message, "\n")
  })
} else {
  cat("✗ getDataset function not found\n")
}

# Check if filterEvents is available
if (exists("filterEvents")) {
  cat("✓ filterEvents function is available\n")
} else {
  cat("✗ filterEvents function not found\n")
}

# List all exported functions
cat("\nAvailable betAS functions:\n")
print(ls("package:betAS"))

cat("\n=== Setup Complete ===\n")
cat("You can now use: library(betAS); getDataset(...)\n")