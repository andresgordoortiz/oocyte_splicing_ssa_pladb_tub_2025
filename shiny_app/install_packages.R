# ============================================================================
# Installation Script for Drug Splicing Effects Explorer
# ============================================================================

cat("\n")
cat("╔════════════════════════════════════════════════╗\n")
cat("║  Drug Splicing Effects Explorer               ║\n")
cat("║  Package Installation                          ║\n")
cat("╚════════════════════════════════════════════════╝\n")
cat("\n")

# ============================================================================
# CRAN Packages
# ============================================================================

cran_packages <- c(
  "shiny",
  "shinydashboard",
  "shinyWidgets",
  "DT",
  "dplyr",
  "ggplot2",
  "plotly",
  "readr",
  "tidyr",
  "ggrepel",
  "scales",
  "viridis",
  "gridExtra",
  "patchwork",
  "writexl"
)

cat("Installing CRAN packages...\n")
cat("─────────────────────────────────────────────────\n")

new_packages <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]

if (length(new_packages) > 0) {
  cat("Installing:", paste(new_packages, collapse = ", "), "\n\n")
  install.packages(new_packages, dependencies = TRUE, repos = "https://cran.rstudio.com/")
} else {
  cat("All CRAN packages already installed!\n\n")
}

# Verify CRAN installations
cat("Verifying CRAN packages:\n")
for (pkg in cran_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("  ✓", pkg, "\n")
  } else {
    cat("  ✗", pkg, "FAILED\n")
  }
}

# ============================================================================
# Bioconductor Packages
# ============================================================================

cat("\n")
cat("Installing Bioconductor packages...\n")
cat("─────────────────────────────────────────────────\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager", repos = "https://cran.rstudio.com/")
}

bioc_packages <- c(
  "ComplexHeatmap",
  "circlize"
)

new_bioc <- bioc_packages[!(bioc_packages %in% installed.packages()[,"Package"])]

if (length(new_bioc) > 0) {
  cat("Installing:", paste(new_bioc, collapse = ", "), "\n\n")
  BiocManager::install(new_bioc, update = FALSE, ask = FALSE)
} else {
  cat("All Bioconductor packages already installed!\n\n")
}

# Verify Bioconductor installations
cat("Verifying Bioconductor packages:\n")
for (pkg in bioc_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("  ✓", pkg, "\n")
  } else {
    cat("  ✗", pkg, "FAILED\n")
  }
}

# ============================================================================
# Summary
# ============================================================================

cat("\n")
cat("═════════════════════════════════════════════════\n")
cat("Installation Complete!\n")
cat("═════════════════════════════════════════════════\n")
cat("\n")

# Check if all packages are installed
all_packages <- c(cran_packages, bioc_packages)
failed_packages <- all_packages[!(all_packages %in% installed.packages()[,"Package"])]

if (length(failed_packages) == 0) {
  cat("✓ All packages successfully installed!\n")
  cat("\n")
  cat("You can now run the application:\n")
  cat("  shiny::runApp('shiny_app')\n")
  cat("\n")
} else {
  cat("✗ The following packages failed to install:\n")
  cat("  ", paste(failed_packages, collapse = ", "), "\n")
  cat("\n")
  cat("Please install them manually:\n")
  for (pkg in failed_packages) {
    if (pkg %in% bioc_packages) {
      cat("  BiocManager::install('", pkg, "')\n", sep = "")
    } else {
      cat("  install.packages('", pkg, "')\n", sep = "")
    }
  }
}

cat("\n")
