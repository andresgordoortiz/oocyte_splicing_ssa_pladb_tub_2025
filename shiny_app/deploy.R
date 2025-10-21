# ============================================================================
# Deployment Script for Drug Splicing Effects Explorer
# ============================================================================

# This script helps deploy the Shiny application to shinyapps.io or other platforms

library(rsconnect)

# ============================================================================
# CONFIGURATION
# ============================================================================

APP_NAME <- "drug-splicing-explorer"
APP_TITLE <- "Drug Splicing Effects Explorer"

# ============================================================================
# DEPLOYMENT TO SHINYAPPS.IO
# ============================================================================

deploy_to_shinyapps <- function(account_name, force_update = TRUE) {
  
  cat("===============================================\n")
  cat("Deploying to shinyapps.io\n")
  cat("===============================================\n\n")
  
  # Check if account is configured
  accounts <- rsconnect::accounts()
  
  if (nrow(accounts) == 0) {
    stop(
      "No shinyapps.io account configured.\n",
      "Please run:\n",
      "  rsconnect::setAccountInfo(name='your-account', token='your-token', secret='your-secret')\n",
      "Get your token from: https://www.shinyapps.io/admin/#/tokens"
    )
  }
  
  cat("Found account(s):", paste(accounts$name, collapse = ", "), "\n\n")
  
  # Deploy
  cat("Deploying application...\n")
  
  tryCatch({
    rsconnect::deployApp(
      appDir = getwd(),
      appName = APP_NAME,
      appTitle = APP_TITLE,
      account = account_name,
      forceUpdate = force_update,
      launch.browser = TRUE
    )
    
    cat("\n✓ Deployment successful!\n")
    cat("Your app should open in a browser shortly.\n")
    
  }, error = function(e) {
    cat("\n✗ Deployment failed!\n")
    cat("Error:", e$message, "\n")
  })
}

# ============================================================================
# PRE-DEPLOYMENT CHECKS
# ============================================================================

check_deployment_readiness <- function() {
  
  cat("===============================================\n")
  cat("Pre-Deployment Checks\n")
  cat("===============================================\n\n")
  
  checks_passed <- TRUE
  
  # Check required files
  required_files <- c(
    "global.R",
    "ui.R",
    "server.R",
    "www/logo_crg.png",
    "www/cfrance_logo.png"
  )
  
  cat("1. Checking required files...\n")
  for (file in required_files) {
    if (file.exists(file)) {
      cat("  ✓", file, "\n")
    } else {
      cat("  ✗", file, "NOT FOUND\n")
      checks_passed <- FALSE
    }
  }
  
  # Check data files in parent directory
  cat("\n2. Checking data files...\n")
  data_files <- c(
    "../tub_fdr.csv",
    "../pladb_fdr.csv",
    "../ssa_fdr.csv",
    "../EVENT_INFO-mm10.tab"
  )
  
  for (file in data_files) {
    if (file.exists(file)) {
      cat("  ✓", basename(file), "\n")
    } else {
      cat("  ⚠", basename(file), "NOT FOUND (optional)\n")
    }
  }
  
  # Check package dependencies
  cat("\n3. Checking R packages...\n")
  required_packages <- c(
    "shiny", "shinydashboard", "shinyWidgets", "DT",
    "dplyr", "ggplot2", "plotly", "readr", "tidyr",
    "ComplexHeatmap", "circlize", "ggrepel", "scales",
    "viridis", "gridExtra", "patchwork", "writexl"
  )
  
  for (pkg in required_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      cat("  ✓", pkg, "\n")
    } else {
      cat("  ✗", pkg, "NOT INSTALLED\n")
      checks_passed <- FALSE
    }
  }
  
  cat("\n")
  if (checks_passed) {
    cat("✓ All critical checks passed!\n")
    cat("Ready for deployment.\n")
  } else {
    cat("✗ Some checks failed. Please fix issues before deploying.\n")
  }
  
  return(checks_passed)
}

# ============================================================================
# PACKAGE INSTALLATION HELPER
# ============================================================================

install_dependencies <- function() {
  
  cat("===============================================\n")
  cat("Installing Dependencies\n")
  cat("===============================================\n\n")
  
  # CRAN packages
  cran_packages <- c(
    "shiny", "shinydashboard", "shinyWidgets", "DT",
    "dplyr", "ggplot2", "plotly", "readr", "tidyr",
    "ggrepel", "scales", "viridis", "gridExtra", 
    "patchwork", "writexl", "rsconnect"
  )
  
  cat("Installing CRAN packages...\n")
  for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("  Installing", pkg, "...\n")
      install.packages(pkg, dependencies = TRUE)
    } else {
      cat("  ✓", pkg, "already installed\n")
    }
  }
  
  # Bioconductor packages
  cat("\nInstalling Bioconductor packages...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  bioc_packages <- c("ComplexHeatmap", "circlize")
  for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("  Installing", pkg, "...\n")
      BiocManager::install(pkg, update = FALSE)
    } else {
      cat("  ✓", pkg, "already installed\n")
    }
  }
  
  cat("\n✓ All dependencies installed!\n")
}

# ============================================================================
# CREATE DEPLOYMENT PACKAGE
# ============================================================================

create_deployment_package <- function(output_file = "drug_splicing_app.zip") {
  
  cat("===============================================\n")
  cat("Creating Deployment Package\n")
  cat("===============================================\n\n")
  
  # Files to include
  files_to_package <- c(
    "global.R",
    "ui.R",
    "server.R",
    "app.R",
    "README.md",
    "www"
  )
  
  # Create zip file
  cat("Creating zip file:", output_file, "\n")
  zip(output_file, files_to_package, flags = "-r9X")
  
  cat("✓ Package created successfully!\n")
  cat("File size:", round(file.size(output_file) / 1024^2, 2), "MB\n")
  
  return(output_file)
}

# ============================================================================
# MAIN DEPLOYMENT WORKFLOW
# ============================================================================

if (interactive()) {
  cat("\n")
  cat("╔════════════════════════════════════════════════╗\n")
  cat("║  Drug Splicing Effects Explorer - Deployment  ║\n")
  cat("╚════════════════════════════════════════════════╝\n")
  cat("\n")
  
  cat("Available commands:\n")
  cat("  1. check_deployment_readiness()   - Run pre-deployment checks\n")
  cat("  2. install_dependencies()         - Install required packages\n")
  cat("  3. deploy_to_shinyapps('account') - Deploy to shinyapps.io\n")
  cat("  4. create_deployment_package()    - Create deployment ZIP\n")
  cat("\n")
  
  cat("Quick deploy example:\n")
  cat("  check_deployment_readiness()\n")
  cat("  deploy_to_shinyapps('your-account-name')\n")
  cat("\n")
}

# Example usage (uncomment to run):
# check_deployment_readiness()
# install_dependencies()
# deploy_to_shinyapps("your-shinyapps-account")
