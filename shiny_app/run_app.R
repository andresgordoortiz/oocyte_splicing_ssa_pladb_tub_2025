#!/usr/bin/env Rscript

# ============================================================================
# Run Drug Splicing Effects Explorer
# ============================================================================

cat("\n")
cat("╔════════════════════════════════════════════════╗\n")
cat("║  Drug Splicing Effects Explorer               ║\n")
cat("║  Starting application...                       ║\n")
cat("╚════════════════════════════════════════════════╝\n")
cat("\n")

# Check if we're in the right directory
if (!file.exists("app.R")) {
  stop(
    "Error: app.R not found.\n",
    "Please run this script from the shiny_app directory:\n",
    "  Rscript run_app.R\n"
  )
}

# Load shiny package
if (!requireNamespace("shiny", quietly = TRUE)) {
  cat("Shiny package not found. Installing...\n")
  install.packages("shiny", repos = "https://cran.rstudio.com/")
}

# Run the application
cat("Loading application...\n")
cat("\n")
cat("The application will open in your default browser.\n")
cat("Press Ctrl+C to stop the server.\n")
cat("\n")

shiny::runApp(
  appDir = ".",
  launch.browser = TRUE,
  host = "0.0.0.0",
  port = 3838
)
