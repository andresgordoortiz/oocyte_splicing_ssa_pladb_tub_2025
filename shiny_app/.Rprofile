# ============================================================================
# .Rprofile for Shiny Application
# ============================================================================

# This file is loaded when the application starts
# Use it to set options and configure the environment

# Set CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Set Bioconductor repositories
options(BioC_mirror = "https://bioconductor.org")

# Increase memory limit (adjust as needed)
if (Sys.info()["sysname"] == "Windows") {
  memory.limit(size = 16000)
}

# Set encoding
options(encoding = "UTF-8")

# Disable scientific notation
options(scipen = 999)

# Set number of digits
options(digits = 4)

# Set timeout for downloads (for package installation)
options(timeout = 600)

# Configure plotting
options(
  bitmapType = "cairo",
  device = "png"
)

# Shiny-specific options
options(
  shiny.maxRequestSize = 100*1024^2,  # 100 MB max upload
  shiny.trace = FALSE,
  shiny.error = function() {
    logging::logerror(sys.calls())
  }
)

# Development vs Production mode
if (Sys.getenv("R_CONFIG_ACTIVE") == "production") {
  options(
    shiny.sanitize.errors = TRUE,
    shiny.fullstacktrace = FALSE
  )
} else {
  options(
    shiny.sanitize.errors = FALSE,
    shiny.fullstacktrace = TRUE
  )
}

# Load functions quietly
options(quietly = TRUE)

# Set working directory if needed
# Uncomment if you need to ensure working directory
# if (interactive()) {
#   setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# }

message("Application environment configured")
