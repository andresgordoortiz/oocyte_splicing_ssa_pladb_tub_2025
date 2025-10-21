# Single-file application launcher
# This file allows running the app with shiny::runApp("shiny_app/app.R")

library(shiny)

# Source the global, UI, and server files
source("global.R", local = TRUE)
source("ui.R", local = TRUE)
source("server.R", local = TRUE)

# Run the application
shinyApp(ui = ui, server = server)
