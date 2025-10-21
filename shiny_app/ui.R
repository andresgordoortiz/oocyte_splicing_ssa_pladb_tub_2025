# ============================================================================
# UI Definition for Oocyte Splicing Effects Explorer
# ============================================================================

ui <- dashboardPage(
  skin = "blue",  # Changed from purple to blue for professional look
  
  # ============================================================================
  # HEADER
  # ============================================================================
  dashboardHeader(
    title = "Oocyte Splicing Explorer",
    titleWidth = 300
  ),
  
  # ============================================================================
  # SIDEBAR
  # ============================================================================
  dashboardSidebar(
    width = 280,
    
    # Add custom CSS
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
      tags$style(HTML("
        .main-header .logo {
          font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
          font-size: 20px;
          font-weight: 600;
        }
      "))
    ),
    
    sidebarMenu(
      id = "sidebar",
      menuItem("Overview", tabName = "overview", icon = icon("home")),
      menuItem("Volcano Plots", tabName = "volcano", icon = icon("mountain")),
      menuItem("Event Analysis", tabName = "events", icon = icon("chart-bar")),
      menuItem("Gene Overlaps", tabName = "upset", icon = icon("project-diagram")),
      menuItem("Heatmaps", tabName = "heatmaps", icon = icon("th")),
      menuItem("GO Enrichment", tabName = "go", icon = icon("dna")),
      menuItem("Sequence Features", tabName = "sequence", icon = icon("microscope")),
      menuItem("RBP Analysis", tabName = "rbp", icon = icon("shapes")),
      menuItem("Data Tables", tabName = "tables", icon = icon("table")),
      menuItem("Downloads", tabName = "downloads", icon = icon("download")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    
    # Sidebar filters
    tags$hr(),
    tags$div(
      style = "padding: 15px;",
      h5("Global Filters", style = "font-weight: bold; margin-bottom: 15px;"),
      
      sliderInput(
        "fdr_threshold",
        "FDR Threshold:",
        min = 0.001,
        max = 0.1,
        value = 0.05,
        step = 0.001
      ),
      
      sliderInput(
        "deltapsi_threshold",
        "Δ PSI Threshold:",
        min = 0.05,
        max = 0.3,
        value = 0.1,
        step = 0.01
      ),
      
      pickerInput(
        "selected_drugs",
        "Select Drugs:",
        choices = c("Tubercidin", "Pladienolide B", "Spliceostatin A"),
        selected = c("Tubercidin", "Pladienolide B", "Spliceostatin A"),
        multiple = TRUE,
        options = list(
          `actions-box` = TRUE,
          `selected-text-format` = "count > 2"
        )
      ),
      
      pickerInput(
        "event_types",
        "Event Types:",
        choices = c("Exon", "Intron", "Alt5", "Alt3"),
        selected = c("Exon", "Intron", "Alt5", "Alt3"),
        multiple = TRUE,
        options = list(
          `actions-box` = TRUE
        )
      )
    )
  ),
  
  # ============================================================================
  # BODY
  # ============================================================================
  dashboardBody(
    # Custom CSS
    tags$head(
      tags$style(HTML("
        /* Main colors */
        .skin-purple .main-header .logo { background-color: #706993; }
        .skin-purple .main-header .navbar { background-color: #706993; }
        .skin-purple .main-sidebar { background-color: #34495e; }
        
        /* Content area */
        .content-wrapper { background-color: #f4f6f9; }
        
        /* Box styling */
        .box { 
          border-top: 3px solid #706993; 
          box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        
        /* Info boxes */
        .info-box { 
          box-shadow: 0 2px 4px rgba(0,0,0,0.1);
          border-radius: 3px;
        }
        
        /* Tabs */
        .nav-tabs-custom > .nav-tabs > li.active { border-top-color: #706993; }
        
        /* Value boxes */
        .small-box { border-radius: 3px; }
        .small-box.bg-aqua { background-color: #70A0AF !important; }
        .small-box.bg-purple { background-color: #706993 !important; }
        .small-box.bg-green { background-color: #A0C1B9 !important; }
        
        /* Custom styling for plots */
        .plot-container { 
          background: white; 
          padding: 15px; 
          border-radius: 5px;
          box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }
        
        /* DataTable styling */
        .dataTables_wrapper { padding: 10px; }
        
        /* Headers */
        h2, h3, h4 { color: #34495e; font-weight: 600; }
        
        /* Buttons */
        .btn-primary { background-color: #706993; border-color: #706993; }
        .btn-primary:hover { background-color: #5a5577; border-color: #5a5577; }
      "))
    ),
    
    tabItems(
      # ========================================================================
      # OVERVIEW TAB
      # ========================================================================
      tabItem(
        tabName = "overview",
        fluidRow(
          column(
            width = 12,
            box(
              width = NULL,
              title = tags$div(
                icon("flask", style = "margin-right: 10px;"),
                "Oocyte Splicing Effects Explorer"
              ),
              status = "primary",
              solidHeader = TRUE,
              tags$div(
                style = "font-size: 16px; line-height: 1.8;",
                p("Welcome to the", strong("Drug Splicing Effects Explorer"), "- an interactive platform 
                  for exploring the transcriptomic effects of three splicing-modulating drugs in mouse oocytes."),
                
                tags$hr(),
                
                h4(icon("pills"), " Study Drugs", style = "color: #706993; margin-top: 20px;"),
                tags$ul(
                  tags$li(strong(tags$span(style = "color: #A0C1B9;", "Tubercidin:")), 
                          " An adenosine analog that affects RNA processing"),
                  tags$li(strong(tags$span(style = "color: #70A0AF;", "Pladienolide B:")), 
                          " A selective inhibitor of SF3B1, a core component of the U2 snRNP"),
                  tags$li(strong(tags$span(style = "color: #706993;", "Spliceostatin A:")), 
                          " A potent splicing inhibitor targeting SF3B1")
                ),
                
                h4(icon("chart-line"), " Available Analyses", style = "color: #706993; margin-top: 20px;"),
                tags$ul(
                  tags$li(strong("Volcano Plots:"), " Visualize differential splicing events (ΔPSI vs FDR)"),
                  tags$li(strong("Event Analysis:"), " Breakdown of splicing event types and quantification"),
                  tags$li(strong("Gene Overlaps:"), " UpSet plots showing shared affected genes across drugs"),
                  tags$li(strong("Heatmaps:"), " PSI value clustering patterns"),
                  tags$li(strong("GO Enrichment:"), " Functional enrichment of affected genes"),
                  tags$li(strong("Sequence Features:"), " GC content and length analysis"),
                  tags$li(strong("RBP Analysis:"), " RNA-binding protein motif enrichment"),
                  tags$li(strong("Data Tables:"), " Interactive exploration of all events")
                ),
                
                h4(icon("filter"), " Using the Filters", style = "color: #706993; margin-top: 20px;"),
                p("Use the sidebar filters to customize your analysis:"),
                tags$ul(
                  tags$li(strong("FDR Threshold:"), " Control statistical stringency (default: 0.05)"),
                  tags$li(strong("ΔPSI Threshold:"), " Set minimum effect size (default: 0.1 = 10%)"),
                  tags$li(strong("Drug Selection:"), " Compare one or multiple drugs"),
                  tags$li(strong("Event Types:"), " Focus on specific splicing events")
                )
              )
            )
          )
        ),
        
        # Summary statistics
        fluidRow(
          valueBoxOutput("vbox_tub", width = 4),
          valueBoxOutput("vbox_pladb", width = 4),
          valueBoxOutput("vbox_ssa", width = 4)
        ),
        
        fluidRow(
          box(
            width = 12,
            title = "Quick Summary Statistics",
            status = "info",
            solidHeader = TRUE,
            DTOutput("summary_table")
          )
        )
      ),
      
      # ========================================================================
      # VOLCANO PLOTS TAB
      # ========================================================================
      tabItem(
        tabName = "volcano",
        fluidRow(
          box(
            width = 12,
            title = "Volcano Plot - Differential Splicing Events",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            
            fluidRow(
              column(
                width = 3,
                selectInput(
                  "volcano_drug",
                  "Select Drug:",
                  choices = c("Tubercidin", "Pladienolide B", "Spliceostatin A"),
                  selected = "Tubercidin"
                ),
                checkboxInput("volcano_labels", "Show Gene Labels", value = TRUE),
                sliderInput(
                  "volcano_top_n",
                  "Number of Labels:",
                  min = 5,
                  max = 50,
                  value = 20,
                  step = 5
                ),
                downloadButton("download_volcano", "Download Plot", class = "btn-primary")
              ),
              column(
                width = 9,
                div(class = "plot-container",
                    plotlyOutput("volcano_plot", height = "600px")
                )
              )
            )
          )
        ),
        
        fluidRow(
          box(
            width = 6,
            title = "Significant Events by Direction",
            status = "info",
            solidHeader = TRUE,
            plotOutput("volcano_direction_plot", height = "300px")
          ),
          box(
            width = 6,
            title = "Event Type Distribution",
            status = "info",
            solidHeader = TRUE,
            plotOutput("volcano_event_types_plot", height = "300px")
          )
        )
      ),
      
      # ========================================================================
      # EVENT ANALYSIS TAB
      # ========================================================================
      tabItem(
        tabName = "events",
        fluidRow(
          box(
            width = 12,
            title = "Splicing Event Quantification",
            status = "primary",
            solidHeader = TRUE,
            
            fluidRow(
              column(
                width = 3,
                radioButtons(
                  "event_view_type",
                  "View Type:",
                  choices = c("Absolute Counts" = "absolute",
                              "Proportions" = "proportion"),
                  selected = "absolute"
                ),
                checkboxInput("event_show_numbers", "Show Values", value = TRUE),
                downloadButton("download_events", "Download Plot", class = "btn-primary")
              ),
              column(
                width = 9,
                div(class = "plot-container",
                    plotOutput("events_barplot", height = "500px")
                )
              )
            )
          )
        ),
        
        fluidRow(
          box(
            width = 6,
            title = "Event Type Summary Table",
            status = "info",
            solidHeader = TRUE,
            DTOutput("events_summary_table")
          ),
          box(
            width = 6,
            title = "Event Distribution by Drug",
            status = "info",
            solidHeader = TRUE,
            selectInput(
              "events_pie_drug",
              "Select Drug:",
              choices = c("Tubercidin", "Pladienolide B", "Spliceostatin A"),
              selected = "Tubercidin"
            ),
            plotOutput("events_pie_chart", height = "400px")
          )
        )
      ),
      
      # ========================================================================
      # UPSET PLOTS TAB
      # ========================================================================
      tabItem(
        tabName = "upset",
        fluidRow(
          box(
            width = 12,
            title = "Gene Overlap Analysis",
            status = "primary",
            solidHeader = TRUE,
            
            fluidRow(
              column(
                width = 3,
                selectInput(
                  "upset_level",
                  "Analysis Level:",
                  choices = c("All Events" = "events",
                              "Genes Only" = "genes"),
                  selected = "genes"
                ),
                selectInput(
                  "upset_event_filter",
                  "Filter by Event Type:",
                  choices = c("All" = "all", "Exon", "Intron", "Alt5", "Alt3"),
                  selected = "all"
                ),
                downloadButton("download_upset", "Download Plot", class = "btn-primary")
              ),
              column(
                width = 9,
                div(class = "plot-container",
                    plotOutput("upset_plot", height = "600px")
                )
              )
            )
          )
        ),
        
        fluidRow(
          box(
            width = 6,
            title = "Overlap Statistics",
            status = "info",
            solidHeader = TRUE,
            DTOutput("upset_stats_table")
          ),
          box(
            width = 6,
            title = "Venn Diagram Representation",
            status = "info",
            solidHeader = TRUE,
            plotOutput("venn_plot", height = "400px")
          )
        )
      ),
      
      # ========================================================================
      # HEATMAPS TAB
      # ========================================================================
      tabItem(
        tabName = "heatmaps",
        fluidRow(
          box(
            width = 12,
            title = "PSI Value Heatmaps",
            status = "primary",
            solidHeader = TRUE,
            
            fluidRow(
              column(
                width = 3,
                selectInput(
                  "heatmap_drug",
                  "Select Drug:",
                  choices = c("Tubercidin", "Pladienolide B", "Spliceostatin A"),
                  selected = "Tubercidin"
                ),
                selectInput(
                  "heatmap_event_type",
                  "Event Type:",
                  choices = c("Exons" = "exon", "Introns" = "intron"),
                  selected = "exon"
                ),
                sliderInput(
                  "heatmap_max_events",
                  "Max Events to Show:",
                  min = 20,
                  max = 200,
                  value = 50,
                  step = 10
                ),
                selectInput(
                  "heatmap_cluster_method",
                  "Clustering Method:",
                  choices = c("Ward" = "ward.D2", "Complete" = "complete", "Average" = "average"),
                  selected = "ward.D2"
                ),
                downloadButton("download_heatmap", "Download Plot", class = "btn-primary")
              ),
              column(
                width = 9,
                div(class = "plot-container",
                    plotlyOutput("heatmap_plot", height = "700px")
                )
              )
            )
          )
        )
      ),
      
      # ========================================================================
      # GO ENRICHMENT TAB
      # ========================================================================
      tabItem(
        tabName = "go",
        fluidRow(
          box(
            width = 12,
            title = "Gene Ontology Enrichment Analysis",
            status = "primary",
            solidHeader = TRUE,
            
            fluidRow(
              column(
                width = 3,
                selectInput(
                  "go_drug",
                  "Select Drug:",
                  choices = c("Tubercidin", "Pladienolide B", "Spliceostatin A"),
                  selected = "Tubercidin"
                ),
                selectInput(
                  "go_ontology",
                  "Ontology:",
                  choices = c(
                    "Biological Process" = "BP",
                    "Molecular Function" = "MF",
                    "Cellular Component" = "CC"
                  ),
                  selected = "BP"
                ),
                sliderInput(
                  "go_top_n",
                  "Number of Terms:",
                  min = 5,
                  max = 30,
                  value = 15,
                  step = 5
                ),
                sliderInput(
                  "go_pval_cutoff",
                  "P-value Cutoff:",
                  min = 0.001,
                  max = 0.1,
                  value = 0.05,
                  step = 0.001
                ),
                downloadButton("download_go", "Download Plot", class = "btn-primary")
              ),
              column(
                width = 9,
                div(class = "plot-container",
                    plotOutput("go_plot", height = "600px")
                )
              )
            )
          )
        ),
        
        fluidRow(
          box(
            width = 12,
            title = "GO Terms Table",
            status = "info",
            solidHeader = TRUE,
            DTOutput("go_table")
          )
        )
      ),
      
      # ========================================================================
      # SEQUENCE FEATURES TAB
      # ========================================================================
      tabItem(
        tabName = "sequence",
        fluidRow(
          box(
            width = 12,
            title = "Sequence Feature Analysis",
            status = "primary",
            solidHeader = TRUE,
            
            tabBox(
              width = NULL,
              
              tabPanel(
                "GC Content",
                fluidRow(
                  column(
                    width = 3,
                    selectInput(
                      "gc_event_type",
                      "Event Type:",
                      choices = c("Exons" = "exon", "Introns" = "intron"),
                      selected = "exon"
                    ),
                    checkboxInput("gc_show_stats", "Show Statistics", value = TRUE),
                    downloadButton("download_gc", "Download", class = "btn-primary")
                  ),
                  column(
                    width = 9,
                    plotOutput("gc_plot", height = "500px")
                  )
                )
              ),
              
              tabPanel(
                "Length Distribution",
                fluidRow(
                  column(
                    width = 3,
                    selectInput(
                      "length_event_type",
                      "Event Type:",
                      choices = c("Exons" = "exon", "Introns" = "intron"),
                      selected = "exon"
                    ),
                    checkboxInput("length_show_stats", "Show Statistics", value = TRUE),
                    downloadButton("download_length", "Download", class = "btn-primary")
                  ),
                  column(
                    width = 9,
                    plotOutput("length_plot", height = "500px")
                  )
                )
              )
            )
          )
        )
      ),
      
      # ========================================================================
      # RBP ANALYSIS TAB
      # ========================================================================
      tabItem(
        tabName = "rbp",
        fluidRow(
          box(
            width = 12,
            title = "RNA-Binding Protein Motif Enrichment",
            status = "primary",
            solidHeader = TRUE,
            
            p(strong("Note:"), " This analysis shows enrichment of RBP binding motifs in regions 
              surrounding differentially spliced events. Analysis based on RNA maps."),
            
            fluidRow(
              column(
                width = 3,
                selectInput(
                  "rbp_drug",
                  "Select Drug:",
                  choices = c("Tubercidin", "Pladienolide B", "Spliceostatin A"),
                  selected = "Tubercidin"
                ),
                selectInput(
                  "rbp_event_type",
                  "Event Type:",
                  choices = c("Exons" = "ex", "Introns" = "int"),
                  selected = "ex"
                ),
                selectInput(
                  "rbp_direction",
                  "Splicing Direction:",
                  choices = c("Both" = "both", "Included/Retained" = "up", "Skipped/Removed" = "down"),
                  selected = "both"
                ),
                sliderInput(
                  "rbp_top_n",
                  "Number of Top RBPs:",
                  min = 10,
                  max = 50,
                  value = 20,
                  step = 5
                ),
                downloadButton("download_rbp", "Download", class = "btn-primary")
              ),
              column(
                width = 9,
                plotOutput("rbp_plot", height = "600px")
              )
            )
          )
        )
      ),
      
      # ========================================================================
      # DATA TABLES TAB
      # ========================================================================
      tabItem(
        tabName = "tables",
        fluidRow(
          box(
            width = 12,
            title = "Interactive Data Tables",
            status = "primary",
            solidHeader = TRUE,
            
            fluidRow(
              column(
                width = 3,
                selectInput(
                  "table_drug",
                  "Select Drug:",
                  choices = c("Tubercidin", "Pladienolide B", "Spliceostatin A"),
                  selected = "Tubercidin"
                ),
                radioButtons(
                  "table_filter",
                  "Show:",
                  choices = c("All Events" = "all",
                              "Significant Only" = "significant",
                              "Non-significant" = "nonsig"),
                  selected = "significant"
                ),
                downloadButton("download_table_csv", "Download CSV", class = "btn-primary"),
                tags$br(), tags$br(),
                downloadButton("download_table_excel", "Download Excel", class = "btn-success")
              ),
              column(
                width = 9,
                DTOutput("data_table")
              )
            )
          )
        )
      ),
      
      # ========================================================================
      # DOWNLOADS TAB
      # ========================================================================
      tabItem(
        tabName = "downloads",
        fluidRow(
          box(
            width = 12,
            title = "Bulk Downloads",
            status = "primary",
            solidHeader = TRUE,
            
            h4("Download Complete Datasets"),
            p("Download all analysis results for offline exploration."),
            
            tags$hr(),
            
            fluidRow(
              column(
                width = 4,
                div(
                  style = "text-align: center; padding: 20px; background: #f9f9f9; border-radius: 5px;",
                  icon("table", style = "font-size: 48px; color: #706993;"),
                  h4("All FDR Tables"),
                  p("Complete results for all three drugs"),
                  downloadButton("download_all_fdr", "Download ZIP", class = "btn-primary")
                )
              ),
              column(
                width = 4,
                div(
                  style = "text-align: center; padding: 20px; background: #f9f9f9; border-radius: 5px;",
                  icon("dna", style = "font-size: 48px; color: #70A0AF;"),
                  h4("GO Enrichment"),
                  p("All GO enrichment results"),
                  downloadButton("download_all_go", "Download ZIP", class = "btn-primary")
                )
              ),
              column(
                width = 4,
                div(
                  style = "text-align: center; padding: 20px; background: #f9f9f9; border-radius: 5px;",
                  icon("database", style = "font-size: 48px; color: #A0C1B9;"),
                  h4("Raw Sequencing Data"),
                  p(style = "color: #E67E22; font-weight: 500;", "Coming Soon"),
                  p(style = "font-size: 12px; color: #7F8C8D;", "RNA-seq FASTQ files will be available via GEO/ENA"),
                  tags$button(
                    class = "btn btn-primary",
                    disabled = "disabled",
                    style = "cursor: not-allowed; opacity: 0.6;",
                    "Access Data (Pending)"
                  )
                )
              )
            )
          )
        )
      ),
      
      # ========================================================================
      # ABOUT TAB
      # ========================================================================
      tabItem(
        tabName = "about",
        fluidRow(
          box(
            width = 12,
            title = "About This Application",
            status = "primary",
            solidHeader = TRUE,
            
            div(
              style = "font-size: 15px; line-height: 1.8;",
              
              h3("Oocyte Splicing Effects Explorer"),
              p("An interactive web application for exploring the transcriptomic effects of 
                splicing-modulating drugs in mouse oocytes."),
              
              tags$hr(),
              
              h4("Citation"),
              p(em("If you use this data or application in your research, please cite:")),
              div(
                style = "background: #f5f5f5; padding: 15px; border-left: 4px solid #706993; margin: 15px 0;",
                p("[Your Paper Citation Here]")
              ),
              
              h4("Data Source"),
              p("This application is based on alternative splicing analysis using the betAS package 
                and VAST-TOOLS for splicing event quantification."),
              
              h4("Methods"),
              tags$ul(
                tags$li(strong("Differential Splicing:"), " Events with FDR ≤ 0.05 and |ΔPSI| ≥ 0.1"),
                tags$li(strong("PSI:"), " Percent Spliced In (0-100 scale)"),
                tags$li(strong("ΔPSI:"), " Change in PSI between control and treated samples"),
                tags$li(strong("FDR:"), " False Discovery Rate (Benjamini-Hochberg correction)")
              ),
              
              h4("Technologies Used"),
              tags$ul(
                tags$li("R Shiny - Interactive web framework"),
                tags$li("betAS - Differential splicing analysis"),
                tags$li("VAST-TOOLS - Splicing event quantification"),
                tags$li("ComplexHeatmap - Advanced heatmap visualization"),
                tags$li("plotly - Interactive plotting")
              ),
              
              h4("Contact"),
              p(
                icon("envelope"), " For questions or issues, please contact: ",
                tags$a(href = "mailto:adel.aljord@crg.eu", "adel.aljord@crg.eu")
              ),
              
              h4("Acknowledgements"),
              p("This work was performed at the Centre for Genomic Regulation (CRG) and Collège de France, and 
                supported by institutional funding."),
              
              tags$hr(),
              
              div(
                style = "text-align: center; padding: 20px;",
                img(src = "logo_crg.png", height = 60, style = "margin: 10px 20px;"),
                img(src = "cfrance_logo.png", height = 60, style = "margin: 10px 20px;")
              ),
              
              hr(),
              p(
                style = "text-align: center; color: #7F8C8D; margin-top: 30px; font-size: 13px;",
                strong("Oocyte Splicing Effects Explorer"), br(),
                "Centre for Genomic Regulation (CRG) & Collège de France", br(),
                "© 2025 | Developed by Andrés Gordo Ortiz", br(),
                tags$a(href = "https://github.com/andresgordoortiz/oocyte_splicing_ssa_pladb_tub_2025", target = "_blank", 
                       icon("github"), " GitHub"),
                " | ",
                tags$a(href = "mailto:adel.aljord@crg.eu", 
                       icon("envelope"), " Contact")
              )
            )
          )
        )
      )
    )
  )
)
