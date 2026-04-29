# DEnote: visually guided scRNA-seq clustering & labeling
# Requires: shiny, bslib, Seurat, plotly, DT, ggplot2, ggrepel, dplyr,
#           shinycssloaders, shinyjs, rmarkdown, scales,
#           RColorBrewer, tidyr, tibble

library(shiny)
library(bslib)
library(Seurat)
library(plotly)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(DT)
library(scales)
library(shinycssloaders)
library(shinyjs)
library(rmarkdown)
library(RColorBrewer)
library(tibble)
library(tidyr)
library(thematic)
thematic_shiny(font = "auto")

# Allow uploads up to 10 GB (adjust as needed)
options(shiny.maxRequestSize = 10 * 1024^3)

# Source helpers & modules
source("helpers/processing.R")
source("helpers/plotting.R")
source("modules/mod_import.R")
source("modules/mod_clustering.R")
source("modules/mod_markers.R")
source("modules/mod_labeling.R")
source("modules/mod_report.R")

# UI
ui <- page_navbar(
  title = tags$span(
    style = "font-weight: 700; letter-spacing: 0.5px;",
    bsicons::bs_icon("diagram-3-fill"), " DEnote: Visually Guided Cell Labeller"
  ),
  theme = bs_theme(
    version    = 5,
    preset     = "darkly",
    bg         = "#222222",
    fg         = "#EEEEEE",
    primary    = "#4FC3F7",
    success    = "#66BB6A",
    info       = "#4DD0E1",
    font_scale = 0.92
  ) |>
    bs_add_rules(
      ".card-header { font-weight: 600; font-size: 0.8rem;
                      text-transform: uppercase; letter-spacing: 0.04em; }
       .navbar-brand { letter-spacing: 0.03em; }
       .accordion-button { font-size: 0.82rem; font-weight: 600; }
       .value-box .value-box-title { font-size: 0.72rem; text-transform: uppercase;
                                     letter-spacing: 0.05em; }"
    ),
  fillable  = TRUE,
  id        = "main_nav",

  # Header right: data loaded indicator
  nav_spacer(),
  nav_item(
    uiOutput("data_indicator")
  ),

  # ── Tab 1: Import ────────────────────────────────────────────────────────
  nav_panel(
    title = tagList(bsicons::bs_icon("upload"), " Import"),
    value = "import",
    mod_import_ui("import")
  ),

  # ── Tab 2: Clustering ────────────────────────────────────────────────────
  nav_panel(
    title = tagList(bsicons::bs_icon("grid-3x3-gap"), " Clustering"),
    value = "clustering",
    mod_clustering_ui("clustering")
  ),

  # Tab 3: Marker Exploration
  nav_panel(
    title = tagList(bsicons::bs_icon("search"), " Markers"),
    value = "markers",
    mod_markers_ui("markers")
  ),

  # Tab 4: Cell Labeling
  nav_panel(
    title = tagList(bsicons::bs_icon("tags"), " Label Cells"),
    value = "labeling",
    mod_labeling_ui("labeling")
  ),

  # Tab 5: Report & Export 
  nav_panel(
    title = tagList(bsicons::bs_icon("file-earmark-text"), " Report"),
    value = "report",
    mod_report_ui("report")
  ),

  # Footer
  nav_item(
    tags$small(class = "text-muted px-2",
               paste0("DEnote  ·  ", format(Sys.Date(), "%Y")))
  )
)

 # Server
server <- function(input, output, session) {

  # reactive state
  rv <- reactiveValues(
    seurat         = NULL,
    markers        = NULL,
    labels         = NULL,
    selected_cells = NULL
  )

  # ── Module servers ─────────────────────────────────────────────────────
  mod_import_server("import",     rv)
  mod_clustering_server("clustering", rv)
  mod_markers_server("markers",   rv)
  mod_labeling_server("labeling", rv)
  mod_report_server("report",     rv)

  # Data loaded indicator (top-right) 
  output$data_indicator <- renderUI({
    if (is.null(rv$seurat)) {
      tags$span(class = "badge bg-secondary", "No data loaded")
    } else {
      cluster_col <- active_cluster_col(rv$seurat)
      n_cl <- length(unique(rv$seurat@meta.data[[cluster_col]]))
      tags$span(
        class = "badge bg-success",
        paste0(scales::comma(ncol(rv$seurat)), " cells · ", n_cl, " clusters")
      )
    }
  })

  # Navigate to Clustering after load
  observeEvent(rv$seurat, {
    req(rv$seurat)
    nav_select("main_nav", "clustering")
  }, once = TRUE, ignoreNULL = TRUE)

}

shinyApp(ui, server)
