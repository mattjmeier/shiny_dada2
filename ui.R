##########################################################################################
## 23-10-2019
## This version was developed by Rose Yan and Matthew J. Meier based on code provided
## by Ben Callahan and Paul J. McMurdie
## Shiny app for DADA2 
##########################################################################################
library("shiny")
library("shinyFiles")
##########################################################################################

# ui filter & trim

ft_side_panel <-  sidebarPanel(
  radioButtons("data_type", label = "Select origin of sequenced files", 
               choices = c("Ion Torrent" = "ion_torrent", "Illumina" = "illumina")),
  tags$hr(),
  fluidRow(column(width=6, h4("Data Selection")
  ),
  # directory button to get the fastq file directory
  column(width = 4, shinyDirButton(id = "fastq_files", 
                                   label = "Choose file directory", 
                                   title = "Select fastq file folder")
  )
  ),
  tags$hr(),
  h5("Input Directory Path"),
  verbatimTextOutput("dir_text", placeholder = TRUE),
  # user input the info file name
  textInput("info_txt", label = "Enter data info file name", value = "info.txt"),
  withSpinner(verbatimTextOutput("info_file_text", placeholder = TRUE),
              type = 7, size = 0.5, proxy.height = "30px"),
  # Updates table once user has entered a file dir
  actionButton(inputId = "update_table", label = "Update Table"),
  tags$hr(),
  # Parameters for plotting the graph and trimming parameters
  ##################################################################
  #NOTE: Following code is taken from ui.R of dada2 docker shiny app
  fluidRow(
    h4("Quality Downsample"),
    column(width = 5,
           sliderInput("nreads", "Reads / Sample (log10)", 
                       value = 3, min = 2, max = 5, step = 1)
    ),
    column(width = 5,
           sliderInput("NSamples", "Number of Samples", 
                       value = 10L, min = 1L, max = 50L, step = 1L, round = TRUE)
    )
  ),
  h5("Max Length"),
  verbatimTextOutput('maxLength', placeholder = TRUE),
  fluidRow(h4("Trimming Prediction"),
           # input$minQual # 15
           column(width = 6,
                  sliderInput("minQual", "Min Smoothed Quality",
                              value = 15, min = 0, max = 40, step = 1)
           ),
           column(width = 6, sliderInput("qquantile", "Quality Quantile",
                                         value = 0.2, min = 0.0, max = 1.0, step = 0.05)
           )
  ),
  fluidRow(h4("Trimming Parameters"),
           column(width = 6, uiOutput("uiForward")
           ),
           column(width = 6, uiOutput("uiReverse")
           )
  ),
  h6("Note: An error will occur if the truncation length is set longer than some reads."),
  fluidRow(
    h4("Parameters for Filter and Trim"),
    column(width = 3, numericInput("maxEE", "maxEE", value = 2L, min = 0L, step = 1L)
    ),
    column(width = 3, sliderInput("multithread", "multithread", 
                                  value = NCores - 1L, min = 1L, max = NCores,
                                  step = 1L, round = TRUE)
    ),
    column(width = 3, sliderInput("n", "n [log10]",
                                  value = 5, min = 2, max = 8, step = 1)
    ),
    column(width = 3, numericInput("truncQ", "truncQ",
                                   value = 2L, min = 0L, step = 1L)
    )
  ),
  ############################################################################
  tags$hr(),
  # table to view input files
  h4("Fastq Info Table"),
  DT::dataTableOutput("files")
)

ft_main_panel <- mainPanel( #added some icons
  actionButton("plot_quality_button", 
               label = "Plot Quality Profile",
               icon("chart-line")),
  
  h5("Samples included in quality profiling:"),
  withSpinner(textOutput("include_samples"), type = 7, size = 0.6, proxy.height = "50px"),
  h5("Quality Plot by Cycle"),
  withSpinner(plotOutput("quality_plot_output"), type = 6, size = 0.8),
  actionButton("actionb_filtertrim", label = "Execute Filter & Trim",
               icon("filter", lib = "glyphicon")),
  withSpinner(textOutput("filter_trim_message"), type = 7,
              size = 0.6, proxy.height = "50px")
)
##########################################################################################
# ui learn errors

le_side_panel <- sidebarPanel(
  textOutput("le_data"),
  tags$hr(), 
  fluidRow(column(width = 8, h5("Input Directory Path"),
                  verbatimTextOutput("le_filt_dir", placeholder = TRUE)),
           column(width = 4, h5("Info file"),
                  verbatimTextOutput("le_filt_info_txt", placeholder = TRUE))
  ),
  actionButton("le_update_filt_table", label = "Update Table"),
  tags$hr(),
  #######################################################################################
  # Parameters are taken from the dada2 docker shiny app
  fluidRow(
    h4("Parameters for Learn Errors"),
    column(width = 3,
           sliderInput("LE_learnSize", "Number of Samples", 
                       value = 6, min = 1, max = 30, step = 1)
    ),
    column(width = 3,
           sliderInput("LE_minSize", "Min. File Size [log10]", 
                       value = 4, min = 3, max = 6, step = 1)
    ),
    column(width = 3,
           sliderInput("LE_nreads", "nreads [log10]", value = 5, 
                       min = 2, max = 8, step = 1)
    ),
    column(width = 3,
           sliderInput("LE_multithread", "multithread", 
            value = NCores - 1L, min = 1L, max = NCores, step = 1L, round = TRUE)
    )
  ),
  tags$hr(),
  # table to view input files
  h4("Filtered Fastq Info Table"),
  DT::dataTableOutput("le_filt_files")
)

le_main_panel <- mainPanel(
  actionButton("learn_errors_button", label = "Learn Errors", icon("book-open")),
  withSpinner(textOutput("le_message"), type = 7, size = 0.6, proxy.height = "50px"),
  tags$hr(),
  uiOutput("plot_error_graphs"),
  h4("Forward Error Rates"),
  fluidRow(column(width = 7, 
                  withSpinner(plotOutput("f_errors_plot"), type = 6, size = 0.8))),
  h4("Reverse Error Rates"),
  fluidRow(column(width = 7,
                  withSpinner(plotOutput("r_errors_plot"), type = 6, size = 0.8)))
)
##########################################################################################
# ui run dada()

dada_side_panel <- sidebarPanel(
  textOutput("rd_data"),
  tags$hr(), 
  fluidRow(column(width = 8, h5("Input Directory Path"),
                  verbatimTextOutput("rd_filt_dir", placeholder = TRUE)),
           column(width = 4, h5("Info file"),
                  verbatimTextOutput("rd_filt_info_txt", placeholder = TRUE))
  ),
  actionButton("rd_update_filt_table", label = "Update Table"),
  tags$hr(),
  #######################################################################################
  # Parameters are taken from the dada2 docker shiny app
  fluidRow(
    h4("Parameters for DADA: Sample Interference"),
    column(width = 3,
           sliderInput("RD_minSize", "Min. File Size [log10]",
                       value = 4, min = 3, max = 6, step = 1)
    ),
    column(width = 3,
           sliderInput("RD_nreads", "nreads [log10]", 
                       value = 5, min = 2, max = 8, step = 1)
    ),
    column(width = 3,
           sliderInput("RD_multithread", "multithread", 
            value = NCores - 1L, min = 1L, max = NCores, step = 1L, round = TRUE)
    )
  ),
  # table to view input files
  h4("Filtered Fastq Info Table"),
  DT::dataTableOutput("dada_filt_files")
)

dada_main_panel <- mainPanel(
  actionButton("run_dada_button", label = "Run DADA", icon("play", lib = "glyphicon")),
  withSpinner(textOutput("dada_message"), type = 7, size = 0.6, proxy.height = "50px"),
  tags$hr(),
  uiOutput("show_table"),
  fluidRow(column(width = 10, DT::dataTableOutput("sequence_table"))
  )
)
##########################################################################################
# ui assign taxonomy 

at_side_panel <- sidebarPanel(
  textOutput("at_data"),
  tags$hr(), 
  fluidRow(column(width = 8, h5("Input Directory Path"),
                  verbatimTextOutput("at_filt_dir", placeholder = TRUE)),
           column(width = 4, h5("Info file"),
                  verbatimTextOutput("at_filt_info_txt", placeholder = TRUE))
  ),
  actionButton("at_update_filt_table", label = "Update Table"),
  tags$hr(),
  selectInput("ref_seq", label = "Select Reference Sequences",
              choices = c("None selected" = "none",
                          "Upload File" = "upload")),
  uiOutput("file_input1"),
  tags$hr(),
  selectInput("ref_species", label = "Select species assignment",
              choices = c("None selected" = "none",
                          "Upload File" = "upload")),
  uiOutput("file_input2"),
  tags$hr(),
  # table to view input files
  h4("Filtered Fastq Info Table"),
  DT::dataTableOutput("at_filt_files")
)

at_main_panel <- mainPanel(
  actionButton("assign_tax", label = "Assign Taxonomy", icon("list-alt",
                                                             lib = "glyphicon")),
  withSpinner(textOutput("tax_message"), type = 7, size = 0.6, proxy.height = "50px"),
  tags$hr(),
  actionButton("assign_species_button", label = "Assign Species",
               icon("tasks", lib = "glyphicon")),
  withSpinner(textOutput("species_message"), type = 7, size = 0.6, proxy.height = "50px"),
  # preview of the assignment worked
  tags$hr(),
  verbatimTextOutput("preview"),
  tags$hr(),
  # Create phyloseq object
  fileInput("metadata_file", label = "Upload metadata file", 
            accept = c(".txt", ".metadata", ".csv", ".tsv")),
  actionButton("create_ps", label = "Create phyloseq data", icon("pen")),
  withSpinner(textOutput("final_message"), type = 7, size = 0.6, proxy.height = "50px"),
  tags$hr(),
  uiOutput("download_ui"),
  tags$hr()
)
##########################################################################################
ui <- navbarPage(theme = shinytheme("yeti"),
                 windowTitle = "Shiny DADA2",
                 title = "Shiny DADA2" %>% a(href="http://benjjneb.github.io/dada2/",
                                             style="color:#ffffff"),
                 tabPanel("Filter and Trim",
                  h3(icon("filter", lib = "glyphicon"), icon("cut"),
                      "Filter & Trim fastq sequences"),
                  {sidebarLayout(ft_side_panel, ft_main_panel)}
                 ),
                 tabPanel("Learn Errors", 
                  h3(icon("book-open"),"Learn Errors"),
                  h5("Learn error matrix from a subset of samples in each direction"),
                  {sidebarLayout(le_side_panel, le_main_panel)}
                 ),
                 tabPanel("Run DADA",
                  h3(icon("pencil", lib = "glyphicon"), "Run DADA"),
                  h5("Sample Interference, Merge Pairs, Create Sequence Table"),
                  {sidebarLayout(dada_side_panel, dada_main_panel)}
                 ),
                 tabPanel("Assign Taxonomy",
                  h3(icon("list-alt", lib = "glyphicon"), "Assign Taxonomy"),
                  h5("Assign taxonomy, add species, handoff to phyloseq"),
                  {sidebarLayout(at_side_panel, at_main_panel)}
  )
)
##########################################################################################

shinyUI(ui)