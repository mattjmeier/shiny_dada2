##########################################################################################
## 23-10-2019
## This version was developed by Rose Yan and Matthew J. Meier based on code provided
## by Ben Callahan and Paul J. McMurdie
## Shiny app for DADA2 
##########################################################################################
library("shiny")
library("shinyFiles") # <- used to access file system within the shiny app
##########################################################################################

shinyServer(function(input, output, session){
  ######################################################################################
  # render filter & trim 
  # render objects for side panel
  # shiny button - using path_home from global.R
  shinyDirChoose(input, "fastq_files", session = session, roots = volumes)
  
  # reactive to get file path directory
  file_path <- reactive({
    if(is.null(input$fastq_files)){
      return(" ")
    }
    parseDirPath(roots = volumes, input$fastq_files)
  }) #print the selected filepath directory
  output$dir_text <- renderText({
    file_path()
  })
  #render the info file in the textbox
  info <- reactive({
    if(!file.exists(paste0(file_path(), "/", 
                           input$info_txt)) || file_ext(input$info_txt) == ""){
      return(paste(input$info_txt, "does not exist."))
    }
    input$info_txt
  })
  output$info_file_text <- renderText({
    info()
  })
  
  # reactive for table
  info_file_path <- reactive({
    if(!file.exists(paste0(file_path(), "/", 
                           input$info_txt)) || file_ext(input$info_txt) == ""){
      return(NULL)
    }
    #IMPORTANT: check info file is valid format!
    else if(file_ext(input$info_txt) != "txt" & file_ext(input$info_txt) != "csv" &
            file_ext(input$info_txt) != "tsv"){
      return(NULL)
    }
    else{
      return(paste0(file_path(), "/", 
                    input$info_txt))
    }
  })
  info_file <- reactive({
    # need to check that info.txt file is present using validate(need)
    validate(need(input$fastq_files, info_file_path(), 
                  message = "Please enter fastq file directory and valid info file"))
    if(is.null(info_file_path())){
      return(NULL)
    }
    #list file is a list, rendertable only accepts matrices or data frames
    fread(input = info_file_path())
  })
  # table output
  output$files <- DT::renderDataTable({
    #reactive from action button
    input$update_table
    isolate(validate(need(info_file_path(), message = "Invalid info file")))
    isolate(info_file())
  })
  
  #render for main panel
  ######################################################################################
  # code from dada2 docker shiny app  
  plot_quality_graph <- quality_plot(input, output, info_file, file_path, info_file_path)
  
  output$quality_plot_output <- renderPlot({input$plot_quality_button
    isolate(validate(need(file_path() != " ", info_file_path(),
                          message = "Fastq file directory required")))
    isolate(validate(need(info_file()$Sample, message = "Invalid info file")))
    isolate(plot_quality_graph())
  })
  # execute filter & trim
  output$filter_trim_message <- renderText({
    # This stops action unless filtertrim action button has been pressed.
    input$actionb_filtertrim
    isolate(validate(need(info_file(), message = "Invalid info file")))
    isolate(validate(need(file_path() != " ", info_file_path(),
        message = "Please enter fastq file directory and valid info file")))
    isolate(filter_trim(input, output, info_file, file_path, info_file_path))
  })
  ######################################################################################
  # render learn errors
  # Rendering for side panel
  output$le_data <- renderText({
    if(input$data_type == "ion_torrent"){
      return("Data origin selected: Ion torrent")
    }
    else{
      return("Data origin selected: Illumina")
    }
  })
  # reactive to get file path directory
  filt_file_path <- reactive({
    input$le_update_filt_table
    input$rd_update_filt_table
    input$at_update_filt_table
    validate(need(file_path(), message = "Filtered fastq files required"))
    validate(need(file.exists(paste0(file_path(), "/FT")), 
                  message = "FT folder not found"))
    return(paste0(file_path(), "/FT"))
  })
  #print the selected filepath directory
  output$le_filt_dir <- renderText({
    filt_file_path()
  })
  
  #get info file
  filt_info_file_path <- reactive({
    if(!file.exists(paste0(filt_file_path(), "/info.txt"))){
      return(NULL)
    }
    return(paste0(filt_file_path(), "/info.txt"))
  }) # read in info file
  filt_info_file <- reactive({
    if(is.null(filt_info_file_path())){
      return(NULL) 
    }
    fread(input = filt_info_file_path())
  })
  #output the info file into the text box
  output$le_filt_info_txt <- renderText({
    if(is.null(filt_info_file_path())){
      return("info.txt does not exist") # <- check correct folder is selected
    }
    "info.txt"
  })
  
  # info table
  output$le_filt_files <- DT::renderDataTable({
    req(input$le_update_filt_table)
    isolate(validate(need(filt_info_file_path(),
                          message = "Filtered fastq files required")))
    isolate(filt_info_file())
  })
  ######################################################################################
  # Rendering for main panel
  output$le_message <- renderText({
    input$learn_errors_button
    isolate(validate(need(filt_file_path(), message  = "Filtered fastq files required")))
    isolate(validate(need(filt_info_file_path(), message = "info.txt file not found")))
    # check filtered files exist
    isolate(file <- filt_info_file())
    isolate(path <- filt_file_path())
    isolate(validate(need(file.exists(paste0(path, "/", file$File)),
                          message = "Filtered files not found")))
    isolate(learn_errors(input, output, session, filt_info_file, filt_file_path))
  })
  
  # error rates graphs
  le_plot_output <- reactive({
    req(filt_info_file, filt_file_path)
    errors_plot(input, output, session, filt_info_file, filt_file_path)
  })
  output$plot_error_graphs <- renderUI({
    input$learn_errors_button
    isolate(validate(need(filt_file_path(), 
                          message  = "Filtered fastq files required")))
    isolate(validate(need(filt_info_file_path(), message = "info.txt file not found")))
    # check filtered files exist
    isolate(file <- filt_info_file())
    isolate(path <- filt_file_path())
    isolate(validate(need(file.exists(paste0(path, "/", file$File)),
                          message = "Filtered files not found")))
    if(!is.null(learn_errors(input, output, session, filt_info_file, filt_file_path))){
      actionButton("plot_error_button", label = "Plot Error Rates", 
                   icon("th", lib = "glyphicon"))
    }
    else{
      NULL
    }
  })
  # plot only forward for ion torrent
  observeEvent(
    input$data_type == "ion_torrent", {
      output$f_errors_plot <- renderPlot({
        req(input$plot_error_button)
        isolate(le_plot_output()$ForwardErrors)
      })
      output$r_errors_plot <- NULL
    })
  observeEvent(
    input$data_type == "illumina", {
      output$f_errors_plot <- renderPlot({
        req(input$plot_error_button)
        isolate(le_plot_output()$ForwardErrors)
      })
      output$r_errors_plot <- renderPlot({
        req(input$plot_error_button)
        isolate(le_plot_output()$ReverseErrors)
      })
    })
  ######################################################################################
  # render for run dada2
  # Rendering for side panel
  output$rd_data <- renderText({
    if(input$data_type == "ion_torrent"){
      return("Data origin selected: Ion torrent")
    }
    else{
      return("Data origin selected: Illumina")
    }
  })
  # print the selected filepath directory
  output$rd_filt_dir <- renderText({
    filt_file_path()
  })
  #output the info file into the text box
  output$rd_filt_info_txt <- renderText({
    if(is.null(filt_info_file_path())){
      return("info.txt does not exist") # <- check correct folder is selected
    }
    "info.txt"
  })
  
  # info table
  output$dada_filt_files <- DT::renderDataTable({
    req(input$rd_update_filt_table)
    isolate(validate(need(filt_info_file_path(),
                          message = "Filtered fastq files required")))
    isolate(filt_info_file())
  })
  ######################################################################################
  # Rendering for main panel
  output$dada_message <- renderText({
    input$run_dada_button
    isolate(validate(need(filt_info_file(), filt_info_file_path(),
                          message = "Filtered fastq files required")))
    isolate(if(isolate(input$data_type == "illumina")){
      isolate(run_dada(input,
                       output, filt_file_path, filt_info_file_path))
    }
    else if(isolate(input$data_type == "ion_torrent")){
      isolate(rd_ion_torrent(input, output, filt_file_path, filt_info_file))
    }
    else{
      NULL
    })
  })
  
  output$show_table <- renderUI({
    if(input$data_type == "illumina" && !is.null(run_dada(input,
            output, filt_file_path, filt_info_file_path))){
      fluidRow(column(width = 4, checkboxGroupInput("select_columns",
        label = "Columns to Show", choices = c("Sample" = "Sample", 
            "Sequence" = "sequence", "Abundance" = "abundance"),
              selected = c("Sample", "abundance"))),
        column(width = 6, actionButton("update_sequence_table",
                                              label = "Update Sequence Table",
                                              icon("th_list", lib = "glyphicon")))
      )
    }
    else if(input$data_type == "ion_torrent" && 
            !is.null(rd_ion_torrent(input, output, filt_file_path, filt_info_file))){
      fluidRow(column(width = 4, checkboxGroupInput("select_columns",
          label = "Columns to Show", choices = c("Sample", "Sequence", "Abundance"),
                                                    selected = c("Sample", "Abundance"))),
               column(width = 6, actionButton("update_sequence_table",
                                              label = "Update Sequence Table",
                                              icon("th_list", lib = "glyphicon")))
      )
    }
    else{
      NULL
    }  
  })
  
  output$sequence_table <- DT::renderDataTable({
    input$update_sequence_table
    req(filt_file_path(), filt_info_file(), filt_info_file_path(), 
        file.exists(paste0(filt_file_path(), "/dadaTabBimeraFilt.RDS")))
    # user needs to select at least one column
    isolate(validate(need(input$select_columns, 
                          message = "Please select a column to show")))
    isolate(
      table <-sequence_table_function(input, output, session, filt_file_path))
    isolate(select(table, input$select_columns))
  })
  ######################################################################################
  # render assign taxonomy
  # Rendering for side panel
  output$at_data <- renderText({
    if(input$data_type == "ion_torrent"){
      return("Data origin selected: Ion torrent")
    }
    else{
      return("Data origin selected: Illumina")
    }
  })
  # print the selected filepath directory
  output$at_filt_dir <- renderText({
    filt_file_path()
  })
  #output the info file into the text box
  output$at_filt_info_txt <- renderText({
    if(is.null(filt_info_file_path())){
      return("info.txt does not exist") # <- check correct folder is selected
    }
    "info.txt"
  })
  
  # info table
  output$at_filt_files <- DT::renderDataTable({
    req(input$at_update_filt_table)
    isolate(validate(need(filt_info_file_path(),
                          message = "Filtered fastq files required")))
    isolate(filt_info_file())
  })
  # render the ui for the file input if file upload is selected
  output$file_input1 <- renderUI({
    if(input$ref_seq == "upload"){
      isolate(fileInput("ref_seq_file", label = "Upload Reference Sequence file"))
    }
    else{
      NULL
    }
  })
  output$file_input2 <- renderUI({
    if(input$ref_species == "upload"){
      isolate(fileInput("ref_species_file", label = "Upload species assignment file"))
    }
    else{
      NULL
    }
  })
  ######################################################################################
  # Rendering for main panel
  # reactive for assigning taxonomy
  taxa <- reactive({NULL}) # initialize
  taxa <- reactive({
    input$assign_tax
    isolate(validate(need(filt_file_path() != " ", 
                          message = "Filtered fastq file folder required")))
    isolate(validate(need(file.exists(paste0(filt_file_path(), "/seqmat.RDS")),
                          message = "seqmat.RDS file not found")))
    # Check that the selected ref file exists or a chosen file is uploaded
    if(input$ref_seq == "upload"){
      isolate(validate(need(input$ref_seq_file, message = "Reference does not exist")))
    }
    isolate(assign_taxonomy(input, output, filt_file_path, file_path))
  })
  output$tax_message <- renderText({
    input$assign_tax
    isolate(validate(need(taxa(), message = "...")))
    isolate("Taxonomy assigned successfully!")
  })
  
  # reactive for assigning species  
  taxa_species <- reactive({
    input$assign_species_button
    isolate(validate(need(filt_file_path(),
                          message = "Filtered fastq file folder required")))
    isolate(validate(need(file.exists(paste0(filt_file_path(), "/seqmat.RDS")),
                          message = "seqmat.RDS file not found")))
    # Check that the selected ref file exists or that a chosen file is uploaded
    if(input$ref_species == "upload"){
      isolate(validate(need(input$ref_species_file,
                            message = "Reference does not exist")))
    }
    if(file.exists(paste0(filt_file_path(), "/taxa.RDS"))){
      taxa <- readRDS(paste0(filt_file_path(), "/taxa.RDS"))
      isolate(assign_species(input, output, filt_file_path, taxa, file_path))
    }
    else{
      isolate(validate(need(taxa(), message = "Taxonomy not yet assigned")))
      isolate(assign_species(input, output, filt_file_path, taxa(), file_path))
    }
  })
  output$species_message <- renderText({
    input$assign_species_button
    isolate(req(taxa_species()))
    isolate(taxa_species())
  })
  # testing to see that the assignment worked
  output$preview <- renderPrint({
    input$assign_species_button
    isolate(req(file.exists(paste0(filt_file_path(), "/taxa_species.RDS"))))
    isolate(taxa_species <- readRDS(paste0(filt_file_path(), "/taxa_species.RDS")))
    isolate(taxa.print <- taxa_species)
    isolate(rownames(taxa.print) <- NULL)
    isolate(head(taxa.print, 10))
  })
  # run the phyloseq function
  ps_object <- reactive({
    input$create_ps
    isolate(validate(need(filt_file_path(),
                          message = "Filtered fastq file folder required")))
    isolate(validate(need(file.exists(paste0(filt_file_path(), "/taxa_species.RDS")), 
                          message = "Taxonomy not yet assigned")))
    isolate(validate(need(file.exists(paste0(filt_file_path(), "/seqmat.RDS")),
                          message = "seqmat.RDS file not found")))
    isolate(validate(need(input$metadata_file, message = "metadata file not found")))
    # check if a file exists already in the filtered files dir
    if(file.exists(paste0(filt_file_path(), "/taxa_species.RDS"))){
      taxa_species <- readRDS(paste0(filt_file_path(), "/taxa_species.RDS"))
      isolate(create_ps_object(input, output, filt_file_path, taxa_species))
    }
    else{
      isolate(validate(need(taxa(), taxa_species(),
                            message = "Taxonomy not yet assigned")))
      isolate(create_ps_object(input, output, filt_file_path, taxa_species()))
    }
  })
  # output message for ps object
  output$final_message <- renderText({
    input$create_ps
    isolate(validate(need(ps_object(), message = "...")))
    # isolate(ps_object())
    isolate("DADA2 pipeline completed!")
  })
  # render download button
  output$download_ui <- renderUI({
    if(!is.null(ps_object())){
      downloadButton("download_button", label = "Download Phyloseq Data")
    }
    else{
      NULL
    }
  })
  # download the ps object
  output$download_button <- downloadHandler(
    filename = paste0("phyloseq_data_", Sys.Date(), ".rda"),
    content = function(save_file){ps_data <- ps_object()
    save(ps_data, file = save_file)}
  )
  ######################################################################################
})