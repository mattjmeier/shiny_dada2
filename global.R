##########################################################################################
## 23-10-2019
## This version was developed by Rose Yan and Matthew J. Meier based on code provided
## by Ben Callahan and Paul J. McMurdie
## Shiny app for DADA2 
##########################################################################################
library("shiny")
library("shinyFiles") # <- used to access file system within the shiny app
library("dada2") # <- NOTE: using dada2 version 1.8. Not using 1.12!
library("DT")
library("data.table")
library("ggplot2")
library("magrittr")
library("tools") 
library("fs") # <- added this library to get the home path from any system
library("dplyr")
library("phyloseq") # <- require this library to create the phyloseq object
library("shinycssloaders")
library("shinythemes")
library("reshape2")
##########################################################################################
theme_set(theme_bw())
home <- fs::path_home() # gets home path from system
root <- "/"
downloads <- "~/Downloads/"
volumes <- c(Home = home, Downloads = downloads, Root = root)
options(shiny.maxRequestSize=100*1024^2) # increase max size to 100MB file upload
NCores <- max(1L, RcppParallel::defaultNumThreads())
options(DT.options = list(scrollX = TRUE, fixedColumns = TRUE))
##########################################################################################
# filter & trim functions
# function to plot the quality graph (taken from dada2 docker shiny app)
plot_quality_by_cycle = function(CycleStats, CycleCounts, TrimTable = NULL){
  p = ggplot(data = CycleStats,
             mapping = aes(Cycle, y = Quality, color = Statistic)) +
    ylim(0, 40) +
    geom_raster(data = CycleCounts,
                mapping = aes(x = Cycle,
                              y = Score,
                              fill = log10(Proportion)),
                inherit.aes = FALSE) +
    # grey0 (black) to grey100 (white)
    # scale_fill_gradient(
    scale_fill_gradient2(
      midpoint = -1,
      low = "grey99",
      mid = "grey90",
      high = "grey10",
      space = "white",
      na.value = "white",
      guide = "colourbar") +
    geom_path(size = 0.25) +
    geom_path(mapping = aes(y = Smooth), size = 1, alpha = 0.35) +
    geom_text(mapping = aes(label = Cycle,
                            y = 1,
                            hjust = ifelse(Side == "Right", yes = 1.1, no = -0.1)),
              data = TrimTable,
              vjust = 0.5,
              color = "black", size = 3) +
    facet_wrap(~Direction, nrow = 2)
  if( !is.null(TrimTable) ){
    p <- p + geom_vline(mapping = aes(xintercept = Cycle), data = TrimTable, size = 0.25)
  }
  return(p)
}

#tabulate quality function required
tabulate_quality = function(fastqFile, nReads = 1e4){
  require("ShortRead")
  require("data.table")
  message("Reading from file: ", fastqFile, "\n")
  # fastqFile = fls[1]
  FQS = FastqSampler(con = fastqFile, n = nReads)
  fq = yield(FQS)
  # Fun borrowed from ShortRead internals
  .qa_perCycleQuality = function(abc, quality){
    if (missing(abc) || dim(abc)[[3]] == 0) {
      df <- data.frame(Cycle=integer(0), Quality=numeric(0),
                       Score=numeric(0), Count=integer(0),
                       lane=character(0))
      return(df)
    }
    abc <- apply(abc, 2:3, sum)
    q <- factor(rownames(abc)[row(abc)], levels=rownames(abc))
    q0 <- as(do.call(class(quality), list(rownames(abc))), "matrix")
    df <- data.frame(Cycle=as.integer(colnames(abc)[col(abc)]),
                     Quality=q, Score=as.integer(q0)[q],
                     Count=as.vector(abc),
                     row.names=NULL)
    df[df$Count != 0, ]
  }
  abc <- alphabetByCycle(fq)
  perCycleQuality <- data.table(.qa_perCycleQuality(abc, quality(fq)))
  perCycleQuality[, fastqFile := fastqFile]
  return(perCycleQuality)
}

#' Write tab delimited table
#' 
#' @inheritParams utils::write.table
#' @param ... additional args passed to [write.table]
#' 
write_table_tab = function(x, file, ...){
  write.table(
    x = x, 
    file = file,
    sep = "\t", 
    row.names = FALSE, 
    col.names = TRUE,
    append = FALSE,
    quote = FALSE,
    ...) 
}

# placed code from dada2 docker shiny app in a function
quality_plot <- function(input, output, info_file, file_path, info_file_path){
  # NOTE: Using code from the dada2 docker shiny app. 
  ########################################
  # Quality Sub-Sampling 
  # Tabulate q-values into a long granular table for each file.
  ########################################
  #need file path to the info file
  fastqFilesTab <- reactive({
    info_file()
  })
  #need this function for the directory path
  inputDirPath <- reactive({
    file_path()
  })
  #need this function 
  # The samples
  includeSamples = reactive({
    validate(need(info_file()$Sample, message = "Invalid info file!"))
    validate(need(fastqFilesTab(), message = "..."))
    validate(need(input$NSamples,
                  message = "Need number of files for downsampled quality evaluation."))
    fastqFilesTab <- copy(fastqFilesTab())
    # Define a subset of sample names
    NSamples = input$NSamples
    # The minimum of either the prescribed number of samples and the number available
    NSamples <- c(NSamples, fastqFilesTab[, uniqueN(Sample)]) %>% min(na.rm = TRUE)
    includeSamples = fastqFilesTab[, unique(Sample)] %>% 
      sample(size = NSamples, replace = FALSE) %>% 
      sort
    return(includeSamples)
  })
  TabulateQuality = reactive({
    # How many reads, at most, to read from each file.
    # validate(need(input$nreads, message = "..."))
    input$nreads %>% need(message = "...") %>% validate
    nReads = 10^(input$nreads)
    # Don't move forward if no sequence table yet.
    validate(need(fastqFilesTab(), message = "Select Info File."))
    fastqFilesTab <- copy(fastqFilesTab())
    includeSamples = includeSamples()
    # Progress bar...
    # progress <- shiny::Progress$new(session, min=1, max=15)
    # on.exit(progress$close())
    qtab <- fastqFilesTab[
      (Sample %chin% includeSamples),
      {message("Tabulating from:\n", File)
        incProgress(
          amount = incProgUnit,
          message = Direction[1],
          detail = File[1])
        tabulate_quality(
          fastqFile = file.path(inputDirPath(), File),
          nReads = nReads)
      }, by = c("File", "Direction")]
    return(qtab)
  })
  maxLength = reactive({
    validate(need(info_file()$Sample, message = "Invalid info file!"))
    qtab = TabulateQuality() %>% copy
    return(
      qtab[, max(Cycle)]
    )
  })
  ########################################
  # Summarize Sub-Sampled Quality, reactively
  ########################################
  QualSummReact = reactive({
    validate(need(input$qquantile, message = "..."))
    # Quantile to show
    desiredQuantile = input$qquantile
    qtab = TabulateQuality() %>% copy
    # `CycleCounts` is collapsed for distributional summary at each cycle
    CycleCounts = qtab[, .(Count = sum(Count, na.rm = TRUE)),
                       by = c("Direction", "Cycle", "Score")]
    setorderv(CycleCounts, c("Cycle", "Score"), order = c(1, 1))
    # Define quantile
    CycleCounts[, Quantile := cumsum(Count)/sum(Count, na.rm = TRUE),
                by = c("Direction", "Cycle")]
    # Define proportion
    CycleCounts[, Proportion := (Count)/sum(Count, na.rm = TRUE),
                by = c("Direction", "Cycle")]
    # Collect summary statistic
    cycsum = CycleCounts[,
                         .(
                           mean = sum(Score * Count, na.rm = TRUE) / sum(Count, 
                                                                         na.rm = TRUE),
                           median = Score[(Quantile >= 0.50)][1],
                           QuantN = Score[(Quantile >= desiredQuantile)][1]
                         ),
                         by = c("Direction", "Cycle")]
    nameQuantN = paste0("Quantile_", desiredQuantile)
    setnames(cycsum, "QuantN", nameQuantN)
    suppressWarnings({
      CycleStats <- melt.data.table(cycsum,
                                    id.vars = c("Direction", "Cycle"),
                                    variable.name = "Statistic",
                                    value.name = "Quality")
    })
    fitTab = CycleStats[, .(fit = list(fit = loess(Quality ~ Cycle, data = .SD))),
                        by = c("Direction", "Statistic")]
    # Define smoothed values at each entry
    # Form: A[B, bb:=i.b, on='a']
    # CycleStats[fitTab, Smooth := predict(i.fit[[1]],
    #newdata = Cycle), on = c("Direction", "Statistic")]
    CycleStats[fitTab, Fit := i.fit, on = c("Direction", "Statistic")]
    CycleStats[, Smooth := predict(Fit[[1]], newdata = Cycle),
               by = c("Direction", "Statistic", "Cycle")]
    CycleStats[, Fit := NULL]
    return(
      list(
        CycleStats = CycleStats,
        CycleCounts = CycleCounts
      )
    )
  })
  # Compute default/predicted trim values
  suggestedTrimTable = reactive({
    # Min-quality at the indicated quantile (smooth)
    # for computing the default trimming parameters
    # e.g. 15L
    validate(need(info_file()$Sample, message = "Invalid info file!"))
    minQual = input$minQual
    LeftTrimDefault = 10L
    CycleStats = QualSummReact()$CycleStats %>% copy
    # Use this to define default Right Trim
    RightTrimTable = CycleStats[, .(Cycle = min(max(Cycle, na.rm = TRUE), 
                                                Cycle[(Smooth <= minQual)],
                                                na.rm = TRUE)), by = "Direction"]
    RightTrimTable[, Side := "Right"]
    # Define default Left Trim
    LeftTrimTable = copy(RightTrimTable)[, c("Side", "Cycle") := list("Left",
                                                                      LeftTrimDefault)]
    TrimTable = list(RightTrimTable, LeftTrimTable) %>% rbindlist
    return(TrimTable)
  })
  output$uiForward <- renderUI({
    validate(need(info_file_path(), message = "..."))
    if(input$data_type == "ion_torrent"){ # change the default trim parame
      ForwardRight = suggestedTrimTable()[(Direction == "F" & Side == "Right")]$Cycle
      return( # default left trimming is 15 rather than 10 for ion torrent
        sliderInput("Forward", "Forward", min = 1, max = maxLength(),
                    step = 5, value = c(15, ForwardRight))
      )
    }
    else{
      ForwardRight = suggestedTrimTable()[(Direction == "F" & Side == "Right")]$Cycle
      return(
        sliderInput("Forward", "Forward", min = 1, max = maxLength(),
                    step = 5, value = c(10, ForwardRight))
      )
    }
  })
  output$uiReverse <- renderUI({
    validate(need(info_file_path(), message = "..."))
    if(input$data_type == "ion_torrent"){ # no need to trim reverse reads
      return(NULL)
    }
    else if(input$data_type == "illumina"){
      ReverseRight = suggestedTrimTable()[(Direction == "R" & Side == "Right")]$Cycle
      return(
        sliderInput("Reverse", "Reverse", min = 1, max = maxLength(),
                    step = 5, value = c(10, ReverseRight))
      )
    }
    else{
      return(NULL)
    }
  })
  # Populate the trim table to be used in the chart and in the executed filtertrimming
  # from the user-input, the dual-slider for each direction.
  TrimTable = reactive({
    if(input$data_type == "ion_torrent"){
      data.table(
        Direction = c("F", "F", "R", "R"), 
        Side = c("Left", "Right", "Left", "Right"),
        Cycle = c(input$Forward[1],
                  input$Forward[2],
                  NA,
                  NA)
      )
    }
    else{
      data.table(
        Direction = c("F", "F", "R", "R"), 
        Side = c("Left", "Right", "Left", "Right"),
        Cycle = c(input$Forward[1],
                  input$Forward[2],
                  input$Reverse[1],
                  input$Reverse[2])
      )
    }
  })
  ########################################
  # FilterTrim Quality Graphic
  ########################################
  # Define the main graphic as a ggplot2 object.
  pQbyC = reactive({
    qualitySummaryList = QualSummReact()
    TrimTable = copy(TrimTable())
    # Modify TrimTable manually/hardcoded as-needed
    # to avoid the weird "cliff" in this particular dataset
    # TrimTable[(Direction == "R1" & Side == "Right"), Cycle := 252]
    pQbyC = plot_quality_by_cycle(CycleStats = qualitySummaryList$CycleStats,
                                  CycleCounts = qualitySummaryList$CycleCounts,
                                  TrimTable = TrimTable)
    return(pQbyC)
  })
  output$maxLength <- renderText({
    validate(need(info_file_path(), message = ""))
    maxLength()})
  
  output$include_samples <- renderText({
    input$plot_quality_button
    input$quality_plot_output
    isolate(validate(need(info_file_path(), message = "")))
    isolate(paste0(includeSamples(), collapse = ", "))
  })
  return(pQbyC)
}

# Execute filter and trim
filter_trim <-function(input, output, info_file, file_path, info_file_path){
  ########################################
  #need file path to the info file
  fastqFilesTab <- reactive({
    info_file()
  })
  #need this function for the directory path
  inputDirPath <- reactive({
    file_path()
  })
  ########################################
  # FilterTrim Execution
  ########################################
  # The filtertrim directory
  ftDir = reactive({
    # Create the filtertrim directory if it doesn't exist yet
    ftDir = "FT" %>% file.path(inputDirPath(), .)
    if(!(ftDir %>% dir.exists)){
      ftDir %>% dir.create()
    }
    return(ftDir)
  })
  
  ftFilesTable = reactive({
    ftDir = ftDir()
    # Define the file paths of the filter-trimmed files
    filterTrimFilesTable = fastqFilesTab() %>% copy
    filterTrimFilesTable[, FileOG := File %>% file.path(inputDirPath(),
                                                        .) %>% normalizePath]
    filterTrimFilesTable[, FileFT := paste0(Sample, "_", Direction,
                                            "-filtrim.fastq") %>% 
                           file.path(ftDir, .)]
    message("FilterTrim Paths:\n\n",
            filterTrimFilesTable$FileFT %>% paste0(collapse = "\n"),
            "\n\n")
    TrimTable = TrimTable() %>% copy
    # The must be ordered such that forward read is always first among pairs of files
    setorderv(filterTrimFilesTable, c("Sample", "Direction"))
    # Before returning, write the info.txt table to the filtertrim dir
    infoTabFT = copy(filterTrimFilesTable)[, .(Sample, Direction, FileFT)]
    infoTabFT[, File := basename(FileFT)]
    infoTabFT[, FileFT := NULL]
    write_table_tab(x = infoTabFT, 
                    file = file.path(ftDir, "info.txt")
    )
    return(filterTrimFilesTable)
  })
  
  prepFT = reactive({
    # should return dummy text, as non-execution workaround
    # to make FT dir and info file
    tab = ftFilesTable()
    return("filtertrim prepared.")
  })
  TrimTable = reactive({
    if(input$data_type == "ion_torrent"){
      data.table(
        Direction = c("F", "F", "R", "R"), 
        Side = c("Left", "Right", "Left", "Right"),
        Cycle = c(input$Forward[1],
                  input$Forward[2],input$Forward[1],
                  input$Forward[2])
      )
    }
    else{
      data.table(
        Direction = c("F", "F", "R", "R"), 
        Side = c("Left", "Right", "Left", "Right"),
        Cycle = c(input$Forward[1],
                  input$Forward[2],
                  input$Reverse[1],
                  input$Reverse[2])
      )
    }
  })
  ExecuteFilterTrim = reactive({
    # This stops action unless filtertrim action button has been pressed.
    # input$actionb_filtertrim %>% need(message = "...") %>% validate
    validate(need(info_file()$Sample, message = "Invalid info file!"))
    message("Filter Trim execution, iteration number:\n",
            input$actionb_filtertrim)
    # The output dir, created if need-be
    ftDir = ftDir()
    # Define the file paths of the filter-trimmed files
    filterTrimFilesTable = ftFilesTable() %>% copy
    # The table that defines the trimming params
    TrimTable = TrimTable() %>% copy
    # The must be ordered such that forward read is always first among pairs of files
    setorderv(filterTrimFilesTable, c("Sample", "Direction"))
    # The must be ordered such that forward read is always first among pairs of files
    setorderv(TrimTable, c("Direction", "Side"))
    (trimLeft <- TrimTable[(Side == "Left")]$Cycle)
    message("trimLeft:\n", trimLeft %>% paste0(collapse = ", "))
    (truncLen <- TrimTable[(Side == "Right")]$Cycle - trimLeft)
    message("truncLen:\n", truncLen %>% paste0(collapse = ", "))
    # write trimming parameters to the top-level folder
    write_table_tab(x = TrimTable, 
                    file = file.path(inputDirPath(), "trimtable.txt")
    )
    # when data = ion torrent, then cannot have reverse in the args list
    if(input$data_type == "ion_torrent"){
      # Define the arguments list that will be passed to filterAndTrim
      ftArgsList = list(
        # File I/O
        fwd = filterTrimFilesTable[(Direction == "F")]$FileOG, 
        filt = filterTrimFilesTable[(Direction == "F")]$FileFT,
        multithread = input$multithread,
        # Filtering and trimming params
        maxLen = Inf,
        minLen = 0,
        trimLeft = trimLeft,
        truncLen = truncLen,
        maxEE = input$maxEE,
        minQ = 0,
        truncQ = input$truncQ,
        maxN = 0,
        n = input$n,
        rm.phix = TRUE,
        compress = FALSE,
        verbose = TRUE)
    }
    else{
      # Define the arguments list that will be passed to filterAndTrim
      ftArgsList = list(
        # File I/O
        fwd = filterTrimFilesTable[(Direction == "F")]$FileOG, 
        filt = filterTrimFilesTable[(Direction == "F")]$FileFT,
        rev = filterTrimFilesTable[(Direction == "R")]$FileOG, 
        filt.rev = filterTrimFilesTable[(Direction == "R")]$FileFT,
        multithread = input$multithread,
        # Filtering and trimming params
        maxLen = Inf,
        minLen = 0,
        trimLeft = trimLeft,
        truncLen = truncLen,
        maxEE = input$maxEE,
        minQ = 0,
        truncQ = input$truncQ,
        maxN = 0,
        n = input$n,
        rm.phix = TRUE,
        compress = FALSE,
        verbose = TRUE)
    }
    # Save these arguments prior to running filterAndTrim
    jsonlite::write_json(x = ftArgsList,
                         path = file.path(inputDirPath(),
                                          "filterAndTrim.json"))
    # Execute the filter-trimming
    do.call(what = "filterAndTrim", args = ftArgsList)
    # return message that filter trim was executed successfully
    return("Filter & trim executed successfully!")
  })
  return(ExecuteFilterTrim())
}

##########################################################################################
# learn errors functions
# This function executes learn errors, when successful it will print out a message
learn_errors <- function(input, output, session, filt_info_file, filt_file_path){
  # paths to the files
  fastqFilesTabFT <- reactive({
    filt_info_file()
  })
  inputDirPathFT <- reactive({
    filt_file_path()
  })
  ########################################
  # Learn Errors Execution
  ########################################
  ExecuteLearnErrors = eventReactive(eventExpr = input$learn_errors_button,
    valueExpr = {
      # The FT dir, from which to read data
      ftDir = inputDirPathFT()
      fastqFilesTab = fastqFilesTabFT() %>% copy
      fastqFilesTab[, FileFull := ftDir %>% normalizePath %>% file.path(., File)]
      fastqFilesTab[, PassFilter := all(file.exists(FileFull)), by = "Sample"]
      fastqFilesTab[, Size := file.size(FileFull)]
      fastqFilesTab$FileFull %>% 
        head %>% 
        paste0(collapse = "\n") %>% 
        message("First five full paths for Learn Errors:\n", .)
                                       
      # Get the user-spec params
      learnSize = input$LE_learnSize
      minSize = 10^(input$LE_minSize)
      nReads = 10^(input$LE_nreads)
      multithread = input$LE_multithread
                                       
      # Define the samples to learn from (random)
      samplesToLearnFrom = fastqFilesTab[(Size > minSize & PassFilter), unique(Sample)]
      stopifnot(length(samplesToLearnFrom) > 1)
      learnSamples = sample(x = samplesToLearnFrom,
        replace = FALSE,
        size = min(learnSize, length(samplesToLearnFrom)))
      setorderv(fastqFilesTab, c("Sample", "Direction"))
      filesLearnForward = fastqFilesTab[(Sample %chin% learnSamples & Direction == "F")]$FileFull
      errF <- learnErrors(
        fls = filesLearnForward,
        nreads = nReads,
        multithread = multithread)
      if(input$data_type == "ion_torrent"){
        # if the data type is ion torrent, do not learn reverse error rates
        message("\n\n Saving error matrices to FT directory...\n\n")
        # Save list of error matrices
        errs = list(forward = errF)
        saveRDS(object = errs, file = file.path(ftDir, "errs.RDS"))
      }
      else{
        # Learn reverse error rates
        filesLearnReverse = fastqFilesTab[(Sample %chin% learnSamples & Direction == "R")]$FileFull
        errR <- learnErrors(
          fls = filesLearnReverse,
          nreads = nReads,
          multithread = multithread)
        message("\n\n Saving error matrices to FT directory...\n\n")
        # Save list of error matrices
        errs = list(forward = errF, reverse = errR)
        saveRDS(object = errs, file = file.path(ftDir, "errs.RDS"))
      }
    return("Learned errors successfully!")
  })
  return(ExecuteLearnErrors())
}

# Function to plot the graphs
errors_plot <- function(input, output, session, filt_info_file, filt_file_path){
  fastqFilesTabFT <- reactive({
    filt_info_file()
  })
  inputDirPathFT <- reactive({
    filt_file_path()
  })
  # Reactive holding the errors
  errs = reactive({
    ftDir = inputDirPathFT()
    errorsFile = file.path(ftDir, "errs.RDS")
    validate(need(file.exists(errorsFile), message = "errs.RDS not found"))
    message("Learned Errors being stored at:\n",
            errorsFile)
    # A reactive object holding the learned errors.
    errs <- reactiveFileReader(intervalMillis = 1000,
                               session = session,
                               filePath = errorsFile,
                               readFunc = readRDS)
    # Expect to at least have $forward errors (even if not paired)
    errs()$forward %>% need(message = "$forward missing or malformed...") %>% validate
    return(errs())
  })
  
  #function plots errors
  pErrors = reactive({
    # plot only forward when ion torrent is used
    if(input$data_type == "ion_torrent"){
      return(list(ForwardErrors = plotErrors(errs()$forward, nominalQ = TRUE)))
    }
    else{
      return(
        list(
          ForwardErrors = plotErrors(errs()$forward, nominalQ=TRUE),
          ReverseErrors = plotErrors(errs()$reverse, nominalQ=TRUE)
        )
      )
    }
  })
  return(pErrors())
}

##########################################################################################
# Run dada() functions

# function for run dada() on ion torrent sequenced fastq files
rd_ion_torrent <- function(input, output, filt_file_path, filt_info_file){
  nReads = 10^(input$RD_nreads)
  multithread = input$RD_multithread
  file_paths <- filt_info_file()$File
  FT_path <- filt_file_path()
  file_paths <- paste0(FT_path, "/", file_paths)
  file_paths <- sort(file_paths)
  
  execute_data <- eventReactive(
    eventExpr = input$run_dada_button,
    valueExpr = {
      # dereplicate
      derepF <- derepFastq(file_paths, n = nReads)
      sample.names <- filt_info_file()$Sample
      names(derepF) <- sample.names
      
      # run dada
      errs <- readRDS(paste0(filt_file_path(), "/errs.RDS"))
      errF <- errs[[1]]
      dadaF <- dada(derepF, err = errF, HOMOPOLYMER_GAP_PENALTY=-1,
                    BAND_SIZE=32, multithread=TRUE)
      # create sequence table
      seqmat <- makeSequenceTable(dadaF)
      seqmat_nochim <- removeBimeraDenovo(seqmat,
        method="consensus", multithread=TRUE, verbose=TRUE)
      
      # save as seqmat.RDS
      saveRDS(seqmat_nochim, paste0(filt_file_path(), "/seqmat.RDS"))
      
      # Create an abundance table to show on app
      abundance_tab <- melt(seqmat_nochim)
      colnames(abundance_tab) <- c("Sample", "Sequence", "Abundance")
      saveRDS(abundance_tab, paste0(filt_file_path(), "/dadaTabBimeraFilt.RDS"))
      
      return("Ran DADA successfully!")
    }
  )
  return(execute_data())
}

#' Wrapper for running DADA2 algorithm from sequence file to sequence result.
#' 
#' Check that this isn't redundant with recent additions. Migrate if so.
#'
wrap_dada2_workflow = function(seqFiles,
                               dadaOutFiles = c("DADA2-Forward.RDS", "DADA2-Reverse.RDS"),
                               err = NULL,
                               selfConsist = TRUE,
                               minOverlap = 20,
                               maxMismatch = 0,
                               # Performance params
                               nReads = 1e6,
                               multithread = TRUE){
  require("dada2")
  merged = NULL
  
  stopifnot(all(file.exists(seqFiles)))
  seqFileF = seqFiles[1]
  seqFileR = seqFiles[2]
  
  stopifnot(
    dir.exists(
      dirname(
        c(dadaOutFiles[1],
          dadaOutFiles[2]))))
  
  dadaOutFileF = dadaOutFiles[1]
  dadaOutFileR = dadaOutFiles[2]
  
  if(length(err) > 2){
    warning("Provided more than two error matrices. Most likely something is wrong.")
  }
  
  if(length(err) == 2){
    # Assume in forward-then-reverse order
    errF = err[[1]]
    errR = err[[2]]
  } else {
    errF = errR = err[[1]]
  }
  dadaF = dadaR = derepF = derepR = NULL
  
  message("Dereplicating forward reads:\n", seqFileF, "\n")
  derepF <- derepFastq(seqFileF, n = nReads)
  message("DADA2-ing:\n", seqFileF, "\n")
  dadaF <- dada(derepF, err=errF, selfConsist = selfConsist, multithread = multithread)
  saveRDS(dadaF, dadaOutFileF)
  
  message("Dereplicating reverse reads:\n", seqFileR, "\n")
  derepR <- derepFastq(seqFileR, n = nReads)
  message("DADA2-ing:\n", seqFileR, "\n")
  dadaR <- dada(derepR, err=errR, selfConsist = selfConsist, multithread = multithread)
  saveRDS(dadaR, dadaOutFileR)
  
  message("Merging DADA2 results for read-pairs:\n",
          paste0(seqFiles, collapse = "\n"), "\n")
  # merger <- mergePairs(ddF, derepF, ddR, derepR)
  trash = try(expr = {
    ## Merge seq directions by ID, return an abundance data.table
    merged = dada2:::mergePairsByID(
      # Forward
      dadaF = dadaF,
      derepF = derepF,
      srF = seqFileF,
      # Reverse
      dadaR = dadaR,
      derepR = derepR,
      srR = seqFileR,
      # Additional params
      minOverlap = minOverlap,
      maxMismatch = maxMismatch,
      returnRejects = FALSE,
      verbose = TRUE)
    message("Sum of read-pairs properly merged after denoising:\n",
            round(100 * merged[, sum(abundance[(accept)])/sum(abundance)],
                  digits = 1), "%")
    merged <- merged[(accept & !is.na(sequence)), .(sequence, abundance)]
    saveRDS(merged, file = file.path(dirname(dadaOutFileF), 
        gsub("\\.RDS$", "merge.RDS", basename(dadaOutFileF))))
  }, silent = TRUE)
  return(merged)
}

# dada function taken from dada2 docker shiny app
run_dada <- function(input, output, filt_file_path, filt_info_file_path){
  # paths to the files
  infoFilePathRD <- reactive({
    filt_info_file_path()
  })
  inputDirPathRD <- reactive({
    filt_file_path()
  })
  ########################################
  # DADA Execution
  ## DADA2::dada()
  ########################################
  # Interpret relevant inputs to update info table.
  fastqFilesTabRD = reactive({
    fastqFilesTab = fread(input = infoFilePathRD())
    # Minimum file size
    minSize = 10^(input$RD_minSize)
    # The FT dir, from which to read data
    ftDir = inputDirPathRD()
    fastqFilesTab = fastqFilesTab %>% copy
    fastqFilesTab[, FileFull := ftDir %>% normalizePath %>% file.path(., File)]
    fastqFilesTab[, PassFilter := all(file.exists(FileFull)), by = "Sample"]
    fastqFilesTab[, Size := file.size(FileFull)]
    fastqFilesTab[, PassFilter := PassFilter & Size > minSize, by = "FileFull"]
    # Filter the files that don't pass the requirements of existence, and minimum size
    fastqFilesTabDADA2 = copy(fastqFilesTab)[(PassFilter)]
    
    samplesLost = nrow(fastqFilesTab) - nrow(fastqFilesTabDADA2)
    samplesLostPerc = 100 * (samplesLost / nrow(fastqFilesTab)) %>% round(digits = 1)
    # Status messages
    message("\n\nDADA2 start:\n Of the ", nrow(fastqFilesTab), " input files\n",
            samplesLost, " (", samplesLostPerc,
            "%) were lost due to existence-check and size filter.\n\n")
    fastqFilesTabDADA2$FileFull %>% 
      head %>% 
      paste0(collapse = "\n") %>% 
      message("First five full paths for Running DADA:\n", .)
    # Send it along to outputs
    return(fastqFilesTabDADA2)
  })
  # RUN DADA2
  ExecuteDADA2 = eventReactive(
    eventExpr = input$run_dada_button,
    valueExpr = {
      message(".\n.\nRunning DADA2...\n.\n.")
      ftDir = inputDirPathRD()
      # Fail early if errors file not present
      errorsFile = file.path(ftDir, "errs.RDS")
      file.exists(errorsFile) %>% 
        need(message = "No learned errors file: `errs.RDS`!\nSee `Learn Errors` tab") %>% 
        validate
      errs = readRDS(errorsFile)
      
      # Run-related user-spec params
      nReads = 10^(input$RD_nreads)
      multithread = input$RD_multithread
      
      fastqFilesTabDADA2 = fastqFilesTabRD() %>% copy
      
      stopifnot(nrow(fastqFilesTabDADA2) > 0)
      
      # Set the dada2 cache directory, helpful for debugging
      dadaCacheDir = file.path(ftDir, "dadaCache")
      if(!dir.exists(dadaCacheDir)){
        dir.create(dadaCacheDir)
      }
      # Names of dada cache files
      fastqFilesTabDADA2[, FileDadaCache := file.path(dadaCacheDir,
            paste0(Sample, "-", Direction, "-dada2.RDS"))]
      # RUN DADA via wrapper function.
      time0 = Sys.time()
      # Run on each sample read-pair
      setorderv(fastqFilesTabDADA2, c("Sample", "Direction"))
      setkeyv(fastqFilesTabDADA2, "Sample")
      # The multi-threaded one-sample-at-a-time approach:
      mergeTab = fastqFilesTabDADA2[(PassFilter),
                                    wrap_dada2_workflow(
                                      seqFiles = FileFull,
                                      dadaOutFiles = FileDadaCache,
                                      err = errs,
                                      # No need for selfConsist if we are confident
                                      # in convergence and consistency of our error model
                                      selfConsist = FALSE,
                                      # performance params
                                      multithread = multithread,
                                      nReads = nReads),
                                    by = c("Sample")]
      
      timeDADA2 = (time0 - Sys.time())
      
      saveRDS(mergeTab,
              file = file.path(ftDir, "dadaTab.RDS"))
      
      ####################
      # Remove Chimeras
      ####################
      seqtab = dcast.data.table(
        data = mergeTab,
        formula = Sample ~ sequence, 
        fun.aggregate = sum, 
        fill = 0,
        value.var = "abundance")
      seqmat <- as(seqtab[, -1L, with = FALSE], "matrix")
      rownames(seqmat) <- seqtab$Sample
      # Chimera filter on the whole table (added parameters from dada2 tutorial)
      seqmat <- removeBimeraDenovo(seqmat, method="consensus",
                                   verbose=TRUE, multithread = NCores)
      dim(seqmat)
      # saveRDS(seqmat, "seqmat.RDS")
      # Filter chimeras from mergeTab
      setkey(mergeTab, sequence)
      mergeTabNoChimera = mergeTab[colnames(seqmat)]
      Nseq = mergeTab[, uniqueN(sequence)]
      NseqNoChimera = mergeTabNoChimera[, uniqueN(sequence)]
      message("\nNumber of sequences determined to be chimeras (de novo):\n",
              (Nseq - NseqNoChimera), " out of ", Nseq, " denoised sequences (", 
              round(100*(Nseq - NseqNoChimera)/Nseq, digits = 1), "%)\n\n"
      )
      saveRDS(seqmat, 
              file = file.path(ftDir, "seqmat.RDS"))
      saveRDS(mergeTabNoChimera, 
              file = file.path(ftDir, "dadaTabBimeraFilt.RDS"))
      # return(
      #   paste0("Time to execute DADA Workflow:\n", timeDADA2)
      # )
      return("Ran DADA successfully!")
    })
  return(ExecuteDADA2())
}

# Function for reading the created sequence table
sequence_table_function <- function(input, output, session, filt_file_path){
  # paths to the files
  inputDirPathRD <- reactive({
    filt_file_path()
  })
  # Reactive holding the dada results
  DADAS = reactive({
    dadasFile = file.path(inputDirPathRD(), "dadaTabBimeraFilt.RDS")
    # Validate that file exists in order to show plot.
    file.exists(dadasFile) %>% need(message = "...") %>% validate
    # message the errors file path
    message("Learned Errors being stored at:\n",
            dadasFile)
    # A reactive object holding the learned errors.
    dadas <- reactiveFileReader(intervalMillis = 1000,
                                session = session,
                                filePath = dadasFile,
                                readFunc = readRDS)
    return(dadas())
  })
  dadas_table <- reactive({
    sequence_table_file <- paste0(inputDirPathRD(), "/dadaTabBimeraFilt.RDS")
    readRDS(sequence_table_file)
  })
  return(dadas_table())
}

##########################################################################################
# functions for assign taxonomy

assign_taxonomy <- function(input, output, filt_file_path, file_path){
  # get the file from dada: /seqmat.RDS. Should be in same FT directory
  seqtab_nochim_path <- paste0(filt_file_path(), "/seqmat.RDS")
  seqtab_nochim <- readRDS(seqtab_nochim_path)
  if(!is.null(input$ref_seq_file) && input$ref_seq == "upload"){
    path <- input$ref_seq_file
    taxa <- assignTaxonomy(seqtab_nochim, path$datapath)
  }
  else{
    validate(need(input$ref_seq != "none", input$ref_seq_file,
                  message = "No reference selected"))
  }
  # save the object as rds in the filtered folder dir
  saveRDS(taxa, paste0(filt_file_path(), "/taxa.RDS"))
  return(taxa)
}

assign_species <- function(input, output, filt_file_path, taxa, file_path){
  seqtab_nochim_path <- paste0(filt_file_path(), "/seqmat.RDS")
  seqtab_nochim <- readRDS(seqtab_nochim_path)
  if(!is.null(input$ref_species_file) && input$ref_species == "upload"){
    path <- input$ref_species_file
    taxa <- addSpecies(taxa, path$datapath)
  }
  else{
    validate(need(input$ref_species != "none", input$ref_species_file,
                  message = "No reference selected"))
  }
  # save the species taxa in the filtered folder dir
  saveRDS(taxa, paste0(filt_file_path(), "/taxa_species.RDS"))
  return("Added species successfully!")
}

# create phyloseq object
create_ps_object <- function(input, output, filt_file_path, taxa){
  # get otu table (need to convert abundance to numeric)
  seqtab_nochim_path <- paste0(filt_file_path(), "/seqmat.RDS")
  seqtab_nochim <- readRDS(seqtab_nochim_path)
  # get sample data
  sample_file <- input$metadata_file
  sample_df <- read.delim(sample_file$datapath)
  # make the row names the sample names
  row.names(sample_df) <- sample_df[,1]
  # get tax table from taxa
  # create ps object
  ps <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows = FALSE), sample_data(sample_df),
                 tax_table(taxa))
  return(ps)
}
##########################################################################################