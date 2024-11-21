library(stringr)
library(ggplot2)
library(gggenes)
library(aplot)
library(ggeasy)
library(RColorBrewer)
library(labeling) #eventually not necessary
library(rebus)    #eventually not necessary
library(stringr)
library(shiny)
library(DT)
library(dplyr)
library(shinyWidgets)
library(shinycssloaders)
library(colourpicker)
library(shinyWidgets)
library(svglite)

source("plot.R")


#setwd("Y:/exchange/Marcus/M_pro/RNAplotter2")

server <- function(input, output) {
  
  output$whichplots <- renderUI({
    if (is.null(input$Name)){
      N <- str_c("Plot_", 1:(ncol(RNA())/2))
    } else{N <- readLines(input$Name$datapath)}
    
    tagList(
      awesomeCheckboxGroup("filter",  label = "Which graphs should be displayed?", 
                           choices = N,
                           status = "danger", selected = N))
  })
  
  output$whichcolor <- renderUI({
    if (is.null(input$Name)){
      N <- str_c("Plot_", 1:(ncol(RNA())/2))
    } else{N <- readLines(input$Name$datapath)}
    
    col <- col1[1:length(N)]    #create colors for RNA-reads
    
    tagList(
      lapply(1:length(N), function(i) {
        colourInput(str_c("P",i), label=N[i], showColour = "both", value = col[i])}))
  })
  
  output$whichplot <- renderUI({
    withSpinner(plotOutput("plot", width = str_c(input$width,"in"), height = str_c(input$height,"in")), type=6, hide.ui = FALSE)
  })
  
  
  RNA <- reactive({
    inFile1 <- input$filefwd
    inFile2 <- input$filerev
    if (is.null(inFile1) & is.null(inFile2)) return(NULL)
    Rf <- read.delim(inFile1$datapath, header=FALSE, comment.char="#")                 #load first grp-file
    if(ncol(Rf)==1){ Rf <- read.table(inFile1$datapath, quote="\"", comment.char="")}
    if(all(Rf[,1]==1:nrow(Rf))){Rf <- Rf[,2:ncol(Rf)]} 
    Rr <- read.delim(inFile2$datapath, header=FALSE, comment.char="#")                 #load second grp-file
    if(ncol(Rr)==1){ Rr <- read.table(inFile2$datapath, quote="\"", comment.char="")}
    if(all(Rr[,1]==1:nrow(Rr))){Rr <- Rr[,2:ncol(Rr)]}
    if (is.null(input$Name)){
      N <- str_c("Plot_", 1:ncol(Rf))
    } else{N <- readLines(input$Name$datapath)}
    colnames(Rf) <- str_c(N,"fwd")
    colnames(Rr) <- str_c(N,"rev")
    cbind(Rf,Rr) 
  })
  
  mapD <- reactive({
    inFile <- input$Map
    if (is.null(inFile)) return(NULL)
    if(str_detect(inFile$name, ".csv"%R%END)){
      data <- load_csv(inFile$datapath)
    }else if(str_detect(inFile$name, ".gff3"%R%END)){
      data <- load_Gff3(inFile$datapath)
    }else if(str_detect(inFile$name, ".gff"%R%END)){
      data <- load_Gff(inFile$datapath)
    }
    data
  })
  
  
  #  observeEvent(input$table_cell_edit, {
  #   mapD <<- editData(mapD(), input$table_cell_edit, 'table')
  #})
  
  
  output$plot <- renderPlot(NULL)
  
  col <- reactive({
    if (is.null(input$Name)){
      N <- str_c("Plot_", 1:(ncol(RNA())/2))
    } else{N <- readLines(input$Name$datapath)}
    c(input$P1, input$P2, input$P3, input$P4, input$P5, input$P6, input$P7, input$P8,
      input$P9, input$P10, input$P11, input$P12, input$P13, input$P14, input$P15, input$P16,
      input$P17, input$P18, input$P19, input$P20, input$P21, input$P22, input$P23, input$P24)[1:length(N)]
  })
  
  R <- reactive({
    RNAplot(RNA(),mapD(), input$start, input$end, 
            alpha = input$alpha,
            graph_size = input$graph_size,
            subgenes = isolate(input$subgenes),
            line_visible = input$line_visible,
            incom_genes = input$incom_genes,
            color = col(),
            max_read = input$max_read,
            arrow_body_height = input$arrow_body_height,
            arrowhead_height = input$arrowhead_height,
            arrowhead_width = input$arrowhead_width,
            filter = input$filter, 
            Gfill= input$Gfill,
            Gsize= input$Gsize,
            msize= input$msize,
            Mwidth = input$width,
            ntlength = input$ntlength)
  })
  
  
  
  observeEvent(input$do, {
    output$plot <- renderPlot({
      
      isolate({R()})
    })
    
  })
  
  
  output$table <- DT::renderDT({
    
    DT::datatable(mapD(), editable = TRUE)
  })
  output$foo = downloadHandler(
    filename = function() {"RNAplot.png"},
    content = function(file){
      
      ggsave(file, plot = R(), width = input$width, height = input$height, units = "in", 
             device = "png")
    })
  output$tabgo <- downloadHandler(
    filename = function(){"Map.csv"}, 
    content = function(fname){
      write.csv2(mapD(), fname)
    }
  )
}

