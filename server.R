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


source("plot.R")




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
    if(str_detect(inFile$name, ".csv")){
      data <- read.csv(inFile$datapath, sep=";")
    }else if(str_detect(inFile$name, ".gff3")){
      data <- load_Gff(inFile$datapath)
    }
    data
  })
  
  observeEvent(input$do, {
    output$plot <- renderPlot({
      
      isolate(RNAplot(isolate(RNA()),isolate(mapD()), isolate(input$start), isolate(input$end), 
                      alpha = isolate(input$alpha),
                      graph_size = isolate(input$graph_size),
                      #color = isolate(input$color),
                      subgenes = isolate(input$subgenes),
                      label = isolate(input$label),
                      line_visible = isolate(input$line_visible),
                      incom_genes = isolate(input$incom_genes),
                      max_read = isolate(input$max_read),
                      arrow_body_height = isolate(input$arrow_body_height),
                      arrowhead_height = isolate(input$arrowhead_height),
                      arrowhead_width = isolate(input$arrowhead_width),
                      filter = isolate(input$filter)))
    })
    
  })
  
  
  output$table <- DT::renderDT({
    mapD2 <- mapD() %>% select(molecule, gene, start, end, orientation)
    DT::datatable(mapD2)#, editable = TRUE)
  })
  output$foo = downloadHandler(
    filename = 'test.png',
    content = function(file) {
      device <- function(..., width, height) {
        grDevices::png(..., width = width, height = height,
                       res = 300, units = "in")
      }
      R <- RNAplot(isolate(RNA()),isolate(mapD()), isolate(input$start), isolate(input$end), 
                   alpha = isolate(input$alpha),
                   graph_size = isolate(input$graph_size),
                   #color = isolate(input$color),
                   subgenes = isolate(input$subgenes),
                   label = isolate(input$label),
                   line_visible = isolate(input$line_visible),
                   incom_genes = isolate(input$incom_genes),
                   max_read = isolate(input$max_read),
                   arrow_body_height = isolate(input$arrow_body_height),
                   arrowhead_height = isolate(input$arrowhead_height),
                   arrowhead_width = isolate(input$arrowhead_width),
                   filter = isolate(input$filter))
      ggsave(file, plot = R, device = device)
    })
}
