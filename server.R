library(stringr)
library(ggplot2)
library(gggenes)
library(aplot)
library(ggeasy)
library(RColorBrewer)
library(labeling) #eventually not necessary
library(rebus)    #eventually not necessary
library(shiny)
library(DT)
library(dplyr)
library(shinyWidgets)
library(shinycssloaders)
library(colourpicker)
library(svglite)
library(tidyr)
library(cowplot)
library(grDevices)
source("plot.R")


#setwd("Y:/exchange/Marcus/M_pro/RNAplotter2")

server <- function(input, output) {
  
  output$whichplots <- renderUI({
    N <- plotName()
    
    tagList(
      awesomeCheckboxGroup("filter",  label = "Which graphs should be displayed?", 
                           choices = N,
                           status = "danger", selected = N))
  })
  
  output$whichcolor <- renderUI({
    
    N <- plotName()
    col <- col1[1:length(N)]    #create colors for RNA-reads
    
    tagList(
      lapply(1:length(N), function(i) {
        colourInput(str_c("P",i), label=N[i], showColour = "both", value = col[i])}))
  })
  
  output$whichplot <- renderUI({
    
    #textOutput(str_c(plotName(), collapse ="_"))
    withSpinner(plotOutput("plot", width = str_c(input$width,"in"), height = str_c(input$height,"in")), type=6, hide.ui = FALSE)
  })
  
  
  RNA <- reactive({
    inFile1 <- input$filefwd
    inFile2 <- input$filerev
    if (is.null(inFile1) & is.null(inFile2)) return(NULL)
    Rf <- read.table(inFile1$datapath, quote="\"", comment.char="#")#load first grp-file
    if(all(Rf[,1]==1:nrow(Rf))){Rf <- Rf[,2:ncol(Rf)]
    q1 <- TRUE}else{ q1<- FALSE} 
    if(is.na(suppressWarnings(as.numeric(Rf[1,1])))){Rf <- as.data.frame(Rf[,2:ncol(Rf)])}
    Rr <- read.table(inFile2$datapath, quote="\"", comment.char="#")#load second grp-file
    if(all(Rr[,1]==1:nrow(Rr))){Rr <- Rr[,2:ncol(Rr)]} 
    if(is.na(suppressWarnings(as.numeric(Rr[1,1])))){Rr <- as.data.frame(Rr[,2:ncol(Rr)])}
    if(min(Rr)>=0){Rr <- Rr*(-1)}
    
    
    if (is.null(input$Name)){
      n1 <- readLines(inFile1$datapath, 10)
      n2 <- readLines(inFile2$datapath, 10)
      
      n1 <- n1[which(str_detect(n1, START%R%"#"))[1]]
      n2 <- n2[which(str_detect(n2, START%R%"#"))[1]]
      n0 <- str_sub(n1, str_locate(n1, WRD)[,1],-1)
      n0 <- str_extract_all(n0, one_or_more(WRD))[[1]]
      if(q1 & n0[1] %in% c("BASE", "base")){n0 <- n0[2:length(n0)]}
      
      if(n1==n2 & length(n0)==ncol(Rf)){
        N <- n0
      }else{N <- str_c("Plot_", 1:ncol(Rf))}
    }else{N <- readLines(input$Name$datapath)
    N <- str_replace(N, or(one_or_more(SPC), "\t")%R%END,"")}
    colnames(Rf) <- str_c(N,"fwd")
    colnames(Rr) <- str_c(N,"rev")
    cbind(Rf,Rr) 
  })
  
  plotName  <-  reactive({
    if(is.null(input$filefwd)){
      return(NULL)
    }else{
      N <- unique(str_sub(colnames(RNA()), 1, -4))
      return(N)}
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
    
    N <- plotName()
    
    c(input$P1, input$P2, input$P3, input$P4, input$P5, input$P6, input$P7, input$P8,
      input$P9, input$P10, input$P11, input$P12, input$P13, input$P14, input$P15, input$P16,
      input$P17, input$P18, input$P19, input$P20, input$P21, input$P22, input$P23, input$P24)[1:length(N)]
  })
  

  observeEvent(input$do, {
    output$plot <- renderPlot({
      isolate({
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
    })
    
  })
  
  
  output$table <- DT::renderDT({
    
    DT::datatable(mapD())#, editable = TRUE)
  })
  
  
  
  
  
  #output$foo <- downloadHandler(
  # filename = function() {"RNAplot.png"},
  #content = function(file){
  # 
  #ggsave(file, plot = R(), width = input$width, height = input$height, units = "in", 
  #      device = "png")
  #})
  output$foo <- downloadHandler(
    filename = function() {
      paste("RNAplot.", input$download_type, sep="")
    },
    content = function(file) {
      data_to_save <- RNAplot(RNA(),mapD(), input$start, input$end, 
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
      
      if(input$download_type=="svg"){
        save_plot(file, plot = data_to_save , base_width = input$width*1.25, base_height = input$height*1.25, units = "in", 
                  device = input$download_type, fix_text_size = FALSE)
      }else{
        save_plot(file, plot = data_to_save , base_width = input$width*1.25, base_height = input$height*1.25,
                  units = "in", device = input$download_type)
      }
    })
  
  
  
  
  
  output$tabgo <- downloadHandler(
    filename = function(){"Map.csv"}, 
    content = function(fname){
      write.csv2(mapD(), fname)
    }
  )
  
  
  G <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    for(i in 1:length(inFile$datapath)){
      if(str_detect(inFile$name[i], DOT%R%or("grp", "txt")%R%END)){
        Rf <- load_grp(inFile$datapath[i])
      }else if(str_detect(inFile$name[i], DOT%R%"bedgraph"%R%END)){
        Rf <- load_bedgraph(inFile$datapath[i], input$Replicon)
      }
      colnames(Rf) <- str_c(str_sub(inFile$name[i], 1, str_locate(inFile$name[i], DOT%R%one_or_more(WRD)%R%END)[,1]-1),
                            "_Plot", 1:ncol(Rf))
      if(i==1){
        R <- Rf
      }else{ R <- cbind(R, Rf)}
    }
    R
  })
  
  output$mod_file <- downloadHandler(
    filename = function() {str_c(input$com_red_choice, "_", str_c(input$red_which, collapse = "_"), ".grp")},
    content = function(file){
      if(input$com_red_choice=="Combine"){Gl <- G()
      }else{Gl <- Red()
      Gl <- Gl[, which(colnames(Gl)%in% input$red_which)]
      }
      G0 <- unlist(unite(Gl, col= "all", colnames(Gl), sep= "\t"))
      G0 <- c(str_c("#", str_c(colnames(Gl), collapse= "\t")) , G0)
      writeLines(G0,file)
    })
  
  output$which_replicon <- renderUI({
    if(is.null(input$file1) | !any(str_detect(input$file1$name, DOT%R%"bedgraph"%R%END))){
      N<- ""
    }else{
      l <- input$file1$datapath[which( str_detect(input$file1$name, DOT%R%"bedgraph"%R%END))]
      H <-  read.delim(l[1], header=FALSE)
      N <- unique(H$V1)
    }
    tagList(
      if(is.null(input$file1)| !any(str_detect(input$file1$name, DOT%R%"bedgraph"%R%END))){
        ""
      }else{
        selectInput("Replicon", label = "Which Replicon should be used?", choices = N)
      })
  })
  
  Red <- reactive({
    inFile1 <- input$file_red
    if (is.null(inFile1)) return(NULL)
    Rf <- load_grp(inFile1$datapath)
    
    r <- readLines(inFile1$datapath, n = 100)
    g <- length(str_extract_all(r[100], one_or_more(DGT))[[1]])
    N0 <- str_extract_all(r[1], one_or_more(WRD))[[1]]
    Nn <- str_extract_all(r[1], one_or_more(DGT))[[1]]
    r0 <- as.integer(str_extract(r[10:100], one_or_more(DGT)))
    base <- all(r0 == r0[1]+0:90)
    name <- g == length(N0) & all(N0 != Nn)
    
    if(base & !name){         N <- str_c("Plot_", 1:(length(N0)))
    }else if(!base & !name){  N <- str_c("Plot_", 1:(length(N0)))
    }else if(base & name){    N <- N0[2:(length(N0))]
    }else if(!base & name){   N <- N0}
    colnames(Rf) <- N
    Rf
  })
  
  
  
  output$which_red_plot <- renderUI({
    if(is.null(input$file_red)){
      N<- ""
    }else{
      N <- colnames(Red())
    }
    
    tagList(
      if(is.null(input$file_red)){
        ""
      }else{
        checkboxGroupButtons(
          inputId = "red_which",
          status = "info", 
          label = "Which plots do you want to keep?",
          choices = N,
          direction = "vertical")
      })
  })
  
  output$Tutorial <- downloadHandler(
    filename = "RNAplotter_manual.pdf",
    content = function(file) {
      file.copy("RNAplotter_manual.pdf", file)
    }
  )
  
}

