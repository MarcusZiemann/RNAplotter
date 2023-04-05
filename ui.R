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
library(shinythemes)
library(shinycssloaders)
library(colourpicker)



source("plot.R")

###shiny
options(shiny.maxRequestSize=200*1024^2)

ui <-fluidPage(
  titlePanel("RNAplotter"),
  theme = shinythemes::shinytheme("spacelab"),
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        tabPanel("Main",
                 fileInput("filefwd", "Choose fwd RNA-file", multiple = FALSE,
                           accept = ".grp"),
                 fileInput("filerev", "Choose rev RNA-file", multiple = FALSE,
                           accept = ".grp"),
                 fileInput("Map", "enter map file [csv/gff3]", multiple = FALSE,
                           accept = c(".csv", ".gff3")),
                 fileInput("Name", "enter names of RNAreads", multiple = FALSE),
                 numericInput("start", "Please enter the start of your read", 1, step = 100),
                 numericInput("end", "Please enter the end of your read", 1e3, step = 100),
                 actionButton("do", "Plot"),
                 downloadButton('foo'),
                 sliderInput("width", "Width:", min = 3, max=30, value =7),
                 sliderInput("height", "Height:", min = 3, max=30, value =7)),
        tabPanel("optional",
                 materialSwitch(inputId = "Gfill", label = "RNAreads Graphs filled?", value= TRUE),
                 conditionalPanel(condition = "input.Gfill",
                                  numericInput("alpha", "Transparency:", 0.8, step = 0.1)),
                 conditionalPanel(condition = "!input.Gfill",
                                  numericInput("Gsize", "Graphsize:", 1.2, step = 0.1)),
                 numericInput("max_read", "max. number of reads:", NA, step = 100),
                 prettySwitch("subgenes", label = "display subgenes", value= TRUE),
                 prettySwitch("label", label = "change labelsize?", value= FALSE),
                 conditionalPanel(condition = "input.label",
                                  numericInput("msize", "label fontsize:", value=4)),
                 conditionalPanel(condition = "input.label",
                                  numericInput("ntlength", "only label genes longer than [nt]:", value=NULL)),
                 prettySwitch("line_visible", label = "line in maps?", value= TRUE),
                 uiOutput("whichplots"),
                 textInput("incom_genes", "ending of incomplete genes:", 
                           value = "_partial", width = NULL, placeholder = NULL),
                 numericInput("graph_size", "map-graph ratio:", 3),
                 numericInput("arrow_body_height", "height of arrowbody:", 7),
                 numericInput("arrowhead_height", "height of arrow:", 10),
                 numericInput("arrowhead_width", "width of arrow:", 8)),
        tabPanel("color",
                 uiOutput("whichcolor")),
        
      )),
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", uiOutput("whichplot")),
        tabPanel("table", DT::DTOutput('table'),
                 downloadButton('tabgo', label = "Download table"))
      )
    )
  )
)

