library(shiny)
library(corrplot)
library(shinyalert)
#install.packages("shinythemes")
#install.packages('rsconnect')
#install.packages('shinyalert')
library(shinythemes)
library(shinyBS)

CoRSIV <- read.csv("./CoRSIV_9926_gencode28.csv")

fluidPage(
  tags$head(includeHTML(("google-analytics.html")),
    tags$style(HTML("
      .shiny-output-error-validation {
        color: red;
        font-size: 200%;
        font-weight: bold;
      }
    "))
  ),
 theme = shinytheme("simplex"),
  titlePanel("CoRSIV Plotter"),
  
  sidebarLayout(
    
    sidebarPanel(

      selectInput("select", h3("Option 1: Select an autosome to download CoRSIV plot/s"),
                   choices = list("Chromosome 1" = 1, "Chromosome 2" = 2,
                                  "Chromosome 3" = 3, "Chromosome 4" = 4,
                                  "Chromosome 5" = 5, "Chromosome 6" = 6,
                                  "Chromosome 7" = 7, "Chromosome 8" = 8,
                                  "Chromosome 9" = 9, "Chromosome 10" = 10,
                                  "Chromosome 11" = 11, "Chromosome 12" = 12,
                                  "Chromosome 13" = 13, "Chromosome 14" = 14,
                                  "Chromosome 15" = 15, "Chromosome 16" = 16,
                                  "Chromosome 17" = 17, "Chromosome 18" = 18,
                                  "Chromosome 19" = 19, "Chromosome 20" = 20,
                                  "Chromosome 21" = 21, "Chromosome 22" = 22),selected = 1),
      
      downloadButton("downloadData", label = "Download"),br(),
      em("Please unzip the downloaded file and open the PDF file/s to visualize the CoRSIVs"),br(),
      br(),
      br(),
      textInput("gene_sym", h3("Option 2: Provide a gene symbol to download CoRSIVs in +/- 1Mb of the gene"),value = "SPATC1L"),
       downloadButton("downloadData3", label = "Download"),br(),
      em("Please open the downloaded PDF file to visualize CoRSIVs"),br(),
      
        br(),br(),
      textInput("ucsc_coord", h3("Option 3: Provide a genomic coordinate(hg38) to download CoRSIVs in +/- 1Mb of the genomic coordinate"),value = "6:32050000"),
      downloadButton("downloadData2", label = "Download"),br(),
      em("Please open the downloaded PDF file to visualize CoRSIVs"),br()
      
    ),
    
    mainPanel(
      h1("A Genomic Atlas of Systemic Interindividual Epigenetic Variation in Humans"),
      br(),
      tabsetPanel(type="tabs",tabPanel("Main",

      div(img(src = "intro.png", height = 225, width = 410), style="text-align: center;"),

      p("CoRSIV* : ",em("Correlated Regions of Systemic Interindividual Variation")),

      span(textOutput("gene_symbol"),style=c("color:black;font-size:200%")),
      tableOutput("coolplot2"),
      span(textOutput("coodinate"),style=c("color:black;font-size:200%")),
      tableOutput("coolplot")
      ),
      tabPanel("Help",div(img(src = "help.png", height = 1000, width = 1000), style="text-align: center;"),textOutput("help")))
       
    )
  )
)