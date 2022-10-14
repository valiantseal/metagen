# map needs on the server: sudo apt install libgdal-dev
# to find capitals coordinates uses free https://simplemaps.com/data/world-cities

#setwd("/home/ubuntu/github/shinyApp/shinyVizApps/uploadDatCov")




library(shiny)
library(DT)

allSamples<-read.csv('/home/ubuntu/ICMC/summaries/combined/allSamples.csv')


linebreaks <- function(n){HTML(strrep(br(), n))}  # introduce multiple breaks in one function



ui <- navbarPage("Summary",
                 tabsetPanel(
                   tabPanel("Summary", fluid=T,
                            fluidPage(
                              mainPanel( 
                                selectInput(inputId = "summaryVar1", label= "Column 1", choices = c(colnames(allSamples)[2:ncol(allSamples)]), selected='Name'), 
                                selectInput(inputId = "summaryVar2", label= "Column 2", choices = c(colnames(allSamples)[2:ncol(allSamples)]), selected='Project'),
                                div(DT::dataTableOutput("NameProjectSummary", width = "100%", height = "auto"), style = "font-size:115%"),
                                linebreaks(3),
                              ))),
                   
                   tabPanel("All Samples", fluid=T,
                            fluidPage(
                              mainPanel(
                                DT::dataTableOutput("AllSamples", width = "150%"),
                                linebreaks(3),
                              ))),
                 ))



server <- function(input, output) {
  
  
  tableSummary<-reactive({
    df<-data.frame(table(allSamples[,input$summaryVar1], allSamples[,input$summaryVar2]))
    sumFiltr<-df[!(df$Freq==0),]
    return(sumFiltr)
  })
  
  output$NameProjectSummary <- DT::renderDataTable(tableSummary(), 
    options = list(scrollX = TRUE), rownames= FALSE)
  
  
  # table with all the samples
  output$AllSamples<- DT::renderDataTable(allSamples, 
    options = list(scrollX = TRUE), rownames= FALSE)
  
}

shinyApp(ui, server)
  
  
