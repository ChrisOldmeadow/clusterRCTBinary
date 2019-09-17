#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
require(ggplot2)
require(DT)
require(lme4)
require(lmerTest)
source("utils.R")

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("Power for cluster RCT, continuous outcomes"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
      sidebarPanel(
        h4("Select group proportions"),
        numericInput(inputId = "p1",
                     label = "Outcome proportion in the intervention group:",
                     value = 0.11),
        numericInput(inputId = "p2",
                     label = "Outcome proportion in the control group:",
                     value = 0.03),
        
        h4 ("Select number of clusters and ICC"),
        # Input: Numeric entry for total sample size ----
        numericInput(inputId = "k1",
                     label = "Number of intervention clusters:",
                     value = 8),
        numericInput(inputId = "k2",
                     label = "Number of control clusters:",
                     value = 10),
        numericInput(inputId = "rho",
                     label = "Intra-class correlation:",
                     value = 0.01, min=0.01, max=0.2),
        h4 ("Select values of average cluster size to range over"),
        numericInput(inputId = "start",
                     label = "min:",
                     value = 5, min=1, max=1000),
        numericInput(inputId = "end",
                     label = "max:",
                     value = 15, min=1, max=1000),
        numericInput(inputId = "by",
                     label = "by:",
                     value = 1, min=1, max=100),
        
        numericInput(inputId = "alpha",
                     label = "Type 1 error rate",
                     value = .05, min=.001, max = 0.1),
        # Input: Numeric entry for number of simulations ----
        numericInput(inputId = "nreps",
                     label = "Number of simulations",
                     value = 500, min=1, max=1500),
        # Input: Select the effects ----
        #radioButtons("type", "Effect type:",
        #             c("Main effect Factor 1" = "mainF1",
        #               "Main effect Factor 2" = "mainF2",
        #               "Interaction" = "F1F2")),
        # run simulation
        h4("Run the Power Simulations"),
        actionButton(inputId = "submit", label = "Submit")
        
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        #tags$p(tags$h3("Selected Data")),
        #plotOutput("meanPlot"),
        tags$p(tags$h3("Power Results")),
        DT::dataTableOutput(outputId = "resultsTable")
      )
    )
)


# Define server logic required to draw a histogram
server <- function(input, output) {
  # SIMULATE 
  observeEvent(input$submit, {
    withProgress(message="Running Power Simulation", value=0, {
      res <- powersims(input$p1,input$p2,
                       input$k1,input$k2,
                       input$rho,
                       input$alpha, 
                       input$nreps, 
                       input$start, 
                       input$end, 
                       input$by)
      
      # OUPUT: build table of power results
      output$resultsTable <-  DT::renderDataTable({ 
        res %>%
          datatable() %>%
          formatRound(columns=colnames(res)[2], digits=c(0,2))
      })
    })
  })
}
# Run the application 
shinyApp(ui = ui, server = server)

