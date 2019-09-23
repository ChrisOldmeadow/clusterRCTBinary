#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
require(plyr)
require(shinyWidgets)
require(shinythemes)
require(ggplot2)
require(DT)
require(lme4)
require(lmerTest)
source("utils.R")

# Define UI for application 
ui <- fluidPage(
    
  theme = shinytheme("cerulean"),
  
    # Application title
    titlePanel("Power for cluster RCT, binary outcomes"),
    
    # Sidebar  
    sidebarLayout(
      sidebarPanel(
        # numeric input for seed 
        h4("Enter seed value"),  
        numericInput(inputId = "seed", 
                     label = NULL, 
                     value = 123),
       
        
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
    
        
        h4("Select the parameters you would like to vary over a range"),
        radioButtons(
          inputId = "vary_range", 
          label = NULL, 
          choices = c("Cluster size", "Effect size", "Intra-class correlation")),
        
        # show if want to vary sample size 
        # i.e. select multiple sample sizes and one effect size
        conditionalPanel(
          condition = "input.vary_range == 'Cluster size'",
          
          numericInput(inputId = "m_start",
                       label = "Avg cluster size (min):",
                       value = 5, min=1, max=1000),
          numericInput(inputId = "m_end",
                       label = "Avg cluster size (max):",
                       value = 15, min=1, max=1000),
          numericInput(inputId = "m_by",
                       label = "by increments:",
                       value = 1, min=1, max=100),
          h4("Effect size parameters"),
          numericInput(inputId = "p1",
                       label = "Outcome proportion in the intervention group:",
                       value = 0.11),
          numericInput(inputId = "p2",
                       label = "Outcome proportion in the control group:",
                       value = 0.03)
          
          ),
        
        conditionalPanel(
          condition = "input.vary_range == 'Effect size'",
          numericInput(inputId = "p2",
                       label = "Outcome proportion in the control group:",
                       value = 0.03),
         
          numericInput(inputId = "p1_start",
                       label = "intervention outcome prevelance (min):",
                       value = .11, min=0.000001, max=0.99999),
          numericInput(inputId = "p1_end",
                       label = "intervention outcome prevelance (max):",
                       value = .21,  min=0.000001, max=0.99999),
          numericInput(inputId = "p1_by",
                       label = "by increments:",
                       value = .02, min=.00001, max=0.99999),
          
          numericInput(inputId = "m",
                       label = "Avg cluster size:",
                       value = 5, min=1, max=1000)
          
        ),
        
        h4("Simulation parameters") ,
        numericInput(inputId = "alpha",
                     label = "Type 1 error rate",
                     value = .05, min=.001, max = 0.1),
          
        # Input: Numeric entry for number of simulations ----
        numericInput(inputId = "nreps",
                     label = "Number of simulations",
                     value = 100, min=1, max=1500),
      
        
        # Copy the line below to make a set of radio buttons
        radioButtons("radio", label = h3("Estimation method"),
                     choices = list("lme4::glmm" = 1, "MASS:glmmPQL" = 2), 
                     selected = 2),
        
        hr(),
        fluidRow(column(3, verbatimTextOutput("value"))),
        
        # run simulation
        h4("Run the Power Simulations"),
        actionButton(inputId = "submit", label = "Submit")
        
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        h3("Power Results"),
        DT::dataTableOutput(outputId = "resultsTable"),
        h3("Power plot"),
        plotOutput("powerPlot")
      )
    )
)


# Define server logic 
server <- function(input, output) {
  # SIMULATE 
  observeEvent(input$submit, {
    withProgress(message="Running Power Simulation", value=0, {
      set.seed(input$seed)
      if (input$vary_range == 'Cluster size'){
        res <- powersims_m(input$p1,input$p2,
                           input$k1,input$k2,
                           input$rho,
                           input$alpha, 
                           input$nreps, 
                           input$m_start, 
                           input$m_end, 
                           input$m_by,
                           input$radio)
      
      }
      if (input$vary_range == 'Effect size'){
        res <- powersims_p(input$p1_start,input$p1_end, input$p1_by,
                           input$p2,
                           input$k1,input$k2,
                           input$rho,
                           input$alpha, 
                           input$nreps, 
                           input$m,
                           input$radio)
       
      }
      
      # OUPUT: build table of power results
      output$resultsTable <-  DT::renderDataTable({ 
        res %>%
          datatable(rownames= FALSE,
                    extensions = c("Buttons", "Select"),
                    options = 
                      list(
                        select = TRUE,
                        dom = "Bfrtip",
                        buttons = list(
                          list(
                            extend = "copy",
                            text = 'Copy to clipboard',
                            exportOptions = list(modifier = list(selected = TRUE))
                          )
                        )
                      )) %>%
          formatRound(columns=colnames(res)[2], digits=c(0,2))
      })
      # OUPUT: build plot of power results 
      #h3(textOutput("Power Plot"))
      output$powerPlot <- renderPlot({
        if(input$vary_range == 'Cluster size'){
          x <- seq(from = input$m_start, to = input$m_end, by=input$m_by)
          xlab <- "Avg cluster Size"
          #title <- paste("Power Simulations for Effect =", input$es)
        }
        else{
          x <- seq(from = input$p1_start, to = input$p1_end, by=input$p1_by)
          xlab <- "Intervention prevalence"
          #title <- paste("Power Simulations for Sample Size =", input$ss)
        }
        ggplot(res, aes(x = x, y = Power)) +
          geom_point() + geom_line() +
          theme(text = element_text(size=20))+
          geom_hline(yintercept = 0.8,
                     colour = "#8b1c3f",
                     linetype = "dashed",
                     size = 2) +
          xlab(xlab)  +
          ylab("Power") 
        #  ggtitle(title)
      })
    })
  })
}
# Run the application 
shinyApp(ui = ui, server = server)

