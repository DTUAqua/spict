library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("SPiCT"),
    
    # Sidebar with a slider input for the number of bins
    sidebarLayout(
        sidebarPanel(
            h2('Inputs'),
            #actionButton("action", label = "Action"),
            textOutput("explogr"),
            selectInput('ycol', 'Y Variable', c(1, 2)),
            fileInput("file", label = h3("File input (csv)"), accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
            h3('Initial parameter values'),
            fluidRow(
                column(3, sliderInput("logr", "logr:", min = -3, max = 3, value = -1, step=1e-3, round=3)),
                column(3, sliderInput("logK", "logK:", min = log(1e2), max = log(1e7), value = log(1e3), step=1e-3, round=3)),
                column(3, sliderInput("logq", "logq:", min = log(1e-5), max = log(1), value = log(1e-3), step=1e-3, round=3))
                ),
            fluidRow(
                column(3, sliderInput("logsdb", "logsdb:", min = log(0.05), max = log(5), value = log(1), step=1e-3, round=3)),
                column(3, sliderInput("logsdf", "logsdf:", min = log(0.05), max = log(5), value = log(1), step=1e-3, round=3))
                ),
            fluidRow(
                column(3, sliderInput("alpha", "alpha:", min = 0.1, max = 5, value = 1, step=1e-3, round=3)),
                column(3, sliderInput("beta", "beta:", min = 0.1, max = 5, value = 1, step=1e-3, round=3))
                ),
                        tabsetPanel(type = "tabs", 
                                    tabPanel("Plot", plotOutput("dataplot")),
                                    tabPanel("Table", tableOutput("contents")),
                                    tabPanel("Input parameters", verbatimTextOutput('inp'))
                        )
            ),

        # Show a plot of the generated distribution
        mainPanel(h2('Ouputs'),
                  plotOutput("resplot", height='800px')
                  )
        )
    ))
