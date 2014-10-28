library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("SPiCT"),
    
    # Sidebar with a slider input for the number of bins
    sidebarLayout(
        sidebarPanel(
            fluidRow(column(9, h3('Inputs')), column(3, actionButton("runspict", 'Run SPiCT'))),
            tabsetPanel(type = "tabs", 
                        tabPanel("File", fileInput("file", label = h3("File input (csv)"), accept=c('text/csv','text/comma-separated-values,text/plain','.csv'))),
                        tabPanel("Read options",
                                 fluidRow(
                column(3, checkboxInput("inphead", label = "Header?", value = TRUE)),
                column(3, numericInput("inpskip", label = "Skip lines", value = 0)),
                column(3, selectInput('inpsep', 'Separator', c(',', ' ', ';')))),
                                 fluidRow(
                column(3, numericInput("timecol", label = "Time column", value = 1)),
                column(3, numericInput("catchcol", label = "Catch column", value = 2)),
                column(3, numericInput("indexcol", label = "Index column", value = 3)))
                                 ),
                        tabPanel("Initial parameters",
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
                )
                                 )
                        ),
            #actionButton("action", label = "Action"),
            #textOutput("parvals"),
            
            tabsetPanel(type = "tabs", 
                        tabPanel("Plot", selectInput('ycol', 'Plot', c('Catch', 'Index')), plotOutput("dataplot")),
                        tabPanel("Table", tableOutput("contents")),
                        tabPanel("Input list", verbatimTextOutput('inp'))
                        )
            ),

        # Show a plot of the generated distribution
        mainPanel(h2('Ouputs'),
                  plotOutput("resplot", height='800px')
                  )
        )
    ))
