library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("SPiCT"),

    #sidebarLayout(
        #sidebarPanel(
    column(4,
            fluidRow(column(8, h3('Inputs')), column(4, actionButton("runspict", 'Run SPiCT', col='green'))),
            tabsetPanel(type = "tabs",
                        tabPanel("Data", wellPanel(
                            radioButtons("radio", label = h5('Source'), choices = list("Own data" = 'own', "Demo data" = 'demo'), selected = 'own'),
                            fileInput("file", label = h5("File input (csv)"), accept=c('text/csv','text/comma-separated-values,text/plain','.csv')))),
                        tabPanel("Read options", wellPanel(
                                 fluidRow(
                column(4, checkboxInput("inphead", label = "Header?", value = TRUE)),
                column(3, selectInput('inpsep', 'Separator', choices=c(',', 'blank', ';'))),
                column(1, span('')),
                column(3, numericInput("inpskip", label = "Skip lines", value = 0))
                                     ),
                                 fluidRow(
                column(3, numericInput("timecol", label = "Time col", value = 1)),
                                     column(1, span('')),
                column(3, numericInput("catchcol", label = "Catch col", value = 2)),
                                     column(1, span('')),
                column(3, numericInput("indexcol", label = "Index col", value = 3)))
                                 )),
                        tabPanel("Initial parameters", wellPanel(
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
                                 ))
                        ),
            #actionButton("action", label = "Action"),
            #textOutput("text2"),
            
            tabsetPanel(type = "tabs", 
                        tabPanel("Table", tableOutput("contents")),
                        tabPanel("Data plot", selectInput('ycol', 'Plot', c('Catch', 'Index')), plotOutput("dataplot")),
                        tabPanel("Input list", verbatimTextOutput('inp'))
                        )
#           wellPanel(verbatimTextOutput('summary'))
            #),
           ),

        # Show a plot of the generated distribution
        #mainPanel(
        column(9,
               h3('Ouputs'),
               tabsetPanel(type='tabs',
                           tabPanel('Plots',
               column(4, plotOutput("bplot")),
               column(4, plotOutput("osarplot")),
               column(4, plotOutput("fbplot")),
               column(4, plotOutput("catchplot")),
               column(4, plotOutput("prodplot")),
               column(4, plotOutput("fplot"))),
                           tabPanel('Summary', verbatimTextOutput('summary'))
               ))
        #)
    ))
