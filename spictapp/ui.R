library(shiny)
options(shiny.trace=TRUE)

shinyUI(fluidPage(
    tags$style(type="text/css", "textarea {width:100%}"),
    # Application title
    titlePanel("SPiCT"),

    #sidebarLayout(
        #sidebarPanel(
    #tags$textarea(id = 'input_text', placeholder = 'Type here', rows = 8, ""),
    #tag$textarea(id="foo", rows=3, cols=40, "Default value"),
    column(8,
            fluidRow(column(8, h3('Inputs')), column(4, actionButton("runspict", 'Run SPiCT', col='green'))),
            tabsetPanel(type = "tabs",
                        tabPanel("Data", wellPanel(
                            radioButtons("radio", label = h5('Source'), choices = list("Own data" = 'own', "Demo data" = 'demo'), selected = 'own'),
                            fileInput("file", label = h5("File input (csv)"), accept=c('text/csv','text/comma-separated-values,text/plain','.csv')))),
                        tabPanel("Read options", wellPanel(
                                 fluidRow(
                column(3, checkboxInput("inphead", label = "Header?", value = TRUE)),
                column(2, selectInput('inpsep', 'Separator', choices=c(',', 'blank', ';'))),
                column(1, span('')),
                column(2, numericInput("inpskip", label = "Skip lines", value = 0))
                                     ),
                                 fluidRow(
                column(2, numericInput("timecol", label = "Time col", value = 1)),
                                     column(1, span('')),
                column(2, numericInput("catchcol", label = "Catch col", value = 2)),
                                     column(1, span('')),
                column(2, numericInput("indexcol", label = "Index col", value = 3)))
                                 )),

                        tabPanel("Set input options", wellPanel(
                                 fluidRow(
                column(4, tags$textarea(id = 'input_text', placeholder = 'Type here', rows = 8, "")),
                column(6, verbatimTextOutput("output_text"))
                #column(4, verbatimTextOutput("inp"))
                                     )
                                 ))
                        ),
            
            tabsetPanel(type = "tabs", 
                        tabPanel("Uploaded data", tableOutput("contents")),
                        tabPanel("Data plot", wellPanel( fluidRow( column(6,plotOutput("dataplot"))))),
                        tabPanel("Input list", verbatimTextOutput('inp')),
                        tabPanel("Hide input")
                        )
           ),

        # Show a plot of the generated distribution
        #mainPanel(
        column(9,
               h3('Ouputs'),
               tabsetPanel(type='tabs',
                           tabPanel('Plots',
                                    column(4, plotOutput("bplot")),
                                    column(4, plotOutput("fplot")),
                                    column(4, plotOutput("catchplot")),
                                    column(4, plotOutput("bbplot")),
                                    column(4, plotOutput("ffplot")),
                                    column(4, plotOutput("fbplot")),
                                    #column(4, plotOutput("osarplot")),
                                    column(4, plotOutput("prodplot")),
                                    column(4, plotOutput("tcplot"))
                                    ),
                           tabPanel('Summary', verbatimTextOutput('summary'))
               ))
        #)
    ))
