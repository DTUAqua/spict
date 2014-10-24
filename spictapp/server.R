library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

  output$explogr <- renderText({ 
      paste('r:', round(exp(input$logr), digits=3),
            'K:', round(exp(input$logK), digits=3),
            'q:', round(exp(input$logq), digits=3),
            'sdb:', round(exp(input$logsdb), digits=3),
            'sdf:', round(exp(input$logsdf), digits=3)
            )
  })

  dataInput <- reactive({
      inFile <- input$file
      if (is.null(inFile)) return(NULL)
      read.csv(inFile$datapath, header=TRUE, sep=',')
  })

  makeinp <- reactive({
      dat <- dataInput()
      inp <- list()
      inp$obsC <- dat[, 2]
      inp$obsI <- dat[, 3]
      inp$ini <- list()
      inp$ini$logr <- input$logr
      inp$ini$logK <- input$logK
      inp$ini$logq <- input$logq
      inp$ini$logsdb <- input$logsdb
      inp$ini$logsdf <- input$logsdf
      inp$ini$alpha <- input$alpha
      inp$ini$beta <- input$beta
      inp
  })

  
  spict <- reactive({
      inp <- makeinp()
      require(spict)
      fit.spict(inp)
  })

  #output$text2 <- renderText({ 
  #    dataInput()
  #})

  # input$file1 will be NULL initially. After the user selects
  # and uploads a file, it will be a data frame with 'name',
  # 'size', 'type', and 'datapath' columns. The 'datapath'
  # column will contain the local filenames where the data can
  # be found.
  output$contents <- renderTable({
      inFile <- input$file
      if (is.null(inFile)) return(NULL)
      read.csv(inFile$datapath, header=TRUE, sep=',')
  })

  
  output$dataplot <- renderPlot({
      dat <- dataInput()
      i <- as.numeric(input$ycol)+1
      if(!is.null(dat)) plot(dat[, 1], dat[, i], typ='l', ylab=names(dat)[i], xlab='Time')
  })

  output$inp <- renderPrint({
      makeinp()
  })
  
  output$resplot <- renderPlot({
      rep <- spict()
      plotspict(rep)
  })

})
