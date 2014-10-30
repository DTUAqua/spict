library(shiny)
library(spict)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

  output$parvals <- renderText({ 
      paste('r:', round(exp(input$logr), digits=3),
            'K:', round(exp(input$logK), digits=3),
            'q:', round(exp(input$logq), digits=3),
            'sdb:', round(exp(input$logsdb), digits=3),
            'sdf:', round(exp(input$logsdf), digits=3)
            )
  })

  dataInput <- reactive({
      inFile <- input$file
      if (is.null(inFile)){
          return(NULL)
      } else {
          sep <- input$inpsep
          if(sep=='blank') sep <- ' '
          return(read.csv(inFile$datapath, header=input$inphead, sep=sep, skip=input$inpskip))
      }
  })

  makeinp <- reactive({
      if(input$radio == 'own'){
          dat <- dataInput()
          if(is.null(dat)){
              return(NULL)
          } else {
              inp <- list()
              timecol <- as.numeric(input$timecol)
              if(timecol>0) inp$timeC <- dat[, timecol]
              inp$obsC <- dat[, as.numeric(input$catchcol)]
              inp$obsI <- dat[, as.numeric(input$indexcol)]
              inp$ini <- list()
              inp$ini$logr <- input$logr
              inp$ini$logK <- input$logK
              inp$ini$logq <- input$logq
              inp$ini$logsdb <- input$logsdb
              inp$ini$logsdf <- input$logsdf
              inp$ini$alpha <- input$alpha
              inp$ini$beta <- input$beta
          }
      }
      if(input$radio == 'demo'){
          data(pol.albacore)
      }
      return(inp)
  })

  
  spict <- reactive({
      inp <- makeinp()
      if(!is.null(inp)){
          rep <- try(fit.spict(inp))
          return(rep)
      } else {
          return(NULL)
      }
  })

  output$text2 <- renderText({
      input$runspict
      isolate(input$file$datapath)
  #    dataInput()
  })

  #output$text2 <- renderText({
  #    input$radio
  #})

  # input$file1 will be NULL initially. After the user selects
  # and uploads a file, it will be a data frame with 'name',
  # 'size', 'type', and 'datapath' columns. The 'datapath'
  # column will contain the local filenames where the data can
  # be found.
  output$contents <- renderTable({
      #inFile <- input$file
      #if (is.null(inFile)) return(NULL)
      #read.csv(inFile$datapath, header=TRUE, sep=',')
      dataInput()
  })

  
  output$dataplot <- renderPlot({
      #dat <- dataInput()
      #i <- as.numeric(input$ycol)+1
      inp <- makeinp()
      #if(!is.null(dat)) plot(dat[, 1], dat[, i], typ='l', ylab=names(dat)[i], xlab='Time')
      if(!is.null(inp)){
          if('timeC' %in% names(inp)){
              time <- inp$timeC
          } else {
              time <- 1:length(inp$obsC)
          }
          if(input$ycol == 'Catch'){
              plot(time, inp$obsC, typ='l', xlab='Time', ylab='Catch')
              points(time, inp$obsC, pch=20, cex=0.7)
          }
          if(input$ycol == 'Index'){
              plot(time, inp$obsI, typ='l', xlab='Time', ylab='Index')
              points(time, inp$obsI, pch=20, cex=0.7)
          }
      }
  })

  output$inp <- renderPrint({
      makeinp()
  })

  output$summary <- renderPrint({
      # Take dependency on action button
      if(input$runspict == 0) return()
      isolate({
          #sink(file='deleteme.txt')
          msg.trap <- capture.output(rep <- spict())
          #sink()
          if(class(rep)=='try-error'){
              cat(paste('There was an Error fitting the model!', rep, '\n'))
          } else {
              if(!is.null(rep)){
                  summaryspict(rep)
              }
          }
      })
  })

  output$bplot <- renderPlot({
      # Take dependency on action button
      if(input$runspict == 0) return()

      # Isolate to ensure spict only runs if button is pressed.
      # For isolate to work it is important to fix the extends of the output such as setting the height and width of a plot.
      isolate({
          rep <- spict()
          if(!is.null(rep)){
              plotspict.biomass(rep)
          }
      })
  })
  output$fbplot <- renderPlot({
      # Take dependency on action button
      if(input$runspict == 0) return()
      isolate({
          rep <- spict()
          if(!is.null(rep)){
              plotspict.fb(rep)
          }
      })
  })
  output$fplot <- renderPlot({
      # Take dependency on action button
      if(input$runspict == 0) return()
      isolate({
          rep <- spict()
          if(!is.null(rep)){
              plotspict.f(rep)
          }
      })
  })
  output$osarplot <- renderPlot({
      # Take dependency on action button
      if(input$runspict == 0) return()
      isolate({
          rep <- spict()
          osar <- calc.osa.resid(rep)
          if(!is.null(osar)){
              plotspict.osar(osar)
          }
      })
  })
  output$catchplot <- renderPlot({
      # Take dependency on action button
      if(input$runspict == 0) return()
      isolate({
          rep <- spict()
          if(!is.null(rep)){
              plotspict.catch(rep)
          }
      })
  })
  output$prodplot <- renderPlot({
      # Take dependency on action button
      if(input$runspict == 0) return()
      isolate({
          rep <- spict()
          if(!is.null(rep)){
              plotspict.production(rep)
          }
      })
  })
  output$tcplot <- renderPlot({
      # Take dependency on action button
      if(input$runspict == 0) return()
      isolate({
          rep <- spict()
          if(!is.null(rep)){
              plotspict.tc(rep)
          }
      })
  })
})
