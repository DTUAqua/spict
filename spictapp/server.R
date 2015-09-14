library(shiny)
library(spict)
data(pol)
polalb <- check.inp(pol$albacore)
anms <- names(polalb)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

    # Upload file
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    dataInput <- reactive({
        inFile <- input$file
        if(is.null(inFile)){
            return(NULL)
        } else {
            sep <- input$inpsep
            if(sep=='blank') sep <- ' '
            return(read.csv(inFile$datapath, header=input$inphead, sep=sep, skip=input$inpskip))
        }
    })

    # Input textarea
    make.inpout <- reactive({
        a <- strsplit(input$input_text, '\n')[[1]]
        inpp <- 'asdkj dtler alskjoie\ndteuler\n1\npoireut'
        inpp <- ''
        inpp <- 'hej'
        #a <- strsplit(inpp, '\n')[[1]]
        b <- NULL
        nanms <- length(anms)
        inpout <- list()
        cls <- character(nanms)
        if(length(a)>0){
            for(i in 1:nanms){
                cls[i] <- class(polalb[[anms[i]]])
                b <- grep(paste0(anms[i], '$'), a)
                if(length(b)==1){
                    inpval <- a[b+1]
                    class(inpval) <- cls[i]
                    inpout[[anms[i]]] <- inpval
                }
            }
        }
        return(inpout)
    })

    # Make inp list
    make.inp <- reactive({
        inp <- make.inpout()
        if(input$radio == 'own'){
            dat <- dataInput()
            if(is.null(dat)){
                return(NULL)
            } else {
                #inp <- list()
                timecol <- as.numeric(input$timecol)
                if(timecol>0){
                    inp$timeC <- dat[, timecol]
                    inp$timeI <- dat[, timecol]
                }
                inp$obsC <- dat[, as.numeric(input$catchcol)]
                inp$obsI <- dat[, as.numeric(input$indexcol)]
            }
        }
        if(input$radio == 'demo'){
            nms <- names(pol$albacore)
            for(i in 1:length(nms)) inp[[nms[i]]] <- pol$albacore[[nms[i]]]
            #inp <- pol$albacore
        }
        return(inp)
    })

    # Translate input inp list to verbatim output
    output$output_text <- renderPrint({
        make.inpout()
    })

    # Translate current inp list to verbatim output
    output$inp <- renderPrint({
        make.inp()
    })

    # Translate input file to a table
    output$contents <- renderTable({
        dataInput()
    })

    # Fit spict
    spict <- reactive({
        inp <- make.inp()
        if(!is.null(inp)){
            rep <- try(fit.spict(inp))
            return(rep)
        } else {
            return(NULL)
        }
    })

    # Capture summary of results
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
                    summary(rep)
                }
            }
        })
    })

    ### --- PLOTS --- ###
    # Data plot
    output$dataplot <- renderPlot({
        inp <- make.inp()
        if(!is.null(inp)){
            # dev.new(width=7, height=9) # This makes the plot a pop-up (may not work if R is not installed)
            plotspict.ci(inp)
        }
    })

    # Result plots
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
    output$catchplot <- renderPlot({
        # Take dependency on action button
        if(input$runspict == 0) return()
        isolate({
            rep <- spict()
            if(!is.null(rep)){
                plotspict.catch(rep, qlegend=FALSE)
            }
        })
    })
    output$bbplot <- renderPlot({
        # Take dependency on action button
        if(input$runspict == 0) return()
        isolate({
            rep <- spict()
            if(!is.null(rep)){
                plotspict.bbmsy(rep, qlegend=FALSE)
            }
        })
    })
    output$ffplot <- renderPlot({
        # Take dependency on action button
        if(input$runspict == 0) return()
        isolate({
            rep <- spict()
            if(!is.null(rep)){
                plotspict.ffmsy(rep)
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
