
#' @name latex.figure
#' @title Generate latex code for including a figure.
#' @param figfile Path to figure file.
#' @param reportfile Path to report file.
#' @param caption This character string will be included as the figure caption.
#' @return Nothing.
latex.figure <- function(figfile, reportfile, caption=''){
    figstr <- paste0('\\begin{figure}\n\\centering\n \\includegraphics[trim=0cm 0cm 0cm 0cm,clip,width=1\\textwidth]{', figfile, '}\n \\caption{', caption, '}\n \\end{figure}\n')
    cat(figstr, file=reportfile, append=TRUE)
}


#' @name make.report
#' @title Creates a pdf file containing the summary output and result plots
#' @details This function probably requires that you are running linux and that you have latex functions installed (pdflatex).
#' @param rep A valid result from fit.spict with OSA residuals.
#' @param reporttitle This character string will be printed as the first line of the report.
#' @param reportfile The generated tex code will be stored in this file.
#' @param summaryoutfile Output of the summary will be stored in this file as plain text.
#' @param keep.figurefiles If TRUE generated figure files will not be cleaned up.
#' @param keep.txtfiles If TRUE generated txt files will not be cleaned up.
#' @return Nothing.
#' @export
make.report <- function(rep, reporttitle='', reportfile='report.tex', summaryoutfile='summaryout.txt', keep.figurefiles=FALSE, keep.txtfiles=FALSE){
    latexstart <- '\\documentclass[12pt]{article}\n\\usepackage{graphicx}\n\\usepackage{verbatim}\n\\begin{document}\n'
    latexend <- '\\end{document}\n'
    cat(reporttitle, file=summaryoutfile)
    summaryout <- capture.output(summary(rep), file=summaryoutfile, append=TRUE)

    # -- Write tex file -- #
    cat(latexstart, file=reportfile)
    # Summary
    cat('\\footnotesize\n', file=reportfile, append=TRUE)
    cat(paste0('\\verbatiminput{', summaryoutfile, '}\n'), file=reportfile, append=TRUE)
    #cat('\\begin{verbatim}\n', file=reportfile, append=TRUE)
    #cat(reporttitle, file=reportfile, append=TRUE)
    #cat(summaryout, sep='\n', file=reportfile, append=TRUE)
    #cat('\\end{verbatim}\n', file=reportfile, append=TRUE)

    # Results plot
    figfile1 <- 'res.pdf'
    pdf(figfile1, width=9, height=10)
    plot(rep)
    dev.off()
    latex.figure(figfile1, reportfile, caption='Results.')

    # Retrospective analysis plot
    if('retro' %in% names(rep)){
        figfile1b <- 'retro.pdf'
        pdf(figfile1b)
        plotspict.retro(rep)
        dev.off()
        latex.figure(figfile1b, reportfile, caption='Retrospective analysis.')
    }
    
    # Diagnostic plot
    figfile2 <- 'diag.pdf'
    pdf(figfile2, width=9, height=10)
    plotspict.diagnostic(rep)
    dev.off()
    latex.figure(figfile2, reportfile, caption='Model diagnostics.')
    
    # Data plot
    figfile3 <- 'data.pdf'
    pdf(figfile3, width=7, height=6+rep$inp$nindex)
    par(mfrow=c(1+rep$inp$nindex,1))
    plot(rep$inp$timeC, rep$inp$obsC, xlab='Time', ylab='Catch', typ='b')
    for(i in 1:rep$inp$nindex) plot(rep$inp$timeI[[i]], rep$inp$obsI[[i]], xlab='Time', ylab=paste('Index', i), typ='b')
    dev.off()
    latex.figure(figfile3, reportfile, caption='Data.')
    
    cat(latexend, file=reportfile, append=TRUE)

    # -- Compile tex file -- #
    #latexcompile <- system(paste('pdflatex -output-directory=../res/', reportfile), intern=TRUE)
    latexcompile <- system(paste('pdflatex', reportfile), intern=TRUE)

    # -- Remove temporary files -- #
    #file.remove(paste0('../res/', substr(reportfile, 1, nchar(reportfile)-4), '.log'))
    #file.remove(paste0('../res/', substr(reportfile, 1, nchar(reportfile)-4), '.aux'))
    file.remove(paste0(substr(reportfile, 1, nchar(reportfile)-4), '.log'))
    file.remove(paste0(substr(reportfile, 1, nchar(reportfile)-4), '.aux'))
    if(!keep.txtfiles) file.remove(summaryoutfile)
    if(!keep.figurefiles){
        file.remove(figfile1)
        file.remove(figfile1b)
        file.remove(figfile2)
        file.remove(figfile3)
    }
}
