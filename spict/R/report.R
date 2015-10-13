# Stochastic surplus Production model in Continuous-Time (SPiCT)
#    Copyright (C) 2015  Martin Waever Pedersen, mawp@dtu.dk or wpsgodd@gmail.com
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


#' @name latex.figure
#' @title Generate latex code for including a figure.
#' @param figfile Path to figure file.
#' @param reportfile Path to report file.
#' @param caption This character string will be included as the figure caption.
#' @return Nothing.
latex.figure <- function(figfile, reportfile, caption=''){
    figstr <- paste0('\\begin{figure}[h!]\n\\centering\n \\includegraphics[trim=0cm 0cm 0cm 0cm,clip,width=1\\textwidth]{', figfile, '}\n \\caption{', caption, '}\n \\end{figure}\n\\newpage')
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

    # Results plot
    figfile1 <- 'res.pdf'
    pdf(figfile1, width=9, height=10)
    plot(rep)
    dev.off()
    latex.figure(figfile1, reportfile, caption='Results.')

    # Summary
    cat('\\scriptsize\n', file=reportfile, append=TRUE)
    cat(paste0('\\verbatiminput{', summaryoutfile, '}\n\\newpage'), file=reportfile, append=TRUE)
    #cat('\\begin{verbatim}\n', file=reportfile, append=TRUE)
    #cat(reporttitle, file=reportfile, append=TRUE)
    #cat(summaryout, sep='\n', file=reportfile, append=TRUE)
    #cat('\\end{verbatim}\n', file=reportfile, append=TRUE)

    # Retrospective analysis plot
    if('retro' %in% names(rep)){
        figfile1b <- 'retro.pdf'
        pdf(figfile1b)
        plotspict.retro(rep)
        dev.off()
        latex.figure(figfile1b, reportfile, caption='Retrospective analysis.')
    }
    
    # Diagnostic plot
    if('osar' %in% names(rep)){
        figfile2 <- 'diag.pdf'
        pdf(figfile2, width=7, height=9)
        plotspict.diagnostic(rep)
        dev.off()
        latex.figure(figfile2, reportfile, caption='Model diagnostics.')
    }    

    # Data plot
    figfile3 <- 'data.pdf'
    pdf(figfile3, width=9, height=10)
    plotspict.ci(rep$inp)
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
        if('retro' %in% names(rep)) file.remove(figfile1b)
        if('osar' %in% names(rep)) file.remove(figfile2)
        file.remove(figfile3)
    }
}
