# Stochastic surplus Production model in Continuous-Time (SPiCT)
#    Copyright (C) 2015-2016  Martin W. Pedersen, mawp@dtu.dk, wpsgodd@gmail.com
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


#' Generate latex code for including a figure.
#' 
#' @param figfile Path to figure file.
#' @param reportfile Path to report file.
#' @param caption This character string will be included as the figure caption.
#' @return Nothing.
latex.figure <- function(figfile, reportfile, caption=''){
    figstr <- paste0('\\begin{figure}[h!]\n\\centering\n \\includegraphics[trim=0cm 0cm 0cm 0cm,clip,width=1\\textwidth]{', figfile, '}\n \\caption{', caption, '}\n \\end{figure}\n\\newpage')
    cat(figstr, file=reportfile, append=TRUE)
}


#' Creates a pdf file containing the summary output and result plots
#' 
#' This function probably works probably only in linux with a LaTeX distribution installed (pdflatex).
#' @param rep A valid result from fit.spict with OSA residuals.
#' @param reporttitle This character string will be printed as the first line of the report.
#' @param reportfile character, filename (ending in '.tex') of the report LaTeX file.
#' @param summaryoutfile character, filename of the summary output file.     
#' @param outdir character, directory where all output files are going to be saved.
#' @param keep.figurefiles If TRUE generated figure files will not be cleaned up.
#' @param keep.txtfiles If TRUE generated txt files will not be cleaned up.
#' @param keep.texfiles If TRUE generated tex file will not be cleaned up.
#' @return Nothing.
#' @export
make.report <- function(rep, reporttitle='', reportfile='report.tex', summaryoutfile='summaryout.txt', outdir = ".", keep.figurefiles=FALSE, keep.txtfiles=FALSE, keep.texfiles=FALSE){
    if(basename(reportfile) == "") reportfile <- 'report.tex'
    if(file.access(outdir, mode=1) != 0)  reportfile <- file.path(getwd(), reportfile)
    
    latexstart <- '\\documentclass[12pt]{article}\n\\usepackage{graphicx}\n\\usepackage{verbatim}\n\\begin{document}\n'
    latexend <- '\\end{document}\n'

    # -- Write tex file -- #
    cat(latexstart, file=reportfile)

    # Results plot
    figfile1 <- file.path(outdir,'res.pdf')
    pdf(figfile1, width=8, height=9)
    plot(rep)
    dev.off()
    latex.figure(figfile1, reportfile, caption='Results.')

    # Summary
    summaryoutfile <- file.path(outdir, summaryoutfile)    
    cat(reporttitle, '- report compiled', as.character(Sys.time()), '\n', file=summaryoutfile)
    summaryout <- capture.output(summary(rep), file=summaryoutfile, append=TRUE)
    cat('\\scriptsize\n', file=reportfile, append=TRUE)
    cat(paste0('\\verbatiminput{', summaryoutfile, '}\n\\newpage'), file=reportfile, append=TRUE)

    # Management summary
    if ('man' %in% names(rep)){
        mansummaryoutfile <- file.path(outdir,'mansummaryout.txt')
        cat('Management results\n\n', file=mansummaryoutfile)
        mansummaryout <- capture.output(mansummary(rep), file=mansummaryoutfile, append=TRUE)
        cat('\\scriptsize\n', file=reportfile, append=TRUE)
        cat(paste0('\\verbatiminput{', mansummaryoutfile, '}\n\\newpage'), file=reportfile, append=TRUE)
    }
    
    # Retrospective analysis plot
    if ('retro' %in% names(rep)){
        figfile1b <- file.path(outdir,'retro.pdf')
        pdf(figfile1b)
        plotspict.retro(rep)
        dev.off()
        latex.figure(figfile1b, reportfile, caption='Retrospective analysis.')
    }
    
    # Diagnostic plot
    if ('osar' %in% names(rep)){
        figfile2 <- file.path(outdir,'diag.pdf')
        pdf(figfile2, width=7, height=9)
        plotspict.diagnostic(rep)
        dev.off()
        latex.figure(figfile2, reportfile, caption='Model diagnostics.')
    }    

    # Data plot
    figfile3 <- file.path(outdir,'data.pdf')
    pdf(figfile3, width=7, height=9)
    plotspict.data(rep$inp)
    dev.off()
    latex.figure(figfile3, reportfile, caption='Data.')

    cat(latexend, file=reportfile, append=TRUE)

    # -- Compile tex file -- 
    #latexcompile <- system(paste('pdflatex -output-directory=../res/', reportfile), intern=TRUE)
    latexcompile <- system(paste(paste0('pdflatex -output-directory=',file.path(outdir)), file.path(reportfile)), intern=TRUE)

    # -- Remove temporary files -- #
    #file.remove(paste0('../res/', substr(reportfile, 1, nchar(reportfile)-4), '.log'))
    #file.remove(paste0('../res/', substr(reportfile, 1, nchar(reportfile)-4), '.aux'))
    file.remove(paste0(substr(reportfile, 1, nchar(reportfile)-4), '.log'))
    file.remove(paste0(substr(reportfile, 1, nchar(reportfile)-4), '.aux'))
    if (!keep.texfiles){
        file.remove(reportfile)
    }
    if (!keep.txtfiles){
        file.remove(summaryoutfile)
        if ('man' %in% names(rep)){
            file.remove(mansummaryoutfile)
        }
    }
    if (!keep.figurefiles){
        file.remove(figfile1)
        if ('retro' %in% names(rep)){
            file.remove(figfile1b)
        }
        if ('osar' %in% names(rep)){
            file.remove(figfile2)
        }
        file.remove(figfile3)
    }
}
