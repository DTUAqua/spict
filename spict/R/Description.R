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


#' Fits a continuous-time surplus production model to data
#' 
#' @aliases spict
#' @name spict
#' @docType package
#' @author Martin W. Pedersen \email{mawp@@dtu.dk}
#' @references
#' \url{https://github.com/mawp/spict/}
#' @keywords production model fisheries assessment
#' @seealso \code{\link{test.spict}}
#' @examples
#' \dontrun{
#' rep <- test.spict()
#' }
#' @importFrom grDevices col2rgb colorRamp dev.cur dev.new dev.off pdf rgb
#' @importFrom graphics abline arrows axis box contour grid legend lines matplot mtext par plot points polygon rect strheight strwidth text title
#' @importFrom methods is
#'@importFrom stats Box.test acf cov2cor dgamma dlnorm dnorm lm na.omit na.pass nlminb optim pchisq pnorm predict qchisq qnorm qqline qqnorm rnorm runif shapiro.test smooth.spline t.test update
#' @importFrom utils capture.output data head packageVersion read.table tail
NULL
