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


#' @useDynLib spict
.onUnload <- function (lib) {
  library.dynam.unload("spict", lib)
}

.onAttach <- function(lib, pkg) {
    #ver <- utils::packageVersion('spict')
    ver <- get.version('spict')
    checkTMBversion()
    packageStartupMessage(paste0('Welcome to ', ver))
}

## Code adapted from TMB package
checkTMBversion <- function() {
  file <- file.path(system.file(package = "spict"), "TMB-version")
  cur.TMB.version <- get.version("TMB")
  if (!file.exists(file)) {
    writeLines(cur.TMB.version, con = file)
  }
  spict.TMB.version <- readLines(file)
  if (!identical(spict.TMB.version, cur.TMB.version)) {
    warning("Package version inconsistency detected.\n",
            "spict was built with TMB version ", spict.TMB.version,
            "\n", "Current TMB version is ", cur.TMB.version,
            "\n", "Please re-install 'spict' from source using remotes::install_github('DTUAqua/spict/spict')")
  }
}
