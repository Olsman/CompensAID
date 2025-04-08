#' @title Identify value closest to visual estimation
#'
#' @param row (numerical): Numerical values of all density-based cut-off detection values.
#' @param closest (numerical): Numerical value determining the preliminary center.
#' 
#' @return (numerical) Returns the value closest to the visual center of the plot.
#' 
#' @seealso \code{\link{CompensAID}, \link{DensityGating}}
#' 
#' @examples 
#' # Import FCS file
#' file <- flowCore::read.FCS("path/to/exampleFCS.fcs")
#' primary.channel <- "PE-Cy7-A"
#' 
#' # Density-based cut-off detection of a single channel.
#' density.detection <- flowDensity::deGate(file, 
#'                                          channel = primary.channel, 
#'                                          verbose = FALSE, 
#'                                          all.cuts = TRUE, 
#'                                          upper = TRUE)
#'                                
#' # Identify cut-off closest to visual estimation
#' center.plot <- 2
#' cut.off <- GetClosestCenter(density.detection, 
#'                             closest = center.plot)                                
#' 
#' @export

GetClosestCenter <- function(row, closest) {
  
  
  # Input validation -----------------------------------------------------------
  checkmate::checkNumeric(row)
  checkmate::checkNumeric(closest)
  
  
  # Identify the index of the estimations closest to the visual estimation -----
  distances <- abs(row - closest)
  closest <- which.min(distances)
  
  
  # Generate output ------------------------------------------------------------
  return(row[closest])
}