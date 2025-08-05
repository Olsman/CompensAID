#' @title Detect if smaller peaks inclusion improved the density-based cut-off detection.
#'
#' @param old.limit (numerical) Density-based cut-off before inclusion of small peaks.
#' @param new.limit (numerical): Density-based cut-off after inclusion of small peaks.
#' @param center.plot (numerical): Numerical value determining the preliminary center.
#'
#' @return (numerical) Returns limit that is closes to the center.plot
#'
#' @seealso \code{\link{GetMarkerCombinations}}, \code{\link{DensityGating}}, \code{\link{EmptyMatrix}}, \code{\link{EmptyMatrixInfo}}, \code{\link{GetPopulations}}, \code{\link{WithinLimit}}
#'
#' @examples
#' # Import FCS file
#' file <- flowCore::read.FCS(system.file("extdata", "68983.fcs", package = "CompensAID"))
#'
#' # Obtain all possible marker combinations
#' mc <- GetMarkerCombinations(og = file)
#'
#' # Density-based cut-off detection
#' center.plot <- 2
#' co <- DensityGating(og = file,
#'                     cp.value = center.plot)
#'
#' # Empty SSI(-info) matrix/dataFrame
#' range.value <- 4
#' separation.distance <- 0.25
#' sm <- EmptyMatrix(og = file)
#' si <- EmptyMatrixInfo(og = file,
#'                       rv.input = range.value,
#'                       mc.input = mc,
#'                       co.input = co,
#'                       sd.input = separation.distance)
#'
#' # Gating positive and negative population
#' primary.marker <- "CD19"
#' secondary.marker <- "CD3"
#' pop <- GetPopulations(og = file,
#'                       primary = primary.marker,
#'                       secondary = secondary.marker,
#'                       co.input = co,
#'                       sd.input = separation.distance)
#'
#' primary.channel <- "PE-Cy5-A" # CD19
#' old.lt <- co$opt[co$channel == primary.channel]
#' 
#' # Adjust limits
#' pop.limit <- WithinLimit(population = pop,
#'                          og = file,
#'                          primary = primary.marker,
#'                          secondary = secondary.marker,
#'                          min = 10,
#'                          max = 90,
#'                          si.input = si,
#'                          sd.input = separation.distance,
#'                          co.input = co,
#'                          cp.value = center.plot)
#' 
#' # Visualize if limit was adjusted                         
#' co.adjusted <- pop.limit[["co.dat"]]                          
#' new.lt <- co.adjusted$opt[co$channel == primary.channel]
#'
#' # Test the chosen
#' limit <- GetClosestLimit(old.limit = old.lt,
#'                          new.limit = new.lt,
#'                          center.plot = center.plot)
#'
#' @export

GetClosestLimit <- function(old.limit, new.limit, center.plot) {
  
  
  # Input validation -----------------------------------------------------------
  checkmate::checkNumeric(old.limit, len = 1)
  checkmate::checkNumeric(new.limit, len = 1)
  checkmate::checkNumeric(center.plot, len = 1)
  
  # Identify distance from center.plot
  distance.old <- abs(old.limit - center.plot)
  distance.new <- abs(new.limit - center.plot)
  
  # Determine closest
  if (distance.old < distance.new) {
  best.limit <- old.limit
  } else {
  best.limit <- new.limit 
  }
  
  
  # Generate output ------------------------------------------------------------
  return(best.limit)
}
