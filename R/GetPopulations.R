#' @title Obtain (segmented) populations
#'
#' @param og (FlowFrame): FlowFrame containing the expression matrix, channel names, and marker names.
#' @param primary (character): Name of the primary marker.
#' @param secondary (character): Name of the secondary marker.
#' @param co.input (dataFrame) dataFrame with the output of the density-based cut-off detection.
#' @param sd.input (numerical): Numerical value determining the distance between the primary negative and positive population.
#' 
#' @return (list) Returns a list with the negative population of the primary and secondary marker, and positive population of the primary marker. 
#' 
#' @seealso \code{\link{CompensAID}, \code{\link{DensityGating}}
#' 
#' @examples 
#' # Import FCS file
#' file <- flowCore::read.FCS("path/to/exampleFCS.fcs")
#' 
#' # Marker names
#' primary.marker <- "CD19"
#' secondary.marker <- "IgL"
#' 
#' # Density-based cut-off detection
#' center.plot <- 2
#' co <- DensityGating(og = file, 
#'                     cp.value = center.plot)
#' 
#' # Parameter for the distance between the primary positive and negative population
#' separation.distance <- 0.25
#' 
#' # Obtain populations
#' pop <- GetPopulations(og = file, 
#'                       primary = primary.marker, 
#'                       secondary = secondary.marker, 
#'                       co.input = co, 
#'                       sd.input = separation.distance)
#' @export

GetPopulations <- function(og, primary, secondary, co.input, sd.input) {
  
  # Input validation -----------------------------------------------------------
  checkmate::assert(is(og, "flowFrame"), "Object is not a flowFrame.")
  checkmate::checkCharacter(primary)
  checkmate::checkCharacter(secondary)
  checkmate::checkDataFrame(co.input)
  checkmate::checkNumeric(sd.input)
  
  # Empty list 
  populations <- list()
  
  # Retrieve channels names and filtered expression matrix
  cp <- names(flowCore::markernames(og))[flowCore::markernames(og) == primary]
  cs <- names(flowCore::markernames(og))[flowCore::markernames(og) == secondary]
  sub <- flowCore::exprs(og)[, c(cp, cs)]
  
  # Retrieve gating
  co.s <- co.input$opt[co.input$channel == cs]
  co.p.n <- co.input$opt[co.input$channel == cp] - sd.input
  co.p.p <- co.input$opt[co.input$channel == cp] + sd.input
  
  
  # Gate populations -----------------------------------------------------------
  populations[["secondary.negative"]] <- sub[sub[, cs] <= co.s, , drop = FALSE]
  sub.s <- populations[["secondary.negative"]]
  populations[["primary.negative"]] <- sub.s[sub.s[, cp] <= co.p.n, , drop = FALSE]
  populations[["primary.positive"]] <- sub.s[sub.s[, cp] >= co.p.p, , drop = FALSE]
  
  
  # Generate output ------------------------------------------------------------
  return(populations)
}