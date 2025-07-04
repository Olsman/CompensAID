#' @title Combine segments that do not meet the requirements
#'
#' @param si.input (dataFrame): dataFrame containing SSI info.
#' @param primary (character): Name of the primary marker.
#' @param secondary (character) Name of the secondary marker.
#' @param ev.input (numerical): Minimum required number of events.
#' @param rv.input (numerical): Number of segments.
#'
#' @return (dataFrame) Returns a dataframe with an update SSI information dataFrame.
#'
#' @seealso \code{\link{CompensAID}}, \code{\link{DensityGating}}, \code{\link{EmptyMatrixInfo}}, \code{\link{UpdateMatrixInfo}}, \code{\link{GetPopulations}}
#'
#' @examples
#' # Import FCS file
#' file <- flowCore::read.FCS(system.file("extdata", "68983.fcs", package = "CompensAID"))
#'
#' # Marker names
#' primary.marker <- "CD19"
#' secondary.marker <- "CD3"
#'
#' # Channel names
#' cp <- "PE-Cy5-A"
#' cs <- "PE-Cy7-A"
#'
#' # Density-based cut-off detection
#' center.plot <- 2
#' co <- DensityGating(og = file,
#'                     cp.value = center.plot)
#'
#' # Get combinations
#' mc <- GetMarkerCombinations(og = file)
#'
#' # Parameter for the distance between the primary positive and negative population
#' separation.distance <- 0.25
#'
#' # Obtain empty SSI matrix information
#' range.value <- 4
#' si <- EmptyMatrixInfo(og = file,
#'                       rv.input = range.value,
#'                       mc.input = mc,
#'                       co.input = co,
#'                       sd.input = separation.distance)
#'
#' # Obtain populations
#' pop <- GetPopulations(og = file,
#'                       primary = primary.marker,
#'                       secondary = secondary.marker,
#'                       co.input = co,
#'                       sd.input = separation.distance)
#'
#' # Select primary positive population
#' population.positive <- pop$primary.positive
#'
#' # Get the events that fall within the fourth segment
#' segment.value <- 4
#' range.value <- (max(population.positive[, cp]) - min(population.positive[, cp]))/segment.value
#'
#' # Adjust the SSI info if there are >20 events in the (segmented) positive population
#' # (output = "PASS"); meaning more events present than threshold
#' si <- UpdateMatrixInfo(si.input = si,
#'                        rv.input = segment.value,
#'                        range.input = range.value,
#'                        primary = primary.marker,
#'                        secondary = secondary.marker,
#'                        output = "PASS",
#'                        population = pop)
#'
#' # Update segments
#' events.value <- 20
#' si <- UpdateSegments(si.input = si,
#'                      primary = primary.marker,
#'                      secondary = secondary.marker,
#'                      ev.input = events.value,
#'                      rv.input = segment.value)
#'
#' @export

UpdateSegments <- function(si.input, primary, secondary, ev.input, rv.input) {


  # Input validation -----------------------------------------------------------
  checkmate::checkDataFrame(si.input)
  checkmate::checkCharacter(primary)
  checkmate::checkCharacter(secondary)
  checkmate::checkNumeric(ev.input)
  checkmate::checkNumeric(rv.input)


  # Merge segments -------------------------------------------------------------
  if (si.input$message[si.input$primary.marker == primary &
                       si.input$secondary.marker == secondary][1] != "No positive/negative population") {

    # Identify which segments need to be merged
    merge <- MergeSegments(si.input$event.count[si.input$primary.marker == primary &
                                                si.input$secondary.marker == secondary],
                           ev.input)
    si.input$mergeGroup <- NA

    # Assign merge IDs
    for (i in seq_along(merge)) {

      si.input$mergeGroup[si.input$primary.marker == primary &
                          si.input$secondary.marker == secondary][merge[[i]]] <- i
    }

    # Adjust segment information
    si.input <- AdjustSegments(si.input,
                               mp = primary,
                               ms = secondary)
  }


  # Generate output ------------------------------------------------------------
  return(si.input)
}
