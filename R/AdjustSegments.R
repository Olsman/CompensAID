#' @title Adjust (merged) segment information
#'
#' @param si.input (dataFrame): DataFrame containing the SSI info.
#' @param mp (character): Name of the primary marker.
#' @param ms (character): Name of the secondary marker.
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
#' # Adjust segment information
#' si <- AdjustSegments(si.input = si,
#'                      mp = primary.marker,
#'                      ms = secondary.marker)
#'
#' @export

AdjustSegments <- function(si.input, mp, ms) {


  # Input validation -----------------------------------------------------------
  checkmate::checkDataFrame(si.input)
  checkmate::checkCharacter(mp)
  checkmate::checkCharacter(ms)


  # Add temporary columns ------------------------------------------------------
  si.input$segment.min.merged <- NA
  si.input$segment.max.merged <- NA
  si.input$event.count.sum <- NA
  si.input$merged.segments <- NA


  # Merge and Adjust segments --------------------------------------------------
  for (grp in unique(si.input$mergeGroup)) {


    # Identify merge
    group.rows <- which(si.input$mergeGroup[si.input$primary.marker == mp &
                                            si.input$secondary.marker == ms] == grp)

    # Skip empty rows
    if (length(group.rows) == 0) { next }


    # Determine new segment ranges and event counts ----------------------------
    min.val <- min(si.input$segment.min[si.input$primary.marker == mp &
                                        si.input$secondary.marker == ms][group.rows])
    max.val <- max(si.input$segment.max[si.input$primary.marker == mp &
                                        si.input$secondary.marker == ms][group.rows])
    sum.val <- sum(si.input$event.count[si.input$primary.marker == mp &
                                        si.input$secondary.marker == ms][group.rows])
    segments.merged <- paste(si.input$segment[si.input$primary.marker == mp &
                                              si.input$secondary.marker == ms][group.rows], collapse = "+")


    # Input values -------------------------------------------------------------
    first.row <- group.rows[1]
    si.input$segment.min.merged[si.input$primary.marker == mp &
                                si.input$secondary.marker == ms][first.row] <- min.val
    si.input$segment.max.merged[si.input$primary.marker == mp &
                                si.input$secondary.marker == ms][first.row] <- max.val
    si.input$event.count.sum[si.input$primary.marker == mp &
                             si.input$secondary.marker == ms][first.row] <- sum.val
    si.input$merged.segments[si.input$primary.marker == mp &
                             si.input$secondary.marker == ms][first.row] <- segments.merged
    }


  # Clean dataFrame ------------------------------------------------------------
  si.input <- CleanSegments(si.input, mp, ms)


  # Generate output ------------------------------------------------------------
  return(si.input)
  }
