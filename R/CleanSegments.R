#' @title Clean adjust segment information
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
#' # Clean segment information
#' si <- CleanSegments(si.input = si,
#'                     mp = primary.marker,
#'                     ms = secondary.marker)
#'
#' @export

CleanSegments <- function(si.input, mp, ms) {


  # Input validation -----------------------------------------------------------
  checkmate::checkDataFrame(si.input)
  checkmate::checkCharacter(mp)
  checkmate::checkCharacter(ms)


  # Move updated segments ------------------------------------------------------
  si.input$segment.min[si.input$primary.marker == mp &
                       si.input$secondary.marker == ms] <- si.input$segment.min.merged[si.input$primary.marker == mp &
                                                                                       si.input$secondary.marker == ms]
  si.input$segment.max[si.input$primary.marker == mp &
                       si.input$secondary.marker == ms] <- si.input$segment.max.merged[si.input$primary.marker == mp &
                                                                                       si.input$secondary.marker == ms]
  si.input$event.count.merge[si.input$primary.marker == mp &
                             si.input$secondary.marker == ms] <- si.input$event.count.sum[si.input$primary.marker == mp &
                                                                                          si.input$secondary.marker == ms]


  # Adjust the number of events after update -----------------------------------
  # No adjustment needed
  if (all(si.input$event.count.merge[si.input$primary.marker == mp &
                                     si.input$secondary.marker == ms] == si.input$event.count[si.input$primary.marker == mp &
                                                                                              si.input$secondary.marker == ms])) {

    # Set merge to NA
    si.input[si.input$primary.marker == mp &
             si.input$secondary.marker == ms &
             is.na(si.input$merged.segments), c("event.count.merge")] <- NA

    } else {

    # Loop per segment
    for (m in si.input$segment[si.input$primary.marker == mp & si.input$secondary.marker == ms]) {

      # Next loop if the number of events remains the same
      if (identical(si.input$event.count.merge[si.input$primary.marker == mp &
                                               si.input$secondary.marker == ms &
                                               si.input$segment == m],
                    si.input$event.count[si.input$primary.marker == mp &
                                         si.input$secondary.marker == ms &
                                         si.input$segment == m])) {

      next

      } else {

      # Add which segments are merged
      si.input$message[si.input$primary.marker == mp &
                       si.input$secondary.marker == ms &
                       si.input$segment == m] <- paste0("Merged segments: ", si.input$merged.segments[si.input$primary.marker == mp &
                                                                                                      si.input$secondary.marker == ms &
                                                                                                      !is.na(si.input$merged.segments) &
                                                                                                      si.input$segment == m])
      }
    }
  }

  # Set merge to NA
  si.input[si.input$primary.marker == mp &
           si.input$secondary.marker == ms &
           is.na(si.input$merged.segments), c("message", "mfi.neg", "mfi.pos", "sd.neg", "ssi")] <- NA

  # Remove temporary columns
  si.input <- subset(si.input, select = -c(merge_group, segment.min.merged, segment.max.merged, event.count.sum, merged.segments))


  # Generate output ------------------------------------------------------------
  return(si.input)
}
