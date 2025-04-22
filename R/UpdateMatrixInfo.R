#' @title Calculate segment information
#'
#' @param si.input (dataFrame): dataFrame containing SSI info.
#' @param range.input (numerical): Width of each segment.
#' @param primary (character): Name of the primary marker.
#' @param secondary (character) Name of the secondary marker.
#' @param output (character): Will be "PASS" or "No positive/negative population" depending on the required events.
#' @param rv.input (numerical): Number of segments.
#' @param segment (numerical): The specific segment for which the information is updated.
#' @param population (matrix): Matrix with the population for which the information is obtained.
#'
#' @return (dataFrame) Returns a dataframe with an update SSI information dataFrame.
#'
#' @seealso \code{\link{CompensAID}}, \code{\link{DensityGating}}, \code{\link{EmptyMatrixInfo}}, \code{\link{GetPopulations}}
#'
#' @examples
#' # Import FCS file
#' file <- flowCore::read.FCS(system.file("extdata", "68983.fcs", package = "compensAID"))
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
#' @export

UpdateMatrixInfo <- function(si.input, range.input = NULL, primary, secondary, output, rv.input, segment = NULL, population = NULL) {


  # Input validation -----------------------------------------------------------
  checkmate::checkDataFrame(si.input)
  checkmate::checkNumeric(range.input)
  checkmate::checkCharacter(primary)
  checkmate::checkCharacter(secondary)
  checkmate::checkCharacter(output)
  checkmate::checkNumeric(rv.input)
  checkmate::checkNumeric(segment)
  checkmate::checkList(population)

  # Channel names
  channel.primary = si.input$primary.channel[si.input$primary.marker == primary & si.input$secondary.marker == secondary][1]
  channel.secondary = si.input$secondary.channel[si.input$primary.marker == primary & si.input$secondary.marker == secondary][1]


  # Update SSI information if there is no positive/negative population ---------
  if (output == "No positive/negative population") {
    col.adjust <- c("segment.min", "segment.max", "event.count", "event.count.merge", "mfi.neg", "sd.neg", "mfi.pos", "ssi")
    si.input[si.input$primary.marker == primary & si.input$secondary.marker == secondary, col.adjust] <- NA
    si.input[si.input$primary.marker == primary & si.input$secondary.marker == secondary, "message"] <- "No positive/negative population"
  }


  # Update SSI information if there is a positive/negative population ----------
  if (output == "PASS") {

    # Information for the positive population
    si.input[si.input$primary.marker == primary & si.input$secondary.marker == secondary, "message"] <- "PASS"
    population.positive <- population$primary.positive


    # Obtain information per segment -------------------------------------------
    for (s in 1:rv.input) {

      # Minimum fluorescence value within the segment
      si.input$segment.min[si.input$primary.marker == primary &
                           si.input$secondary.marker == secondary &
                           si.input$segment == s] <- GetSegment(population = population.positive,
                                                                primary.channel = channel.primary,
                                                                secondary.channel = channel.secondary,
                                                                segment = s,
                                                                range.input = range.input)[,"min"]

      # Maximum fluorescence value within the segment
      si.input$segment.max[si.input$primary.marker == primary &
                           si.input$secondary.marker == secondary &
                           si.input$segment == s] <- GetSegment(population = population.positive,
                                                                primary.channel = channel.primary,
                                                                secondary.channel = channel.secondary,
                                                                segment = s,
                                                                range.input = range.input)[,"max"]

      # Number of events within the segment
      si.input$event.count[si.input$primary.marker == primary &
                           si.input$secondary.marker == secondary &
                           si.input$segment == s] <- GetSegment(population = population.positive,
                                                                primary.channel = channel.primary,
                                                                secondary.channel = channel.secondary,
                                                                segment = s,
                                                                range.input = range.input)[,"count"]

      # Mean fluorescence intensity of the segment
      si.input$mfi.pos[si.input$primary.marker == primary &
                       si.input$secondary.marker == secondary &
                       si.input$segment == s] <- GetSegment(population = population.positive,
                                                            primary.channel = channel.primary,
                                                            secondary.channel = channel.secondary,
                                                            segment = s,
                                                            range.input = range.input)[,"mfi.pos"]
    }


    # Obtain information for the negative population ---------------------------
    population.negative <- population$primary.negative

    # Mean fluorescence intensity of the negative population
    si.input$mfi.neg[si.input$primary.marker == primary &
                     si.input$secondary.marker == secondary] <- mean(population.negative[, channel.secondary])

    # Standard deviation of the negative population
    si.input$sd.neg[si.input$primary.marker == primary &
                    si.input$secondary.marker == secondary] <- stats::sd(population.negative[, channel.secondary])


    # Calculate the Secondary Stain Index --------------------------------------
    for (s in 1:rv.input) {

      # SSI per segment
      si.input$ssi[si.input$primary.marker == primary &
                   si.input$secondary.marker == secondary &
                   si.input$segment == s] <- CalculateSSI(si.input = si.input,
                                                          primary.channel = channel.primary,
                                                          secondary.channel = channel.secondary,
                                                          segment = s)
    }
  }


  # Generate output ------------------------------------------------------------
  return(si.input)
}
