#' @title Update Secondary Stain Index
#'
#' @param si.input (dataFrame): DataFrame containing the SSI info.
#' @param rv.input (numerical): Numerical value for the number of segments.
#' @param primary (character): Name of the primary marker
#' @param secondary (character): Name of the secondary marker
#' @param population (list): List with the population matrices
#'
#' @return (numerical) Returns the Secondary Stain Index score.
#'
#' @seealso \code{\link{CompensAID}}, \code{\link{EmptyMatrixInfo}}, \code{\link{DensityGating}}
#'
#' @examples
#' # Import FCS file
#' file <- flowCore::read.FCS(system.file("extdata", "68983.fcs", package = "compensAID"))
#'
#' # Parameter for the number of segments
#' range.value <- 4
#'
#' # Density-based cut-off detection
#' center.plot <- 2
#' co <- DensityGating(og = file,
#'                     cp.value = center.plot)
#'
#' # All possible marker combinations
#' mc <- GetMarkerCombinations(og = file)
#'
#' # Parameter for the distance between the primary positive and negative population
#' separation.distance <- 0.25
#'
#' # Obtain SSI info dataFrame
#' si <- EmptyMatrixInfo(og = file,
#'                       rv.input = range.value,
#'                       mc.input = mc,
#'                       co.input = co,
#'                       sd.input = separation.distance)
#'
#' # Channel names
#' primary.marker <- "CD19"
#' secondary.marker <- "CD3"
#'
#' # Get populations
#' pop <- GetPopulations(og = file,
#'                       primary = primary.marker,
#'                       secondary = secondary.marker,
#'                       co.input = co,
#'                       sd.input = separation.distance)
#'
#' #' # Channel names
#' cp <- "PE-Cy5-A"
#' cs <- "PE-Cy7-A"
#'
#' # Update matrix information
#' range <- (max(pop$primary.positive[, cp]) - min(pop$primary.positive[, cp]))/range.value
#' si <- UpdateMatrixInfo(si.input = si,
#'                        rv.input = range.value,
#'                        range.input = range,
#'                        primary = primary.marker,
#'                        secondary = secondary.marker,
#'                        output = "PASS",
#'                        population = pop)
#'
#' # Calculate SSI for the last segment
#' si <- UpdateSSI(si.input = si,
#'                 rv.input = range.value,
#'                 primary = primary.marker,
#'                 secondary = secondary.marker,
#'                 population = pop)
#'
#' @export

UpdateSSI <- function(si.input, rv.input, primary, secondary, population) {


  # Input validation -----------------------------------------------------------
  checkmate::checkDataFrame(si.input)
  checkmate::checkNumeric(rv.input)
  checkmate::checkCharacter(primary)
  checkmate::checkCharacter(secondary)
  checkmate::checkList(population)

  # Obtain information
  population.positive <- population$primary.positive
  primary.channel = si.input$primary.channel[si.input$primary.marker == primary & si.input$secondary.marker == secondary][1]
  secondary.channel = si.input$secondary.channel[si.input$primary.marker == primary & si.input$secondary.marker == secondary][1]

  # Update SSI per segment -----------------------------------------------------
  for (s in 1:rv.input) {

    # Skip iteration if segment is ok
    if (si.input$message[si.input$primary.channel == primary.channel &
                         si.input$secondary.channel == secondary.channel &
                         si.input$segment == s] == "PASS") { next }

    # Skip iteration if segment is merged
    if (is.na(si.input$mfi.neg[si.input$primary.channel == primary.channel &
                               si.input$secondary.channel == secondary.channel &
                               si.input$segment == s])) { next }

    # Get segment minimum
    new.min <- si.input$segment.min[si.input$primary.channel == primary.channel &
                                    si.input$secondary.channel == secondary.channel &
                                    si.input$segment == s]

    # Get segment maximum
    new.max <- si.input$segment.max[si.input$primary.channel == primary.channel &
                                    si.input$secondary.channel == secondary.channel &
                                    si.input$segment == s]


    # Obtain new segments ------------------------------------------------------
    positive.sub <- population.positive[population.positive[, primary.channel] >= new.min & population.positive[, primary.channel] <= new.max, , drop = FALSE]

    # Get mean fluorescence intensity value
    si.input$mfi.pos[si.input$primary.channel == primary.channel &
                     si.input$secondary.channel == secondary.channel &
                     si.input$segment == s] <- mean(positive.sub[, secondary.channel])

    # Get secondary stain index
    si.input$ssi[si.input$primary.channel == primary.channel &
                 si.input$secondary.channel == secondary.channel &
                 si.input$segment == s] <- CalculateSSI(si.input = si.input,
                                                        primary.channel = primary.channel,
                                                        secondary.channel = secondary.channel,
                                                        segment = s)
    }


  # Generate output ------------------------------------------------------------
  return(si.input)
}
