#' @title Calculate segment information
#'
#' @param population (matrix): Matrix of the population that is segmented.
#' @param primary.channel (character): Name of the primary channel.
#' @param secondary.channel (character): Name of the secondary channel.
#' @param segment (numerical) Segment number for which the information is obtained.
#' @param range.input (numerical): Numerical value defining the width of each segment.
#'
#' @return (dataFrame) Returns a dataframe with the minimum segment value, maximum segment value, events within the segment, and mean fluorescence intensity.
#'
#' @seealso \code{\link{CompensAID}}, \code{\link{DensityGating}}, \code{\link{UpdateMatrixInfo}}, \code{\link{GetPopulations}}
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
#' # Parameter for the distance between the primary positive and negative population
#' separation.distance <- 0.25
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
#' info <- GetSegment(population = population.positive,
#'                    primary.channel = cp,
#'                    secondary.channel = cs,
#'                    segment = segment.value,
#'                    range.input = range.value)
#' @export

GetSegment <- function(population, primary.channel, secondary.channel, segment, range.input) {


  # Input validation -----------------------------------------------------------
  checkmate::checkMatrix(population)
  checkmate::checkCharacter(primary.channel)
  checkmate::checkCharacter(secondary.channel)
  checkmate::checkNumeric(range.input)
  checkmate::checkNumeric(range.input)


  # Calculate segment information of the first segment -------------------------
  if (segment == 1) {

    # Filter segment
    min.div <- min(population[, primary.channel])
    max.div <- max(population[population[, primary.channel] <= (min.div + (segment*range.input)), , drop = FALSE][,primary.channel])
    div <- population[population[, primary.channel] >= min.div & population[, primary.channel] <= max.div, , drop = FALSE]

    # Determine the number of events and mean fluorescence intensity
    div.count <- nrow(div)
    div.mfi.pos <- mean(div[, secondary.channel])


    # Calculate segment information of the other segment(s) --------------------
    } else {

      # Filter segment
      min.div <- min(population[, primary.channel])
      max.div <- max(population[population[, primary.channel] <= (min.div + (segment*range.input)), , drop = FALSE][,primary.channel])
      min.div <- max(population[population[, primary.channel] <= (min.div + ((segment-1)*range.input)), , drop = FALSE][,primary.channel])
      div <- population[population[, primary.channel] >= min.div & population[, primary.channel] <= max.div, , drop = FALSE]

      # Determine the number of events and mean fluorescence intensity
      div.count <- nrow(div)
      div.mfi.pos <- mean(div[, secondary.channel])
      }


  # Combine segment information ------------------------------------------------
  segment.info <- data.frame("min" = min.div,
                             "max" = max.div,
                             "count" = div.count,
                             "mfi.pos" = div.mfi.pos)


  # Generate output ------------------------------------------------------------
  return(segment.info)
}
