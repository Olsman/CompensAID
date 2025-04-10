#' @title Calculate Secondary Stain Index
#'
#' @param si.input (dataFrame): DataFrame containing the SSI info.
#' @param primary.channel (character): Name of the primary channel.
#' @param secondary.channel (character): Name of the secondary channel.
#' @param segment (numerical): Numerical value determining for which segment the SSI should be calculated.
#'
#' @return (numerical) Returns the Secondary Stain Index score.
#'
#' @seealso \code{\link{CompensAID}}, \code{\link{EmptyMatrixInfo}}, \code{\link{DensityGating}}
#'
#' @examples
#' # Import FCS file
#' file <- flowCore::read.FCS("path/to/exampleFCS.fcs")
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
#'
#' # Parameter for the distance between the primary positve and negative population
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
#' cp <- "PE-Cy7-A"
#' cs <- "APC-H7-A"
#'
#' # Calculate SSI for the last segment
#' CalculateSSI(si.input = si,
#'              primary.channel = cp,
#'              secondary.channel = cs
#'              secondary,segment = range.value)
#'
#' @export

CalculateSSI <- function(si.input, primary.channel, secondary.channel, segment) {


  # Input validation -----------------------------------------------------------
  checkmate::checkDataFrame(si.input)
  checkmate::checkCharacter(primary.channel)
  checkmate::checkCharacter(secondary.channel)
  checkmate::checkNumeric(segment)


  # Calculate Secondary Stain Index --------------------------------------------
  ssi <- si.input[si.input$primary.channel == primary.channel &
                  si.input$secondary.channel == secondary.channel &
                  si.input$segment == segment,] %>%
    dplyr::mutate(ssi = round((mfi.pos - mfi.neg)/(2*sd.neg), digits = 2)) %>%
    dplyr::pull(ssi)


  # Generate output ------------------------------------------------------------
  return(ssi)
}
