#' @title Obtain empty SSI dataFrame
#'
#' @param og (FlowFrame): FlowFrame containing the expression matrix, channel names, and marker names.
#' @param rv.input (numerical): Numerical value of the number of segments in the primary positive population.
#' @param mc.input (dataFrame): dataFrame containing all possible marker combinations.
#' @param co.input (dataFrame) dataFrame with the output of the density-based cut-off detection.
#' @param sd.input (numerical): Numerical value determining the distance between the primary negative and positive population.
#'
#' @return (dataFrame) Returns an empty SSI dataFrame
#'
#' @seealso \code{\link{CompensAID}}, \code{\link{GetMarkerCombinations}}, \code{\link{DensityGating}}
#'
#' @examples
#' # Import FCS file
#' file <- flowCore::read.FCS(system.file("extdata", "68983.fcs", package = "CompensAID"))
#'
#' # Parameter for the number of segments
#' range.value <- 4
#'
#' # All possible marker combinations
#' mc <- GetMarkerCombinations(og = file)
#'
#' # Density-based cut-off detection
#' center.plot <- 2
#' co <- DensityGating(og = file,
#'                     cp.value = center.plot)
#'
#' # Parameter for the distance between the primary positve and negative population
#' separation.distance <- 0.25
#'
#' # Empty SSI dataFrame
#' si <- EmptyMatrixInfo(og = file,
#'                       rv.input = range.value,
#'                       mc.input = mc,
#'                       co.input = co,
#'                       sd.input = separation.distance)
#'
#' @export

EmptyMatrixInfo <- function(og = ff,
                            rv.input = range.value,
                            mc.input = mc,
                            co.input = co,
                            sd.input = separation.distance) {


  # Input validation -----------------------------------------------------------
  checkmate::assert(methods::is(og, "flowFrame"), "Object is not a flowFrame.")
  checkmate::checkNumeric(rv.input)
  checkmate::checkDataFrame(mc.input)
  checkmate::checkDataFrame(co.input)
  checkmate::checkNumeric(sd.input)

  # Obtain sample name
  fn <- gsub(".*/", "", og@description[["FILENAME"]])


  # Obtain empty SSI dataFrame -------------------------------------------------
  mi <- mc.input %>%
    dplyr::mutate(primary.channel = names(flowCore::markernames(og))[match(primary.marker, flowCore::markernames(og))],
                  secondary.channel = names(flowCore::markernames(og))[base::match(secondary.marker, flowCore::markernames(og))]) %>%
    dplyr::mutate(pretty.primary = paste0(primary.channel, ": ", primary.marker),
                  pretty.secondary = paste0(secondary.channel, ": ", secondary.marker)) %>%
    dplyr::mutate(file = fn,
                  primary.cutoff.neg = (co.input$opt[match(primary.channel, co.input$channel)] - sd.input),
                  primary.cutoff.pos = (co.input$opt[match(primary.channel, co.input$channel)] + sd.input),
                  secondary.cutoff = co.input$opt[match(secondary.channel, co.input$channel)],
                  mfi.neg = NA,
                  mfi.pos = NA,
                  sd.neg = NA,
                  ssi = NA,
                  message = NA,
                  segment.min = NA,
                  segment.max = NA,
                  event.count = NA,
                  event.count.merge = NA) %>%
    dplyr::slice(rep(seq_len(nrow(.)), each = rv.input)) %>%
    dplyr::group_by(primary.marker, secondary.marker) %>%
    dplyr::mutate(segment = rep(seq_len(rv.input))) %>%
    dplyr::select(file, primary.marker, primary.channel, pretty.primary, secondary.marker, secondary.channel, pretty.secondary,
                  segment, secondary.cutoff, primary.cutoff.neg, primary.cutoff.pos, segment.min, segment.max, event.count, event.count.merge,
                  mfi.neg, sd.neg, mfi.pos, ssi, message)


  # Generate output ------------------------------------------------------------
  return(mi)
}
