#' @title Perform density-based cut-off detection
#'
#' @param og (FlowFrame): FlowFrame containing the expression matrix, channel names, and marker names.
#' @param cp.value (numerical): Numerical value determining the preliminary center.
#'
#' @return (dataFrame) Returns a dataFrame with the output of the density-based cut-off detection
#'
#' @seealso \code{\link{CompensAID}, \link{GetClosestCenter}}
#'
#' @examples
#' # Import FCS file
#' file <- flowCore::read.FCS(system.file("extdata", "68983.fcs", package = "CompensAID"))
#'
#' # Perform density-based cut-off detection
#' center.plot <- 2
#' co <- DensityGating(og = file,
#'                     cp.value = center.plot)
#'
#' @export

DensityGating <- function(og, cp.value) {


  # Input validation -----------------------------------------------------------
  checkmate::assert(methods::is(og, "flowFrame"), "Object is not a flowFrame.")
  checkmate::checkNumeric(cp.value)


  # Perform density-based cut-off detection ------------------------------------
  d <- lapply(names(flowCore::markernames(og)), function(channel) {
    cutoffs <- flowDensity::deGate(og, channel = channel, verbose = FALSE, all.cuts = TRUE, upper = TRUE)
    tibble::tibble(channel = channel, cutoffs = list(cutoffs))})

  # Alter format
  d <- dplyr::bind_rows(d) %>%
    tidyr::unnest_wider(cutoffs, names_sep = "_") %>%
    as.data.frame()


  # Retrieve value closest to the visual determined cut-off --------------------
  if (ncol(d) != 2) {
    d$opt <- apply(d[, grepl("cutoffs", names(d))], 1, function(row) GetClosestCenter(row, closest = cp.value))
    } else {
      d$opt <- d[, grepl("cutoffs", names(d))]
      }


  # Obtain best gating ---------------------------------------------------------
  d <- d[, c("channel", "opt")]


  # Generate output ------------------------------------------------------------
  return(d)
}
