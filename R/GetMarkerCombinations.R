#' @title Obtain all possible marker combinations
#'
#' @param og (FlowFrame): FlowFrame containing the expression matrix, channel names, and marker names.
#'
#' @return (dataFrame) Returns a dataFrame containing all possible marker combinations
#'
#' @seealso \code{\link{CompensAID}}
#'
#' @examples
#' # Import FCS file
#' file <- flowCore::read.FCS(system.file("extdata", "68983.fcs", package = "CompensAID"))
#'
#' # Get combinations
#' mc <- GetMarkerCombinations(og = file)
#'
#' @export

GetMarkerCombinations <- function(og) {


  # Input validation -----------------------------------------------------------
  checkmate::assert(methods::is(og, "flowFrame"), "Object is not a flowFrame.")


  # Get all possible marker combinations ---------------------------------------
  m <- expand.grid("primary.marker" = flowCore::markernames(og),
                   "secondary.marker" = flowCore::markernames(og)) %>%
    dplyr::filter(!primary.marker == secondary.marker) %>%
    lapply(., function(x) utils::type.convert(as.character(x), as.is = TRUE)) %>%
    as.data.frame()


  # Generate output ------------------------------------------------------------
  return(m)
}
