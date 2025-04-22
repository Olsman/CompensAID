#' @title Obtain empty SSI matrix
#'
#' @param og (FlowFrame): FlowFrame containing the expression matrix, channel names, and marker names.
#'
#' @return (matrix) Returns an empty SSI matrix
#'
#' @seealso \code{\link{CompensAID}}
#'
#' @examples
#' # Import FCS file
#' file <- flowCore::read.FCS(system.file("extdata", "68983.fcs", package = "compensAID"))
#'
#' # Empty matrix
#' sm <- EmptyMatrix(og = file)
#'
#' @export

EmptyMatrix <- function(og) {


  # Input validation -----------------------------------------------------------
  checkmate::assert(methods::is(og, "flowFrame"), "Object is not a flowFrame.")


  # Generate empty matrix ------------------------------------------------------
  m <- matrix(data = NA,
              nrow = length(names(flowCore::markernames(og))),
              ncol = length(names(flowCore::markernames(og))),
              dimnames = list(names(flowCore::markernames(og)), names(flowCore::markernames(og))))


  # Generate output ------------------------------------------------------------
  return(m)
}
