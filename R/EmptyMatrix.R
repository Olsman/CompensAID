# Internal function - Obtain empty SSI matrix
.EmptyMatrix <- function(og) {
  
  
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
