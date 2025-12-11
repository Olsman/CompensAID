# Internal function - Obtain all possible marker combinations
.GetMarkerCombinations <- function(og) {
  
  
  # Input validation -----------------------------------------------------------
  checkmate::assert(methods::is(og, "flowFrame"), "Object is not a flowFrame.")
  
  
  # Get all possible marker combinations ---------------------------------------
  m <- expand.grid("primary.marker" = flowCore::markernames(og),
                   "secondary.marker" = flowCore::markernames(og)) |>
    dplyr::filter(primary.marker != secondary.marker) |>
    lapply(function(x) utils::type.convert(as.character(x), as.is = TRUE)) |>
    as.data.frame()
  
  
  # Generate output ------------------------------------------------------------
  return(m)
}
