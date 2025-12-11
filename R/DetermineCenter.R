# Internal function - Estimate the center of the FCS file for robust cutoff estimation.
.DetermineCenter <- function(og) {

  
  # Input validation -----------------------------------------------------------
  checkmate::assert(methods::is(og, "flowFrame"), "Object is not a flowFrame.")
  
  
  # Obtain marker and assess distribution --------------------------------------
  useCols <- names(flowCore::markernames(og))
  center <- graphics::hist(og@exprs[, useCols],
                           breaks = seq(round(min(og@exprs[, useCols]))-1,
                                        round(max(og@exprs[, useCols]))+1, by = 1), plot = FALSE)
  
  # Calculate center
  counts <- data.frame(counts = center$counts,
                       breaks = ceiling(center[["mids"]]))
  threshold <- 0.05 * max(counts)
  max.bin <- max(counts$breaks[counts$counts >= threshold])
  
  
  # Obtain values --------------------------------------------------------------
  center.plot <- max.bin/2
  separation.distance <- center.plot*0.10
  
  # Combine
  center <- list(center.plot = center.plot,
                 separation.distance = separation.distance)
  
  
  # Generate output ------------------------------------------------------------
  return(center)
}