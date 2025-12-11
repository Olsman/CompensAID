# Internal function - Detect if smaller peaks inclusion improved the density-based cut-off detection.
.GetClosestLimit <- function(old.limit, new.limit, center.plot) {
  
  
  # Input validation -----------------------------------------------------------
  checkmate::checkNumeric(old.limit, len = 1)
  checkmate::checkNumeric(new.limit, len = 1)
  checkmate::checkNumeric(center.plot, len = 1)
  
  # Identify distance from center.plot
  distance.old <- abs(old.limit - center.plot)
  distance.new <- abs(new.limit - center.plot)
  
  # Determine closest
  if (distance.old < distance.new) {
    best.limit <- old.limit
  } else {
    best.limit <- new.limit 
  }
  
  
  # Generate output ------------------------------------------------------------
  return(best.limit)
}