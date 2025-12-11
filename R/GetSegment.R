# Internal function - Calculate segment information
.GetSegment <- function(population, primary.channel, secondary.channel, segment, range.input) {
  
  
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
    div.mfi.pos <- stats::median(div[, secondary.channel])
    
    
    # Calculate segment information of the other segment(s) --------------------
  } else {
    
    # Filter segment
    min.div <- min(population[, primary.channel])
    max.div <- max(population[population[, primary.channel] <= (min.div + (segment*range.input)), , drop = FALSE][,primary.channel])
    min.div <- max(population[population[, primary.channel] <= (min.div + ((segment-1)*range.input)), , drop = FALSE][,primary.channel])
    div <- population[population[, primary.channel] >= min.div & population[, primary.channel] <= max.div, , drop = FALSE]
    
    # Determine the number of events and mean fluorescence intensity
    div.count <- nrow(div)
    div.mfi.pos <- stats::median(div[, secondary.channel])
  }
  
  
  # Combine segment information ------------------------------------------------
  segment.info <- data.frame("min" = min.div,
                             "max" = max.div,
                             "count" = div.count,
                             "mfi.pos" = div.mfi.pos)
  
  
  # Generate output ------------------------------------------------------------
  return(segment.info)
}
