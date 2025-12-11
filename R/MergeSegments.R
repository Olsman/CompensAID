# Internal function - Asssess which segments do not meet the requirements
.MergeSegments <- function(segment.event.count, ev.input) {
  
  
  # Input validation -----------------------------------------------------------
  checkmate::checkNumeric(segment.event.count)
  checkmate::checkNumeric(ev.input)
  
  
  # Temporary empty lists
  merged <- list()
  buffer <- 0
  temp <- c()
  indices <- c()
  merged.indices <- list()
  
  
  # Iterate over segments ------------------------------------------------------
  for (i in rev(seq_along(segment.event.count))) {
    
    # Identify event count and segment numbers
    buffer <- buffer + segment.event.count[i]
    temp <- c(temp, segment.event.count[i])
    indices <- c(indices, i)
    
    # Assess event requirement
    if (buffer >= ev.input) {
      merged <- append(merged, list(temp))
      merged.indices <- append(merged.indices, list(indices))
      
      # Reset for iteration
      buffer <- 0
      temp <- c()
      indices <- c()
    }
  }
  
  
  # Append segments not meeting the requirement --------------------------------
  if (length(temp) > 0) {
    
    if (length(merged.indices) > 0) {
      
      # Add remaining segments to the last merged group
      merged.indices[[length(merged.indices)]] <- c(merged.indices[[length(merged.indices)]], indices)
      
    } else {
      
      # Create groups
      merged.indices <- list(indices)
    }
  }
  
  
  # Generate output ------------------------------------------------------------
  return(merged.indices)
}