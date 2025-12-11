# Internal function - Assess event requirement
.EventRequirement <- function(negative, positive, events.value) {
  
  
  # Input validation -----------------------------------------------------------
  checkmate::checkMatrix(negative)
  checkmate::checkMatrix(positive)
  checkmate::checkNumeric(events.value)
  
  # Obtain numer of events within each population
  eventPos <- nrow(positive)
  eventNeg <- nrow(negative)
  
  
  # Assess event requirement ---------------------------------------------------
  outcome <- ifelse(eventPos < events.value | eventNeg < events.value, TRUE, FALSE)
  
  
  # Generate output ------------------------------------------------------------
  return(outcome)
}
