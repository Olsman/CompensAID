#' @title Assess event requirement
#'
#' @param negative (matrix): All events that fall within the negative primary population.
#' @param positive (matrix): All events that fall within the positive primary population.
#' @param events.value (numerical): Numerical value defining the minimum requirement of events.
#' 
#' @return (boolean) Returns TRUE/FALSE depending on the required number of events.
#' 
#' @seealso \code{\link{CompensAID}, \link{GetPopulations}, \link{}}
#' 
#' @examples 
#' # Import FCS file
#' file <- flowCore::read.FCS("path/to/exampleFCS.fcs")
#' 
#' # Marker names
#' primary.marker <- "CD19"
#' secondary.marker <- "IgL"
#' 
#' # Obtain density-based cut-off detection
#' center.plot <- 2
#' co <- DensityGating(og = file, 
#'                     cp.value = center.plot)
#'               
#' # Obtain populations
#' separation.distance <- 0.25
#' pop <- GetPopulations(og = file,
#'                       primary = primary.marker,
#'                       secondary = secondary.marker,
#'                       co.input = co,
#'                       sd.input = separation.distance)
#'                       
#' # Select populations                       
#' primary.negative <- pop$primary.negative
#' primary.positive <- pop$primary.positive
#'
#' # Assess event requirement
#' events.value <- 20
#' EventRequirement(primary.negative, primary.positive, events.value))
#' 
#' @export

EventRequirement <- function(negative, positive, events.value) {
  
  
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