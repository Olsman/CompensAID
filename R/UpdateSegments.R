#' @title Combine segments that do not meet the requirements
#'
#' @param si.input (dataFrame): dataFrame containing SSI info.
#' @param primary (character): Name of the primary marker.
#' @param secondary (character) Name of the secondary marker.
#' @param ev.input (numerical): Minimum required number of events.
#' @param rv.input (numerical): Number of segments.
#' 
#' @return (dataFrame) Returns a dataframe with an update SSI information dataFrame.
#' 
#' @seealso \code{\link{CompensAID}}, \code{\link{DensityGating}}, \code{\link{EmptyMatrixInfo}}, \code{\link{UpdateMatrixInfo}}, \code{\link{GetPopulations}}
#' 
#' @examples 
#' # Import FCS file
#' file <- flowCore::read.FCS("path/to/exampleFCS.fcs")
#' 
#' # Marker names
#' primary.marker <- "CD19"
#' secondary.marker <- "IgL"
#' 
#' # Channel names
#' cp <- "PE-Cy7-A"
#' cs <- "APC-H7-A"
#' 
#' # Density-based cut-off detection
#' center.plot <- 2
#' co <- DensityGating(og = file, 
#'                     cp.value = center.plot)
#' 
#' # Parameter for the distance between the primary positive and negative population
#' separation.distance <- 0.25
#' 
#' # Obtain empty SSI matrix information
#' range.value <- 4
#' si <- EmptyMatrixInfo(og = file, 
#'                       rv.input = range.value, 
#'                       mc.input = mc, co.input = co, sd.input = separation.distance)
#' 
#' # Obtain populations
#' pop <- GetPopulations(og = file, 
#'                       primary = primary.marker, 
#'                       secondary = secondary.marker, 
#'                       co.input = co, 
#'                       sd.input = separation.distance)
#' 
#' # Select primary positive population                      
#' population.positive <- pop$primary.positive
#'
#' # Get the events that fall within the fourth segment
#' segment.value <- 4
#' range.value <- (max(population.positive[, cp]) - min(population.positive[, cp])/segment.value                       
#' 
#' # Adjust the SSI information if there are more than 20 events in the (segmented) positive population (output = "PASS")
#' si <- UpdateMatrixInfo(si.input = si,
#'                        rv.input = segment.value,
#'                        range.input = range.value,
#'                        primary = primary.marker,
#'                        secondary = secondary.marker,
#'                        output = "PASS",
#'                        population = pop)
#' 
#' # Update segments
#' events.value <- 20                       
#' si <- updateSegments(si.input = si,
#'                      primary = primary.marker, 
#'                      secondary = secondary.marker, 
#'                      ev.input = events.value, 
#'                      rv.input = segment.value)        
#'         
#' @export

UpdateSegments <- function(si.input, primary, secondary, ev.input, rv.input) {
  
  
  # Input validation -----------------------------------------------------------
  checkmate::checkDataFrame(si.input)
  checkmate::checkCharacter(primary)
  checkmate::checkCharacter(secondary)
  checkmate::checkNumeric(ev.input)
  checkmate::checkNumeric(rv.input)
  
  
  # Update information per segment ---------------------------------------------
  for (s in c(1, rv.input, seq(rv.input)[c(-1, -rv.input)])) {
    
    # Skip segment if it does not have enough events
    if (si.input$message[si.input$primary.marker == primary &
                         si.input$secondary.marker == secondary &
                         si.input$segment == s] == "No positive/negative population") { 
      next 
      }
    
    
    # Adjust first segment -----------------------------------------------------
    if (s == 1 & si.input$event.count[si.input$primary.marker == primary &
                                    si.input$secondary.marker == secondary &
                                    si.input$segment == s] < ev.input) {
      
      # Adjust segment ranges
      si.input$segment.min[si.input$primary.marker == primary &
                         si.input$secondary.marker == secondary &
                         si.input$segment == 1 & 
                         si.input$event.count < ev.input] <- si.input$segment.max[si.input$primary.marker == primary &
                                                                                  si.input$secondary.marker == secondary &
                                                                                  si.input$segment == 2]
      
      # Update other information
      si.input[si.input$primary.marker == primary &
               si.input$secondary.marker == secondary &
               si.input$segment == 2,
               c("mfi.neg", "sd.neg", "mfi.pos", "ssi")] <- NA
      
      # Update other information
      si.input[si.input$primary.marker == primary &
               si.input$secondary.marker == secondary &
               si.input$segment == 2,
               c("message")] <- "Combined with next segment"
      
      # Calculate number of events
      si.input$event.count.merge[si.input$primary.marker == primary &
                                 si.input$secondary.marker == secondary &
                                 si.input$segment == 1] <- sum(si.input$event.count[si.input$primary.marker == primary &
                                                                                    si.input$secondary.marker == secondary &
                                                                                    si.input$segment == 1],
                                                               si.input$event.count[si.input$primary.marker == primary &
                                                                                    si.input$secondary.marker == secondary &
                                                                                    si.input$segment == 2])
    }
    
    
    # Adjust last segment ------------------------------------------------------
    if (s == rv.input & s != 1 & si.input$event.count[si.input$primary.marker == primary &
                                                      si.input$secondary.marker == secondary &
                                                      si.input$segment == s] < ev.input) {
      
      # Adjust segment ranges
      si.input$segment.min[si.input$primary.marker == primary &
                           si.input$secondary.marker == secondary &
                           si.input$segment == s & 
                           si.input$event.count < ev.input] <- si.input$segment.min[si.input$primary.marker == primary &
                                                                                    si.input$secondary.marker == secondary &
                                                                                    si.input$segment == s-1]
      
      # Update other information
      si.input[si.input$primary.marker == primary &
               si.input$secondary.marker == secondary &
               si.input$segment == s-1,
               c("mfi.neg", "sd.neg", "mfi.pos", "ssi")] <- NA
      
      # Update other information
      si.input[si.input$primary.marker == primary &
               si.input$secondary.marker == secondary &
               si.input$segment == s-1,
               c("message")] <- "Combined with next segment"
    
      # Update other information
      si.input[si.input$primary.marker == primary &
               si.input$secondary.marker == secondary &
               si.input$segment == s,
               c("message")] <- "Combined with previous segment"
      
      # Calculate number of events
      si.input$event.count.merge[si.input$primary.marker == primary &
                                 si.input$secondary.marker == secondary &
                                 si.input$segment == s] <- sum(si.input$event.count[si.input$primary.marker == primary &
                                                                                    si.input$secondary.marker == secondary &
                                                                                    si.input$segment == s],
                                                               si.input$event.count[si.input$primary.marker == primary &
                                                                                    si.input$secondary.marker == secondary &
                                                                                    si.input$segment == s-1])
      
      
      # Combine more segments if the merge was not enough ----------------------
      if (si.input$event.count.merge[si.input$primary.marker == primary &
                                     si.input$secondary.marker == secondary &
                                     si.input$segment == s] < ev.input & s != 2) {
        
        # Adjust segment ranges
        si.input$segment.min[si.input$primary.marker == primary &
                             si.input$secondary.marker == secondary &
                             si.input$segment == s & 
                             si.input$event.count < ev.input] <- si.input$segment.min[si.input$primary.marker == primary &
                                                                                      si.input$secondary.marker == secondary &
                                                                                      si.input$segment == s-2]
        
        # Adjust other information
        si.input[si.input$primary.marker == primary &
                 si.input$secondary.marker == secondary &
                 si.input$segment == s-2,
                 c("mfi.neg", "sd.neg", "mfi.pos", "ssi")] <- NA
        
        # Adjust other information
        si.input[si.input$primary.marker == primary &
                 si.input$secondary.marker == secondary &
                 si.input$segment == s-2,
                 c("message")] <- "Combined with next segment"
        
        # Calculate number of events
        si.input$event.count.merge[si.input$primary.marker == primary &
                                   si.input$secondary.marker == secondary &
                                   si.input$segment == s] <- sum(si.input$event.count[si.input$primary.marker == primary &
                                                                                      si.input$secondary.marker == secondary &
                                                                                      si.input$segment == s],
                                                                 si.input$event.count[si.input$primary.marker == primary &
                                                                                      si.input$secondary.marker == secondary &
                                                                                      si.input$segment == s-1],
                                                                 si.input$event.count[si.input$primary.marker == primary &
                                                                                      si.input$secondary.marker == secondary &
                                                                                      si.input$segment == s-2])
    }
    
      
    # Combine more segments if the merge was not enough ------------------------
    if (si.input$event.count.merge[si.input$primary.marker == primary &
                                   si.input$secondary.marker == secondary &
                                   si.input$segment == s] < ev.input & s != 2) {
      
      # Adjust segment ranges
      si.input$segment.min[si.input$primary.marker == primary &
                           si.input$secondary.marker == secondary &
                           si.input$segment == s & 
                           si.input$event.count < ev.input] <- si.input$segment.min[si.input$primary.marker == primary &
                                                                                    si.input$secondary.marker == secondary &
                                                                                    si.input$segment == s-3]
      
      # Adjust other information
      si.input[si.input$primary.marker == primary &
               si.input$secondary.marker == secondary &
               si.input$segment == s-3,
               c("mfi.neg", "sd.neg", "mfi.pos", "ssi")] <- NA
      
      # Adjust other information
      si.input[si.input$primary.marker == primary &
               si.input$secondary.marker == secondary &
               si.input$segment == s-3,
               c("message")] <- "Combined with next segment"
      
      # Calculate number of events
      si.input$event.count.merge[si.input$primary.marker == primary &
                                 si.input$secondary.marker == secondary &
                                 si.input$segment == s] <- sum(si.input$event.count[si.input$primary.marker == primary &
                                                                                    si.input$secondary.marker == secondary &
                                                                                    si.input$segment == s],
                                                               si.input$event.count[si.input$primary.marker == primary &
                                                                                    si.input$secondary.marker == secondary &
                                                                                    si.input$segment == s-1],
                                                               si.input$event.count[si.input$primary.marker == primary &
                                                                                    si.input$secondary.marker == secondary &
                                                                                    si.input$segment == s-2],
                                                               si.input$event.count[si.input$primary.marker == primary &
                                                                                    si.input$secondary.marker == secondary &
                                                                                    si.input$segment == s-3])
      }
    }
    
    
    # Assess if segments in the middle require merging -------------------------
    if (s != 1 & 
        s != rv.input & 
        si.input$message[si.input$primary.marker == primary &
                         si.input$secondary.marker == secondary &
                         si.input$segment == s ] == "PASS" &
        si.input$event.count[si.input$primary.marker == primary &
                             si.input$secondary.marker == secondary &
                             si.input$segment == s ] < ev.input) {
      
      # Adjust segment ranges
      si.input$segment.min[si.input$primary.marker == primary &
                           si.input$secondary.marker == secondary &
                           si.input$segment == s+1 & 
                           si.input$event.count < ev.input] <- si.input$segment.min[si.input$primary.marker == primary &
                                                                                    si.input$secondary.marker == secondary &
                                                                                    si.input$segment == s]
      
      # Adjust other information
      si.input[si.input$primary.marker == primary &
               si.input$secondary.marker == secondary &
               si.input$segment == s,
               c("mfi.neg", "sd.neg", "mfi.pos", "ssi")] <- NA
      
      # Adjust other information
      si.input[si.input$primary.marker == primary &
               si.input$secondary.marker == secondary &
               si.input$segment == s,
               c("message")] <- "Combined with next segment"
      
      # Adjust other information
      si.input[si.input$primary.marker == primary &
               si.input$secondary.marker == secondary &
               si.input$segment == s+1,
               c("message")] <- "Combined with previous segment"
      
      # Calculate number of events
      si.input$event.count.merge[si.input$primary.marker == primary &
                                 si.input$secondary.marker == secondary &
                                 si.input$segment == s+1] <- sum(si.input$event.count[si.input$primary.marker == primary &
                                                                                      si.input$secondary.marker == secondary &
                                                                                      si.input$segment == s],
                                                                 si.input$event.count[si.input$primary.marker == primary &
                                                                                      si.input$secondary.marker == secondary &
                                                                                      si.input$segment == s+1])
    }
  }
  
  
  # Generate output ------------------------------------------------------------
  return(si.input)
}