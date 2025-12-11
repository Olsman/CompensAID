# Internal function - Combine segments that do not meet the requirements
.UpdateSegments <- function(si.input, primary, secondary, ev.input, rv.input) {
  
  
  # Input validation -----------------------------------------------------------
  checkmate::checkDataFrame(si.input)
  checkmate::checkCharacter(primary)
  checkmate::checkCharacter(secondary)
  checkmate::checkNumeric(ev.input)
  checkmate::checkNumeric(rv.input)
  
  
  # Merge segments -------------------------------------------------------------
  if (si.input$message[si.input$primary.marker == primary &
                       si.input$secondary.marker == secondary][1] != "No positive/negative population") {
    
    # Identify which segments need to be merged
    merge <- .MergeSegments(si.input$event.count[si.input$primary.marker == primary &
                                                  si.input$secondary.marker == secondary],
                           ev.input)
    si.input$mergeGroup <- NA
    
    # Assign merge IDs
    for (i in seq_along(merge)) {
      
      si.input$mergeGroup[si.input$primary.marker == primary &
                            si.input$secondary.marker == secondary][merge[[i]]] <- i
    }
    
    # Adjust segment information
    si.input <- .AdjustSegments(si.input,
                               mp = primary,
                               ms = secondary)
  }
  
  
  # Generate output ------------------------------------------------------------
  return(si.input)
}