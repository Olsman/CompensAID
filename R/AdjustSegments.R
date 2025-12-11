# Internal function - Adjust (merged) segment information
.AdjustSegments <- function(si.input, mp, ms) {
  
  
  # Input validation -----------------------------------------------------------
  checkmate::checkDataFrame(si.input)
  checkmate::checkCharacter(mp)
  checkmate::checkCharacter(ms)
  
  
  # Add temporary columns ------------------------------------------------------
  si.input$segment.min.merged <- NA
  si.input$segment.max.merged <- NA
  si.input$event.count.sum <- NA
  si.input$merged.segments <- NA
  
  
  # Merge and Adjust segments --------------------------------------------------
  for (grp in unique(si.input$mergeGroup)) {
    
    
    # Identify merge
    group.rows <- which(si.input$mergeGroup[si.input$primary.marker == mp &
                                            si.input$secondary.marker == ms] == grp)
    
    # Skip empty rows
    if (length(group.rows) == 0) { next }
    
    
    # Determine new segment ranges and event counts ----------------------------
    min.val <- min(si.input$segment.min[si.input$primary.marker == mp &
                                        si.input$secondary.marker == ms][group.rows])
    max.val <- max(si.input$segment.max[si.input$primary.marker == mp &
                                        si.input$secondary.marker == ms][group.rows])
    sum.val <- sum(si.input$event.count[si.input$primary.marker == mp &
                                        si.input$secondary.marker == ms][group.rows])
    segments.merged <- paste(si.input$segment[si.input$primary.marker == mp &
                                              si.input$secondary.marker == ms][group.rows], collapse = "+")
    
    
    # Input values -------------------------------------------------------------
    first.row <- group.rows[1]
    si.input$segment.min.merged[si.input$primary.marker == mp &
                                si.input$secondary.marker == ms][first.row] <- min.val
    si.input$segment.max.merged[si.input$primary.marker == mp &
                                si.input$secondary.marker == ms][first.row] <- max.val
    si.input$event.count.sum[si.input$primary.marker == mp &
                             si.input$secondary.marker == ms][first.row] <- sum.val
    si.input$merged.segments[si.input$primary.marker == mp &
                             si.input$secondary.marker == ms][first.row] <- segments.merged
  }
  
  
  # Clean dataFrame ------------------------------------------------------------
  si.input <- .CleanSegments(si.input, mp, ms)
  
  
  # Generate output ------------------------------------------------------------
  return(si.input)
}
