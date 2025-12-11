# Internal function - Clean adjust segment information
.CleanSegments <- function(si.input, mp, ms) {
  
  
  # Input validation -----------------------------------------------------------
  checkmate::checkDataFrame(si.input)
  checkmate::checkCharacter(mp)
  checkmate::checkCharacter(ms)
  
  
  # Move updated segments ------------------------------------------------------
  si.input$segment.min[si.input$primary.marker == mp &
                         si.input$secondary.marker == ms] <- si.input$segment.min.merged[si.input$primary.marker == mp &
                                                                                           si.input$secondary.marker == ms]
  si.input$segment.max[si.input$primary.marker == mp &
                         si.input$secondary.marker == ms] <- si.input$segment.max.merged[si.input$primary.marker == mp &
                                                                                           si.input$secondary.marker == ms]
  si.input$event.count.merge[si.input$primary.marker == mp &
                               si.input$secondary.marker == ms] <- si.input$event.count.sum[si.input$primary.marker == mp &
                                                                                              si.input$secondary.marker == ms]
  
  
  # Adjust the number of events after update -----------------------------------
  # No adjustment needed
  if (all(si.input$event.count.merge[si.input$primary.marker == mp &
                                     si.input$secondary.marker == ms] == si.input$event.count[si.input$primary.marker == mp &
                                                                                              si.input$secondary.marker == ms])) {
    
    # Set merge to NA
    si.input[si.input$primary.marker == mp &
               si.input$secondary.marker == ms &
               is.na(si.input$merged.segments), c("event.count.merge")] <- NA
    
  } else {
    
    # Loop per segment
    for (m in si.input$segment[si.input$primary.marker == mp & si.input$secondary.marker == ms]) {
      
      # Next loop if the number of events remains the same
      if (identical(si.input$event.count.merge[si.input$primary.marker == mp &
                                               si.input$secondary.marker == ms &
                                               si.input$segment == m],
                    si.input$event.count[si.input$primary.marker == mp &
                                         si.input$secondary.marker == ms &
                                         si.input$segment == m])) {
        
        next
        
      } else {
        
        # Add which segments are merged
        si.input$message[si.input$primary.marker == mp &
                           si.input$secondary.marker == ms &
                           si.input$segment == m] <- paste0("Merged segments: ", si.input$merged.segments[si.input$primary.marker == mp &
                                                                                                            si.input$secondary.marker == ms &
                                                                                                            !is.na(si.input$merged.segments) &
                                                                                                            si.input$segment == m])
      }
    }
  }
  
  # Set merge to NA
  si.input[si.input$primary.marker == mp &
             si.input$secondary.marker == ms &
             is.na(si.input$merged.segments), c("message", "mfi.neg", "mfi.pos", "sd.neg", "ssi")] <- NA
  
  # Remove temporary columns
  si.input <- subset(si.input, select = -c(mergeGroup, segment.min.merged, segment.max.merged, event.count.sum, merged.segments))
  
  
  # Generate output ------------------------------------------------------------
  return(si.input)
}