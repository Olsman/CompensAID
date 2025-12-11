#' @title CompensAID
#' @description Run CompensAID to assess the potential presence of reference errors within the FCS file. The outcome is a list which contains an SSI matrix. The matrix quickly visualizes the marker combinations for each marker combinations. Additional information regarding the number of events per segments, which segments were merged, and the corresponding SSI values can be found in the Info dataframe. 
#' 
#' @param ff (FlowFrame): FlowFrame containing the expression matrix, channel names, and marker names.
#' @param segment.value (numerical): Numerical value of the number of segments in the primary positive population.
#' @param events.value (numerical): Numerical value of the minimum number of events per population/segment.
#'
#' @return (list) Returns a list containing the full SSI output and SSI matrix
#'
#' @importFrom checkmate assert checkNumeric
#' @importFrom methods is
#' @importFrom ParallelLogger logInfo
#' @importFrom flowCore markernames
#' @importFrom dplyr filter pull
#' 
#' @examples
#' # Import FCS file
#' file <- flowCore::read.FCS(system.file("extdata", "68983.fcs", package = "CompensAID"))
#'
#' # Run CompensAID tool
#' compensAID.res <- CompensAID(ff = file)
#'
#' @export

CompensAID <- function(ff, segment.value = 4, events.value = 50) {
  
  
  # Input validation -----------------------------------------------------------
  checkmate::assert(methods::is(ff, "flowFrame"), "Object is not a flowFrame.")
  checkmate::checkNumeric(segment.value)
  checkmate::checkNumeric(events.value)
  sprintf('Importing sample: %s', ff@description[["FILENAME"]]) |> ParallelLogger::logInfo()
  
  
  # Modify expression matrix ---------------------------------------------------
  # Add a rows of zero's to prevent errors during the density-based cut-off detection
  ff@exprs <- rbind(ff@exprs, rep(0, ncol(ff@exprs)))
  
  
  # Obtain all possible marker combinations ------------------------------------
  mc <- .GetMarkerCombinations(og = ff)
  
  
  # Determine center plot and separation distance
  center <- .DetermineCenter(og = ff)
  center.plot <- center$center.plot
  separation.distance <- center$separation.distance
  
  # Perform density-based cut-off detection ------------------------------------
  co <- .DensityGating(og = ff,
                      cp.value = center.plot)
  
  
  # Empty SSI(-info) matrix/dataFrame ------------------------------------------
  sm <- .EmptyMatrix(og = ff)
  si <- .EmptyMatrixInfo(og = ff,
                        rv.input = segment.value,
                        mc.input = mc,
                        co.input = co,
                        sd.input = separation.distance)
  
  
  # Perform density-based cut-off detection ------------------------------------
  for (i in seq_len(nrow(mc))) {
    
    
    # Perform density-based cut-off detection
    co.renew <- co
    
    
    # Gating positive and negative population ----------------------------------
    pop <- .GetPopulations(og = ff,
                          primary = mc$primary.marker[i],
                          secondary = mc$secondary.marker[i],
                          co.input = co.renew,
                          sd.input = separation.distance)
    
    
    # Adjust limits ------------------------------------------------------------
    pop.limit <- .WithinLimit(population = pop,
                             og = ff,
                             primary = mc$primary.marker[i],
                             secondary = mc$secondary.marker[i],
                             min = 10, max = 90,
                             si.input = si,
                             sd.input = separation.distance,
                             co.input = co.renew,
                             cp.value = center.plot)
    
    # Update dataFrame
    co.renew <- pop.limit$co.dat
    si <- pop.limit$si.dat
    
    
    # Obtain dataFrame with adjusted limits ------------------------------------
    pop <- .GetPopulations(og = ff,
                          primary = mc$primary.marker[i],
                          secondary = mc$secondary.marker[i],
                          co.input = co.renew,
                          sd.input = separation.distance)
    
    
    # Assess presence of positive/negative population ----------------------------
    if (.EventRequirement(pop$primary.negative, pop$primary.positive, events.value)) {
      
      # Update SSI dataFrame
      si <- .UpdateMatrixInfo(si.input = si,
                             rv.input = segment.value,
                             range.input = range,
                             primary = mc$primary.marker[i],
                             secondary = mc$secondary.marker[i],
                             output = "No positive/negative population")
      
    } else {
      
      
      # Segment the positive population ----------------------------------------
      range <- (max(pop$primary.positive) - si$primary.cutoff.pos[si$primary.marker == mc$primary.marker[i] & si$secondary.marker == mc$secondary.marker[i]][1])/segment.value
      si <- .UpdateMatrixInfo(si.input = si,
                             rv.input = segment.value,
                             range.input = range,
                             primary = mc$primary.marker[i],
                             secondary = mc$secondary.marker[i],
                             output = "PASS",
                             population = pop)
    }
    
    
    # Update the segment of the positive population ----------------------------
    si <- .UpdateSegments(si.input = si,
                         primary = mc$primary.marker[i],
                         secondary = mc$secondary.marker[i],
                         ev.input = events.value,
                         rv.input = segment.value)
    
    
    # Update SSI score based on adjusted segments ------------------------------
    si <- .UpdateSSI(si.input = si,
                    rv.input = segment.value,
                    primary = mc$primary.marker[i],
                    secondary = mc$secondary.marker[i],
                    population = pop)
    
    
    # Fill out the SSI matrix ---------------------------------------------------
    # Obtain full channel names
    pc <- names(flowCore::markernames(ff))[flowCore::markernames(ff) == mc$primary.marker[i]]
    sc <- names(flowCore::markernames(ff))[flowCore::markernames(ff) == mc$secondary.marker[i]]
    
    # Adjust SSI matrix values
    if (all(is.na(si$ssi[si$primary.channel == pc & si$secondary.channel == sc]))) {
      sm[sc, pc] <- NA
    } else {
      sm[sc, pc] <- si |>
        dplyr::filter(primary.channel == pc, secondary.channel == sc, !is.na(ssi)) |>
        dplyr::filter(segment == max(segment)) |>
        dplyr::pull(ssi)
    }
  }
  
  
  # Generate output ------------------------------------------------------------
  result <- list("matrix" = sm,
                 "matrixInfo" = si)
  
  
  # Generate output ------------------------------------------------------------
  return(result)
}