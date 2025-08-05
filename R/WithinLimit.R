#' @title Detect if smaller peaks need included for density-based cut-off detection.
#'
#' @param population (list) List with the negative population of the primary and secondary marker, and positive population of the primary marker.
#' @param og (FlowFrame): FlowFrame containing the expression matrix, channel names, and marker names.
#' @param primary (character): Name of the primary channel.
#' @param secondary (character): Name of the secondary channel.
#' @param min (numerical): Minimum percentage of events required.
#' @param max (numerical): Maximum percentage of events required.
#' @param si.input (dataFrame): dataFrame containing SSI info.
#' @param co.input (dataFrame) dataFrame with the output of the density-based cut-off detection.
#' @param sd.input (numerical): Numerical value determining the distance between the primary negative and positive population.
#' @param cp.value (numerical): Numerical value determining the preliminary center.
#'
#' @return (list) Returns a list with the adjusted cutoffs and update secondary stain index information dataFrame.
#'
#' @seealso \code{\link{CompensAID}}, \code{\link{EmptyMatrixInfo}}, \code{\link{DensityGating}}
#'
#' @examples
#' # Import FCS file
#' file <- flowCore::read.FCS(system.file("extdata", "68983.fcs", package = "CompensAID"))
#'
#' # Obtain all possible marker combinations
#' mc <- GetMarkerCombinations(og = file)
#'
#' # Density-based cut-off detection
#' center.plot <- 2
#' co <- DensityGating(og = file,
#'                     cp.value = center.plot)
#'
#' # Empty SSI(-info) matrix/dataFrame
#' range.value <- 4
#' separation.distance <- 0.25
#' sm <- EmptyMatrix(og = file)
#' si <- EmptyMatrixInfo(og = file,
#'                       rv.input = range.value,
#'                       mc.input = mc,
#'                       co.input = co,
#'                       sd.input = separation.distance)
#'
#' # Gating positive and negative population
#' primary.marker <- "CD19"
#' secondary.marker <- "CD3"
#' pop <- GetPopulations(og = file,
#'                       primary = primary.marker,
#'                       secondary = secondary.marker,
#'                       co.input = co,
#'                       sd.input = separation.distance)
#'
#' # Adjust limits
#' pop.limit <- WithinLimit(population = pop,
#'                          og = file,
#'                          primary = primary.marker,
#'                          secondary = secondary.marker,
#'                          min = 10,
#'                          max = 90,
#'                          si.input = si,
#'                          sd.input = separation.distance,
#'                          co.input = co,
#'                          cp.value = center.plot)
#'
#' @export

WithinLimit <- function(population, og, primary, secondary, min = 10, max = 90, si.input, sd.input, co.input, cp.value) {


  # Input validation -----------------------------------------------------------
  checkmate::checkList(population)
  checkmate::assert(methods::is(og, "flowFrame"), "Object is not a flowFrame.")
  checkmate::checkCharacter(primary)
  checkmate::checkCharacter(secondary)
  checkmate::checkNumeric(min)
  checkmate::checkNumeric(max)
  checkmate::checkDataFrame(si.input)
  checkmate::checkNumeric(sd.input)
  checkmate::checkDataFrame(co.input)
  checkmate::checkNumeric(cp.value)


  # Obtain percentages from secondary marker -----------------------------------
  percentage <- round((nrow(population[["secondary.negative"]])*100)/nrow(og))

  # Get channel names
  cp <- names(flowCore::markernames(og))[flowCore::markernames(og) == primary]
  cs <- names(flowCore::markernames(og))[flowCore::markernames(og) == secondary]


  # Adjust density-based cut-off detection parameters --------------------------
  if (percentage <= min | percentage >= max) {

    # Perform density-based cut-off detection
    co.adjust <- flowDensity::deGate(og, channel = cs, all.cuts = TRUE, tinypeak.removal = 0.0001, verbose = FALSE, upper = TRUE)
    
    # Assess if new cut-off improves gating
    co.adjust <- GetClosestLimit(old.limit = co.input$opt[co.input$channel == cs], 
                                 new.limit = GetClosestCenter(co.adjust, cp.value),
                                 center.plot = cp.value)
    
    # Adjust values
    si.input$secondary.cutoff[si.input$primary.marker == primary & si.input$secondary.marker == secondary] <- GetClosestCenter(co.adjust, cp.value)
    co.input$opt[co.input$channel == cs] <- GetClosestCenter(co.adjust, cp.value)

    # Obtain populations with new cut-off
    pop <- GetPopulations(og = og,
                          primary = primary,
                          secondary = secondary,
                          co.input = co.input,
                          sd.input = sd.input)
  }


  # Obtain percentages from primary marker -----------------------------------
  percentage.positive <- round((nrow(population[["primary.positive"]])*100)/nrow(population[["secondary.negative"]]))
  percentage.negative <- round((nrow(population[["primary.negative"]])*100)/nrow(population[["secondary.negative"]]))


  # Adjust density-based cut-off detection parameters --------------------------
  if (percentage.positive <= min | percentage.positive >= max | percentage.negative <= min | percentage.negative >= max) {

    # Perform density-based cut-off detection
    co.adjust <- flowDensity::deGate(og, channel = cp, all.cuts = TRUE, tinypeak.removal = 0.0001, verbose = FALSE, upper = TRUE)
    
    # Asses
    # Assess if new cut-off improves gating
    co.adjust <- GetClosestLimit(old.limit = co.input$opt[co.input$channel == cp], 
                                 new.limit = GetClosestCenter(co.adjust, cp.value),
                                 center.plot = cp.value)
    
    # Adjust
    si.input$primary.cutoff.neg[si.input$primary.marker == primary & si.input$secondary.marker == secondary] <- GetClosestCenter(co.adjust, cp.value) - sd.input
    si.input$primary.cutoff.pos[si.input$primary.marker == primary & si.input$secondary.marker == secondary] <- GetClosestCenter(co.adjust, cp.value) + sd.input
    co.input$opt[co.input$channel == cp] <- GetClosestCenter(co.adjust, cp.value)

    # Obtain populations with new cut-off
    pop <- GetPopulations(og = og,
                          primary = primary,
                          secondary = secondary,
                          co.input = co.input,
                          sd.input = sd.input)
  }


  # Combine output -------------------------------------------------------------
  output <- list(co.dat = co.input,
                 si.dat = si.input)


  # Generate output ------------------------------------------------------------
  return(output)
}
