#' @title Calculate Secondary Stain Index
#'
#' @param si.input (dataFrame): DataFrame containing the SSI info.
#' @param primary.channel (character): Name of the primary channel.
#' @param secondary.channel (character): Name of the secondary channel.
#' @param segment (numerical): Numerical value determining for which segment the SSI should be calculated.
#'
#' @return (numerical) Returns the Secondary Stain Index score.
#'
#' @seealso \code{\link{CompensAID}}, \code{\link{EmptyMatrixInfo}}, \code{\link{DensityGating}}
#'
#' @examples
#' # Import FCS file
#' file <- flowCore::read.FCS("path/to/exampleFCS.fcs")
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
#' secondary.marker <- "IgL"
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
#'                          co.input = co)
#'
#' @export

WithinLimit <- function(population, og, primary, secondary, min = 10, max = 90, si.input, sd.input, co.input) {


  # Input validation -----------------------------------------------------------
  checkmate::checkList(population)
  checkmate::assert(is(og, "flowFrame"), "Object is not a flowFrame.")
  checkmate::checkCharacter(primary)
  checkmate::checkCharacter(secondary)
  checkmate::checkNumeric(min)
  checkmate::checkNumeric(max)
  checkmate::checkDataFrame(si.input)
  checkmate::checkNumeric(sd.input)
  checkmate::checkDataFrame(co.input)

  # Obtain percentages
  percentage <- round((nrow(population[["secondary.negative"]])*100)/nrow(og), 2)
  si.input$secondary.perc.before[si.input$primary.marker == primary & si.input$secondary.marker == secondary] <- percentage

  # Get channel names
  cp <- names(flowCore::markernames(og))[flowCore::markernames(og) == primary]
  cs <- names(flowCore::markernames(og))[flowCore::markernames(og) == secondary]


  # Adjust density-based cut-off detection parameters --------------------------
  if (percentage < min | percentage > max) {

    # Perform density-based cut-off detection
    co.adjust <- flowDensity::deGate(og, channel = cs, all.cuts = TRUE, tinypeak.removal = 0.0001, verbose = FALSE, upper = TRUE)
    si.input$secondary.cutoff[si.input$primary.marker == primary & si.input$secondary.marker == secondary] <- GetClosestCenter(co.adjust, 2)
    si.input$secondary.cutoff.adjusted[si.input$primary.marker == primary & si.input$secondary.marker == secondary] <- TRUE
    co.input$opt[co.input$channel == cs] <- GetClosestCenter(co.adjust, 2)

    # Obtain populations with new cut-off
    pop <- GetPopulations(og = ff,
                          primary = primary,
                          secondary = secondary,
                          co.input = co.input,
                          sd.input = sd.input)

    # Obtain percentages
    percentage <- round((nrow(pop[["secondary.negative"]])*100)/nrow(og), 2)
    si.input$secondary.perc.after[si.input$primary.marker == primary & si.input$secondary.marker == secondary] <- percentage

    # Adjust information
    si.input$secondary.cutoff.adjusted <- ifelse(si.input$secondary.perc.before == si.input$secondary.perc.after, FALSE, TRUE)


    # Density-based cut-off detection parameters remained the same -------------
    } else {

      # Adjust information
      si.input$secondary.cutoff.adjusted[si.input$primary.marker == primary & si.input$secondary.marker == secondary] <- FALSE
      si.input$secondary.perc.after[si.input$primary.marker == primary & si.input$secondary.marker == secondary] <- percentage
      }


  # Combine output -------------------------------------------------------------
  output <- list(co.dat = co.input,
                 si.dat = si.input)


  # Generate output ------------------------------------------------------------
  return(output)
}
