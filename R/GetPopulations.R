# Internal function - Obtain (segmented) populations
.GetPopulations <- function(og, primary, secondary, co.input, sd.input) {
  
  # Input validation -----------------------------------------------------------
  checkmate::assert(methods::is(og, "flowFrame"), "Object is not a flowFrame.")
  checkmate::checkCharacter(primary)
  checkmate::checkCharacter(secondary)
  checkmate::checkDataFrame(co.input)
  checkmate::checkNumeric(sd.input)
  
  # Empty list
  populations <- list()
  
  # Retrieve channels names and filtered expression matrix
  cp <- names(flowCore::markernames(og))[flowCore::markernames(og) == primary]
  cs <- names(flowCore::markernames(og))[flowCore::markernames(og) == secondary]
  sub <- flowCore::exprs(og)[, c(cp, cs)]
  
  # Retrieve gating
  co.s <- co.input$opt[co.input$channel == cs]
  co.p.n <- co.input$opt[co.input$channel == cp] - sd.input
  co.p.p <- co.input$opt[co.input$channel == cp] + sd.input
  
  
  # Gate populations -----------------------------------------------------------
  populations[["secondary.negative"]] <- sub[sub[, cs] <= co.s, , drop = FALSE]
  sub.s <- populations[["secondary.negative"]]
  populations[["primary.negative"]] <- sub.s[sub.s[, cp] <= co.p.n, , drop = FALSE]
  populations[["primary.positive"]] <- sub.s[sub.s[, cp] >= co.p.p, , drop = FALSE]
  
  
  # Generate output ------------------------------------------------------------
  return(populations)
}