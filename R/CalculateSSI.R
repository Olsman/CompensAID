# Internal function - Calculate Secondary Stain Index
.CalculateSSI <- function(si.input, primary.channel, secondary.channel, segment) {
  
  
  # Input validation -----------------------------------------------------------
  checkmate::checkDataFrame(si.input)
  checkmate::checkCharacter(primary.channel)
  checkmate::checkCharacter(secondary.channel)
  checkmate::checkNumeric(segment)
  
  
  # Calculate Secondary Stain Index --------------------------------------------
  ssi <- si.input[si.input$primary.channel == primary.channel &
                    si.input$secondary.channel == secondary.channel &
                    si.input$segment == segment,] |>
    dplyr::mutate(ssi = round((mfi.pos - mfi.neg)/(2*sd.neg), digits = 2)) |>
    dplyr::pull(ssi)
  
  
  # Generate output ------------------------------------------------------------
  return(ssi)
}