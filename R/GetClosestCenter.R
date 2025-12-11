# Internal function - Identify value closest to visual estimation
.GetClosestCenter <- function(row, closest) {
  
  
  # Input validation -----------------------------------------------------------
  checkmate::checkNumeric(row)
  checkmate::checkNumeric(closest)
  
  
  # Remove NAs from the array --------------------------------------------------
  valid.index <- which(!is.na(row))
  valid.row <- row[valid.index]
  
  
  # Identify the index of the estimations closest to the visual estimation -----
  distances <- abs(row - closest)
  closest <- which.min(distances)
  
  # Get index
  index <- valid.index[closest]
  
  
  # Generate output ------------------------------------------------------------
  return(row[index])
}