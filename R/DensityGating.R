# Internal function - Perform density-based cut-off detection
.DensityGating <- function(og, cp.value) {
  
  
  # Input validation -----------------------------------------------------------
  checkmate::assert(methods::is(og, "flowFrame"), "Object is not a flowFrame.")
  checkmate::checkNumeric(cp.value)
  
  
  # Perform density-based cut-off detection ------------------------------------
  d <- lapply(names(flowCore::markernames(og)), function(channel) {
    cutoffs <- flowDensity::deGate(og, channel = channel, verbose = FALSE, all.cuts = TRUE, upper = TRUE)
    tibble::tibble(channel = channel, cutoffs = list(cutoffs))})
  
  # Alter format
  d <- dplyr::bind_rows(d) |>
    tidyr::unnest_wider(cutoffs, names_sep = "_") |>
    as.data.frame()
  
  
  # Retrieve value closest to the visual determined cut-off --------------------
  if (ncol(d) != 2) {
    d$opt <- apply(d[, grepl("cutoffs", names(d))], 1, function(row) .GetClosestCenter(row, closest = cp.value))
  } else {
    d$opt <- d[, grepl("cutoffs", names(d))]
  }
  
  
  # Obtain best gating ---------------------------------------------------------
  d <- d[, c("channel", "opt")]
  
  
  # Generate output ------------------------------------------------------------
  return(d)
}