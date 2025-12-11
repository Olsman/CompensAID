# Internal function - Obtain empty SSI dataFrame
.EmptyMatrixInfo <- function(og = ff,
                            rv.input = range.value,
                            mc.input = mc,
                            co.input = co,
                            sd.input = separation.distance) {
  
  
  # Input validation -----------------------------------------------------------
  checkmate::assert(methods::is(og, "flowFrame"), "Object is not a flowFrame.")
  checkmate::checkNumeric(rv.input)
  checkmate::checkDataFrame(mc.input)
  checkmate::checkDataFrame(co.input)
  checkmate::checkNumeric(sd.input)
  
  # Obtain sample name
  fn <- gsub(".*/", "", og@description[["FILENAME"]])
  
  
  # Obtain empty SSI dataFrame -------------------------------------------------
  mi <- mc.input |>
    dplyr::mutate(primary.channel = names(flowCore::markernames(og))[match(primary.marker, flowCore::markernames(og))],
                  secondary.channel = names(flowCore::markernames(og))[base::match(secondary.marker, flowCore::markernames(og))]) |>
    dplyr::mutate(pretty.primary = paste0(primary.channel, ": ", primary.marker),
                  pretty.secondary = paste0(secondary.channel, ": ", secondary.marker)) |>
    dplyr::mutate(file = fn,
                  primary.cutoff.neg = (co.input$opt[match(primary.channel, co.input$channel)] - sd.input),
                  primary.cutoff.pos = (co.input$opt[match(primary.channel, co.input$channel)] + sd.input),
                  secondary.cutoff = co.input$opt[match(secondary.channel, co.input$channel)],
                  mfi.neg = NA,
                  mfi.pos = NA,
                  sd.neg = NA,
                  ssi = NA,
                  message = NA,
                  segment.min = NA,
                  segment.max = NA,
                  event.count = NA,
                  event.count.merge = NA) |>
    (\(x) dplyr::slice(x, rep(seq_len(nrow(x)), each = rv.input)))() |>
    dplyr::group_by(primary.marker, secondary.marker) |>
    dplyr::mutate(segment = rep(seq_len(rv.input))) |>
    dplyr::select(file, primary.marker, primary.channel, pretty.primary, secondary.marker, secondary.channel, pretty.secondary,
                  segment, secondary.cutoff, primary.cutoff.neg, primary.cutoff.pos, segment.min, segment.max, event.count, event.count.merge,
                  mfi.neg, sd.neg, mfi.pos, ssi, message)
  
  
  # Generate output ------------------------------------------------------------
  return(mi)
}