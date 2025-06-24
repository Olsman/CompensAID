#' @title Plot Secondary Stain Index matrix
#'
#' @param output.compensAID (matrix): Matrix containing the SSI output.
#'
#' @return (ggplot2) Returns figure containing a matrix of the total compensAID output.
#'
#' @seealso \code{\link{CompensAID}}
#'
#' @examples
#' # Import FCS file
#' file <- flowCore::read.FCS(system.file("extdata", "68983.fcs", package = "compensAID"))
#'
#' # Run compensAID tool
#' segments <- 4
#' min.events <- 20
#' center <- 2
#' distance.populations <- 0.25
#' compensAID.res <- CompensAID(ff = file,
#'                              range.value = segments,
#'                              events.value = min.events,
#'                              center.plot = center,
#'                              separation.distance = distance.populations)
#'
#' # Plot matrix
#' figure <- PlotMatrix(output.compensAID = compensAID.res)
#'
#' @export

PlotMatrix <- function(output.compensAID) {


  # Input validation -----------------------------------------------------------
  checkmate::checkList(output.compensAID)


  # Obtain output --------------------------------------------------------------
  comp.temp <- output.compensAID[["matrix"]]
  comp.temp[is.na(comp.temp)] <- 0
  colnames(comp.temp) <- output.compensAID[["matrixInfo"]]$pretty.primary[match(colnames(comp.temp), output.compensAID[["matrixInfo"]]$primary.channel)]
  rownames(comp.temp) <- output.compensAID[["matrixInfo"]]$pretty.secondary[match(rownames(comp.temp), output.compensAID[["matrixInfo"]]$secondary.channel)]

  dfm <- reshape2::melt(comp.temp) %>%
    dplyr::mutate(value = as.numeric(value)) %>%
    dplyr::mutate(value2 = ifelse(value < 1 & value > -1, NA, value))

  order.figure <- colnames(comp.temp)

  limit.max <- ceiling(max(dfm$value))
  limit.min <- floor(min(dfm$value))

  # Adjust if all are within threshold -----------------------------------------
  if (all(is.na(dfm$value2))) {

    dfm <- dfm %>%
      dplyr::mutate(value2 = 0)
    }


  # Plot SSI matrix ------------------------------------------------------------
  p <- dfm %>%
    dplyr::mutate(Var1 = factor(Var1, levels = rev(order.figure)),
                  Var2 = factor(Var2, levels = order.figure),
                  color = "#FFFFFF") %>%
    dplyr::mutate(color = ifelse(value <= -1, "#ef8a62", ifelse(value >= 1, "#67a9cf", color))) %>%
    ggplot2::ggplot(., ggplot2::aes(x = Var2, y = Var1, fill = color)) +
    ggplot2::geom_tile(ggplot2::aes(fill = color)) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(value == 0.00, 0, format(value, nsmall = 2))), size = 2, fontface = "bold") +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(x = "Primary marker",
                  y = "Secondary marker") +
    ggplot2::theme(legend.position = 'none',
                  legend.direction = 'vertical',
                  text = ggplot2::element_text(size = 12, family = 'Helvetica', face = 'bold'),
                  panel.grid.minor = ggplot2::element_blank(),
                  panel.grid.major = ggplot2::element_blank(),
                  panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 2),
                  axis.ticks.x = ggplot2::element_blank(),
                  axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
                  axis.ticks.y = ggplot2::element_blank(),
                  strip.background = ggplot2::element_blank(),
                  strip.text = ggplot2::element_blank(),
                  panel.background = ggplot2::element_rect(fill = 'white'))


  # Generate output ------------------------------------------------------------
  return(p)
}
