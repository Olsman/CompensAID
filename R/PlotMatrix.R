#' @title Plot Secondary Stain Index matrix
#' @description This function plot the SSI matrix, quickly allowing the identification of potential reference errors if a SSI < -1 is obtained. 
#' 
#' @param output (matrix): Matrix containing the SSI output.
#'
#' @return (ggplot2) Returns figure containing a matrix of the total compensAID output.
#'
#' @importFrom checkmate checkList
#' @importFrom reshape2 melt
#' @importFrom dplyr mutate case_when
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_identity labs theme element_text element_blank element_rect
#' 
#' @seealso \code{\link{CompensAID}}
#'
#' @examples
#' # Import FCS file
#' file <- flowCore::read.FCS(system.file("extdata", "68983.fcs", package = "CompensAID"))
#'
#' # Run compensAID tool
#' compensAID.res <- CompensAID(ff = file)
#'
#' # Plot matrix
#' figure <- PlotMatrix(output = compensAID.res)
#' plot(figure)
#' 
#' @export

PlotMatrix <- function(output) {
  
  
  # Input validation -----------------------------------------------------------
  checkmate::assertList(output)
  
  
  # Obtain output --------------------------------------------------------------
  comp.temp <- output[["matrix"]]
  comp.temp[is.na(comp.temp)] <- 0
  colnames(comp.temp) <- output[["matrixInfo"]]$pretty.primary[match(colnames(comp.temp), output[["matrixInfo"]]$primary.channel)]
  rownames(comp.temp) <- output[["matrixInfo"]]$pretty.secondary[match(rownames(comp.temp), output[["matrixInfo"]]$secondary.channel)]
  
  dfm <- reshape2::melt(comp.temp) |>
    dplyr::mutate(value = as.numeric(value)) |>
    dplyr::mutate(value2 = ifelse(value < 1 & value > -1, NA, value))
  
  order.figure <- colnames(comp.temp)
  
  limit.max <- ceiling(max(dfm$value))
  limit.min <- floor(min(dfm$value))
  
  # Adjust if all are within threshold -----------------------------------------
  if (all(is.na(dfm$value2))) {
    
    dfm <- dfm |>
      dplyr::mutate(value2 = 0)
  }
  
  
  # Plot SSI matrix ------------------------------------------------------------
  max_alpha_value <- max(limit.max, limit.min)
  
  p.df <- dfm |>
    dplyr::mutate(Var1 = factor(Var1, levels = rev(order.figure)),
                  Var2 = factor(Var2, levels = order.figure),
                  color = "#FFFFFF") |>
    dplyr::mutate(color = ifelse(value <= -1, "#ef8a62", ifelse(value >= 1, "#67a9cf", color)),
                  alpha.val = dplyr::case_when(abs(value) < 1 ~ 1, TRUE ~ pmin(1, (abs(value) - 1) / (max_alpha_value - 1))))
  
  
  p <- ggplot2::ggplot(p.df, ggplot2::aes(x = Var2, y = Var1, fill = color)) +
    ggplot2::geom_tile(ggplot2::aes(fill = color, alpha = alpha.val)) +
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