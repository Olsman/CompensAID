#' @title Plot dot plot of Secondary Stain Index scores
#'
#' @param output.compensAID (matrix): Matrix containing the SSI output.
#' @param og (flowFrame): FlowFrame containing the expression matrix, channel names, and marker names.
#' @param primary (character): Minimum SSI value shown in matrix.
#' @param secondary (character): Maximum SSI value shown in matrix.
#' @param showScores (boolean): Boolean variable determining if scores and gating should be visualized.
#'
#' @return (ggplot2) Returns figure containing a dot plot of the compensAID output.
#'
#' @seealso \code{\link{CompensAID}}
#'
#' @examples
#' # Import FCS file
#' file <- flowCore::read.FCS(system.file("extdata", "68983.fcs", package = "CompensAID"))
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
#' # Marker names
#' primary.marker <- "CD19"
#' secondary.marker <- "CD3"
#'
#' # Plot matrix
#' figure <- PlotDotSSI(output.compensAID = compensAID.res,
#'                      og = file,
#'                      primary = primary.marker,
#'                      secondary = secondary.marker,
#'                      showScores = TRUE)
#'
#' @export

PlotDotSSI <- function(output.compensAID, og, primary, secondary, showScores = FALSE) {


  # Input validation -----------------------------------------------------------
  checkmate::checkMatrix(output.compensAID)
  checkmate::assert(methods::is(og, "flowFrame"), "Object is not a flowFrame.")
  checkmate::checkCharacter(primary)
  checkmate::checkCharacter(secondary)


  # Obtain output --------------------------------------------------------------
  extraInformation <- output.compensAID[["matrixInfo"]] %>%
    dplyr::filter(primary.marker == primary & secondary.marker == secondary)


  if(isFALSE(showScores)) {


    # Visualize dot plot -------------------------------------------------------
    p <- ggcyto::autoplot(og,
                          x = primary,
                          y = secondary,
                          bins = 100) +
      ggplot2::guides(fill = "none")+
      ggplot2::theme_minimal() +
      ggplot2::xlab(extraInformation$pretty.primary[1]) +
      ggplot2::ylab(extraInformation$pretty.secondary[1]) +
      ggplot2::theme(legend.position = 'right',
                     legend.direction = 'horizontal',
                     text = ggplot2::element_text(size = 12, family = 'Helvetica', face = 'bold'),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 2),
                     axis.ticks.x = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     strip.background = ggplot2::element_blank(),
                     strip.text = ggplot2::element_blank()) +
      ggplot2::theme(aspect.ratio=1)
    p <- ggcyto::as.ggplot(p)


  } else {


    # Visualize dot plot -------------------------------------------------------
    p <- ggcyto::autoplot(og,
                          x = primary,
                          y = secondary,
                          bins = 100) +
      ggplot2::guides(fill = "none") +
      ggplot2::theme_minimal() +
      ggplot2::xlab(extraInformation$pretty.primary[1]) +
      ggplot2::ylab(extraInformation$pretty.secondary[1]) +
      ggplot2::theme(legend.position = 'right',
                     legend.direction = 'horizontal',
                     text = ggplot2::element_text(size = 12, family = 'Helvetica', face = 'bold'),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 2),
                     axis.ticks.x = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     strip.background = ggplot2::element_blank(),
                     strip.text = ggplot2::element_blank()) +
      ggplot2::theme(aspect.ratio=1)
    p <- ggcyto::as.ggplot(p)

    # Add annotation: adjust gating because of bins
    raw.annotation <- p$data %>%
      dplyr::select(extraInformation$primary.channel[1], extraInformation$secondary.channel[1]) %>%
      dplyr::filter(!!rlang::sym(extraInformation$secondary.channel[1]) < extraInformation$secondary.cutoff[1])

    min.neg <- raw.annotation %>% dplyr::select(extraInformation$primary.channel[1]) %>% min()
    max.pos <- raw.annotation %>% dplyr::select(extraInformation$primary.channel[1]) %>% max()

    min.sec <- raw.annotation %>% dplyr::select(extraInformation$secondary.channel[1]) %>% min()
    max.sec <- raw.annotation %>% dplyr::select(extraInformation$secondary.channel[1]) %>% max()

    range <- raw.annotation %>%
      dplyr::filter(!!rlang::sym(extraInformation$primary.channel[1]) > extraInformation$primary.cutoff.pos[1]) %>%
      dplyr::pull(extraInformation$primary.channel[1])

    range.value <- (max(range)-min(range))/max(extraInformation$segment)


    # Show scores --------------------------------------------------------------
    # Visualize gating
    p <- p +
      ggplot2::annotate("rect",
                        xmin = min.neg,
                        xmax = extraInformation$primary.cutoff.neg[1],
                        ymin = min.sec,
                        ymax = max.sec,
                        alpha = 0,
                        color = "red",
                        lwd = 0.5) +
      ggplot2::annotate("rect",
                        xmin = extraInformation$primary.cutoff.pos[1],
                        xmax = max.pos,
                        ymin = min.sec,
                        ymax = max.sec,
                        alpha = 0,
                        color = "red",
                        lwd = 0.5) +
      ggplot2::annotate("rect",
                        xmin = extraInformation$primary.cutoff.pos[1],
                        xmax = max.pos,
                        ymin = max.sec,
                        ymax = max.sec + 0.5,
                        alpha = 0.8,
                        fill = "white",
                        color = "black",
                        lwd = 0.5)

    # Visualize segments and SSI
    for (range.visual in seq_len(length(!is.na(extraInformation$ssi)))) {

      if (is.na(extraInformation$ssi)[range.visual]) { next }
      p <- p +
        ggplot2::annotate("rect",
                          xmin = extraInformation$primary.cutoff.pos[1],
                          xmax = extraInformation$primary.cutoff.pos[1] + range.visual*range.value,
                          ymin = min.sec,
                          ymax = max.sec,
                          alpha = 0,
                          color = "red",
                          lwd = 0.5) +
        ggplot2::annotate("text",
                          label = ifelse(extraInformation$ssi[range.visual] == 0.00, 0, format(extraInformation$ssi[range.visual], nsmall = 2)),
                          x = extraInformation$primary.cutoff.pos[1] + (range.visual*range.value) - (range.value/2),
                          y = max.sec + 0.25,
                          size = 2,
                          colour = "red",
                          fontface = "bold",
                          angle = 90)

      }
  }


  # Generate output ------------------------------------------------------------
  return(p)
}
