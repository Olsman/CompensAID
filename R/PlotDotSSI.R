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

    set.seed(20)
    p <- ggcyto::autoplot(og,
                          x = primary,
                          y = secondary,
                          bins = 100) +
      guides(fill = "none")+
      theme_minimal() +
      xlab(extraInformation$pretty.primary[1]) +
      ylab(extraInformation$pretty.secondary[1]) +
      theme(legend.position = 'right',
            legend.direction = 'horizontal',
            text = element_text(size = 12, family = 'Helvetica', face = 'bold'),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=2),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            strip.background = element_blank(),
            strip.text = element_blank()) +
      ggplot2::theme(aspect.ratio=1)
    p <- ggcyto::as.ggplot(p)

  } else {

    set.seed(20)
    p <- ggcyto::autoplot(og,
                          x = primary,
                          y = secondary,
                          bins = 100) +
      guides(fill = "none")+
      theme_minimal() +
      xlab(extraInformation$pretty.primary[1]) +
      ylab(extraInformation$pretty.secondary[1]) +
      theme(legend.position = 'right',
            legend.direction = 'horizontal',
            text = element_text(size = 12, family = 'Helvetica', face = 'bold'),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=2),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            strip.background = element_blank(),
            strip.text = element_blank()) +
      ggplot2::theme(aspect.ratio=1)
    p <- ggcyto::as.ggplot(p)

    # add annotation: adjust gating because of bins
    raw.annotation <- p$data %>%
      dplyr::select(extraInformation$primary.channel[1], extraInformation$secondary.channel[1]) %>%
      dplyr::filter(!!sym(extraInformation$secondary.channel[1]) < extraInformation$secondary.cutoff[1])

    min.neg <- min(raw.annotation[, extraInformation$primary.channel[1]])
    max.pos <- max(raw.annotation[, extraInformation$primary.channel[1]])

    min.sec <- min(raw.annotation[[extraInformation$primary.channel[2]]])
    max.sec <- max(raw.annotation[[extraInformation$primary.channel[2]]])

    # add annotation to figure
    p <- p + ggplot2::annotate("rect",
                               xmin = min.neg,
                               xmax = extraInformation$primary.cutoff.neg[1],
                               ymin = min.sec,
                               ymax = max.sec)
  }




    p <- p + ggplot2::annotate("rect",
                               xmin = p[["data"]] %>% as.data.frame() %>%
                                 dplyr::select(marker1, marker2) %>%
                                 dplyr::filter(.data[[marker2]] < extra$secondary.cutoff[1]) %>%
                                 dplyr::select(all_of(marker1)) %>%
                                 min(),
                               xmax = extra$primary.cutoff.neg,
                               ymin = p[["data"]] %>% as.data.frame() %>%
                                 dplyr::select(marker2) %>%
                                 base::min(.),
                               ymax = extra$secondary.cutoff[1],
                               alpha = 0,
                               color = "red",
                               lwd = 0.5) +
      ggplot2::geom_rect(xmin = min(extra$segment.min, na.rm = TRUE), xmax = max(extra$segment.max, na.rm = TRUE),
                         ymin = extra$secondary.cutoff[1]+ 0.1, ymax = extra$secondary.cutoff[1] + 1.3,
                         alpha = 0.01, fill = "white", colour = "black")

    for (range.visual in 1:length(!is.na(extra$ssi))) {
      p <- p +
        # borders of the ranges to calculate the ssi
        ggplot2::annotate("rect",
                          xmin = extra$segment.min[range.visual],
                          xmax = extra$segment.max[range.visual],
                          ymin = p[["data"]] %>%
                            as.data.frame() %>%
                            dplyr::select(marker2) %>%
                            base::min(.),
                          ymax = extra$secondary.cutoff[1],
                          alpha = 0,
                          color = "red",
                          lwd = 0.5) +
        # values of the ssi and mfi of the ranges
        ggplot2::annotate("text",
                          label = round(extra$ssi[range.visual], 2),
                          x = (extra$segment.min[range.visual] +
                                 (extra$segment.max[range.visual] - extra$segment.min[range.visual])/2),
                          # (extra.info.figure.sub$primary.cutoff.neg[range.visual]/2),
                          # y = min(ff@exprs[, c(marker1, marker2)]) + 0.1*min(ff@exprs[, c(marker1, marker2)]),
                          y = extra$secondary.cutoff[1]+ 0.7,
                          size = 3,
                          colour = "red",
                          fontface = "bold",
                          angle = 90)}

  }


  # Generate output ------------------------------------------------------------
  return(p)
}
