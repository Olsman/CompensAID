#' @title Plot dot plot of Secondary Stain Index scores
#' @description This function plots the marker combinations of interest. Additionally, the gating, segmentation, and SSI values can be visualized by setting showScore to TRUE. 
#' 
#' @param output (matrix): Matrix containing the SSI output from the CompensAID tool.
#' @param og (flowFrame): FlowFrame containing the expression matrix, channel names, and marker names.
#' @param primary (character): Primary marker (x-axis)
#' @param secondary (character): Seconday marker (y-axis)
#' @param showScores (boolean): Boolean variable determining if scores and gating should be visualized.
#'
#' @return (ggplot2) Returns figure containing a dot plot of the CompensAID output.
#'
#' @importFrom checkmate checkList assert checkCharacter
#' @importFrom methods is
#' @importFrom ggcyto autoplot as.ggplot
#' @importFrom ggplot2 guides theme_minimal xlab ylab theme element_text element_blank element_rect annotate
#' @importFrom dplyr filter select pull
#' @importFrom rlang sym
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
#' # Marker names
#' primary.marker <- "CD19"
#' secondary.marker <- "CD3"
#'
#' # Plot matrix
#' figure <- PlotDotSSI(output = compensAID.res,
#'                      og = file,
#'                      primary = primary.marker,
#'                      secondary = secondary.marker,
#'                      showScores = TRUE)
#' plot(figure)
#'
#' @export

PlotDotSSI <- function(output, og, primary, secondary, showScores = TRUE) {
  
  
  # Input validation -----------------------------------------------------------
  checkmate::assertList(output)
  checkmate::assert(methods::is(og, "flowFrame"), "Object is not a flowFrame.")
  checkmate::assertCharacter(primary)
  checkmate::assertCharacter(secondary)
  
  
  # Obtain output --------------------------------------------------------------
  extraInformation <- output[["matrixInfo"]] |>
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
    
    # Add annotation: adjust visualization because of bins from ggcyto
    raw.annotation <- p$data |>
      dplyr::select(extraInformation$primary.channel[1], extraInformation$secondary.channel[1]) |>
      dplyr::filter(!!rlang::sym(extraInformation$secondary.channel[1]) < extraInformation$secondary.cutoff[1])
    
    # Retrieve min and max values per 
    min.neg <- raw.annotation |> dplyr::select(extraInformation$primary.channel[1]) |> min()
    max.pos <- raw.annotation |> dplyr::select(extraInformation$primary.channel[1]) |> max()
    min.sec <- raw.annotation |> dplyr::select(extraInformation$secondary.channel[1]) |> min()
    max.sec <- raw.annotation |> dplyr::select(extraInformation$secondary.channel[1]) |> max()
    
    # Retrieve the range for visualization
    range <- raw.annotation |> 
      dplyr::filter(!!rlang::sym(extraInformation$primary.channel[1]) > extraInformation$primary.cutoff.pos[1]) |> 
      dplyr::pull(extraInformation$primary.channel[1])
    
    if (length(range) > 0) {
      range.value <- (max(range) - min(range)) / max(extraInformation$segment)
    } else {
      range.value <- Inf
    }
    
    
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
                        lwd = 0.5)
    
    if (!all(extraInformation$message == "No positive/negative population")) {
      # Visualize segments and SSI
      
      p <- p +
        ggplot2::annotate("rect",
                          xmin = extraInformation$primary.cutoff.pos[1],
                          xmax = max.pos,
                          ymin = max.sec,
                          ymax = max.sec + 0.5,
                          alpha = 0.8,
                          fill = "white",
                          color = "black",
                          lwd = 0.5)
      
      for (range.visual in extraInformation$segment[!is.na(extraInformation$ssi)]) {
        
        
        # Visualization distances for text and boxes
        if (extraInformation$message[range.visual] == "PASS") {
          dist <- range.visual
          dist.text <- 1
        } else {
          dist <- range.visual + lengths(gregexpr("\\+", extraInformation$message[range.visual])) 
          dist.text <- length(strsplit(extraInformation$message[range.visual], "\\+")[[1]])
        }
        
        # Visualize scores
        p <- p +
          ggplot2::annotate("rect",
                            xmin = extraInformation$primary.cutoff.pos[1],
                            xmax = extraInformation$primary.cutoff.pos[1] + (dist * range.value),
                            ymin = min.sec,
                            ymax = max.sec,
                            alpha = 0,
                            color = "red",
                            lwd = 0.5) +
          ggplot2::annotate("text",
                            label = ifelse(extraInformation$ssi[range.visual] == 0.00, 0, format(extraInformation$ssi[range.visual], nsmall = 2)),
                            x = extraInformation$primary.cutoff.pos[1] + ((range.visual-1) * range.value) + (dist.text * range.value)/2,
                            y = max.sec + 0.25,
                            size = 2,
                            colour = "red",
                            fontface = "bold",
                            angle = 90)
        
      }
    }
  }
  
  
  # Generate output ------------------------------------------------------------
  return(p)
}