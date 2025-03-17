#' @title Determine the density-based gating closest to the visual estimation of the center.plot
#'
#' @param row (double): Number of the row that should be checked
#' @param closest (double): Visual estimation of the center of the plot.
#'
#' @return (numerical) Returns value closest to the closest parameter
#' @export

closest_to_center <- function(row, closest) {

  # Input validation --------------------------------------------------------
  checkmate::assertDouble(row)
  checkmate::assertDouble(closest)


  # Calculate the absolute differences between each value in the row and the target
  distances <- abs(row - closest)

  # Find the index of the minimum distance
  closest_index <- which.min(distances)

  # Return the value in the row that is closest to the target
  return(row[closest_index])
}


#' @title Run compensAID tool; automated detection of compensation issues
#'
#' @param path.QC (character): FlowFrame containing the expression matrix, channel names, and marker names.
#' @param data.QC (character): Name of the marker presented on the x-axis.
#' @param range.value (numerical): Numerical value determining the number of segments in the positive population.
#' @param events.value (numerical): Numerical value determining the minimum number of events per segment.
#' @param center.plot (numerical): Numerical value determined by visual inspection the preliminary center.
#' @param seperation.distance (numerical): Numerical value determining the extra margin added between the postive and negative population.
#'
#' @return (dataFrame) Returns a dataframe containing the full SSI output and SSI matrix
#' @export
#'
#' ####ADJUST STILL
compensationError.range.gating.spectral <- function(path.QC, data.QC, range.value, events.value, center.plot, seperation.distance) {


  # Input validation --------------------------------------------------------
  checkmate::assertCharacter(data.QC)
  checkmate::assertAccess(path.QC, access = 'r')
  checkmate::checkNumeric(range.value)
  checkmate::checkNumeric(events.value)
  checkmate::checkNumeric(center.plot)
  checkmate::checkNumeric(seperation.distance)
  sprintf('Importing sample ID: %s', data.QC) %>% ParallelLogger::logInfo()




  # Import data --------------------------------------------------------
  ff.select <- flowCore::read.FCS(paste0(path.QC, data.QC), emptyValue = FALSE)
  ff.select@exprs <- rbind(ff.select@exprs, rep(0, ncol(ff.select@exprs)))




  # Empty SSI matrix --------------------------------------------------------
  final.SSI <- base::matrix(data = NA,
                            nrow = base::length(base::names(flowCore::markernames(ff.select))),
                            ncol = base::length(base::names(flowCore::markernames(ff.select))),
                            dimnames = base::list(base::names(flowCore::markernames(ff.select)), base::names(flowCore::markernames(ff.select))))




  # All marker combinations --------------------------------------------------------
  all.combination <- base::expand.grid("Primary" = base::names(flowCore::markernames(ff.select)),
                                       "Secondary" = base::names(flowCore::markernames(ff.select))) %>%
    dplyr::filter(!Primary == Secondary) %>%
    base::lapply(., function(x) utils::type.convert(base::as.character(x), as.is = TRUE)) %>%
    base::as.data.frame()

  # Density-based gating
  list.gating <- list()
  for (marker in base::names(flowCore::markernames(ff.select))) {
    list.gating[["primary"]][[marker]] <- flowDensity::deGate(ff.select, channel = marker, verbose = FALSE, all.cuts = TRUE, upper = TRUE)
  }

  # Bind list together
  list.gating <- tibble::enframe(list.gating[["primary"]], name = "markers", value = "V") %>%
    tidyr::unnest_wider(V, names_sep = "") %>%
    dplyr::mutate(across(everything(), ~replace(., is.na(.), 0)))

  # Assess gating closest to.center.plot
  if (ncol(list.gating) != 2) {
    list.gating$new <- apply(list.gating[, grepl("V", names(list.gating))], 1, function(row) closest_to_center(row, closest = center.plot))
  } else {
    list.gating$new <- list.gating$V1
  }

  # Obtain best gating
  list.gating <- list.gating[,c("new", "markers")]
  base::names(list.gating) <- c(base::paste0("Value.", data.QC),
                                base::paste0("Channel.", data.QC))




  # Gating and SSI calculation --------------------------------------------------------
  for (k in 1:nrow(all.combination)) {

    # Gating ----------
    # Select markers
    ff.df.total <- ff.select@exprs[, c(all.combination$Primary[k], all.combination$Secondary[k])]

    # # get the automated detected gates
    # marker.gates.total.sub <- list.gating.total %>% dplyr::select(paste0("Value.", data.QC),
    #                                                               paste0("Channel.", data.QC))

    # Obtain cut-off
    cutoff <- as.numeric(c(list.gating[,1][list.gating[,2] == all.combination$Primary[k]][1],
                           list.gating[,1][list.gating[,2] == all.combination$Secondary[k]][1]))


    # Gate negative population for the secondary marker (y-axis)
    # Additional distance
    n <- seperation.distance
    ff.df.secondary <- ff.df.total[ff.df.total[,all.combination$Secondary[k]] < (cutoff[2]), ,drop = FALSE]

    # Lower tiny.peak removal threshold if the gating results in <10% or more than 90% of the data within the gate (primary marker)
    ff.df.primary <- ff.df.total[ff.df.total[,all.combination$Primary[k]] < (cutoff[1]), ,drop = FALSE]
    if (nrow(ff.df.primary) > nrow(ff.df.total)*0.9 | nrow(ff.df.primary) < nrow(ff.df.total)*0.1) {
      temp.vals <- flowDensity::deGate(ff.select, channel = all.combination$Primary[k], all.cuts = TRUE, tinypeak.removal = 0.0001, verbose = FALSE, upper = TRUE)
      new.cutoff <- ifelse(cutoff[1] == closest_to_center(temp.vals, closest = center.plot), FALSE, TRUE)
      cutoff[1] <- closest_to_center(temp.vals, closest = center.plot)
      new.cutoff <- TRUE
    }

    # Lower tiny.peak removal threshold if the gating results in <10% or more than 90% of the data within the gate (secondary marker)
    if (nrow(ff.df.secondary) > nrow(ff.df.total)*0.9 | nrow(ff.df.secondary) < nrow(ff.df.total)*0.1) {
      temp.vals <- flowDensity::deGate(ff.select, channel = all.combination$Secondary[k], all.cuts = TRUE, tinypeak.removal = 0.0001, verbose = FALSE, upper = TRUE)
      new.cutoff <- ifelse(cutoff[2] == closest_to_center(temp.vals, closest = center.plot), FALSE, TRUE)
      cutoff[2] <- closest_to_center(temp.vals, closest = center.plot)
      ff.df.secondary = ff.df.total[ff.df.total[,all.combination$Secondary[k]] < (cutoff[2]), ,drop = FALSE]
    } else {
      new.cutoff <- FALSE
    }

    # Obtain positive and negative population of the primary marker
    ff.df.secondary.neg <- ff.df.secondary[ff.df.secondary[,all.combination$Primary[k]] < (cutoff[1] - n), ,drop = FALSE]
    ff.df.secondary.pos <- ff.df.secondary[ff.df.secondary[,all.combination$Primary[k]] > (cutoff[1] + n), ,drop = FALSE]



    # SSI calculation ----------
    # Assess presence of positive/negative population
    if (nrow(ff.df.secondary.neg) < events.value | nrow(ff.df.secondary.pos) < events.value) {

      # Construct additional information compensAID
      range <- NA
      extra.data <- data.frame(fileId = rep(data.QC, range.value),
                               primary.channel = rep(all.combination$Primary[k], range.value),
                               primary.chosen.cutoff.negative = rep((cutoff[1] - n), range.value),
                               primary.chosen.cutoff.positive = rep((cutoff[1] + n), range.value),
                               secondary.channel = rep(all.combination$Secondary[k], range.value),
                               secondary.chosen.cutoff = rep((cutoff[2]), range.value),
                               mfi.neg = NA,
                               sd.neg = NA,
                               message = rep("No positive/negative population", range.value),
                               range.n = paste0("range", seq(range.value)),
                               range.min = NA,
                               range.max = NA,
                               count = NA,
                               count.new = NA,
                               mfi.pos = NA,
                               ssi = NA,
                               cutoff.adjusted = new.cutoff) %>%
        dplyr::mutate(primary.marker = flowCore::markernames(ff.select)[base::match(primary.channel, base::names(flowCore::markernames(ff.select)))],
                      secondary.marker = flowCore::markernames(ff.select)[base::match(secondary.channel, base::names(flowCore::markernames(ff.select)))]) %>%
        dplyr::mutate(pretty.primary = base::paste0(primary.channel, ": ", primary.marker),
                      pretty.secondary = base::paste0(secondary.channel, ": ", secondary.marker))

    } else {

      # Determine the range - positive and negative population are present
      range <- (max(ff.df.secondary.pos[,1]) - (cutoff[1] + n))/range.value
      extra.data <- data.frame(fileId = rep(data.QC, range.value),
                               primary.channel = rep(all.combination$Primary[k], range.value),
                               primary.chosen.cutoff.negative = rep((cutoff[1] - n), range.value),
                               primary.chosen.cutoff.positive = rep((cutoff[1] + n), range.value),
                               secondary.channel = rep(all.combination$Secondary[k], range.value),
                               secondary.chosen.cutoff = rep((cutoff[2]), range.value),
                               mfi.neg = rep(base::mean(ff.df.secondary.neg[,all.combination$Secondary[k]]), range.value),
                               sd.neg = stats::sd(ff.df.secondary.neg[,all.combination$Secondary[k]]),
                               message = rep("No positive/negative population", range.value),
                               range.n = paste0("range", seq(range.value)),
                               range.min = NA,
                               range.max = NA,
                               count = NA,
                               count.new = NA,
                               mfi.pos = NA,
                               ssi = NA,
                               cutoff.adjusted = new.cutoff) %>%
        dplyr::mutate(primary.marker = flowCore::markernames(ff.select)[base::match(primary.channel, base::names(flowCore::markernames(ff.select)))],
                      secondary.marker = flowCore::markernames(ff.select)[base::match(secondary.channel, base::names(flowCore::markernames(ff.select)))]) %>%
        dplyr::mutate(pretty.primary = base::paste0(primary.channel, ": ", primary.marker),
                      pretty.secondary = base::paste0(secondary.channel, ": ", secondary.marker))

      # SSI calculation per range ----------
      for (test.val in seq(range.value)) {

        # First segment
        if (test.val == 1) {

          # Filter data
          a <- ff.df.secondary.pos[ff.df.secondary.pos[,1] >= (cutoff[1] + n) + ((test.val - 1) * range) &
                                     ff.df.secondary.pos[,1] <= (cutoff[1] + n) + (test.val * range),
                                   ,drop = FALSE]

          # Asses criteria of min num of events
          if (nrow(a) >= events.value) {

            # Add details
            extra.data$count[test.val] <- nrow(a)
            extra.data$range.min[test.val] <- (cutoff[1] + n) + ((test.val - 1) * range)
            extra.data$range.max[test.val] <- (cutoff[1] + n) + (test.val * range)
            extra.data$mfi.pos[test.val] <- ff.df.secondary.pos[,2][ff.df.secondary.pos[,1] > extra.data$range.min[test.val] &
                                                                      ff.df.secondary.pos[,1] < extra.data$range.max[test.val]] %>% base::mean()
            extra.data$ssi[test.val] <- round((extra.data$mfi.pos[test.val] - extra.data$mfi.neg[test.val])/(2*extra.data$sd.neg[test.val]), digits = 2)
            extra.data$message[extra.data$range.n == paste0("range", test.val)] <- "PASS"

            # If first segments has < min number of events
          } else {

            # Add number of events
            extra.data$count[test.val] <- nrow(a)

            # Combine with second segments
            a <- ff.df.secondary.pos[ff.df.secondary.pos[,1] >= (cutoff[1] + n) + ((test.val - 1) * range) &
                                       ff.df.secondary.pos[,1] <= (cutoff[1] + n) + ((test.val + 1) * range),
                                     ,drop = FALSE]

            # Add details
            extra.data$range.min[test.val] <- (cutoff[1] + n) + ((test.val - 1) * range)
            extra.data$range.max[test.val] <- (cutoff[1] + n) + ((test.val + 1) * range)
            extra.data$mfi.pos[test.val] <- ff.df.secondary.pos[,2][ff.df.secondary.pos[,1] > extra.data$range.min[test.val] &
                                                                      ff.df.secondary.pos[,1] < extra.data$range.max[test.val]] %>% base::mean()
            extra.data$ssi[test.val] <- round((extra.data$mfi.pos[test.val] - extra.data$mfi.neg[test.val])/(2*extra.data$sd.neg[test.val]), digits = 2)
            extra.data$message[extra.data$range.n == paste0("range", test.val)] <- "Combined with next segment"


          }

          # All other segments segments
        } else {


          # Check: was the segment combined with the previous segment
          if (extra.data$message[extra.data$range.n == paste0("range", test.val-1)] == "Combined with next segment") {

            # # Filter data
            # a <- ff.df.secondary.pos[ff.df.secondary.pos[,1] >= (cutoff[1] + n) + ((test.val - 1) * range) &
            #                          ff.df.secondary.pos[,1] <= (cutoff[1] + n) + (test.val * range),
            #                          ,drop = FALSE]

            # Add number of events
            extra.data$count[test.val] <- nrow(a)

            # Filter data
            a <- ff.df.secondary.pos[ff.df.secondary.pos[,1] >= (cutoff[1] + n) + ((test.val - 1) * range) &
                                       ff.df.secondary.pos[,1] <= (cutoff[1] + n) + (test.val * range),
                                     ,drop = FALSE]

            # Add details - NA because combined with previous segment
            extra.data$count.new[test.val] <- NA
            extra.data$range.min[test.val] <- NA
            extra.data$range.max[test.val] <- NA
            extra.data$mfi.pos[test.val] <- NA
            extra.data$ssi[test.val] <- NA
            extra.data$message[extra.data$range.n == paste0("range", test.val)] <- "Combined with previous segment"


            # Next: segment not combined with any other segment
          } else {

            # Filter data
            a <- ff.df.secondary.pos[ff.df.secondary.pos[,1] >= (cutoff[1] + n) + ((test.val - 1) * range) &
                                       ff.df.secondary.pos[,1] <= (cutoff[1] + n) + (test.val * range),
                                     ,drop = FALSE]

            # Check: does the last segment have enough events?
            if (test.val == max(range.value) & nrow(a) < events.value) {

              # Add details - last segment does not have enough events
              extra.data$count[test.val] <- nrow(a)
              extra.data$range.min[test.val] <- NA
              extra.data$range.max[test.val] <- NA
              extra.data$mfi.pos[test.val] <- NA
              extra.data$ssi[test.val] <- NA
              extra.data$message[extra.data$range.n == paste0("range", test.val)] <- "Last segment: too few events"
              extra.data$range.max[extra.data$range.n == paste0(na.omit(extra.data$range.n[extra.data$range.max == max(extra.data$range.max, na.rm = TRUE)]))] <- (cutoff[1] + n) + (test.val * range)
            }

            # Segment has enough events?
            if (nrow(a) >= events.value) {

              # Add details
              extra.data$count[test.val] <- nrow(a)
              extra.data$range.min[test.val] <- (cutoff[1] + n) + ((test.val - 1) * range)
              extra.data$range.max[test.val] <- (cutoff[1] + n) + (test.val * range)
              extra.data$mfi.pos[test.val] <- ff.df.secondary.pos[,2][ff.df.secondary.pos[,1] > extra.data$range.min[test.val] &
                                                                        ff.df.secondary.pos[,1] < extra.data$range.max[test.val]] %>% base::mean()
              extra.data$ssi[test.val] <- round((extra.data$mfi.pos[test.val] - extra.data$mfi.neg[test.val])/(2*extra.data$sd.neg[test.val]), digits = 2)
              extra.data$message[extra.data$range.n == paste0("range", test.val)] <- "PASS"

              # Segment (that is not the first/last segment) has enough events?
            } else if (test.val != max(range.value) && nrow(a) < events.value) {

              # Add number of events
              extra.data$count[test.val] <- nrow(a)

              # Filter data
              a <- ff.df.secondary.pos[ff.df.secondary.pos[,1] >= (cutoff[1] + n) + ((test.val - 1) * range) &
                                         ff.df.secondary.pos[,1] <= (cutoff[1] + n) + ((test.val + 1) * range),
                                       ,drop = FALSE]

              # Add details
              extra.data$count.new[test.val] <- nrow(a)
              extra.data$range.min[test.val] <- (cutoff[1] + n) + ((test.val - 1) * range)
              extra.data$range.max[test.val] <- (cutoff[1] + n) + ((test.val + 1) * range)
              extra.data$mfi.pos[test.val] <- ff.df.secondary.pos[,2][ff.df.secondary.pos[,1] > extra.data$range.min[test.val] &
                                                                        ff.df.secondary.pos[,1] < extra.data$range.max[test.val]] %>% base::mean()
              extra.data$ssi[test.val] <- round((extra.data$mfi.pos[test.val] - extra.data$mfi.neg[test.val])/(2*extra.data$sd.neg[test.val]), digits = 2)
              extra.data$message[extra.data$range.n == paste0("range", test.val)] <- "Combined with next segment"
            }
          }
        }
      }

      # Segments after combining still too small?
      if (!is.na(any(extra.data$count.new < events.value))) {

        # Is one segment still too small?
        if (length(which(extra.data$count.new < events.value)) == 1) {

          # Get info from segments
          extra.data.sub <- slice_head(extra.data, n = min(which(extra.data$count.new < events.value)))
          min.sub <- extra.data.sub$range.min[!is.na(extra.data.sub$range.min)][length(extra.data.sub$range.min[!is.na(extra.data.sub$range.min)])-1]
          max.sub <- extra.data.sub$range.max[which(extra.data$count.new < events.value)]

          # Filter data
          a <- ff.df.secondary.pos[ff.df.secondary.pos[,1] >= min.sub &
                                     ff.df.secondary.pos[,1] <= max.sub,
                                   ,drop = FALSE]

          # Combine segment ranges
          range.end <- min(which(extra.data$count.new < events.value))
          range.start <- which(!is.na(extra.data.sub$range.min))
          range.start <- max(range.start[!range.start %in% which(extra.data$count.new < events.value)])

          # Adjust information
          extra.data$range.max[extra.data$range.n == paste0("range", range.start)] <- max.sub
          extra.data$count.new[extra.data$range.n == paste0("range", range.start)] <- nrow(a)
          extra.data$mfi.pos[extra.data$range.n == paste0("range", range.start)] <- ff.df.secondary.pos[,2][ff.df.secondary.pos[,1] > extra.data$range.min[range.start] & ff.df.secondary.pos[,1] < extra.data$range.max[range.start]] %>% base::mean()
          extra.data$ssi[extra.data$range.n == paste0("range", range.start)] <- round((extra.data$mfi.pos[range.start] -
                                                                                         extra.data$mfi.neg[range.start])/(2*extra.data$sd.neg[range.start]), digits = 2)
          extra.data$message[extra.data$range.n == paste0("range", range.start)] <- "Multiple segments"

          # Remove old data
          extra.data$range.min[extra.data$range.n == paste0("range", range.end)] <- NA
          extra.data$range.max[extra.data$range.n == paste0("range", range.end)] <- NA
          extra.data$count.new[extra.data$range.n == paste0("range", range.end)] <- NA
          extra.data$mfi.pos[extra.data$range.n == paste0("range", range.end)] <- NA
          extra.data$ssi[extra.data$range.n == paste0("range", range.end)] <- NA
          extra.data$message[extra.data$range.n == paste0("range", range.end)] <- "Combined with previous segment"

          # Are multiple segments still too small?
        } else {

          # Get info from segments
          extra.data.sub.max <- slice_head(extra.data, n = max(which(extra.data$count.new < events.value)))
          extra.data.sub.min <- slice_head(extra.data, n = min(which(extra.data$count.new < events.value)))
          min.sub <- extra.data.sub.max$range.min[!is.na(extra.data.sub.max$range.min)][length(extra.data.sub.min$range.min[!is.na(extra.data.sub.min$range.min)])-1]
          max.sub <- extra.data.sub.max$range.max[max(which(extra.data$count.new < events.value))]

          # Filter data
          a <- ff.df.secondary.pos[ff.df.secondary.pos[,1] >= min.sub &
                                     ff.df.secondary.pos[,1] <= max.sub,
                                   ,drop = FALSE]

          # Combine segment ranges
          range.end <- max(which(extra.data$count.new < events.value))
          range.start <- which(!is.na(extra.data.sub.max$range.min))
          range.start <- max(range.start[!range.start %in% which(extra.data$count.new < events.value)])

          # Adjust information
          extra.data$range.max[extra.data$range.n == paste0("range", range.start)] <- max.sub
          extra.data$count.new[extra.data$range.n == paste0("range", range.start)] <- nrow(a)
          extra.data$mfi.pos[extra.data$range.n == paste0("range", range.start)] <- ff.df.secondary.pos[,2][ff.df.secondary.pos[,1] > extra.data$range.min[range.start] & ff.df.secondary.pos[,1] < extra.data$range.max[range.start]] %>% base::mean()
          extra.data$ssi[extra.data$range.n == paste0("range", range.start)] <- round((extra.data$mfi.pos[range.start] -
                                                                                         extra.data$mfi.neg[range.start])/(2*extra.data$sd.neg[range.start]), digits = 2)
          extra.data$message[extra.data$range.n == paste0("range", range.start)] <- "Multiple segments"

          # Remove old data
          for (j in which(extra.data$count.new < events.value)) {
            extra.data$range.min[extra.data$range.n == paste0("range", j)] <- NA
            extra.data$range.max[extra.data$range.n == paste0("range", j)] <- NA
            extra.data$count.new[extra.data$range.n == paste0("range", j)] <- NA
            extra.data$mfi.pos[extra.data$range.n == paste0("range", j)] <- NA
            extra.data$ssi[extra.data$range.n == paste0("range", j)] <- NA
            extra.data$message[extra.data$range.n == paste0("range", j)] <- "Combined with previous segment"
          }
        }
      }
    }




    # Add outcome to SSI matrix --------------------------------------------------------
    if (all(is.na(extra.data$ssi))) {
      final.SSI[all.combination$Primary[k], all.combination$Secondary[k]] <- NA
    } else {
      final.SSI[all.combination$Primary[k], all.combination$Secondary[k]] <- extra.data %>%
        dplyr::filter(ssi == min(ssi, na.rm = TRUE)) %>%
        dplyr::distinct(ssi, .keep_all = TRUE) %>%
        dplyr::pull(ssi)
    }




    # Combine all outcomes --------------------------------------------------------
    if (k == 1) {
      extra.data.combined <- extra.data
    } else {
      extra.data.combined <- base::rbind(extra.data.combined, extra.data)
    }
  }




  # Adjust names --------------------------------------------------------
  base::colnames(final.SSI) <- base::paste0(base::colnames(final.SSI), ": ",
                                            flowCore::markernames(ff.select)[base::match(base::colnames(final.SSI),
                                                                                         base::names(flowCore::markernames(ff.select)))])
  base::rownames(final.SSI) <- base::paste0(base::rownames(final.SSI), ": ",
                                            flowCore::markernames(ff.select)[base::match(base::rownames(final.SSI),
                                                                                         base::names(flowCore::markernames(ff.select)))])



  # Output --------------------------------------------------------
  result <- list("SSI.desc" = final.SSI,
                 "Extra" = extra.data.combined)

  return(result)
}
