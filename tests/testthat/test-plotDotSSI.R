testthat::test_that("Plot Dot Plot Errors", {
  
  flowFrame <- flowCore::read.FCS(system.file("extdata", "68983.fcs", package = "CompensAID"))
  compensAID.res <- CompensAID::CompensAID(flowFrame)
  matrix.SSI <- compensAID.res$matrix
  matrix.info <- compensAID.res$matrixInfo
  
  # Correct
  testthat::expect_no_error(PlotDotSSI(compensAID.res,
                                                 flowFrame,
                                                 "GFP",
                                                 "CD8",
                                                 showScores = TRUE))
  
  # No flowFrame added to function
  testthat::expect_error(PlotDotSSI(compensAID.res,
                                              "GFP",
                                              "CD8",
                                              showScores = TRUE))
  
  # Channel names instead of marker names
  testthat::expect_error(PlotDotSSI(compensAID.res,
                                              flowFrame,
                                              "FITC-A",
                                              "Pacific Blue-A",
                                              showScores = TRUE))
  
  # Partial CompensAID output used
  testthat::expect_error(PlotDotSSI(matrix.SSI,
                                              flowFrame,
                                              "GFP",
                                              "CD8",
                                              showScores = TRUE))
  
  # Partial CompensAID output used
  testthat::expect_error(PlotDotSSI(matrix.info,
                                              flowFrame,
                                              "GFP",
                                              "CD8",
                                              showScores = TRUE))
})