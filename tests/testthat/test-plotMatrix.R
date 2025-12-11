testthat::test_that("Plot Matrix Errors", {
  
  flowFrame <- flowCore::read.FCS(system.file("extdata", "68983.fcs", package = "CompensAID"))
  compensAID.res <- CompensAID(flowFrame)
  matrix.SSI <- compensAID.res$matrix
  matrix.info <- compensAID.res$matrixInfo
  
  # Correct
  testthat::expect_no_error(PlotMatrix(compensAID.res))
  
  # Partial CompensAID output used
  testthat::expect_error(PlotMatrix(matrix.SSI))
  
  # Partial CompensAID output used
  testthat::expect_error(PlotMatrix(matrix.info))
})