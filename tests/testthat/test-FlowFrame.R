testthat::test_that("Accept valid FlowFrame input", {
  ff <- flowCore::read.FCS(system.file("extdata", "68983.fcs", package = "CompensAID"))
  
  # Correct
  testthat::expect_silent(checkmate::assert(methods::is(ff, "flowFrame"), "Object is not a flowFrame."))
})

testthat::test_that("Reject non-FlowFrame input", {
  
  # dataFrame not a flowFrame
  testthat::expect_error(checkmate::assert(methods::is(data.frame(a = 1), "flowFrame"), "Object is not a flowFrame."))
  
  # No data
  testthat::expect_error(checkmate::assert(methods::is(NULL, "flowFrame"), "Object is not a flowFrame."))
  
  # list not a dataFrame
  testthat::expect_error(checkmate::assert(methods::is(list(), "flowFrame"), "Object is not a flowFrame."))
  
  # Expression matrix not a flowFrame
  testthat::expect_error(checkmate::assert(methods::is(ff@exprs, "flowFrame"), "Object is not a flowFrame."))
})
