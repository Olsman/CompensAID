testthat::test_that("CompensAID", {
  
  flowSet <- flowCore::read.flowSet(system.file("extdata", "68983.fcs", package = "CompensAID"))
  flowFrame <- flowCore::read.FCS(system.file("extdata", "68983.fcs", package = "CompensAID"))
  
  # Correct
  testthat::expect_output(CompensAID(flowFrame))
  testthat::expect_output(CompensAID(flowFrame, segment.value = 3, events.value = 70))
  
  # FlowSet not a FlowFrame
  testthat::expect_error(CompensAID(flowSet))
  
  # Segment.value not numerical
  testthat::expect_error(CompensAID(flowFrame, segment.value = "3", events.value = 70))
})