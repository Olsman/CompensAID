test_that("Accept valid FlowFrame input", {
  ff <- flowCore::read.FCS(system.file("extdata", "68983.fcs", package = "CompensAID"))
  expect_silent(checkmate::assert(methods::is(ff, "flowFrame"), "Object is not a flowFrame."))
})

test_that("Reject non-FlowFrame input", {
  expect_error(checkmate::assert(methods::is(data.frame(a = 1), "flowFrame"), "Object is not a flowFrame."))
  expect_error(checkmate::assert(methods::is(NULL, "flowFrame"), "Object is not a flowFrame."))
  expect_error(checkmate::assert(methods::is(list(), "flowFrame"), "Object is not a flowFrame."))
})
