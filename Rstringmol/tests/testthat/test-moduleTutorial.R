

context("Rcpp modules load correctly")


test_that("doReaction fails gracefully",{
  result <- doReaction("")
  #expect_error(doReaction(""))
  expect_match(result$status,"bad number of input strings")
  rm(result)
})
