

context("Smith-Waterman algorthms align correctly")


test_that("doSWAlign errors gracefully",{
  result <- doSWAlign("")
  expect_match(result$status,"bad number of input strings")
  rm(result)
})



test_that("Standard C++ SW works",{
  result <- doSWAlign(c("OYHOB","BLUBO"))
  expect_match(result$status,"Aligned")
  rm(result)
})
