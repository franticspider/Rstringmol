
context("conf file data loads")
library(stringr)


test_that("rconf data reads",{

  datafile <- system.file("extdata", "out1_880000.conf", package = "Rstringmol")
  data <- rconf_rdata(datafile)
  #expect_
  #expect_match(,"bad number of input strings")

})

test_that("rconf parameters read",{

  datafile <- system.file("extdata", "out1_880000.conf", package = "Rstringmol")
  params <- rconf_params(datafile)

  expect_equal(params$gridx,125)
  expect_equal(params$gridy,100)



})
