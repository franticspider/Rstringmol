
context("conf file data loads")


test_that("rconf data reads",{

  datafile <- system.file("extdata", "out1_880000.conf", package = "Rstringmol")
  data <- rconf_rdata(datafile)
  #expect_
  #expect_match(,"bad number of input strings")

})

