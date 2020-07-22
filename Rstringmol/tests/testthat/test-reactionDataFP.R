

context("molecule-molecule interactions; FP method")


test_that("runReactionFP fails gracefully",{
  result <- runReactionFP("")
  expect_match(result$status,"bad number of input strings")

})



test_that("FP replicator replicates",{

  result <- runReactionFP(c("OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB","OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB"))
  expect_equal(result$bprob, 0.576381, tolerance=1e-6)
  expect_true(result$product == "OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB")

  result <- runReactionFP(c("$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYH","$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYH"))
  expect_equal(result$bprob, 0.695410, tolerance=1e-6)
  expect_true(result$product == "$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYH")

  rm(result)
})



test_that("FP inexact execution is inexact",{

  result <- runReactionFP(c("OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BBBBBBBBBB^B>C$=?>$$BBBBBBBBB%}OOOONNOOOO",
                            "OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BBBBBBBBBB^B>C$=?>$$BBBBBBBBB%}OOOONNOOOO"))
  expect_equal(result$bprob, 0.576381, tolerance=1e-6)

  rm(result)
})

test_that("Copy Operator counters count",{
  result <- runReactionFP(c("OOOO$BLUB>C^B=========OYHOXXX","BBBBYY"))
  #TODO: need to count calls to '=' that don't fire at all (both R and W pointers are off-string)
  expect_equal(result$ccopy, 6)
  expect_equal(result$cover, 3)


})

test_that("FPnon-deterministic bind is detected",{


  result <- runReactionFP(c("$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYH","$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYH"))
  expect_false(result$deterministicBind)

  rm(result)
})


test_that("FP inexact execution paths can be detected",{

  rr <- runReactionFP(c("AAA$BBBBBBXXXOOPPO","NNN"))

  expect_false(rr$deterministicExec)

})

test_that("inexact execution paths can be counted",{

  results = data.frame()
  for(i in 1:200){
    rr <- runReactionFP(c("OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BBBBBBBBBB^B>C$=?>$$BBBBBBBBB%}OOOONNOOOO",
                       "OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BBBBBBBBBB^B>C$=?>$$BBBBBBBBB%}OOOONNOOOO"))
    results[i,1] <- rr$product
  }

  outputs <- unique(results)

  expect_equal(nrow(outputs), 3)

  #rm(results)
  #rm(outputs)
})


test_that("TOGGLE '^' toggles and COPY '=' copies",{

  result <- runReactionFP(c("AAA^B=","NNN"))

  #NB - using 'expect_true' here because 'expect_match' doesn't handle '^' vs '\^' properly
  expect_true(result$mActive == "NAA^B=")
  expect_true(result$mPassive == "NNN")
  expect_equal(result$count,6)

  rm(result)
})


test_that("count of products works",{
  ms <- "OOOOOOOOBBBBBBBBE$BLU^B>C=====%=====%OYH"

  rs <- runReactionFP(c(ms,ms))

  expect_equal(rs$nprod,2)

})

