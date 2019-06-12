

context("molecule-molecule interactions")


test_that("doReaction fails gracefully",{
  result <- doReaction("")
  expect_match(result$status,"bad number of input strings")

})



test_that("replicator replicates",{

  result <- doReaction(c("OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB","OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB"))
  expect_equal(result$bprob, 0.576381, tolerance=1e-6)
  expect_true(result$product == "OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB")

  #Here's a different replicator that doesn't appear to work:

  result <- doReaction(c("$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYH","$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYH"))
  expect_equal(result$bprob, 0.695410, tolerance=1e-6)
  expect_true(result$product == "$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYH")

})



test_that("inexact execution is inexact",{


  result <- doReaction(c("OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BBBBBBBBBB^B>C$=?>$$BBBBBBBBB%}OOOONNOOOO","OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BBBBBBBBBB^B>C$=?>$$BBBBBBBBB%}OOOONNOOOO"))
  expect_equal(result$bprob, 0.576381, tolerance=1e-6)

})



test_that("non-deterministic bind is detected",{


  result <- doReaction(c("$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYH","$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYH"))
  expect_false(result$deterministicBind)

})


test_that("inexact execution paths can be detected",{

  rr <- doReaction(c("AAA$BBBBBBXXXOOPPO","NNN"))

  expect_false(rr$deterministicExec)

})

test_that("inexact execution paths can be counted",{

  results = data.frame()
  for(i in 1:200){
    rr <- doReaction(c("OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BBBBBBBBBB^B>C$=?>$$BBBBBBBBB%}OOOONNOOOO",
                       "OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BBBBBBBBBB^B>C$=?>$$BBBBBBBBB%}OOOONNOOOO"))
    results[i,1] <- rr$product
  }

  outputs <- unique(results)

  expect_equal(nrow(outputs), 3)

})


test_that("TOGGLE '^' toggles and COPY '=' copies",{

  result <- doReaction(c("AAA^B=","NNN"))

  #NB - using 'expect_true' here because 'expect_match' doesn't handle '^' vs '\^' properly
  expect_true(result$mActive == "NAA^B=")
  expect_true(result$mPassive == "NNN")
  expect_equal(result$count,6)

})





