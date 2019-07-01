

context("molecule-molecule interactions; FP method")


test_that("runReactionFP fails gracefully",{
  result <- runReactionFP(c("AAA","NNN"),tempfn)
  expect_match(result$status,"bad number of input strings")
  rm(result)
})



test_that("FP replicator replicates",{

  result <- runReactionFP(c("OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB","OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB"),tempfn)
  expect_equal(result$bprob, 0.576381, tolerance=1e-6)
  expect_true(result$product == "OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB")

  #Here's a different replicator that doesn't appear to work:

  result <- runReactionFP(c("$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYH","$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYH"))
  expect_equal(result$bprob, 0.695410, tolerance=1e-6)
  expect_true(result$product == "$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYH")

  rm(result)
})



test_that("FP inexact execution is inexact",{


  result <- runReactionFP(c("OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BBBBBBBBBB^B>C$=?>$$BBBBBBBBB%}OOOONNOOOO","OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BBBBBBBBBB^B>C$=?>$$BBBBBBBBB%}OOOONNOOOO"),tempfn)
  expect_equal(result$bprob, 0.576381, tolerance=1e-6)

  rm(result)
})



test_that("FPnon-deterministic bind is detected",{


  result <- runReactionFP(c("$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYH","$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYH"),tempfn)
  expect_false(result$deterministicBind)

  rm(result)
})


test_that("FP inexact execution paths can be detected",{

  rr <- runReactionFP(c("AAA$BBBBBBXXXOOPPO","NNN"),tempfn)

  expect_false(rr$deterministicExec)

  rm(result)
})

test_that("inexact execution paths can be counted",{

  results = data.frame()
  for(i in 1:200){
    rr <- runReactionFP(c("OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BBBBBBBBBB^B>C$=?>$$BBBBBBBBB%}OOOONNOOOO",
                       "OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BBBBBBBBBB^B>C$=?>$$BBBBBBBBB%}OOOONNOOOO"),tempfn)
    results[i,1] <- rr$product
  }

  outputs <- unique(results)

  expect_equal(nrow(outputs), 3)

  rm(result)
  rm(outputs)
})


test_that("TOGGLE '^' toggles and COPY '=' copies",{

  result <- runReactionFP(c("AAA^B=","NNN"))

  #NB - using 'expect_true' here because 'expect_match' doesn't handle '^' vs '\^' properly
  expect_true(result$mActive[1] == "NAA^B=")
  expect_true(result$mPassive[1] == "NNN")
  expect_equal(result$count[1],6)

  rm(result)
})
