

context("molecule-molecule interactions")


test_that("doReaction fails gracefully",{
  result <- doReaction("")
  expect_match(result$status,"bad number of input strings")

})



test_that("replicator replicates",{

  result <- doReaction(c("OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB","OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB"))
  expect_equal(result$bprob, 0.576381, tolerance=1e-6)
  expect_true(result$product == "OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB")

})



test_that("inexact execution is inexact",{


  result <- doReaction(c("OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BBBBBBBBBB^B>C$=?>$$BBBBBBBBB%}OOOONNOOOO","OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BBBBBBBBBB^B>C$=?>$$BBBBBBBBB%}OOOONNOOOO"))
  expect_equal(result$bprob, 0.576381, tolerance=1e-6)

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
  expect_true(result$m0 == "NAA^B=")
  expect_true(result$m1 == "NNN")
  expect_equal(result$count,6)

})

