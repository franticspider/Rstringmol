

context("molecule-molecule interactions")


test_that("doReaction fails gracefully",{
  result <- doReaction("")
  expect_match(result$status,"bad number of input strings")

})



test_that("replicator replicates",{

  result <- doReaction(c("OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB","OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB"))
  expect_equal(result$bprob, 0.576381, tolerance=1e-6)
  expect_match(result$product,"OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB")

})
