

context("molecule-molecule interactions")


test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})



test_that("doReaction fails gracefully",{
  result <- doReaction("")
  expect_match(result$status,"bad number of input strings")

})



test_that("replicator replicates",{

  result <- doReaction(c("OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB","OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB"))

  expect_match(result$product,"OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB")

})
