



context("Reaction classification")


test_that("reactions are classified correctly",{

  mA <- "ABAAABAA"
  m1 <- "NMMMNNNN"

  result <- doReaction(c(mA,m1))

  expect_equal(result$bprob, 0.458261, tolerance=1e-6)
  expect_true(result$product == "empty")

  typ <- reaction_type(result$mActive,result$mPassive)

})
