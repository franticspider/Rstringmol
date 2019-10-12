



context("Reaction classification")


test_that("reactions are classified correctly",{

  mA <- "ABAAABAA"
  m1 <- "NMMMNNNN"

  result <- doReaction(c(mA,m1))
  conv <- doReaction(c(m1,mA))

  expect_equal(result$bprob, 0.458261, tolerance=1e-6)
  expect_true(result$product == "empty")

  # Test no product
  typ <- reaction_type(result,conv)
  expect_true(typ$product == "NoProduct",info = sprintf("no bind: typ$product is \"%s\"",typ$product))

  # Test no bind
  dummy <- result
  dummy$bprob <- 0
  typ <- reaction_type(result,conv)
  expect_true(typ$product == "Macromutation",info = sprintf("no bind: typ$product is \"%s\"",typ$product))

  # Test self-replicator
  dummy <- result
  dummy$mActive <- mA
  dummy$mPassive <- mA
  dummy$product <- mA
  conv <- dummy
  typ <- reaction_type(result,conv)
  expect_true(typ$product == "SelfSelfReplicator",info = sprintf("no bind: typ$product is \"%s\"",typ$product))


  rm(result)
})
