



context("Reaction classification")


test_that("self-self reactions are classified correctly",{

  mA <- "ABAAABAA"
  m1 <- "NMMMNNNN"

  result <- doReaction(c(mA,m1))
  conv <- doReaction(c(m1,mA))

  expect_equal(result$bprob, 0.458261, tolerance=1e-6)
  expect_true(result$product == "empty")

  # SELF-SELF Reactions
  # Test no product
  typ <- reaction_type(mA,m1,result,NA)
  expect_true(typ$type == "NonSelfNoProduct",info = sprintf("no bind: output is \"%s\"",typ$type))

  # Test no bind
  dummy <- result
  dummy$bprob <- 0
  dummy$mActive = "MACROMUTATION"
  typ <- reaction_type(mA,m1,dummy,NA)
  expect_true(typ$type == "Macromutation",info = sprintf("no bind: output is \"%s\"",typ$type))

  # Test self-self-replicator
  rep1 <- "$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYH"
  result <- runReactionFP(c(rep1,rep1))
  typ <- reaction_type(rep1,rep1,result,NA)
  expect_true(typ$type == "SelfSelfReplicator",info = sprintf("no bind: output is \"%s\"",typ$type))

})

test_that("non-self reactions are classified correctly",{

  mA <- "ABAAABAA"
  m1 <- "NMMMNNNN"

  result <- doReaction(c(mA,m1))
  conv <- doReaction(c(m1,mA))

  expect_equal(result$bprob, 0.458261, tolerance=1e-6)
  expect_true(result$product == "empty")

  # Test no bind
  # Simple way is to set bprob to zero
  dummy <- result
  dummy$bprob <- 0
  typ <- reaction_type(mA,m1,dummy,NA)
  expect_true(typ$type == "NonSelfNoProduct",info = sprintf("no bind: output is \"%s\"",typ$type))

  # Test replicator
  rep1 <- "$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYH"
  rep2 <- "$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYHPYH"
  result <- runReactionFP(c(rep1,rep2))
  typ <- reaction_type(rep1,rep2,result,NA)
  expect_true(typ$type == "NonSelfReplicator",info = sprintf("non-self replicator: output is \"%s\"",typ$type))

  # Test no product
  result <- runReactionFP(c(mA,m1))
  typ <- reaction_type(mA,m1,result,NA)
  expect_true(typ$type == "NonSelfNoProduct",info = sprintf("no bind: output is \"%s\"",typ$type))

  # Test different product (faking it)
  rep1 <- "$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYH"
  rep2 <- "$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYHPYH"
  result <- runReactionFP(c(rep1,rep2))
  result$product <- "DIFFERENT^PRODUCT"
  typ <- reaction_type(rep1,rep2,result,NA)
  expect_true(typ$type == "NonSelfdiffp",info = sprintf("non-self replicator: output is \"%s\"",typ$type))



})


