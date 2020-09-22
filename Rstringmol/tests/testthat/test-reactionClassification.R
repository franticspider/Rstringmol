



context("Reaction classification")






# Helper function to run pairwise props
pwc <- function(ain,pin,d.r){
  fdata <- data.frame(actseq = ain, passeq = pin, stringsAsFactors = F)
  ur<-pairwise.properties(fdata,dummy.result = d.r)#,dummy.result=dummy.result)
  return(ur)
}


test_that("Self-preserving and self-modifying are correctly identified",{

  # NOTE! We might not have to use genuine reactions to test all this - so we can put sequences in for hypotheticals

  ain <- "AAAA"
  pin <- "NNNN"
  d.r <- runReactionFP(c(ain,pin))

  ur <- pwc(ain,pin,d.r)

  expect_true(ur$pp_NoProduct)
  expect_false(ur$pp_SelfMod)

})



# OBSOLETE! use pairwise.properties instead
test_that("self-self reactions are classified correctly",{

  skip("obsolete method to be deleted when pairwise.properties and network.properties are running")

  mA <- "ABAAABAA"
  m1 <- "NMMMNNNN"
  np1 <- "$=?>G^AQC"
  rep1 <- "$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYH"


  # SELF-SELF Reactions
  # Test no bind
  result <- doReaction(c(mA,mA))
  typ <- reaction_type(mA,mA)
  expect_true(typ$type == "SelfSelfNoProduct",info = sprintf("no bind: output is \"%s\"",typ$type))

  # Test replicator
  result <- doReaction(c(rep1,rep1))
  typ <- reaction_type(rep1,rep1)
  expect_true(typ$type == "SelfSelfReplicator",info = sprintf("replictor: output is \"%s\"",typ$type))

  # Test no product
  result <- doReaction(c(np1,np1))
  typ <- reaction_type(np1,np1,result)
  expect_true(typ$type == "SelfSelfNoProduct",info = sprintf("no product: output is \"%s\"",typ$type))

  # Test replicator different product
  result <- doReaction(c(rep1,rep1))
  result$product = "DIFFERENTPRODUCT"
  typ <- reaction_type(rep1,rep1,result)
  expect_true(typ$type == "SelfSelfDifferentProduct",info = sprintf("different product: output is \"%s\"",typ$type))

  # Test Macromutation
  result <- doReaction(c(rep1,rep1))
  result$mActive = "MACROMUTATION"
  typ <- reaction_type(rep1,rep1,result)
  expect_true(typ$type == "Macromutation",info = sprintf("macromutation: output is \"%s\"",typ$type))

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
  typ <- reaction_type(mA,m1,dummy)
  expect_true(typ$type == "NonSelfNoProduct",info = sprintf("no bind: output is \"%s\"",typ$type))

  # Test replicator
  rep1 <- "$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYH"
  rep2 <- "$=?>G^AQD$=?>G^BQC$=?>E$BLUO%}PYH"
  result <- runReactionFP(c(rep1,rep2))
  conv <- runReactionFP(c(rep2,rep1))
  typ <- reaction_type(rep1,rep2,result,conv)
  expect_true(typ$type == "NonSelfReplicator",info = sprintf("non-self replicator: output is \"%s\"",typ$type))


  # Test parasite
  rep1 <- "$=?>$EBUB^B$=?>$AVBO%}HOB"
  rep2 <- "$=?>$B"
  result <- runReactionFP(c(rep1,rep2))
  conv <- runReactionFP(c(rep2,rep1))
  typ <- reaction_type(rep1,rep2,result,conv)
  expect_true(typ$type == "Parasite",info = sprintf("non-self replicator: output is \"%s\"",typ$type))



  # Test no product
  result <- runReactionFP(c(mA,m1))
  typ <- reaction_type(mA,m1,result)
  expect_true(typ$type == "NonSelfNoProduct",info = sprintf("no bind: output is \"%s\"",typ$type))

  # Test different product (faking it)
  rep1 <- "$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYH"
  rep2 <- "$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}PYHPYH"
  result <- runReactionFP(c(rep1,rep2))
  result$product <- "DIFFERENT^PRODUCT"
  typ <- reaction_type(rep1,rep2,result)
  expect_true(typ$type == "NonSelfDifferentProduct",info = sprintf("non-self replicator: output is \"%s\", product is \"%s\"",typ$type,result$product))



})


