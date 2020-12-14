



context("Reaction classification")






# Helper function to run pairwise props
pwc <- function(ain,pin,d.r){
  fdata <- data.frame(actseq = ain, passeq = pin, stringsAsFactors = F)
  ur<-pairwise.properties(fdata,dummy.result = d.r)#,dummy.result=dummy.result)
  return(ur)
}


# Helper function to run network props
nwc <- function(ain,pin,d.r,rd.r){
  fdata <- data.frame(actseq = ain, passeq = pin, stringsAsFactors = F)
  ur<-pairwise.properties(fdata,dummy.result = d.r)#,dummy.result=dummy.result)

  nr<-network.properties(ur,rd.r)

  return(nr)
}


test_that("Self-preserving and self-modifying are correctly identified",{

  # NOTE! We might not have to use genuine reactions to test all this - so we can put sequences in for hypotheticals

  ain <- "AAAA"
  pin <- "NNNN"
  d.r <- runReactionFP(c(ain,pin))

  ur <- pwc(ain,pin,d.r)
  expect_false(ur$pp_SelfMod)

  # Change output 1 - pp_SelfMod should now be true
  d.r$mActive <- "AAAAX"
  ur <- pwc(ain,pin,d.r)
  expect_true(ur$pp_SelfMod)

  # Change output 2 - pp_SelfMod should now be true
  d.r$mActive <- ain
  d.r$mPassive <- "NNNNX"
  ur <- pwc(ain,pin,d.r)
  expect_true(ur$pp_SelfMod)

})


test_that("No-product and New-product properties are correctly identified",{

  # NOTE! We might not have to use genuine reactions to test all this - so we can put sequences in for hypotheticals

  ain <- "AAAA"
  pin <- "NNNN"
  d.r <- runReactionFP(c(ain,pin))

  ur <- pwc(ain,pin,d.r)

  expect_true(ur$pp_NoProduct)


  # Generate a copy of mol1:
  d.r$product <- ain
  ur <- pwc(ain,pin,d.r)
  expect_false(ur$pp_NewProduct)

  # Generate a copy of mol2:
  d.r$product <- pin
  ur <- pwc(ain,pin,d.r)
  expect_false(ur$pp_NewProduct)

  # generate a new product
  d.r$product <- "PRODUCT"
  ur <- pwc(ain,pin,d.r)
  expect_true(ur$pp_NewProduct)

})



test_that("Replicator properties are correctly identified",{

  # NOTE! We might not have to use genuine reactions to test all this - so we can put sequences in for hypotheticals

  ain <- "AAAA"
  pin <- "NNNN"
  d.r <- runReactionFP(c(ain,pin))

  ur <- pwc(ain,pin,d.r)

  expect_true(ur$pp_NoProduct)

  # Generate a copy of mol2:
  d.r$product <- pin
  ur <- pwc(ain,pin,d.r)
  expect_true(ur$pp_Repl2)

  # modify the outputs so that there's no increase in number - repl is F
  d.r$product <- pin
  d.r$mPassive <- "XXXX"
  ur <- pwc(ain,pin,d.r)
  expect_false(ur$pp_Repl2)

  # outputs are swapped inputs, but a new copy is produced:
  d.r$product <- pin
  d.r$mActive <- pin
  d.r$mPassive <- ain
  ur <- pwc(ain,pin,d.r)
  expect_true(ur$pp_Repl2)

  # no new products but both outputs are mol2: repl2 is T
  d.r$product <- "empty"
  d.r$mActive <- pin
  d.r$mPassive <- pin
  ur <- pwc(ain,pin,d.r)
  expect_true(ur$pp_Repl2)


  # Generate a copy of mol1:
  d.r$product <- ain
  d.r$mActive <- ain
  d.r$mPassive <- pin
  ur <- pwc(ain,pin,d.r)
  expect_true(ur$pp_Repl1)

  # modify the outputs so that there's no increase in number - repl is F
  d.r$product <- ain
  d.r$mActive <- "XXXX"
  d.r$mPassive <- pin
  ur <- pwc(ain,pin,d.r)
  expect_false(ur$pp_Repl1)

  # outputs are swapped inputs, but a new copy is produced:
  d.r$product <- ain
  d.r$mActive <- pin
  d.r$mPassive <- ain
  ur <- pwc(ain,pin,d.r)
  expect_true(ur$pp_Repl1)

  # no new products but both outputs are mol2: repl2 is T
  d.r$product <- "empty"
  d.r$mActive <- ain
  d.r$mPassive <- ain
  ur <- pwc(ain,pin,d.r)
  expect_true(ur$pp_Repl1)

})


test_that("Auto-Replicator properties are correctly identified",{
  ain <- "AAAA"
  d.r <- runReactionFP(c(ain,ain))

  # 'standard' version
  d.r$product <- ain
  d.r$mActive <- ain
  d.r$mPassive <- ain
  ur <- pwc(ain,ain,d.r)
  expect_true(ur$pp_SelfReplicator)

  # fails:
  d.r$product  <- ain
  d.r$mActive  <- "XXXX"
  d.r$mPassive <- ain
  ur <- pwc(ain,ain,d.r)
  expect_false(ur$pp_SelfReplicator)

  # fails:
  d.r$product  <- "XXXX"
  d.r$mActive  <- ain
  d.r$mPassive <- ain
  ur <- pwc(ain,ain,d.r)
  expect_false(ur$pp_SelfReplicator)

  # fails:
  d.r$product  <- ain
  d.r$mActive  <- ain
  d.r$mPassive <- "XXXX"
  ur <- pwc(ain,ain,d.r)
  expect_false(ur$pp_SelfReplicator)

})

# TODO: Test jumper reactions here to keep order with the paper


test_that("Parasitic properties are correctly identified",{

  # NOTE: this one is a *genuine* parasitic reaction, not a 'dummy'
  ain <- "$=?>$EBUB^B$=?>$AVBO%}HOB"
  pin <- "$=?>$B"
  d.r <- runReactionFP(c(ain,pin))
  rd.r <- runReactionFP(c(pin,ain))
  rd.p <- pwc(pin,ain,rd.r)

  ur <- nwc(ain,pin,d.r,rd.p)

  expect_true(ur$np_Parasite2)
  expect_true(ur$np_Parasite)
})



test_that("Mutual replicator properties are correctly identified",{

  # NOTE: this one is a *genuine* mutual reaction, not a 'dummy'
  ain <- "$=?>$EBUB^B$=?>$AVBO%}HOB"
  pin <- "$=?>$EBUB^B$=?>$AVBOX%}HOB"
  d.r <- runReactionFP(c(ain,pin))
  rd.r <- runReactionFP(c(pin,ain))
  rd.p <- pwc(pin,ain,rd.r)

  ur <- nwc(ain,pin,d.r,rd.p)

  expect_true(ur$np_MutualRepl)
})


test_that("Hypercycle properties are correctly identified",{

  # NOTE: this one is a *genuine* hypercycle reaction, not a 'dummy'
  ain <- "UUZSNSHJBLRWVR$BLUBO^B>C$=?>$DBLUP%}OYH"
  pin <- "HHMYURFUVJKJRJLUWR$BLUCO^B>C$=?>$DBLU%}OYH"
  d.r <- runReactionFP(c(ain,pin))
  rd.r <- runReactionFP(c(pin,ain))
  rd.p <- pwc(pin,ain,rd.r)

  ur <- nwc(ain,pin,d.r,rd.p)

  expect_true(ur$np_MutualRepl)
  expect_true(ur$np_Hypercycle)
})


