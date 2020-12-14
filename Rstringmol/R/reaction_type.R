



# TODO: do we need these?
# OBSOLETE!!!
# List of reaction types - not sure this is the best place for them...
rtypes = c(
  "SelfSelfNoProduct",
  "SelfSelfReplicator",
  "SelfSelfDifferentProduct",
  "Parasite",
  "NonSelfNoProduct",
  "NonSelfReplicator",
  "NonSelfDifferentProduct",
  "Macromutation"
)

# TODO: do we need these?
rcols = c("grey",
          "blue",
          "orange",
          "red",
          "darkgrey",
          "deepskyblue",
          "yellow",
          "white")


print_mols <- function(type,
                       act,
                       pas,
                       pro,
                       conv,
                       actout = NA,
                       pasout = NA) {
  message(sprintf("%s:", type))
  message(sprintf("ACT          = %s", act))
  message(sprintf("PAS          = %s", pas))
  message(sprintf("PRO          = %s", pro))
  message(sprintf("CONV         = %s", conv))
  if (!is.na(actout)) {
    message(sprintf("ACTOUT       = %s", actout))
    message(sprintf("PASOUT       = %s", pasout))
  }
  message("")
}










############################################
## DETECTION REACTION PROPERTIES



#' Identify the network-level properties in a reaction
#'
#' @param fdata data frame containing actseq and passeq fields
#' @param dummy.result an optional result of runReactionFP passed in (for testing)
#' @export
pairwise.properties <- function(fdata,dummy.result=NULL){

  #create a data frame of unique reactions for this file:
  ur <- unique(data.frame(actseq = fdata$actseq, passeq = fdata$passeq, stringsAsFactors = F))

  #ur$type <- "unknown"
  ur$product <- ""

  #These are the properties agreed with Susan
  ur$pp_ActiveMod <- F
  ur$pp_PassiveMod <- F
  ur$pp_SelfMod <- F

  ur$pp_NoProduct <-F
  ur$pp_NewProduct <-F

  ur$pp_biolRep <-F
  ur$pp_SelfReplicator <-F
  ur$pp_Repl1 <- F # formerly ActiveReplicator <-F
  ur$pp_Repl2 <- F # formerly PassiveReplicator <- F

  ur$pp_Jumper <-F

  ur$actlen <- str_length(ur$actseq)
  ur$paslen <- str_length(ur$passeq)
  ur$outputA <- ""
  ur$outputP <- ""

  ur$ccopy <-0
  ur$cmove <-0
  ur$cover <-0
  ur$ctogg <-0

  ur$bprob <-0
  ur$nsteps <- 0

  ur$nprod <- 0

  ur$oneway <-F

  for (rr in 1:nrow(ur)){
    #################################################
    # Reaction properties

    #rt <- reaction_type(ur$actseq[rr],ur$passeq[rr])
    # Run the reaction(s):

    if(is.null(dummy.result))
      result <- runReactionFP(c(ur$actseq[rr],ur$passeq[rr]))
    else{
      message("Calculating properties based on external result")
      result <- dummy.result
    }

    ur$oneway[rr] <- result$deterministicBind

    if(ur$actseq[rr] != result$mActive)
      ur$pp_ActiveMod[rr] <- T

    if(ur$passeq[rr] != result$mPassive)
      ur$pp_PassiveMod[rr] <- T

    ur$pp_SelfMod[rr] <- ur$pp_PassiveMod[rr] | ur$pp_ActiveMod[rr]

    if(result$product == "empty")
      ur$pp_NoProduct[rr] <- T
    else{
      if(result$product != ur$actseq[rr] & result$product != ur$passeq[rr])
        ur$pp_NewProduct[rr] <- T
    }

    ur$ccopy[rr] <- result$ccopy
    ur$cmove[rr] <- result$cmove
    ur$cover[rr] <- result$cover
    ur$ctogg[rr] <- result$ctogg

    ur$bprob[rr] <- result$bprob
    ur$nsteps[rr] <- result$count

    ur$nprod[rr] <- result$nprod
    ur$outputA[rr] <- result$mActive
    ur$outputP[rr] <- result$mPassive

    #####################################################
    #Replicator types

    # Count the inputs copies of active and passive as described in the tech report
    actin<-1
    pasin<-1
    if(ur$actseq[rr]==ur$passeq[rr]){
      actin <- actin+1
      pasin <- pasin+1
    }

    act <- 0
    if(ur$actseq[rr]==result$mActive)  act <- act + 1
    if(ur$actseq[rr]==result$mPassive) act <- act + 1
    if(ur$actseq[rr]==result$product)  act <- act + 1
    if(act>actin)ur$pp_Repl1[rr] <-T

    pct <- 0
    if(ur$passeq[rr]==result$mActive)  pct <- pct + 1
    if(ur$passeq[rr]==result$mPassive) pct <- pct + 1
    if(ur$passeq[rr]==result$product)  pct <- pct + 1
    if(pct>pasin)ur$pp_Repl2[rr] <-T


    # The technical repl property is now "networked": see tech report
    # Here we are using pp_biolRep as a naive term - NOT TO BE USED TO BUILD OTHER PROPERTIES
    ur$pp_biolRep[rr] <- ur$pp_Repl1[rr] | ur$pp_Repl2[rr]

    # We can still spot self-replicators here:
    if((ur$pp_Repl1[rr] | ur$pp_Repl2[rr]) & (ur$actseq[rr]==ur$passeq[rr])){
      ur$pp_SelfReplicator[rr] <- T
    }


    #####################################################
    #Jumper types


    Jumper1 <- F
    Jumper2 <- F
    if(result$mActive != ur$actseq[rr] & ur$actseq[rr] == result$product)
      Jumper1 <- T
    if(result$mPassive != ur$passeq[rr] & ur$passeq[rr] == result$product)
      Jumper2 <- T
    if(Jumper1 | Jumper2)
      ur$pp_Jumper[rr] <- T


    #####################################################
    # Other useful data from the reaction:

    ur$product[rr] <- result$product
    ur$dexec[rr] <- result$deterministicExec
    ur$dbind[rr] <- result$deterministicBind
    ur$nobs[rr] <- nrow(fdata[fdata$actseq == ur$actseq[rr]
                              & fdata$passeq == ur$passeq[rr],])
  }

  return(ur)
}


#' Identify the network-level properties in a reaction
#'
#' @param fn the file name
#' @export
network.properties <- function(nwp,dummy.reverse.pwp = NULL){

  nwp$np_Parasite2 <- F
  nwp$np_Parasite1 <- F
  nwp$np_Parasite <- F

  nwp$np_MutualRepl <- F
  nwp$np_Hypercycle <- F

  for(rr in 1:nrow(nwp)){
     #Run the 'flip' reaction
     if(nwp$actseq[rr] == nwp$passeq[rr]){
       #Can't be a parasite - no need to change
       #Can't be a hypercycle - no need to change
       if(nwp$pp_SelfReplicator[rr]){
         nwp$np_MutualRepl[rr] <- T
       }
     }else{
        if(is.null(dummy.reverse.pwp)){
          flipdata <- data.frame(actseq=nwp$passeq[rr],passeq=nwp$actseq[rr],stringsAsFactors = F)
          flip <- pairwise.properties(flipdata)
        }else{
          flip <- dummy.reverse.pwp
        }

       #PARASITE
       if(nwp$pp_Repl2[rr] & ! (nwp$pp_Repl1[rr] | flip$pp_Repl2)){
           nwp$np_Parasite2[rr] <- T
       }

       if(nwp$pp_Repl1[rr] & ! (nwp$pp_Repl2[rr] | flip$pp_Repl1)){
         nwp$np_Parasite1[rr] <- T
       }

       nwp$np_Parasite[rr] <- nwp$np_Parasite2[rr] | nwp$np_Parasite1[rr]



       #MUTUAL

       if(( nwp$pp_Repl2[rr] | flip$pp_Repl1) & (nwp$pp_Repl1[rr] | flip$pp_Repl2)){
         nwp$np_MutualRepl[rr]<-T
       }



       #HYPERCYCLE
       if(nwp$np_MutualRepl[rr]){
         ina <- data.frame(actseq=nwp$actseq[rr],passeq=nwp$actseq[rr],stringsAsFactors = F)
         pwpa <- pairwise.properties(ina)
         inp <- data.frame(actseq=nwp$passeq[rr],passeq=nwp$passeq[rr],stringsAsFactors = F)
         pwpp <- pairwise.properties(inp)

         if(nwp$np_MutualRepl[rr] & !pwpa$pp_SelfReplicator & !pwpp$pp_SelfReplicator)
           nwp$np_Hypercycle[rr] <- T

       }
     }
  }
  return(nwp)
}




############################################
## PARSING BIG DATA FILES


#' Create a data structure from a *.conf* file for passing into the standard *pairwise.properties* function
#'
#' @param fn the file name
#' @export
pairwise.properties.from.conf <- function(fn){

  fdata <- rconf_rdata(fn)
  fdata$actseq <- as.character(fdata$actseq)
  fdata$passeq <- as.character(fdata$passeq)

  ur <- pairwise.properties(fdata)

  return(ur)
}



#' Obtain the reaction properties from an active and passive sequence
#'
#' @param actseq the active sequence
#' @param passeq the passive sequence
#' @export
propfromseqs <- function(actseq,passeq){

  fd <- data.frame(actseq=actseq,passeq=passeq,stringsAsFactors = F)
  pwp <- pairwise.properties(fd)
  ng <- network.properties(pwp)
  return(ng)
}


#' Obtain the reaction properties of all the reactions in a stringmol run
#'
#' @param froot the pathway to the data files
#' @param from the first timestep
#' @param to the last timestep
#' @param step the timestep increment
#' @param outfn the name of an R data object to save it to
#' @param verbose verbose output
#' @export
runproplist <- function(froot,from=20000,to=2000000,step=20000,outfn=NULL,verbose = F){

  rundata <- list()
  dp <- 1
  for(rr in seq(from,to,step)){

    #testinfn <- sprintf("~/Desktop/paulien/smsp/1705smsp/out3/out1_%d.conf",rr)
    testinfn <- sprintf("%s%d.conf",froot,rr)
    if(verbose)message(sprintf("%d: file is %s",rr,testinfn))

    ug <- pairwise.properties.from.conf(testinfn)

    ng <- network.properties(ug)

    rundata[[dp]] <- ng

    dp <- dp + 1

  }

  if(!is.null(outfn)){
    save(rundata,file=outfn)
  }

  return(rundata)

}


#############################################
## DEPRECIATED FUNCTIONS


## DEPRECIATED FUNCTION
#' Determine the type of stringmol reaction
#'
#' @param act the active stringmol
#' @param pas the passive stringmol
#' @keywords "Mass spectrum",
#' @export
#' @examples
#' p <- reaction_type("AAA","NNN")
reaction_typeFP_old <- function(act,
                                pas,
                                verbose = F,
                                detail = F) {


  .Deprecated(msg = "'reaction_typeFP_old' will be removed in the next version")

  #pro <- doReaction(c(act,pas))
  pro <- runReactionFP(c(act, pas))
  conv <- runReactionFP(c(pas, act))

  catalytic <- T

  # Check for changed parents!
  if (pro$mActive != act) {
    catalytic <- F
  }

  if (pro$mPassive != pas) {
    catalytic <- F
  }


  pro$rtype = "Macromutation"

  if (catalytic) {
    if (pro$mActive == pro$mPassive) {
      if (pro$product == "empty") {
        pro$rtype = "SelfSelfNoProduct"
      }
      else{
        if (pro$product == pro$mPassive) {
          pro$rtype = "SelfSelfReplicator"
          #print_mols(pro$rtype,act,pas,pro$product,conv$product)
        }
        else{
          pro$rtype = "SelfSelfDifferentProduct"
        }
      }
    }
    else{
      #conv <- doReaction(c(pas,act))
      #conv <- runReactionFP(c(pas,act))

      if (pro$product == "empty") {
        if (conv$product == "empty") {
          pro$rtype = "NonSelfNoProduct"
        }
        else{
          if (conv$product == conv$mPassive)
            pro$rtype = "Parasite"
          else{
            pro$rtype = "NonSelfDifferentProduct"
          }
        }
      }
      else{
        if (pro$product == pro$mPassive) {
          #TODO: We need to know if the assignment to active/passive is a coin toss!
          if (conv$product == "empty")
            pro$rtype = "Parasite"
          else{
            if (conv$product == pro$mPassive) {
              #TODO: extend this if detail = T
              pro$rtype = "Parasite"
            } else{
              if (conv$product == conv$mPassive) {
                #TODO: extend this if detail = T
                pro$rtype = "NonSelfReplicator"
                #print_mols(pro$rtype,act,pas,pro$product,conv$product)
              } else{
                pro$rtype = "NonSelfDifferentProduct"
                #print_mols(pro$rtype,act,pas,pro$product,conv$product)
              }
            }
          }
        }
        else{
          pro$rtype = "NonSelfDifferentProduct"
        }
      }
    }
  }
  else{
    pro$rtype = "Macromutation"
    if (verbose)
      print_mols(pro$rtype,
                 act,
                 pas,
                 pro$product,
                 conv$product,
                 pro$mActive,
                 pro$mPassive)
  }
  return(pro)
}



#' Determine the type of stringmol reaction
#'
#' @param act the active stringmol
#' @param pas the passive stringmol
#' @keywords "Mass spectrum",
#' @export
#' @examples
#' p <- reaction_type("AAA","NNN")
reaction_type <- function(act,pas,
                          result=NULL,
                          conv=NULL,
                          verbose = F,
                          detail = F) {


  .Deprecated(msg = "'reaction_typ' will be removed in the next version")

  if(is.null(result))
    result <- runReactionFP(c(act,pas))


  result$type = "undefined"

  #TODO: it should be possible to remove the "if(act == pas)" statement from the outer loop and nest it more deeply when looking for parasites
  if(act == pas){
    if((result$bprob < 0.0000001) || (result$product == "empty"))
      result$type <- "SelfSelfNoProduct"
    else{
      if(result$mActive != act || result$mPassive !=pas)
        result$type <- "Macromutation"
      else{
        if(result$product == pas)
          result$type <- "SelfSelfReplicator"
        if(result$product != pas && result$product != act)
          result$type <- "SelfSelfDifferentProduct"
      }
    }
  }else{
    #generate the converse reaction if it isn't supplied:
    if(is.null(conv))
      conv <- runReactionFP(c(pas,act))

    if((result$bprob < 0.0000001) || (result$product == "empty"))
      result$type <- "NonSelfNoProduct"
    else{
      if(result$mActive != act || result$mPassive !=pas)
        result$type <- "Macromutation"
      else{
        if(result$product == pas){
          if(conv$product == act)
            result$type <- "NonSelfReplicator"
          else
            result$type <- "Parasite"
        }
        if(result$product != pas && result$product != act)
          result$type <- "NonSelfDifferentProduct"
      }
    }
  }

  # for debugging:
  #pro$product = "BBBBBBB"
  return(result)
}


