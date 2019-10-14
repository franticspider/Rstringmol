




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
                          result,
                          conv,
                          verbose = F,
                          detail = F) {


  result$type = "undefined"

  if(act == pas){
    if((result$bprob < 0.0000001) || (result$product == "empty"))
      result$type <- "SelfSelfNoProduct"
    if(result$product == pas)
      result$type <- "SelfSelfReplicator"
  }else{
    if((result$bprob < 0.0000001) || (result$product == "empty"))
      result$type <- "NonSelfNoProduct"
    if(result$product == pas)
      result$type <- "NonSelfReplicator"
    if((result$mActive != act ) && (result$mActive != pas))
      result$type <- "Macromutation"
    if((result$mPassive != act ) && (result$mPassive != pas))
      result$type <- "Macromutation"
  }

  # for debugging:
  #pro$product = "BBBBBBB"
  return(result)
}
