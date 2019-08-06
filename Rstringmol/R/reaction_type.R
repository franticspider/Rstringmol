

# List of reaction types - not sure this is the best place for them...
rtypes = c("SelfSelfNoProduct",
           "SelfSelfReplicator",
           "SelfSelfDifferentProduct",
           "Parasite",
           "NonSelfNoProduct",
           "NonSelfReplicator",
           "NonSelfDifferentProduct")

rcols = c("grey",
          "blue",
          "orange",
          "red",
          "darkgrey",
          "deepskyblue",
          "yellow")


#' Determine the type of stringmol reaction
#'
#' @param act the active stringmol
#' @param pas the passive stringmol
#' @keywords "Mass spectrum",
#' @export
#' @examples
#' p <- reaction_type("AAA","NNN")
reaction_type <- function(act,pas){
  pro <- doReaction(c(act,pas))

  if(pro$mActive == pro$mPassive){
    if(pro$product == "empty"){
      pro$rtype = "SelfSelfNoProduct"
    }
    else{
      if(pro$product == pro$mPassive){
        pro$rtype = "SelfSelfReplicator"
      }
      else{
        pro$rtype = "SelfSelfDifferentProduct"
      }
    }
  }
  else{
    conv <- doReaction(c(pas,act))

    if(pro$product == "empty"){
      pro$rtype = "NonSelfNoProduct"
    }
    else{
      if(pro$product == pro$mPassive){
        #TODO: We need to know if the assignment to active/passive is a coin toss!
        if(conv$product == "empty")
          pro$rtype = "Parasite"
        else
          pro$rtype = "NonSelfReplicator"
      }
      else{
        pro$rtype = "NonSelfDifferentProduct"
      }
    }

  }

  return(pro)

}





#' Determine the type of stringmol reacion
#'
#' @param act the active stringmol
#' @param pas the passive stringmol
#' @keywords "Mass spectrum",
#' @export
#' @examples
#' p <- reaction_type("AAA","NNN")
reaction_typeFP <- function(act,pas){
  #pro <- doReaction(c(act,pas))
  pro <- runReactionFP(c(act,pas))

  if(pro$mActive == pro$mPassive){
    if(pro$product == "empty"){
      pro$rtype = "SelfSelfNoProduct"
    }
    else{
      if(pro$product == pro$mPassive){
        pro$rtype = "SelfSelfReplicator"
      }
      else{
        pro$rtype = "SelfSelfDifferentProduct"
      }
    }
  }
  else{
    conv <- doReaction(c(pas,act))

    if(pro$product == "empty"){
      pro$rtype = "NonSelfNoProduct"
    }
    else{
      if(pro$product == pro$mPassive){
        #TODO: We need to know if the assignment to active/passive is a coin toss!
        if(conv$product == "empty")
          pro$rtype = "Parasite"
        else
          pro$rtype = "NonSelfReplicator"
      }
      else{
        pro$rtype = "NonSelfDifferentProduct"
      }
    }

  }
  return(pro)
}
