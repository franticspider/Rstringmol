


#' Determine the type of stringmol reacion
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
    if(pro$product == "empty"){
      pro$rtype = "NonSelfNoProduct"
    }
    else{
      if(pro$product == pro$mPassive){
        pro$rtype = "NonSelfReplicator"
      }
      else{
        pro$rtype = "NonSelfDifferentProduct"
      }
    }

  }

  return(pro)

}
