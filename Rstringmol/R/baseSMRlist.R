




#' Create a list of fields for a stringmol reaction
#'
#' @export
#' @examples
#' p <- baseSMRlist()
baseSMRlist <- function(){


  result <- list()

  result$bprob <- 0
  result$count <- 0
  result$mActive <- " "
  result$mPassive <- " "
  result$product <- " "
  result$deterministicBind <- T
  result$deterministicExec <- T

  return(result)

}
