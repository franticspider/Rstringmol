
#' Calculate Simpson's Diversity Index D1
#'
#' @param pdata popdy data file loaded from a stringmol output file
#' @param time the time for which the index will be calculated
#' @return The Simpson's Diversity Index D1
#' @export
sdi <- function(pdata, time){

  slice <- pdata[pdata$time == time,]

  tot <- sum(slice$count)
  p<- slice$count/tot
  return(sum(p*p))


}
