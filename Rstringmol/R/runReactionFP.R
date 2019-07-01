




#' run a Stringmol reaction and get results via a file.
#' @export
runReactionFP <- function(reactants,verbose = F){

  fn <- tempfile(fileext = ".csv")

  doReactionFP(reactants,fn)

  if(file.exists(fn)){
    data <- read.table(tempfn,sep = ",",stringsAsFactors = F)
    #c <- as.data.frame(c)
    c <- t(data)
    colnames(c) <- c[1,]
    c = c[2,]


    #delete the file
    file.remove(fn)

    return (c)
  }
  else
    return (NA)



}
