




#' run a Stringmol reaction and get results via a file.
#' @export
runReactionFP <- function(reactants,verbose = F){

  fn <- tempfile(fileext = ".csv")

  doReactionFP(reactants,fn)

  if(file.exists(fn)){
    data <- read.table(fn,sep = ",",stringsAsFactors = F)

    result <- list()
    #c <- data
    #c <- t(data)
    #c <- as.data.frame(c)
    #colnames(c) <- c[1,]
    #c = c[2,]

    result$bprob <- as.numeric(data[data[,1]=="bprob",2])
    result$count <- as.numeric(data[data[,1]=="count",2])


    result$mActive <- data[data[,1]=="mActive",2]
    result$mPassive <- data[data[,1]=="mPassive",2]
    result$product <- data[data[,1]=="product",2]


    dtb <- data[data[,1]=="deterministicBind",2]
    if(dtb=="TRUE")
      result$deterministicBind <- TRUE
    else
      result$deterministicBind <- FALSE


    dtb <- data[data[,1]=="deterministicExec",2]
    if(dtb=="TRUE")
      result$deterministicExec <- TRUE
    else
      result$deterministicExec <- FALSE



    #delete the file
    file.remove(fn)

    return (result)
  }
  else
    return (NA)



}
