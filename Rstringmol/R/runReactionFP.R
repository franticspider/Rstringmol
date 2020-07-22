




#' run a Stringmol reaction and get results via a file.
#' @export
runReactionFP <- function(reactants,verbose = F,remove=T){


  result <- list()

  if(length(reactants)!=2){
    message(sprintf("Wrong number of arguments to runReactionFP. Expecting 2, got %d",length(reactants)))
    result$status<-"bad number of input strings"
    return(result)
  }




  fn <- tempfile(fileext = ".csv")

  doReactionFP(reactants,fn,verbose)

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

    result$ccopy <- as.numeric(data[data[1]=="ccopy",2])
    result$cmove <- as.numeric(data[data[1]=="cmove",2])
    result$cover <- as.numeric(data[data[1]=="cover",2])
    result$ctogg <- as.numeric(data[data[1]=="ctogg",2])


    result$nprod <- as.numeric(data[data[1]=="nprod",2])


    #delete the file
    if(remove)
      file.remove(fn)
    else
      message(sprintf("Created file %s",fn))

    return (result)
  }
  else
    return (NA)



}
