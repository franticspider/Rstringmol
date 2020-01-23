
##TODO: Rename this function - it simply loads the file...!
##      - see rconf_stats for better analysis

#' Load a splist file and format the data
#'
#' @param fn the file name of the species list
#' @param tmin the start time
#' @param verbose messages or not
#' @return species list
#' @export
splist_stats <- function(fn,tmin=0,verbose=F){

  #read the splist file
  xx <- read.table(fn,stringsAsFactors = F,fill=T,sep=",")

  if(ncol(xx)==5){
    #TODO: Handle this more elegantly
    message(sprintf("RESTART DETECTED - check file %s",fn))

    colnames(xx) <- c("spp","act","pass","n3","seq")
    xx$spp <- as.numeric(xx$spp)
    xx$act <- as.numeric(xx$act)
    xx$pass <- as.numeric(xx$pass)
    xx$n3 <- as.numeric(xx$n3)

    return(xx)
  }
  else{

    #get the origin molecule(s)
    origin <- xx[xx[,2]==-1,]

    #message(sprintf("%d origin molecules found",nrow(origin)))

    #get the descendent reactions
    xx <- xx[xx[,2]>-1,]

    colnames(xx) <- c("spp","act","pass","float","obsn","obst1","n3","seq")

    xx$spp <- as.numeric(xx$spp)
    xx$act <- as.numeric(xx$act)
    xx$pass <- as.numeric(xx$pass)
    xx$obsn <- as.numeric(xx$obsn)
    xx$obst1 <- as.numeric(xx$obst1)
    xx$n3 <- as.numeric(xx$n3)

    if(verbose)message(sprintf("%d species found",length(unique(xx$spp))))

    if(verbose)message("Commonest reaction:")
    #com1 <- xx[xx$obsn == max(xx$obsn) && xx$obst1 > tmin,]

    recent <- xx[xx$obst1 > tmin,]
    com1 <- recent[recent$obsn == max(recent$obsn),]

    if(verbose)message(sprintf("%d commonest reactions, seen %d times",nrow(com1),com1$obsn[1]))
    for (cc in 1:nrow(com1)){
      #Get the active molecule
      ac <- xx[xx$spp == com1$act[1],]
      if(verbose)message(sprintf("First seen at t = %d",com1$obst1[1]))
      if(verbose)message(sprintf("Active  molecule is %s",ac$seq[1]))
      pa <- xx[xx$spp == com1$act[1],]
      if(verbose)message(sprintf("Passive molecule is %s",pa$seq[1]))
      if(verbose)message(sprintf("Product is          %s",com1$seq[1]))

    }

  }
  return(xx)
}
