


splist_stats <- function(fn){

  #read the splist file
  xx <- read.table(fn,stringsAsFactors = F,fill=T,sep=",")

  #get the origin molecule(s)
  origin <- xx[xx[,2]==-1,]

  message(sprintf("%d origin molecules found",nrow(origin)))

  #get the descendent reactions
  xx <- xx[xx[,2]>-1,]

  colnames(xx) <- c("spp","act","pass","float","obsn","obst1","n3","seq")

  xx$spp <- as.numeric(xx$spp)
  xx$act <- as.numeric(xx$act)
  xx$pass <- as.numeric(xx$act)
  xx$obsn <- as.numeric(xx$obsn)
  xx$obst1 <- as.numeric(xx$obst1)
  xx$n3 <- as.numeric(xx$n3)

  message(sprintf("%d species found",length(unique(xx$spp))))

  message("Commonest reaction:")
  com1 <- xx[xx$obsn == max(xx$obsn),]

  message(sprintf("%d commonest reactions, seen ",nrow(com1)))
  for (cc in 1:nrow(com1)){


  }

  return(xx)
}
