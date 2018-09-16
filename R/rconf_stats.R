

####Function 'arguments'
#fn <- sprintf("%sout1_200000.conf",froot)
#verbose <- T

rconf_stats <- function(fn,spfn,verbose=F){

  spldata <-  splist_stats(spfn)

  ####Function 'body'

  cf <- read.table(fn,as.is = T, fill = T, sep="\n",comment.char="",stringsAsFactors = F)

  #get the row number of the start of each reaction entry
  rpos <- which(str_sub(cf[,1],1,12)=="%%% REACTION")

  nr <- length(rpos)
  actno <- vector(length = nr)
  actseq <- vector(length = nr)
  pasno <- vector(length = nr)
  passeq <- vector(length = nr)

  #build the data for each reaction - similar to the splist analysis
  for(rr in 1:nr){
    words <-strsplit(cf[rpos[rr]+2,1],split=" ")
    actno[rr] <- as.numeric(words[[1]][3])
    actseq[rr] <- words[[1]][4]
    words <-strsplit(cf[rpos[rr]+6,1],split=" ")
    pasno[rr] <- as.numeric(words[[1]][3])
    passeq[rr] <- words[[1]][4]
    #message(sprintf("Active  molecule %s seq %s",actno[rr],actseq[rr]))
    #message(sprintf("Passive molecule %s seq %s",actno[rr],actseq[rr]))
  }

  eachreaction<-data.frame(actno,actseq,pasno,passeq)

  reactions <- unique(eachreaction)
  reactions$type = "<unknown>"

  for(rr in 1:nrow(reactions)){
    ca <- eachreaction[eachreaction$actno == reactions$actno[rr] & eachreaction$pasno == reactions$pasno[rr],]
    reactions$count[rr] = nrow(ca)
  }
  reactions <- reactions[order(reactions$count,decreasing = T),]

  #Get the product molecule from the splist:
  for(rr in 1:nrow(reactions)){

    if(is.na(reactions$actno[rr]) || is.na(reactions$pasno[rr])){

      if(verbose)message(sprintf("Reaction %d has NAs  (NA)",rr))
      reactions$type[rr]= "NA"
      next
    }


    dd <- spldata[(spldata$act == reactions$actno[rr]) & (spldata$pass == reactions$pasno[rr]), ]
    if(nrow(dd)==0){
      if(verbose)message(sprintf("Reaction %d has no product (NP)",rr))
      reactions$type[rr]= "NP"
    }
    if(reactions$actno[rr] == reactions$pasno[rr]){
      if(nrow(dd)==1){
        if(dd$spp == reactions$actno[rr]){
          if(verbose)message(sprintf("Reaction %d is an exact self:self replicator (ER)*",rr))
          reactions$type[rr]= "ER"
        }
        else{
          if(verbose)message(sprintf("Self:self reaction with different product (DP)*",rr))
          reactions$type[rr]= "DP"
        }
      }
      if(nrow(dd)>1){
        if(verbose)message(sprintf("Reaction %d is a pathological replicator (PR)",rr))
        reactions$type[rr]= "PR"
      }
    }
    else{
      if(nrow(dd)==1){
        if(dd$spp == reactions$pasno[rr]){
          if(verbose)message(sprintf("Reaction %d is parsitic (PA)*",rr))
          reactions$type[rr]= "PA"
        }
        else{
          if(verbose)message(sprintf("none-self reaction with different product (DP)*",rr))
          reactions$type[rr]= "DN"
        }
      }
      if(nrow(dd)>1){
        if(verbose)message(sprintf("Reaction %d is a pathological non-self replicator (PN)",rr))
        reactions$type[rr]= "PN"
      }
    }
  }

  return (reactions)

}
