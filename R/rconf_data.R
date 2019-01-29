


rconf_agents<- function(fn){

  cf <- read.table(fn,as.is = T, fill = T, sep="\n",comment.char="",stringsAsFactors = F)

  agents <- cf[str_sub(cf[,1],1,5)=="AGENT",]
  agstr <- str_sub(agents,7,str_length(agents)-4)
  #ag_strlen <- str_length(agstr)

  return(agstr)

}



rconf_rdata <- function(fn,verbose = F,summarize=T){

  pdebug <- T

  cf <- read.table(fn,as.is = T, fill = T, sep="\n",comment.char="",stringsAsFactors = F)

  #get the row number of the start of each reaction entry
  rpos <- which(str_sub(cf[,1],1,12)=="%%% REACTION")


  #create some data structures to hold the reaction data
  nr <- length(rpos)
  actno <- vector(length = nr)
  actseq <- vector(length = nr)
  pasno <- vector(length = nr)
  passeq <- vector(length = nr)

  if(verbose)message(sprintf("Found %d reactions", nr))
  count_bad_format = 0

  #build the data for each reaction - similar to the splist analysis
  for(rr in 1:nr){
    words <-strsplit(cf[rpos[rr]+2,1],split=" ")
    actno[rr] <- as.numeric(words[[1]][3])
    actseq[rr] <- words[[1]][4]
    words <-strsplit(cf[rpos[rr]+6,1],split=" ")

    ##if(is.na(as.integer(words[[1]][3])))message("BOOM")
    ##if(pdebug)message(sprintf("words[[1]][3] = %s",words[[1]][3]))
    ##TODO: we get 'NAs introduced by coercion' when this is e.g. "$OYHOB"
    pasno[rr] <- as.numeric(words[[1]][3])
    passeq[rr] <- words[[1]][4]
    #TODO handle bad format better where passive mol is above rpos!!
    if(is.na(pasno[rr])){
      words <-strsplit(cf[rpos[rr]-1,1],split=" ")
      pasno[rr] <- as.numeric(words[[1]][3])
      passeq[rr] <- words[[1]][4]
      count_bad_format <- count_bad_format + 1
      if(verbose)message(sprintf("Bad format number %d encountered", count_bad_format))
    }
    #message(sprintf("Active  molecule %s seq %s",actno[rr],actseq[rr]))
    #message(sprintf("Passive molecule %s seq %s",actno[rr],actseq[rr]))
  }

  eachreaction<-data.frame(actno,actseq,pasno,passeq)

  if(summarize){
    reactions <- unique(eachreaction)
    reactions$type = "<unknown>"

    for(rr in 1:nrow(reactions)){
      ca <- eachreaction[eachreaction$actno == reactions$actno[rr] & eachreaction$pasno == reactions$pasno[rr],]
      reactions$count[rr] = nrow(ca)
    }
    reactions <- reactions[order(reactions$count,decreasing = T),]

    #return(eachreaction)
    return(reactions)
  }
  else{
    return(eachreaction)
  }
}
