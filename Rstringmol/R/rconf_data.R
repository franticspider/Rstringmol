


rconf_agents<- function(fn){

  cf <- read.table(fn,as.is = T, fill = T, sep="\n",comment.char="",stringsAsFactors = F)

  agents <- cf[str_sub(cf[,1],1,5)=="AGENT",]
  agstr <- str_sub(agents,8,str_length(agents)-4)
  #ag_strlen <- str_length(agstr)


  return(agstr)

}



#' Load the data from a stringmol .conf file
#' fn is the file name
#' verbose whether to run quietly or not
#' summarize whether to compress the data into a summary of each species
#' @export
rconf_rdata <- function(fn,verbose = F,summarize=F){

  pdebug <- T

  cf <- read.table(fn,as.is = T, fill = T, sep="\n",comment.char="",stringsAsFactors = F,blank.lines.skip = F)

  #get the row number of the start of each reaction entry
  rpos <- which(str_sub(cf[,1],1,12)=="%%% REACTION")


  #create some data structures to hold the reaction data
  nr <- length(rpos)
  actno <- vector(length = nr)
  actseq <- vector(length = nr)
  actx <- vector(length = nr)
  acty <- vector(length = nr)
  pasno <- vector(length = nr)
  passeq <- vector(length = nr)
  passx <- vector(length = nr)
  passy <- vector(length = nr)

  if(verbose)message(sprintf("Found %d reactions", nr))
  count_bad_format = 0

  #build the data for each reaction - similar to the splist analysis
  if(nr>0){
  	if(verbose)message("Checking Reactions")
	  for(rr in 1:nr){
	    words <-strsplit(cf[rpos[rr]+2,1],split=" ")
	    actno[rr] <- as.numeric(words[[1]][3])
	    actseq[rr] <- words[[1]][4]

	    words <-strsplit(cf[rpos[rr]+3,1],split=" ")
	    actx[rr] <- as.numeric(words[[1]][19])
	    acty[rr] <- as.numeric(words[[1]][20])

	    foundpassive <- F




	    if(rpos[rr]+7 < nrow(cf)){
		    if(str_length(cf[rpos[rr]+7,1])>14){

		      words <-strsplit(cf[rpos[rr]+7,1],split=" ")

		      if(verbose)message(sprintf("Act: rr = %d Searching line %d:\n\t\t%s",rr,rpos[rr]+7, cf[rpos[rr]+7,1]))
		      ##if(is.na(as.integer(words[[1]][3])))message("BOOM")

		      if(words[[1]][1] == "###############"){

			##if(pdebug)message(sprintf("words[[1]][3] = %s",words[[1]][3]))
			##TODO: we get 'NAs introduced by coercion' when this is e.g. "$OYHOB"
			pasno[rr] <- suppressWarnings(  as.numeric(words[[1]][3]) )
			passeq[rr] <- words[[1]][4]

			words <-strsplit(cf[rpos[rr]+6,1],split=" ")
			passx[rr] <- as.numeric(words[[1]][5])
			passy[rr] <- as.numeric(words[[1]][6])

			foundpassive <- T
		      }
		    }
		}

	    if(!foundpassive){

	      count_bad_format <- count_bad_format + 1
	      if(verbose)#message(sprintf("Bad format number %d encountered, rpos[rr] is %d, rr is %d, no is %d ax,ay = %d,%d", count_bad_format, rpos[rr], rr, pasno[rr],actx[rr],acty[rr]))
		message(sprintf("Bad format number %d encountered",count_bad_format))

	      #foundpassive <- F
	      loff <- -1
	      while(!foundpassive){
		lno <- rpos[rr]+loff
		if(str_length(cf[lno,1])>14){
		  if(verbose)message(sprintf("rr = %d Searching line %d:\n\t\t%s",rr,lno, cf[lno,1]))

		  words <-strsplit(cf[lno,1],split=" ")
		  if(verbose)for(ww in 1:length(words[[1]]))message(sprintf("words[[1]][%d] is %s",ww,words[[1]][ww]))
		  if(words[[1]][1] == "###############"){
		    if(verbose)message("Found passive line...")
		    foundpassive <- T
		  }
		}
		loff <- loff-1
		if(loff < -10){
		  message(sprintf("Unable to find passive partner; returning from rconf_data"))
		  return(NA)
		}
	      }

	      #TODO handle bad format better where passive mol is above rpos!! (and find out how this happens...)
	      #if(is.na(pasno[rr])){
	      pasno[rr] <- as.numeric(words[[1]][3])
	      passeq[rr] <- words[[1]][4]
	      #go up one more line to get gridx and gridy
	      words <-strsplit(cf[lno-1,1],split=" ")
	      passx[rr] <- as.numeric(words[[1]][5])
	      passy[rr] <- as.numeric(words[[1]][6])

	      if(verbose)message(sprintf("Found passeq:%s, passx:%s, passy:%s",passeq[rr],passx[rr],passy[rr]))
	    }
	    #message(sprintf("Active  molecule %s seq %s",actno[rr],actseq[rr]))
	    #message(sprintf("Passive molecule %s seq %s",actno[rr],actseq[rr]))
	  }
  }

  eachreaction<-data.frame(actno,actseq,actx,acty,pasno,passeq,passx,passy,stringsAsFactors = F)

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
