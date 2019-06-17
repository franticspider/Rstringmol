

rconf_parse <- function(data,pname,spregx="\\s+",argno=2){

  pos<-which(str_sub(data[,1],1,str_length(pname))==pname)
  if(length(pos)==1){
    ww <- strsplit(data[pos,1],split = spregx)

    return(ww[[1]][argno])
  }
  else{

    return(NA)
  }

}

#' Load the data from a stringmol .conf file
#' fn is the file name
#' verbose whether to run quietly or not
#' summarize whether to compress the data into a summary of each species
#' @export
rconf_params <- function(fn,verbose = F){

  cf <- read.table(fn,as.is = T, fill = T, sep="\n",comment.char="",stringsAsFactors = F,blank.lines.skip = F)

  smpar <-list()

  #gridx
  smpar$gridx <- as.numeric(rconf_parse(cf,"GRIDX"))
  #gridy
  smpar$gridy <- as.numeric(rconf_parse(cf,"GRIDY"))

  return(smpar)
}
