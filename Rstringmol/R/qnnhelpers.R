



#' Get the unbound molecular species counts from a stringmol config file
#'
#' @param froot the root of the filename
#' @param tt the timestep of the config file
#' @param verbose whether to write messages whilst processing
#' @return a dataframe holding the following fields:
#'   \item{t}{the timestep}
#'   \item{n}{The number of individuals}
#'   \item{seq}{The species sequence}
unboundpopfromconf <- function(froot,tt,verbose=F){
  
  fn <- sprintf("%s%05d.conf",froot,tt)
  
  if(!file.exists(fn)){
    message(sprintf("file %s doesn't exist- CRASH!",fn))
  }
    
  if(verbose)message(sprintf("Loading file %s",fn))
  raw <- read.table(fn,stringsAsFactors = F,fill=T,sep=",")
  lx <- data.frame(raw[str_detect(str_sub(raw$V1,1,8),"SPECIES"),]
                   ,stringsAsFactors = F)
  
  #get the unbound agent count - create object 'ubx'
  ubx <- data.frame(raw[str_detect(str_sub(raw$V1,1,5),"AGENT"),]
                    ,stringsAsFactors = F)
  
  if(nrow(ubx)<1){
    lx = data.frame(t=integer(),n=integer(),seq=character(),stringsAsFactors = F)
    return(lx)
  }
  
  ubx$seq = ""
  ubx$n = 0
  for(rr in 1:nrow(ubx)){
    split <- strsplit(ubx[rr,1], " +")[[1]]
    ubx$seq[rr]=split[2]
  }
  ubx <- ubx[,2:3]
  
  lx$t = tt
  lx$n = 0
  lx$seq = ""
  for(rr in 1:nrow(lx)){
    split <- strsplit(lx[rr,1], " +")[[1]]
    lx$seq[rr] <- split[3]
    #use ubx to count the unbound mols: 
    lx$n[rr] <- nrow(ubx[ubx$seq == lx$seq[rr],])
  }
  
  lx <- lx[,2:4]
  return(lx)
}





#' Get all the species in a run
#'
#' @param rundata data object of 
#' @param froot the root of the filename
#' @return a list of dataframes holding the following fields:
#'   \item{t}{the timestep}
#'   \item{n}{The number of individuals}
#'   \item{seq}{The species sequence}
#' @export
getallspp <- function(rundata,froot,maxtime=2000000,timestep=20000){
  
  pdl <- list()
  idx = 1
  for(tt in seq(0,maxtime,timestep)){
    #message(tt)
    ubdata <- unboundpopfromconf(froot,tt) 
    
    if(nrow(ubdata)<1){
      ubdata$unb <- integer()
      ubdata$act <- integer()
      ubdata$pas <- integer()
    }
    else{
      ubdata$unb <- ubdata$n
      ubdata$act <- 0
      ubdata$pas <- 0
      
      # NB! nothing is bound at t=0!
      if(idx>1){
        rdata <- rundata[[(idx-1)]]
        
        for(rr in 1:nrow(ubdata)){
          ad <- rdata[rdata$actseq==ubdata$seq[rr],]
          aobs <- sum(ad$nobs)
          pd <- rdata[rdata$passeq==ubdata$seq[rr],]
          pobs <- sum(pd$nobs)
          
          ubdata$act[rr] <- aobs
          ubdata$pas[rr] <- pobs
          ubdata$n[rr] <- ubdata$n[rr] + aobs + pobs
        }
      }
    }
    pdl[[idx]]<-ubdata
    idx = idx+1
  }
  
  return(pdl)
  
}



#' QNN analysis of a popdy object from the QNN package
#'
#' @param x the popdy object
#' @param scale boolean, whether to divide the score by the current population size, default FALSE
#' @return qnn score for the popdy object
#' @export
qa2 <- function (x, scale = FALSE) 
{
  population.per.time <- apply(x@tracks, 1, sum)
  species.proportions <- sweep(x@tracks, 1, population.per.time, 
                               FUN = "/")
  expected.species.proportions.per.time <- rbind(rep(NA, nspecies(x)), 
                                                 species.proportions[1:(nrow(species.proportions) - 1), 
                                                 ])
  diff.expected <- (species.proportions - expected.species.proportions.per.time)
  diff.expected[diff.expected < 0] <- 0
  diff.expected[1, ] <- 0
  diff.expected <- diff.expected^2
  if (identical(scale, TRUE)) 
    diff.expected <- diff.expected * population.per.time
  output <- as(x, "popdy")
  output@tracks <- as(diff.expected, "dgCMatrix")
  return(output)
}



mkrawsp <- function(pd,t){
  pd <- pd[,1:3]
  pd$start <- t
  pd$end <-t
  return(pd)
}



makespptable <- function(pdl){
  
  spp <- mkrawsp(pdl[[1]],0)
  
  # Make the raw data block
  for(ll in 2:length(pdl)){
    tdata <- pdl[[ll]]
    if(nrow(tdata)>1){
      tdata <- mkrawsp(tdata,tdata$t[1])
      spp <- rbind(spp,tdata)
    }
  }
  
  # now set the count, start and end of each identical spp to the same number - then we can do it 'unique'
  seqlist <- unique(spp$seq)
  
  for(ss in 1:length(seqlist)){
    seq <- seqlist[ss]
    spp$start[spp$seq == seq]<- min(spp$start[spp$seq == seq])
    spp$end[spp$seq == seq]<- max(spp$end[spp$seq == seq])
    spp$n[spp$seq == seq]<- sum(spp$n[spp$seq == seq])
  }
  
  spp <- spp[,2:5]
  nonunique <- nrow(spp)
  spp <- unique(spp)
  
  spp = spp[with(spp, order(start,end,-n)),]
  
  # finally, give a species ID column
  spp$ID = ""
  for(ii in 1:nrow(spp)){
    spp$ID[ii] = sprintf("sp%d",ii)
  }
  
  # Re-order the columns so the read nicely
  spp = data.frame(ID = spp$ID,
                   start = spp$start,
                   end = spp$end,
                   n = spp$n,
                   seq = spp$seq,
                   stringsAsFactors = F)
}



getIDfromseq <- function(spp,seq){
  spids <- spp$ID[spp$seq==seq]  
  
  if(length(spids)!=1){
    message(sprintf("ERROR finding species ID for sequence %s\nreturning NA",seq))
    return(NA)
  }
  else{
    return(spids[1])
  }
}

makepopdy <- function(pdl,spp,rno,verbose = T){
  for(ll in 1:length(pdl)){
    if(verbose)message(sprintf("Calcing row %d",ll))
    pdata <- pdl[[ll]]
    if(nrow(pdata)<1){
      pdata$ID <- character()
    }
    else{
      pdata$ID = ""
      for(rr in 1:nrow(pdata)){
        pdata$ID[rr] = getIDfromseq(spp,pdata$seq[rr])
      }
    }
    
    
    pdata <- pdata[,c(1,7,2,4,5,6,3)]
    write.csv(pdata,file = sprintf("runs/run%d/sppcounts%06d.csv",rno,pdata$t[1]),quote = F,row.names = F)
    
    if(ll==1){
      popdy <- pdata[,1:3]
    }
    else{
      popdy <- rbind(popdy,pdata[,1:3])
    }
  }
  
  return(popdy)
}




popdyplot2 <- function (x, xlab = "Time", ylab = "Count", plot.zero = FALSE, 
                        col.fun = rainbow, colours = NA, minpop = 1, logrey = F, 
                        xlim = NA, ylim = NA, smooth = 0, verbose = F) 
{
  if (identical(xlim, NA)) 
    xlim = range(x@times)
  
  if (identical(ylim,NA))
    ylim = c(0, max(x@tracks))

  if (identical(colours, NA)) 
    colours <- col.fun(ncol(x@tracks))
  plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)
  if (logrey) {
    for (j in 1:ncol(x@tracks)) {
      if (max(x@tracks[, j] < minpop)) {
        lines(x = x@times, y = x@tracks[, j], col = "grey")
      }
    }
  }
  for (j in 1:ncol(x@tracks)) {
    if (max(x@tracks[, j] >= minpop)) {
      if (verbose) 
        message(sprintf("Plotting run %d with colour %s", 
                        j, colours[j]))
      if (smooth > 0) {
        lx <- x@times
        ly <- x@tracks[, j]
        lo <- loess(ly ~ lx, span = smooth)
        lines(lx, lo$fitted, col = colours[j], lwd = 2)
      }
      else {
        lines(x = x@times, y = x@tracks[, j], col = colours[j])
      }
    }
  }
}



