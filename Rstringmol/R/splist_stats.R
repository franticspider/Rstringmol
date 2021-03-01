
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
    #This happens after restart
    if(verbose)message(sprintf("RESTART DETECTED - check file %s",fn))

    colnames(xx) <- c("spp","act","pass","n3","seq")
    xx$spp <- as.numeric(xx$spp)
    xx$act <- as.numeric(xx$act)
    xx$pass <- as.numeric(xx$pass)
    xx$n3 <- as.numeric(xx$n3)

    yy=data.frame(spp=xx$spp,act=xx$act,pass=xx$pass,float=rep(0,nrow(xx)),
                  obsn=xx$n3,obst1=rep(-1,nrow(xx)),n3=rep(0,nrow(xx)),
                  seq=xx$seq,restart=rep(T,nrow(xx)))

    return(yy)
  }
  else{

    #get the origin molecule(s)
    origin <- xx[xx[,2]==-1,]
    colnames(origin) <- c("spp","act","pass","n3","seq")
    origin=data.frame(spp=origin$spp,
                      act=origin$act,
                      pass=origin$pass,
                      float=rep(0,nrow(origin)),
                      obsn=origin$n3,
                      obst1=rep(-1,nrow(origin)),
                      n3=rep(0,nrow(origin)),
                      seq=origin$seq,stringsAsFactors = F)

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

    xx<-rbind(xx,origin)

    xx$restart=F

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




#' Make the master ancestry data object
#'
#' @param froot path to the "splist" data files output from a stringmol run
#' @return a master ancestry data object
#' @export
makemanc <- function(froot){

  manc <- list()
  ll=1
  for(tt in seq(2000000,20000,-20000)){
    df <- splist_stats(fn=sprintf("%ssplist%d.dat",froot,tt),verbose = F)
    df$rectime <- tt
    manc[[ll]] <- df
    ll <- ll + 1
  }

  require(data.table)
  manc <- rbindlist(manc)
  #unfortunately rbindlist converts strings to factors, so we have to fix that:
  manc <- data.frame(lapply(manc, as.character), stringsAsFactors=FALSE)
  manc$spp <- as.numeric(manc$spp)
  manc$act <- as.numeric(manc$act)
  manc$pass <- as.numeric(manc$pass)
  manc$obsn <- as.numeric(manc$obsn)
  manc$obst1 <- as.numeric(manc$obst1)
  manc$rectime <- as.numeric(manc$rectime)
  manc$restart <- as.logical(manc$restart)

  return(manc)
}







#' get the molecules from the top 10 reactions at t=2m:
#'
#' @param rundata the reaction data object for this run
#' @param rdmanc the ancestry data object for this run
#' @param verbose verbose output
#' @return a data frame listing the details of the molecules in the top 10 reactions
#' @export
makeanchead <- function(rundata,rdmanc,verbose = F){
  r2m <- rundata[[2000000/20000]]
  r2m <- r2m[order(r2m$nobs,decreasing = T),]

  seqs <- unique(c(r2m$actseq[1:10],r2m$passeq[1:10]))

  rdhead <- data.frame(idx=NA,seq=seqs,t=2000000,birtht=0,actparent=NA,pasparent=NA,stringsAsFactors = F)
  rdhead$gen <- 0
  rdhead$trace <- T

  for(pp in 1:nrow(rdhead)){
    time = 2000000
    idx<-0
    found<-F
    while(!found){

      anc<-rdmanc[rdmanc$rectime == time,]
      parents <- anc[anc$seq == rdhead$seq[pp],]

      if(nrow(parents)>1){
        if(max(parents$obst1)>0){
          found <- T
          if(verbose)message(sprintf("Found parents of seq %s at time %d: %d and %d in splist %d",rdhead$seq[pp], parents$obst1[1], parents$act[1], parents$pass[1],time))

          rdhead$idx[pp]<- parents$spp[1]
          rdhead$birtht[pp]<-parents$obst1[1]
          rdhead$actparent[pp] <- parents$act[1]
          rdhead$pasparent[pp] <- parents$pass[1]

          break
        }
      }

      time = time - 20000
      if(time<1)break
    }
  }
  return(rdhead)

}




#' get the molecules from the top 10 reactions at t=2m:
#'
#' @param phyl the phylogeny
#' @param colbylen whether to colour the nodes by the length of the molecule
#' @param xlim the xlim for the plot
#' @return a data frame listing the details of the molecules in the top 10 reactions
#' @export
plotphyl <- function(phyl,colbylen=T,xlim=NULL){

  if(is.null(xlim))
    xlim=range(c(phyl$idx,phyl$actparent,phyl$pasparent))

  plot(x=phyl$actparent,y=phyl$birtht,
       ylim = c(min(phyl$birtht),2000000),
       xlim=xlim,
       pch=19,col="red",cex=1.5,xlab="Species number",ylab="time")


  if(colbylen){
    #rect(xleft = min(c(phyl$idx,phyl$actparent,phyl$pasparent)),xright = max(c(phyl$idx,phyl$actparent,phyl$pasparent)),ybottom = min(phyl$birtht),ytop=2000000,col="black")

    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black")
    #par(bg = "black")
  }

  points(x=phyl$actparent,y=phyl$birtht,pch=19,col="red")
  points(x=phyl$pasparent,y=phyl$birtht,pch=19,col="blue")

  segments(x0=phyl$idx,x1=phyl$actparent,y0=phyl$birtht,lwd=2,col="red")
  segments(x0=phyl$idx,x1=phyl$pasparent,y0=phyl$birtht,col="blue")

  if(colbylen){
    phyl$len <- str_length(phyl$seq)
    lr <- max(phyl$len)

    #cols<-heat.colors(n = 65)
    #FROM setupSM.cpp in strimgmol library
    cols <- matlab.like(70)

    for(pp in 1:nrow(phyl)){
      if(phyl$len[pp]<70)
        cpp <- cols[phyl$len[pp]]
      else
        cpp = "white"
      segments(x0=phyl$idx[pp],y0=phyl$birtht[pp],y1=phyl$t[pp],col=cpp)
      points(x=phyl$idx[pp],y=phyl$t[pp],      pch=19,col=cpp)
      points(x=phyl$idx[pp],y=phyl$birtht[pp], pch=19,col=cpp)
    }
  }
  else{
    points(x=phyl$idx,y=phyl$t,pch=19,col="green")
    points(x=phyl$idx,y=phyl$birtht,pch=19,col="green")
    segments(x0=phyl$idx,y0=phyl$birtht,y1=phyl$t,col="green")
  }

}





#TODO: Check that this works generally
getidxdata <- function(idx,data){
  spent <- unique(data[data$spp == idx,])

  rtimes <- unique(spent$rectime)

  spent <- spent[spent$rectime == min(rtimes),]
  return(spent[1,])

}




anctable <- function(head,time,data,ngen=3,verbose=F){

  head$gen <- 0

  for(gg in 1:ngen){
    gidx <- head$idx

    pidx <- unique(c(head$actparent,head$pasparent))
    if(verbose)message(sprintf("Found %d unique parents for generation %d",length(pidx),gg))
    #Now build the generation data:
    for(ii in 1:length(pidx)){

      idd<-getidxdata(pidx[ii],data)
      if(verbose)message(sprintf("%d",ii))
      if(verbose)message(sprintf("Recording details for %0.0f: %s",idd$spp,idd$seq))

      # idx, seq, t, birtht, actparent, pasparent
      pd <- data.frame(stringsAsFactors = F,
                       idx=pidx[ii],
                       seq=idd$seq,
                       t=max(unique(data$rectime[data$spp == pidx[ii]])),#idd$rectime,
                       birtht=idd$obst1,
                       actparent=idd$act,
                       pasparent=idd$pass
      )
      if(ii==1)
        newgen<-pd
      else
        newgen<-rbind(newgen,pd)
    }
  }
  return (newgen)
}




nextgen<- function(bigspeciestable,childspp,verbose=F,doplot=F){

  next.gen <- anctable(childspp[childspp$gen == max(childspp$gen) & childspp$trace,],2000000,bigspeciestable,ngen=1)
  genno <- max(childspp$gen) + 1
  next.gen$gen <- genno

  #remove entries in this gen that were also in previous gens:
  #p1[!(p1$idx %in% anchead$idx),]
  next.gen <- next.gen[!(next.gen$idx %in% childspp$idx),]

  if(nrow(next.gen)>0){

    next.gen$trace <- T


    #Don't trace rows with -1 birtht
    next.gen$trace[next.gen$birtht == -1] <- F
    #Don't trace rows where idx == actparent or pasparent
    next.gen$trace[next.gen$idx == next.gen$actparent | next.gen$idx == next.gen$pasparent] <- F

    #TODO - this line causes warnings - probably better to go through each individually (if we care about what next.gen$t is for bad traces)
    #next.gen$t[!next.gen$trace] <- min(bigspeciestable$rectime[bigspeciestable$spp==next.gen$idx[!next.gen$trace]])


    #childspp<-childspp[childspp$gen<(genno-1),]
    childspp<-rbind(childspp,next.gen)



    for(ii in 1:nrow(childspp))childspp$t[ii] <- max(bigspeciestable$rectime[bigspeciestable$spp == childspp$idx[ii]])
  }

  if(doplot){
    plotphyl(childspp[childspp$trace,])
    title(main=sprintf("%d Speciations",genno))

    points(x=childspp$idx[!childspp$trace],y=childspp$t[!childspp$trace],pch=19,cex=2)
  }

  return(childspp)
}


makephyl <- function(fmanc,fanchead,verbose = T){

  phyl <- nextgen(fmanc,fanchead)
  newrows<- nrow(phyl)
  depth<-0
  while(newrows>0){
    oldr<- nrow(phyl)
    phyl<-nextgen(fmanc,phyl,verbose=F,doplot=F)

    newr<- nrow(phyl)

    newrows <- newr-oldr
    depth <- depth + 1
    if(verbose)message(sprintf("\nDepth %d: Found %d new species",depth,newrows))


  }

  return(phyl)
}

