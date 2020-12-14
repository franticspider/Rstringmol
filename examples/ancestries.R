



anctable <- function(head,time,data,ngen=1,verbose=F){

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

  #Handle the case where the input string has birtht -1 because it was created after a log but before a restart:
  childspp$trace[childspp$birtht == -1]<-F

  mychild <- childspp[childspp$gen == max(childspp$gen) & childspp$trace,]

  if(verbose)message(sprintf("mychild has %d rows",nrow(mychild)))

  if(nrow(mychild)>0){

    next.gen <- anctable(mychild,2000000,bigspeciestable,ngen=1,verbose)
    genno <- max(childspp$gen) + 1
    next.gen$gen <- genno

    #remove entries in this gen that were also in previous gens:
    #p1[!(p1$idx %in% anchead$idx),]
    next.gen <- next.gen[!(next.gen$idx %in% childspp$idx),]

    if(nrow(next.gen)>0){

      next.gen$trace <- T


      #Don't trace rows with -1 birtht
      next.gen$trace[next.gen$birtht == -1] <- F
      #Don't trace rows where actparent or pasparent is NA
      next.gen$trace[is.na(next.gen$actparent)]<-F
      next.gen$trace[is.na(next.gen$pasparent)]<-F
      #Don't trace rows where idx == actparent or pasparent
      next.gen$trace[next.gen$idx == next.gen$actparent | next.gen$idx == next.gen$pasparent] <- F

      #TODO - this line causes warnings - probably better to go through each individually (if we care about what next.gen$t is for bad traces)
      #next.gen$t[!next.gen$trace] <- min(bigspeciestable$rectime[bigspeciestable$spp==next.gen$idx[!next.gen$trace]])


      #childspp<-childspp[childspp$gen<(genno-1),]
      childspp<-rbind(childspp,next.gen)

      for(ii in 1:nrow(childspp)){
        childspp$t[ii] <- max(bigspeciestable$rectime[bigspeciestable$spp == childspp$idx[ii]])
      }
    }

    if(doplot){
      plotphyl(childspp[childspp$trace,])
      title(main=sprintf("%d Speciations",genno))
      points(x=childspp$idx[!childspp$trace],y=childspp$t[!childspp$trace],pch=19,cex=2)
    }

    return(childspp)
  }
  else{
    return(NULL)
  }
}



makephyl <- function(fmanc,fanchead,verbose = T){

  phyl <- nextgen(fmanc,fanchead)

  if(!is.null(phyl)){

    newrows<- nrow(phyl)
    depth<-0

      while(newrows > 0){
      oldr<- nrow(phyl)
      newphyl<-nextgen(fmanc,phyl,verbose=verbose,doplot=F)


      if(!is.null(newphyl)){
        newr<- nrow(newphyl)

        newrows <- newr-oldr
        depth <- depth + 1
        if(verbose)message(sprintf("\nDepth %d: Found %d new species",depth,newrows))
        phyl <- newphyl
      }else{
        break
      }
    }
  }
  return(phyl)
}


