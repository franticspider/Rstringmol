

print_mols <- function(type,
                       act,
                       pas,
                       pro,
                       conv,
                       actout = NA,
                       pasout = NA) {
  message(sprintf("%s:", type))
  message(sprintf("ACT          = %s", act))
  message(sprintf("PAS          = %s", pas))
  message(sprintf("PRO          = %s", pro))
  message(sprintf("CONV         = %s", conv))
  if (!is.na(actout)) {
    message(sprintf("ACTOUT       = %s", actout))
    message(sprintf("PASOUT       = %s", pasout))
  }
  message("")
}










############################################
## DETECTION REACTION PROPERTIES



#' Identify the network-level properties in a reaction
#'
#' @param fdata data frame containing actseq and passeq fields
#' @param dummy.result an optional result of runReactionFP passed in (for testing)
#' @export
pairwise.properties <- function(fdata,dummy.result=NULL){

  #create a data frame of unique reactions for this file:
  ur <- unique(data.frame(actseq = fdata$actseq, passeq = fdata$passeq, stringsAsFactors = F))

  #ur$type <- "unknown"
  ur$product <- ""

  #These are the properties agreed with Susan
  ur$pp_ActiveMod <- F
  ur$pp_PassiveMod <- F
  ur$pp_SelfMod <- F

  ur$pp_NoProduct <-F
  ur$pp_NewProduct <-F

  ur$pp_biolRep <-F
  ur$pp_SelfReplicator <-F
  ur$pp_Repl1 <- F # formerly ActiveReplicator <-F
  ur$pp_Repl2 <- F # formerly PassiveReplicator <- F

  ur$pp_Jumper <-F

  ur$actlen <- str_length(ur$actseq)
  ur$paslen <- str_length(ur$passeq)
  ur$outputA <- ""
  ur$outputP <- ""

  ur$ccopy <-0
  ur$cmove <-0
  ur$cover <-0
  ur$ctogg <-0

  ur$bprob <-0
  ur$nsteps <- 0

  ur$nprod <- 0

  ur$oneway <-F

  for (rr in 1:nrow(ur)){
    #################################################
    # Reaction properties

    #rt <- reaction_type(ur$actseq[rr],ur$passeq[rr])
    # Run the reaction(s):

    if(is.null(dummy.result))
      result <- runReactionFP(c(ur$actseq[rr],ur$passeq[rr]))
    else{
      message("Calculating properties based on external result")
      result <- dummy.result
    }

    ur$oneway[rr] <- result$deterministicBind

    if(ur$actseq[rr] != result$mActive)
      ur$pp_ActiveMod[rr] <- T

    if(ur$passeq[rr] != result$mPassive)
      ur$pp_PassiveMod[rr] <- T

    ur$pp_SelfMod[rr] <- ur$pp_PassiveMod[rr] | ur$pp_ActiveMod[rr]

    if(result$product == "empty")
      ur$pp_NoProduct[rr] <- T
    else{
      if(result$product != ur$actseq[rr] & result$product != ur$passeq[rr])
        ur$pp_NewProduct[rr] <- T
    }

    ur$ccopy[rr] <- result$ccopy
    ur$cmove[rr] <- result$cmove
    ur$cover[rr] <- result$cover
    ur$ctogg[rr] <- result$ctogg

    ur$bprob[rr] <- result$bprob
    ur$nsteps[rr] <- result$count

    ur$nprod[rr] <- result$nprod
    ur$outputA[rr] <- result$mActive
    ur$outputP[rr] <- result$mPassive

    #####################################################
    #Replicator types

    # Count the inputs copies of active and passive as described in the tech report
    actin<-1
    pasin<-1
    if(ur$actseq[rr]==ur$passeq[rr]){
      actin <- actin+1
      pasin <- pasin+1
    }

    act <- 0
    if(ur$actseq[rr]==result$mActive)  act <- act + 1
    if(ur$actseq[rr]==result$mPassive) act <- act + 1
    if(ur$actseq[rr]==result$product)  act <- act + 1
    if(act>actin)ur$pp_Repl1[rr] <-T

    pct <- 0
    if(ur$passeq[rr]==result$mActive)  pct <- pct + 1
    if(ur$passeq[rr]==result$mPassive) pct <- pct + 1
    if(ur$passeq[rr]==result$product)  pct <- pct + 1
    if(pct>pasin)ur$pp_Repl2[rr] <-T


    # The technical repl property is now "networked": see tech report
    # Here we are using pp_biolRep as a naive term - NOT TO BE USED TO BUILD OTHER PROPERTIES
    ur$pp_biolRep[rr] <- ur$pp_Repl1[rr] | ur$pp_Repl2[rr]

    # We can still spot self-replicators here:
    if((ur$pp_Repl1[rr] | ur$pp_Repl2[rr]) & (ur$actseq[rr]==ur$passeq[rr])){
      ur$pp_SelfReplicator[rr] <- T
    }


    #####################################################
    #Jumper types


    Jumper1 <- F
    Jumper2 <- F
    if(result$mActive != ur$actseq[rr] & ur$actseq[rr] == result$product)
      Jumper1 <- T
    if(result$mPassive != ur$passeq[rr] & ur$passeq[rr] == result$product)
      Jumper2 <- T
    if(Jumper1 | Jumper2)
      ur$pp_Jumper[rr] <- T


    #####################################################
    # Other useful data from the reaction:

    ur$product[rr] <- result$product
    ur$dexec[rr] <- result$deterministicExec
    ur$dbind[rr] <- result$deterministicBind
    ur$nobs[rr] <- nrow(fdata[fdata$actseq == ur$actseq[rr]
                              & fdata$passeq == ur$passeq[rr],])
  }

  return(ur)
}




#' Identify the network-level properties in a reaction
#'
#' @param fn the file name
#' @export
network.properties <- function(nwp,dummy.reverse.pwp = NULL){

  nwp$np_Parasite2 <- F
  nwp$np_Parasite1 <- F
  nwp$np_Parasite <- F

  nwp$np_MutualRepl <- F
  nwp$np_Hypercycle <- F

  for(rr in 1:nrow(nwp)){
     #Run the 'flip' reaction
     if(nwp$actseq[rr] == nwp$passeq[rr]){
       #Can't be a parasite - no need to change
       #Can't be a hypercycle - no need to change
       if(nwp$pp_SelfReplicator[rr]){
         nwp$np_MutualRepl[rr] <- T
       }
     }else{
        if(is.null(dummy.reverse.pwp)){
          flipdata <- data.frame(actseq=nwp$passeq[rr],passeq=nwp$actseq[rr],stringsAsFactors = F)
          flip <- pairwise.properties(flipdata)
        }else{
          flip <- dummy.reverse.pwp
        }

       #PARASITE
       if(nwp$pp_Repl2[rr] & ! (nwp$pp_Repl1[rr] | flip$pp_Repl2)){
           nwp$np_Parasite2[rr] <- T
       }

       if(nwp$pp_Repl1[rr] & ! (nwp$pp_Repl2[rr] | flip$pp_Repl1)){
         nwp$np_Parasite1[rr] <- T
       }

       nwp$np_Parasite[rr] <- nwp$np_Parasite2[rr] | nwp$np_Parasite1[rr]



       #MUTUAL

       if(( nwp$pp_Repl2[rr] | flip$pp_Repl1) & (nwp$pp_Repl1[rr] | flip$pp_Repl2)){
         nwp$np_MutualRepl[rr]<-T
       }



       #HYPERCYCLE
       if(nwp$np_MutualRepl[rr]){
         ina <- data.frame(actseq=nwp$actseq[rr],passeq=nwp$actseq[rr],stringsAsFactors = F)
         pwpa <- pairwise.properties(ina)
         inp <- data.frame(actseq=nwp$passeq[rr],passeq=nwp$passeq[rr],stringsAsFactors = F)
         pwpp <- pairwise.properties(inp)

         if(nwp$np_MutualRepl[rr] & !pwpa$pp_SelfReplicator & !pwpp$pp_SelfReplicator)
           nwp$np_Hypercycle[rr] <- T

       }
     }
  }
  return(nwp)
}




############################################
## PARSING BIG DATA FILES

#' Create a data structure from a *.conf* file for passing into the standard *pairwise.properties* function
#'
#' @param fn the file name
#' @export
pairwise.properties.from.conf <- function(fn){

  fdata <- rconf_rdata(fn)
  fdata$actseq <- as.character(fdata$actseq)
  fdata$passeq <- as.character(fdata$passeq)

  ur <- pairwise.properties(fdata)

  return(ur)
}




#' Obtain the reaction properties from an active and passive sequence
#'
#' @param actseq the active sequence
#' @param passeq the passive sequence
#' @export
propfromseqs <- function(actseq,passeq){

  fd <- data.frame(actseq=actseq,passeq=passeq,stringsAsFactors = F)
  pwp <- pairwise.properties(fd)
  ng <- network.properties(pwp)
  return(ng)
}




#' Obtain the reaction properties of all the reactions in a stringmol run
#'
#' @param froot the pathway to the data files
#' @param from the first timestep
#' @param to the last timestep
#' @param step the timestep increment
#' @param outfn the name of an R data object to save it to
#' @param verbose verbose output
#' @export
runproplist <- function(froot,from=20000,to=2000000,step=20000,outfn=NULL,verbose = F){

  rundata <- list()
  dp <- 1
  for(rr in seq(from,to,step)){

    #testinfn <- sprintf("~/Desktop/paulien/smsp/1705smsp/out3/out1_%d.conf",rr)
    testinfn <- sprintf("%s%d.conf",froot,rr)
    if(verbose)message(sprintf("%d: file is %s",rr,testinfn))

    ug <- pairwise.properties.from.conf(testinfn)

    ng <- network.properties(ug)

    rundata[[dp]] <- ng

    dp <- dp + 1

  }

  if(!is.null(outfn)){
    save(rundata,file=outfn)
  }

  return(rundata)

}




#' Obtain the species count and population sizes from a stringmol datafile
#'
#' @param froot the path and file stem of the data files, e.g. "~/Desktop/smsp/out3/out1_%d.conf"
#' @param trange the time range and step size of the log files
#' @export
spdata <- function(froot="",trange=seq(20000,2000000,20000)){

  nspp <- vector(length=length(trange))
  nmols <- vector(length=length(trange))

  idx <- 1
  for(tt in trange){
    fn <- sprintf("%s%d.conf",froot,tt)
    #message(sprintf("Loading file %s",fn))
    raw <- read.table(fn,stringsAsFactors = F,fill=T,sep=",")

    ##lx <- raw[str_detect("TOTSPPCT",str_sub(raw$V1,1,8)),]

    lx <- data.frame(raw[str_detect(str_sub(raw$V1,1,40),"extant"),],stringsAsFactors = F)
    if(nrow(lx) != 1){
      message(sprintf("%d: problem finding extant species"))
    }
    else{
      nspp[idx] <- as.numeric(str_sub(lx,17,str_length(lx)))
    }
    #message(sprintf("%d: Found %d species",tt,nspp[tt]))

    lx <- data.frame(raw[str_detect(str_sub(raw$V1,1,40),"NUMAGENTS"),],stringsAsFactors = F)
    if(nrow(lx) != 1){
      message(sprintf("%d: problem finding number of agents"))
    }
    else{
      nmols[idx] <- as.numeric(str_sub(lx,11,str_length(lx)))
    }
    #message(sprintf("%d: Found %d mols and %d species",tt,nmols[idx],nspp[idx]))
    idx <- idx + 1
  }

  rdf <- data.frame(time=trange,nspp=nspp,nmols=nmols,stringsAsFactors = F)

  return (rdf)
}




#' Obtain the species count and population sizes from a stringmol datafile
#' From https://www.dataanalytics.org.uk/make-transparent-colors-in-r/
#'
#' @param col the colour
#' @param percent the percentage transparency
#' @export
opcol <- function(col,percent = 70){
  rgbc<-col2rgb(col)
  t.col <- rgb(rgbc["red",1], rgbc["green",1], rgbc["blue",1],
               max = 255,
               alpha = (100 - percent) * 255 / 100)
  return(t.col)
}




#' Deterministic calculation of the carrying capacity of a self-replicator
#'
#' @param bprob the bind probability
#' @param nsteps the number of timesteps to replicate
#' @param initpop the initial population
#' @param maxt the time limit on the calculation
#' @param death the death rate
#' @param lim the arena size
#' @export
well.mixed.dynamics <- function(bprob,nsteps,initpop=100,maxt=20000,death=0.0005,lim=(125*100)){

  t <- 1

  pop <- vector(length=maxt)
  bpop <- vector(length=maxt)
  pop[] <- 0
  bpop[] <- 0
  pop[1] <- initpop


  for(tt in 2:maxt){

    #decay
    pop[tt]<- pop[tt-1]*(1-death)
    bpop[tt]<- bpop[tt-1]*(1-death)

    #bind - count a pair of bound reactants as 1 in the bpop array - keeps the maths easier
    bpop[tt] <- bpop[tt]+(0.5*(pop[tt]*bprob))
    pop[tt]  <-  pop[tt]-((pop[tt]*bprob))

    #occupancy
    occ <- pop[tt]+(2*bpop[tt])

    #birthrate is dependant upon occupancy, number of bound mols, and rep.rate
    pop[tt] <- pop[tt] +(((lim-occ)/lim)  * (bpop[tt]) * (1/nsteps) )

    #unbinding after replication (assume replication and unbinding are the same rate)
    pop[tt]  <- pop[tt]  + (2 * bpop[tt] * 1/nsteps)
    bpop[tt] <- bpop[tt] - (bpop[tt] * 1/nsteps)
  }

  return(data.frame(pop,bpop))
}




#' Get reaction data from a time, and sort it by frequency of occurrence
#' NB: it is currently assumed that the data is in timesteps of 20000
#'
#' @param time the time
#' @param rd the data
#' @param plax whether to plot the left axis
#' @param plax whether to plot the right axis
#' @param plax whether to plot the bottom axis
#' @param mymar the figure margins
#' @param cols the plot colours (defaults to orange and blue)
#' @param rno the run number
#' @export
dts <- function(time,rd,sr = T){
  dt <- rd[[time/20000]]
  dt <- dt[order(dt$nobs,decreasing = T),]
  if(sr)
    dt <- dt[dt$pp_SelfReplicator,]
  return(dt)
}




#' Figure of population dynamics in a stringmol run
#'
#' @param rundata the reaction data from the run
#' @param popd the population data from the run (including unbound molecules)
#' @param plax whether to plot the left axis
#' @param plax whether to plot the right axis
#' @param plax whether to plot the bottom axis
#' @param mymar the figure margins
#' @param cols the plot colours (defaults to orange and blue)
#' @param rno the run number
#' @export
figpopdy <- function(rundata,popd,plax=F,prax=F,pbax=F,mymar = c(0.5,0.5,0.5,0.5),cols = c("orange","blue"),rno=NULL){
  time <- seq(1:length(rundata))
  val <- time
  var <- time
  par(mar=mymar)

  for(tt in 1:length(time)){
    val[tt] <- sum(rundata[[tt]]$nobs)*2
    var[tt] <- length(unique(c(rundata[[tt]]$actseq,rundata[[tt]]$passeq)))
  }
  plot(x=time*20000,axes=F,y=val,type="l",xlab="time",ylab="population",col=cols[1],ylim=c(0,12500),lty=2)
  box()
  lines(x=popd$time,y=popd$nmols,lty=1,col=cols[1])
  if(pbax)axis(1)
  if(plax)axis(2)
  par(new = T)
  plot(x=time*20000,y=var,col=cols[2],axes=F,type="l",ylim=c(0,2000),ylab="",xlab="",lty=2)
  lines(x=popd$time,y=popd$nspp,lty=1,col=cols[2])
  if(prax)axis(4)
  if(!is.null(rno))text(labels = sprintf("Run %d",rno),x=0,y=1900,pos = 4)

}




#' Figure of molecule length and reaction time
#'
#' @param rundata the reaction data from the run
#' @param prop the proportion of all reacitons to be included (defaults to 0.1, i.e the top 10 percent)
#' @param ylim axis limits for the length (left) axis
#' @param ylim2 axis limits for the reaction time (right) axis
#' @param plax whether to plot the left axis
#' @param prax whether to plot the right axis
#' @param pbax whether to plot the bottom axis
#' @param mymar the figure margins
#' @param cols the plot colours (defaults to orange and blue)
#' @param rno the run number
#' @param plotlen whether to plot the length
#' @param quartiles whether to plot the quartiles
#' @export
figlentime <- function(rundata,prop=0.1,ylim=NULL,ylim2=NULL,plax,prax,pbax,mymar = c(0.5,0.5,0.5,0.5),cols = c("orange","blue"),rno=NULL,plotlen=T,quartiles=T){

  time <- seq(1:length(rundata))
  #val <- time
  #val2 <- time

  # these will be arrays to hold the median and quartile values:
  valm <- time
  valql <- time
  valqh <- time

  lvm <-time
  lvl <-time
  lvu <-time

  par(mar=mymar)

  if(is.null(ylim))
    ylim=c(0,120)

  if(is.null(ylim2))
    ylim2=c(0,400)

  repvals <- list()

  #for each time
  for(tt in 1:length(time)){
    # Get the data
    rd <- rundata[[tt]]
    rd <- rd[order(rd$nobs,decreasing = T),]

    #work out the end point - this is to get the average of the top *n* reactions
    totrp <- sum(rd$nobs)*prop
    count<- 0
    rsum <- 0
    rsum2 <- 0
    idx <- 1
    firstentry <- T
    while(count<totrp){
      rsum <- rsum + (rd$nobs[idx]*0.5*(rd$actlen[idx] + rd$paslen[idx]))
      rsum2 <- rsum2 + (rd$nobs[idx]*(rd$nsteps[idx]))
      count <- count + rd$nobs[idx]

      #gather the data for medians
      if(firstentry){
        mdata <- rep( (rd$actlen[idx]), rd$nobs[idx])
        ldata <- rep(  (rd$nsteps[idx]) , rd$nobs[idx])
        firstentry <- F
      }
      else{
        mdata <- c(mdata,rep((rd$actlen[idx]), rd$nobs[idx]))
        ldata <- c(ldata,rep(  (rd$nsteps[idx]) , rd$nobs[idx]))
      }
      idx<-idx+1
    }

    #val[tt] <-rsum/count
    #val2[tt] <-rsum2/count
    q <- quantile(mdata)
    valm[tt] <- q[3]
    valql[tt] <- q[2] #25th
    valqh[tt] <- q[4] #75th
    #q <- quantile(mdata)

    q <- quantile(ldata)
    lvm[tt]<-q[3]
    lvl[tt]<-q[2]
    lvu[tt]<-q[4]
  }

  if(plotlen){
    # Plot string length in orange on the firxt axis:
    plot(x=time*20000,y=valm,type="l",xlab="time",ylab="string length",col=cols[1],axes=F,ylim=ylim)

    #plot quartiles first so line is overlaid
    if(quartiles){
      shf <- data.frame(t=time*20000,v=valql)
      shr <- data.frame(t = rev(time*20000), v = rev(valqh))
      sh <- rbind(shf,shr)
      polygon(x=sh$t,y=sh$v,col=opcol(cols[1],percent = 40),border=F)
    }

    lines(x = time*20000,y=valm,col=cols[1])
    #plot the lhs y axis now before we get a new ylim
    if(plax)axis(2)

    par(new=T)
  }

  plot(x=time*20000,y=lvm,type="l",xlab="time",ylab="",col=cols[2],ylim=ylim2,axes=F)

  if(quartiles){
    shf <- data.frame(t=time*20000,v=lvl)
    shr <- data.frame(t = rev(time*20000), v = rev(lvu))
    sh <- rbind(shf,shr)
    polygon(x=sh$t,y=sh$v,col=opcol(cols[2]),border=F)
  }
  lines(x = time*20000,y=lvm,col=cols[2])

  box()
  if(pbax)axis(1)
  if(prax)axis(4)
  if(!is.null(rno))
     text(labels = sprintf("Run %d",rno),x=0,y=380,pos = 4)

  #emphasise the orange line:
  if(plotlen){
    par(new=T)
    plot(x=time*20000,y=valm,type="l",xlab="time",ylab="string length",col=cols[1],axes=F,ylim=ylim)
  }
}



#' Figure of number of mutual-replicating vs parasitic reactions
#'
#' @param rundata the reaction data from the run
#' @param ylim y axis limits
#' @param plax whether to plot the left axis
#' @param prax whether to plot the right axis
#' @param pbax whether to plot the bottom axis
#' @param mymar the figure margins
#' @param cols the plot colours (defaults to orange and blue)
#' @param ratio whether to plot the ratio
#' @param v2 whether to count all replicators (v2=F) or just mututal replicators (v2=T)
#' @export
figrepvpar <- function(rundata,ylim=NULL,nsteps=F,plax,prax,pbax,mymar = c(0.5,0.5,0.5,0.5),cols = c("orange","blue"),ratio=F,v2=T){
  time <- seq(1:length(rundata))
  nrep <- time
  nrep <- time
  npar <- time
  par(mar=mymar)

  for(tt in 1:length(time)){
    if(nsteps){
      if(v2){
        nrep[tt] <- mean(rundata[[tt]]$nsteps[  rundata[[tt]]$pp_Replicator | rundata[[tt]]$pp_SelfReplicator  ])
      }else{
        nrep[tt] <- mean(rundata[[tt]]$nsteps[  rundata[[tt]]$np_MutualRepl | rundata[[tt]]$pp_SelfReplicator  ])
      }
      npar[tt] <- mean(rundata[[tt]]$nsteps[  rundata[[tt]]$np_Parasite])
    }
    else{
      if(v2){
        nrep[tt] <- sum(rundata[[tt]]$nobs[rundata[[tt]]$np_MutualRepl |rundata[[tt]]$pp_SelfReplicator])
      }else{
        nrep[tt] <- sum(rundata[[tt]]$nobs[rundata[[tt]]$pp_Replicator |rundata[[tt]]$pp_SelfReplicator])
      }

      npar[tt] <- sum(rundata[[tt]]$nobs[rundata[[tt]]$np_Parasite])
    }
  }

  if(!ratio){
    if(is.null(ylim))
      ylim=c(0,5000)
    plot(x=time*20000,y=nrep,type="l",xlab="time",ylab="population",col=cols[1],ylim=ylim,axes=F)
    box()
    if(pbax)axis(1)
    if(plax)axis(2)
    par(new = T)
    plot(x=time*20000,y=npar,col=cols[2],axes=F,type="l",ylim=ylim,ylab="",xlab="")
    if(prax)axis(4)
  }
  else{
    nrep <- 100*npar/(nrep+npar)
    if(is.null(ylim))
      ylim = range(nrep)
    plot(x=time*20000,y=nrep,type="l",xlab="time",ylab="population",col=cols[2],ylim=ylim,axes=F)
    box()
    if(pbax)axis(1)
    if(plax)axis(2)

  }
}





#' Figure of number of mutual-replicating vs parasitic reactions
#'
#' @param rates datframe of the two reaction rates: bind probability (bprob) and
#' @param title plot title
#' @param dox whether to plot x axis
#' @param doy whether to plot y axis
#' @export
ccplot <- function(rates,title,dox=F,doy=F){
  maxt = 10000
  lim = 125*100

  #drr  <-  d2r
  rates <- rates[rates$nobs > 1,]
  rates$carrycap <- 0

  #todo: figure out how to avoid redoing this from scratch!
  rp <- runReactionFP(c("WWGEWLHHHRLUEUWJJJRJXUUUDYGRHJLRWWRE$BLUBO^B>C$=?>$$BLUBO%}OYHOB","WWGEWLHHHRLUEUWJJJRJXUUUDYGRHJLRWWRE$BLUBO^B>C$=?>$$BLUBO%}OYHOB"))
  d.seed <- well.mixed.dynamics(rp$bprob,rp$count,maxt = maxt)

  plot(x=seq(1:maxt), y=d.seed$pop+(2*d.seed$bpop),type="l",xlim=c(0,10000),ylim=c(0,lim*1.1),axes=F)#xlab="timesteps",ylab="population")
  for(rr in 1:min(10,nrow(rates))){#nrow(d5rr)){
    d.rr <- well.mixed.dynamics(rates$bprob[rr],rates$nsteps[rr],maxt = maxt)

    lines(x=seq(1:maxt),y=d.rr$pop+(2*d.rr$bpop),col="green")

    rates$carrycap[rr]<-max(d.rr$bpop)
  }
  title(main = title,line=-1.5)
  lines(x=seq(1:maxt),y=d.seed$pop+(2*d.seed$bpop))
  segments(x0=0,x1=20000,y0=lim,lty=2,col="red")
  segments(x0=0,x1=20000,y0=(5671*2),lty=2,col="blue")

  if(dox){
    axis(side=1)
  }else{
    axis(side=1,labels=F)
  }
  if(doy){
    axis(side=2,labels=F)
    axis(side=2)
  }else{
  }
}



#############################################
## DEPRECIATED FUNCTIONS


## DEPRECIATED FUNCTION
#' Determine the type of stringmol reaction
#'
#' @param act the active stringmol
#' @param pas the passive stringmol
#' @keywords "Mass spectrum",
#' @export
#' @examples
#' p <- reaction_type("AAA","NNN")
reaction_typeFP_old <- function(act,
                                pas,
                                verbose = F,
                                detail = F) {


  .Deprecated(msg = "'reaction_typeFP_old' will be removed in the next version")

  #pro <- doReaction(c(act,pas))
  pro <- runReactionFP(c(act, pas))
  conv <- runReactionFP(c(pas, act))

  catalytic <- T

  # Check for changed parents!
  if (pro$mActive != act) {
    catalytic <- F
  }

  if (pro$mPassive != pas) {
    catalytic <- F
  }


  pro$rtype = "Macromutation"

  if (catalytic) {
    if (pro$mActive == pro$mPassive) {
      if (pro$product == "empty") {
        pro$rtype = "SelfSelfNoProduct"
      }
      else{
        if (pro$product == pro$mPassive) {
          pro$rtype = "SelfSelfReplicator"
          #print_mols(pro$rtype,act,pas,pro$product,conv$product)
        }
        else{
          pro$rtype = "SelfSelfDifferentProduct"
        }
      }
    }
    else{
      #conv <- doReaction(c(pas,act))
      #conv <- runReactionFP(c(pas,act))

      if (pro$product == "empty") {
        if (conv$product == "empty") {
          pro$rtype = "NonSelfNoProduct"
        }
        else{
          if (conv$product == conv$mPassive)
            pro$rtype = "Parasite"
          else{
            pro$rtype = "NonSelfDifferentProduct"
          }
        }
      }
      else{
        if (pro$product == pro$mPassive) {
          #TODO: We need to know if the assignment to active/passive is a coin toss!
          if (conv$product == "empty")
            pro$rtype = "Parasite"
          else{
            if (conv$product == pro$mPassive) {
              #TODO: extend this if detail = T
              pro$rtype = "Parasite"
            } else{
              if (conv$product == conv$mPassive) {
                #TODO: extend this if detail = T
                pro$rtype = "NonSelfReplicator"
                #print_mols(pro$rtype,act,pas,pro$product,conv$product)
              } else{
                pro$rtype = "NonSelfDifferentProduct"
                #print_mols(pro$rtype,act,pas,pro$product,conv$product)
              }
            }
          }
        }
        else{
          pro$rtype = "NonSelfDifferentProduct"
        }
      }
    }
  }
  else{
    pro$rtype = "Macromutation"
    if (verbose)
      print_mols(pro$rtype,
                 act,
                 pas,
                 pro$product,
                 conv$product,
                 pro$mActive,
                 pro$mPassive)
  }
  return(pro)
}



#' Determine the type of stringmol reaction
#'
#' @param act the active stringmol
#' @param pas the passive stringmol
#' @keywords "Mass spectrum",
#' @export
#' @examples
#' p <- reaction_type("AAA","NNN")
reaction_type <- function(act,pas,
                          result=NULL,
                          conv=NULL,
                          verbose = F,
                          detail = F) {


  .Deprecated(msg = "'reaction_typ' will be removed in the next version")

  if(is.null(result))
    result <- runReactionFP(c(act,pas))


  result$type = "undefined"

  #TODO: it should be possible to remove the "if(act == pas)" statement from the outer loop and nest it more deeply when looking for parasites
  if(act == pas){
    if((result$bprob < 0.0000001) || (result$product == "empty"))
      result$type <- "SelfSelfNoProduct"
    else{
      if(result$mActive != act || result$mPassive !=pas)
        result$type <- "Macromutation"
      else{
        if(result$product == pas)
          result$type <- "SelfSelfReplicator"
        if(result$product != pas && result$product != act)
          result$type <- "SelfSelfDifferentProduct"
      }
    }
  }else{
    #generate the converse reaction if it isn't supplied:
    if(is.null(conv))
      conv <- runReactionFP(c(pas,act))

    if((result$bprob < 0.0000001) || (result$product == "empty"))
      result$type <- "NonSelfNoProduct"
    else{
      if(result$mActive != act || result$mPassive !=pas)
        result$type <- "Macromutation"
      else{
        if(result$product == pas){
          if(conv$product == act)
            result$type <- "NonSelfReplicator"
          else
            result$type <- "Parasite"
        }
        if(result$product != pas && result$product != act)
          result$type <- "NonSelfDifferentProduct"
      }
    }
  }

  # for debugging:
  #pro$product = "BBBBBBB"
  return(result)
}


