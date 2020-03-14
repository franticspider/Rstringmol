



require(vioplot)
require(stringr)

#Here's the basic plotting function:
# ptype is the switch btween 'beta', 'k' and 'length'#' Load mass spec data from a 2-column space-delimited text file
#'
#' @param data is the datafile - an Rstringmol data frame of properties
#' @param cname the name of property column that is to be plotted. Set to NULL to plot all
#' @param text the name of the plot
#' @param start the lower x axis limit
#' @param from the first sample time
#' @param to the last sample time
#' @param by the stepsize between samples
#' @param ptype the plot type, can be "beta", "k", "len" and "plen", and "nos" (for raw numbers)
#' @param wex the scalar for the vioplot width. defaults to 15
#' @param col the colour of the vioplot. defaults to "black"
#' @export
#' @examples
onevio <- function(data,cname="",text="",start=0,from=0,to=2000000,by=20000,ptype,wex=15,col="black",doplot=F){

  dlen <- length(data)

  #message("inside onevio()!")
  if(doplot)plot(NA,xlim=c(start,to),ylim=myylim,xaxt='n')
  dp <- 1
  times <- seq(from,to,by)

  summary <- data.frame(time = times, val = 0)

  for(rr in seq(from,to,by)){
    if(dp<=dlen){

      ng <- data[[dp]]

      if(is.null(cname)){
        reps <- ng
      }else{
        reps <- ng[ ng[,cname] ,]
      }

      if(ptype == "beta"){
        rct <- rep(reps$nsteps,reps$nobs)
      }

      if(ptype == "k"){
        rct <- rep(1/reps$bprob,reps$nobs)
        rct[rct>100]<-5
      }

      if(ptype == "len"){
        rct <- rep(str_length(reps$actseq),reps$nobs)
      }

      if(ptype == "plen"){
        rct <- rep(str_length(reps$passeq),reps$nobs)
      }

      if(ptype == "nos"){
        rct <- reps$nobs
      }

      if(length(unique(rct))>1){
        #vioplot( rct,wex= 15 * sum(reps$nobs),range = 1, width = 100000, at= rr-0,add=T,border=NA,col="black",drawRect = F)
        if(doplot){
          vioplot( rct,wex= wex * sum(reps$nobs),range = 1, width = 100000, at= rr-0,add=T,border=NA,col=col,drawRect = F)
          points(x=rr,y=mean(rct),col="red",pch=19,cex = 1)
        }
        summary$val[dp]<- mean(rct)

      }else{
        if(length(unique(rct))==1){
          if(doplot)points(x=rr,y=mean(rct),col="red",pch=19,cex = 1)
          summary$val[dp]<- mean(rct)
        }else{
          if(doplot)points(x=rr,y=0,col="black",pch=4,cex = 1)
          summary$val[dp]<- 0
        }
      }

      #nos is a special case
      if(ptype == "nos"){
        summary$val[dp] <- sum(reps$nobs)
      }

    }else{
      break
    }
    dp <- dp + 1


    #to emphasise timesteps by point size, do this:
    #if(!(rr%%200000))
    #  points(x=rr,y=mean(rct),col="red",pch=19,cex = 1.5)

  }
  if(doplot)text(text,x=xtitlepos,y=ytitlepos,srt=90,cex = 2)

  return(summary)
}






sumvp <- function(rundata,pltyp,myylim=c(0,500)){


  #ALL:
  sumall <- onevio(rundata,NULL,"All Reactions",0,20000,2000000,20000,pltyp,wex=7.5,col="blue")

  # REPS
  sumreps <- onevio(rundata,"pp_Replicator","Replicator",0,20000,2000000,20000,pltyp)

  # PARASITES
  sumpars <- onevio(rundata,"np_Parasite","Parasites",0,20000,2000000,20000,pltyp)

  ## Hypercycles
  #sumhyp <- onevio(rundata,"np_Hypercycle","Hypercycles",0,20000,2000000,20000,pltyp)
  # MUTUAL REPLICATORS
  sumhyp <- onevio(rundata,"np_MutualRepl","Mutual Repls",0,20000,2000000,20000,pltyp)


  # New Product
  sumnew <- onevio(rundata,"pp_NewProduct","New Product",0,20000,2000000,20000,pltyp)

  # No Product
  sumno <- onevio(rundata,"pp_NoProduct","No Product",0,20000,2000000,20000,pltyp)

  # Jumper
  sumjum <- onevio(rundata,"pp_Jumper","Jumper",0,20000,2000000,20000,pltyp)

  plot(sumall,type="l",lwd = 4,ylim=myylim,xlab="Timesteps",ylab=pltyp)
  lines(sumreps, col="red",    lwd=2)
  lines(sumpars, col="green",  lwd=2)
  lines(sumhyp,  col="purple", lwd=2)
  lines(sumnew,  col="blue",   lwd=2)
  lines(sumno,   col="orange", lwd=2)
  lines(sumjum,  col="brown",  lwd=2)

  legend("topleft",ncol=2,
         legend = c("All Reactions", "Replicators", "Parasites", "Mutual Reps", "New Product", "No Product", "Jumpers"),
         col    = c("black",         "red",         "green",     "purple",      "blue",        "orange",     "brown"),
         lwd    = c(4,2,2,2,2,2,2)
  )


}



















