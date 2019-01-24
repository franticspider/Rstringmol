plot.onepop <- function(infn,outfn,w,h, col.fun=rainbow, totpop = F, plotnspp = F, ploteachspp = T, log = "", minpop=0,title=NA,time.lim=c(0,2000000),ylim=c(0,12000)){
  x	<- read.table(infn,sep=",", col.names=c('time', 'species', 'count'))

  #diagnostics
  #message(sprintf("Reading file %s",infn))
  #message(sprintf("Read %d lines",nrow(x)))


  if(is.na(time.lim))
    time.lim <- range(x$time)

  valid.species <- unique(x$species)
  colours <- col.fun(length(valid.species))

  if(totpop){
    times <- as.data.frame(unique(x$time))
    times$pop <- 0
    colnames(times) <- c("time","pop")
    for(tt in 1:nrow(times)){
      data <- x[x$time == times$time[tt],]
      times$pop[tt] = sum(data$count)
    }

    tpmax =max(times$pop)
  }
  else{
    tpmax=0
  }



  if(plotnspp){
    nspp <- as.data.frame(unique(x$time))
    nspp$nspp <- 0
    colnames(nspp) <- c("time","nspp")
    for(tt in 1:nrow(times)){
      data <- x[x$time == times$time[tt],]
      nspp$nspp[tt] <- nrow(data)
    }
    #return(nspp)
  }

  #y.lim <- c(1, max(x$count,tpmax))#max(x$count))
  # pdf(file=outfn, height=h, width=w, title=outfn)
  #  	par(mar=c(5,5, 0.1, 0.1))
  plot(NA, xlim=time.lim,ylim=ylim, ann=FALSE, axes=FALSE, log = log )
  if(ploteachspp){
    sapply(1:length(valid.species), function(i){
      data <- x[x$species==valid.species[i],]
      if(max(data$count > minpop))
        lines(x=data$time, y=data$count + 0.1, col=colours[i])
    })
  }

  if(totpop)
    lines(x=times$time,y=times$pop,lwd=2,col="red")


  axis(2, las=2)
  title(ylab='population')
  x.positions <- seq(from=time.lim[1], to=time.lim[2], length.out=20)
  axis(1, at=x.positions, labels=sprintf('%d', as.integer(x.positions / 1e3)))


  if(plotnspp){
    par(new = T)
    plot(x=nspp$time,y=nspp$nspp,axes=F,xlab=NA,ylab=NA,type="l",lty=2,col="green", xlim=time.lim, ylim = c(0,2000))
    axis(side=4)
  }

  title(xlab=expression(time%*%10^3))
  if(!is.na(title))
  #else
    title(title)
  #box()
  # dev.off()

  #legend()

}
