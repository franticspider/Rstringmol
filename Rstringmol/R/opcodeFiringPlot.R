



opcode.firing.plot <- function(rundata,norm=F,title = "Plot",ylim=c(0,250000)){

  cvals <- vector(length=length(rundata))
  ovals <- vector(length=length(rundata))
  mvals <- vector(length=length(rundata))
  tvals <- vector(length=length(rundata))

  for(rr in 1:length(rundata)){
    md<-rundata[[rr]]
    #md <- d[d$nprod>1,]
    #if(nrow(md)>0){
    #md <- md[order(md$nobs,decreasing = T),]
    #message(sprintf("\nMax nprod for run %d is %d",rr*20000,max(d$nprod)))
    #for(kk in 1:nrow(md)){
    #  message(sprintf("row %d: nobs = %d, nprod = %d\nact=%s\npas=%s\n",kk,md$nobs[kk],md$nprod[kk],md$actseq[kk],md$passeq[kk]))
    #}
    cvals[rr]<- sum(md$ccopy)#/sum(md$ccopy)#-sum(md$cover)
    ovals[rr]<- sum(md$cover)#/sum(md$ccopy)
    mvals[rr]<- sum(md$cmove)#/sum(md$ccopy)
    tvals[rr]<- sum(md$ctogg)#/sum(md$ccopy)
    #}
  }

  xvals =seq(20000,2000000,20000)
  plot(x=xvals,y=cvals,type="l",ylim=ylim#max(cvals,ovals,mvals,tvals))
       ,xlim=c(0,2000000),xlab="time",ylab="Number of muti-producut reactions",main=title)
  lines(x=xvals,y=ovals,col="red")

  #axis(4,ylim=c(0,4))
  lines(x=xvals,y=mvals,col="green")
  lines(x=xvals,y=tvals*25,col="blue")

  #message(sprintf("max mvals = %d, max tvals = %d",max(mvals),max(tvals)))
  legend("topleft",legend=c("copy","move","toggle * 25","overwrite"),lty=1,col=c("black","green","blue","red"))
}
