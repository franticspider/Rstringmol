




rp.ratio.plot <- function(rundata,norm=F,title = "Plot"){

  rpr <- vector(length=length(rundata))

  for(rr in 1:length(rundata)){
    md<-rundata[[rr]]
    rpr[rr] = sum(md$nobs[md$np_Parasite])/sum(md$nobs[md$pp_Replicator | md$pp_SelfReplicator])
  }

  xvals =seq(20000,2000000,20000)
  plot(x=xvals,y=rpr,type="l"#,ylim=c(0,250000)#max(cvals,ovals,mvals,tvals))
       ,xlim=c(0,2000000),xlab="time",ylab="R-P raio",main=title)
}
