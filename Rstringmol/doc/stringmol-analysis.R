## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----demoraw,eval=T, include=T,message=F--------------------------------------
require("stringr")
require("Rstringmol")

## ----load_data, eval=FALSE, include=TRUE--------------------------------------
#  smsp <- list()
#  
#  for(rr in 1:8){
#    fn <- sprintf("rdata/run%d/rprops%02d.RData",rr,rr)
#    load(fn)
#    smsp[[rr]]<-rundata
#  }
#  

## ----config, eval=FALSE, include=TRUE-----------------------------------------
#  
#  require(stringr)
#  myoma = c(5.5,4.5,4,4.5)
#  mymar = c(0.5,0.5,0.5,0.5)
#  lax=c(T,F,F,F,T,F,F,F)
#  rax=c(F,F,F,T,F,F,F,T)
#  bax=c(F,F,F,F,T,T,T,T)
#  cols = c("orange","blue")

## ----popfig, fig.height=6, fig.width=10, message = F, warning=FALSE, eval=FALSE, include=TRUE----
#  
#  par(mfrow=c(2,4),oma=myoma)
#  
#  for(rr in 1:length(smsp))
#    figpopdy(smsp[[rr]],popdata[[rr]],lax[rr],rax[rr],bax[rr],mymar,cols,rr)
#  title("Population and number of species",outer=T)
#  
#  # To make a legend add a single empty figure over the whole graphics device,
#  # and then put the legend on *that* figure
#  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
#  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
#  legend('bottom',legend=c("population (left axis)",
#                           "bound population (left axis)",
#                           "n spp (right axis)",
#                           "bound n spp (right axis)"),
#         col=c(cols[1],cols[1],cols[2],cols[2]),
#         xpd = T, horiz = TRUE,  seg.len=3,
#         bty = 'n',inset=c(0,0),lty=c(1,2,1,2))

## ----figlentime,fig.height=6, fig.width=10, eval=FALSE, include=TRUE----------
#  
#  par(mfrow=c(2,4),oma=myoma)
#  for(rr in 1:length(smsp))
#    figlentime(smsp[[rr]],prop=1,ylim=c(10,110),ylim2 = c(0,400),
#               lax[rr],rax[rr],bax[rr],mymar,cols)
#  
#  title("Median string length (orange, l axis) and reaction time (blue, r axis)",outer=T)
#  
#  # Legend
#  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
#  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
#  legend('bottom',legend=c("string length (left axis)",
#                           "reaction time (right axis)"),
#         col=c(cols[1],cols[2]),lty=1, xpd = T, horiz = TRUE,  seg.len=3,
#         bty = 'n',inset=c(0,0))
#  
#  mtext(text="time",side=1,line=-3,outer=TRUE)
#  mtext(text="string length",side=2,line=-1.5,outer=TRUE)
#  mtext(text="reaction time",side=4,line=-1.5,outer=TRUE)

## ----fig.height=6, fig.width=10, eval=FALSE, include=TRUE---------------------
#  par(mfrow=c(2,4),oma=myoma,mar=c(3,3,0,0))
#  for(rr in 1:length(smsp))
#    figrepvpar(smsp[[rr]],ylim=c(0,70),nsteps=F,lax[rr],rax[rr],bax[rr]
#               ,mymar,cols,ratio=T,v2=T)
#  title(sprintf("Parasitic reactions as a percentage of all reactions"), outer = T)

## ----fig.height=2.5, fig.width=10,echo=T,warning=F, eval=FALSE, include=TRUE----
#  
#  t1 <-  90000
#  t2 <-  600000
#  t3 <-  680000
#  t4 <- 1080000
#  t4a<- 1240000
#  t5 <- 1500000
#  t6 <- 2000000
#  
#  par(mfrow=c(1,2),oma=c(2.5,2.5,0,2.5))#,oma=c(0,0,0,0))
#  rr<-2
#  
#  figrepvpar(smsp[[rr]],ylim=c(0,70),nsteps=F,lax[rr],rax[rr],bax[rr],ratio=T)
#  axis(1)
#  axis(2)
#  yl <- 66
#  yt <- 68
#  segments(x0=t1,y0=0,y1=yl,lty=3);text(labels = "t1",x=t1,y=yt)
#  segments(x0=t2,y0=0,y1=yl,lty=3);text(labels = "t2",x=t2,y=yt)
#  segments(x0=t3,y0=0,y1=yl,lty=3);text(labels = "t3",x=t3,y=yt)
#  segments(x0=t4,y0=0,y1=yl,lty=3);text(labels = "t4",x=t4,y=yt)
#  segments(x0=t5,y0=0,y1=yl,lty=3);text(labels = "t5",x=t5,y=yt)
#  segments(x0=t6,y0=0,y1=yl,lty=3);text(labels = "t6",x=t6,y=yt)
#  
#  mtext("Time", side = 1, line = 2, cex = 1)
#  mtext("Parasitic Reactions (%)", side = 2, line = 2, cex = 1)
#  
#  figlentime(smsp[[rr]],prop=1,ylim2=c(70,270),lax[rr],rax[rr],bax[rr],mymar,cols,plotlen = F,quartiles = F)
#  axis(1)
#  axis(4,)
#  yl <- 264
#  yt <- 267
#  segments(x0=t1,y0=0,y1=yl,lty=3);text(labels = "t1",x=t1,y=yt)
#  segments(x0=t2,y0=0,y1=yl,lty=3);text(labels = "t2",x=t2,y=yt)
#  segments(x0=t3,y0=0,y1=yl,lty=3);text(labels = "t3",x=t3,y=yt)
#  segments(x0=t4,y0=0,y1=yl,lty=3);text(labels = "t4",x=t4,y=yt)
#  segments(x0=t5,y0=0,y1=yl,lty=3);text(labels = "t5",x=t5,y=yt)
#  segments(x0=t6,y0=0,y1=yl,lty=3);text(labels = "t6",x=t6,y=yt)
#  
#  mtext("Time", side = 1, line = 2, cex = 1)
#  mtext("Reaction Time", side = 4, line = 2, cex = 1)

## ----ccplots,fig.height=4.8, fig.width=8, results='hide',message=FALSE,warning=FALSE, eval=FALSE, include=TRUE----
#  rundata = smsp[[2]]
#  par(mfrow=c(2,3),oma=c(7,4,0,0)+0.1, mar = c(0,0,1,1) + 0.1)
#  
#  ccplot(dts(90000,rundata),title="t1",doy=T)
#  ccplot(dts(600000,rundata),title="t2")
#  ccplot(dts(680000,rundata),title="t3")
#  ccplot(dts(1080000,rundata),title="t4",dox=T,doy=T)
#  ccplot(dts(1500000,rundata),title="t5",dox=T)
#  ccplot(dts(2000000,rundata),title="t6",dox=T)
#  
#  title(xlab = "Timesteps",
#        ylab = "Population",
#        outer = TRUE, line = 3)
#  
#  # Legend
#  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
#  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
#  legend("bottom",horiz=T,border = 'n',
#         legend=c("arena size","carrying capacity","replicators","seed replicator")
#         ,lty=c(2,2,1,1),col=c("red","blue","green","black"),bty="n")

## ----plotpyl, fig.height=6, fig.width=12,echo=T,eval=T,message=F, eval=FALSE, include=TRUE----
#  par(mfrow=c(2,4))
#  library(png)
#  imnos<- c(t1,t2,t3,t4,t5,t6)
#  par(mar=c(0.2,0.2,0.2,0.2),oma=c(2,2,0,0))
#  
#  mxl <- c(0,650000)
#  
#  for(rr in 1:8){
#    phyl<-readRDS(sprintf("rdata/run%d/phyldata%02d.RDS",rr,rr))
#    if(rr==5){
#      # SPECIAL CASE: We have to set trace<-F for one extra molecule
#      # because it wasn't recorded properly in the logfiles:
#      phyl$trace[phyl$idx==285937]<-F
#    }
#    plotphyl(phyl[phyl$trace & !is.na(phyl$idx),],xlim=mxl)
#  }
#  

