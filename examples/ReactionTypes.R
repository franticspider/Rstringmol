




#froot="D:/sjh/stringmol/out3smsp/"
#froot="D:/sjh/stringmol/smsp/1705smsp/out5/"
froot="~/Desktop/paulien/smsp/1705smsp/out3/"

pdf(file = "~/smsp5types.pdf",width = 12,height = 10)

tstep <- 20000 #sometimes its 10000
tt <-   1500000 + tstep
sps <- list()
nsps <- 1
gotdata<-T
while(gotdata){
  fn <- sprintf("%sout1_%d.conf",froot,tt)

  message(sprintf("Working on file %s",fn))

  if(file.exists(fn)){

    #Get the data...
    data <- rconf_rdata(fn,summarize = F)

    actset <- unique(data$actseq)
    for(aa in 1:length(actset)){
      passet <- unique(data$passeq[data$actseq == actset[aa]])
      for(pp in 1:length(passet)){
        #message(sprintf("reacting %s with %s",actset[aa],passet[pp]))
        ty <- reaction_typeFP(actset[aa],passet[pp])
        data$type[data$actseq == actset[aa] & data$passeq == passet[pp]] = ty$rtype
        data$detBind[data$actseq == actset[aa] & data$passeq == passet[pp]] = ty$deterministicBind
        data$detExec[data$actseq == actset[aa] & data$passeq == passet[pp]] = ty$deterministicExec

      }
    }

  # Now do the plot...
  ps <- rconf_params(fn)

  plot(NA,xlim=c(0,ps$gridx+50),ylim=c(0,ps$gridy),asp=1,main=sprintf("T= %d",tt))
  rect(xleft = 0,xright=ps$gridx,ybottom = 0,ytop=ps$gridy,col="black",lty = "blank")

  for(mm in 1:nrow(data)){
  xx <- data$actx[mm]
  yy <- data$acty[mm]
  rc <- rcols[which(rtypes == data$type[mm])]
  rect(xleft = xx,xright = xx+1,ybottom = yy,ytop=yy+1,col=rc,lty = "blank")
  xx <- data$passx[mm]
  yy <- data$passy[mm]
  rect(xleft = xx,xright = xx+1,ybottom = yy,ytop=yy+1,col=rc,lty = "blank")
  legend("topright",pch=15,col = rcols,legend=rtypes)

  }
  }
  else{
  gotdata<-F
  }

  tt <- tt+tstep
}

dev.off()
