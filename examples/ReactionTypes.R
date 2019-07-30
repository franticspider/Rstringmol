

rtype_plot <- function(fn){

  #fn <- sprintf("%sout1_%d.conf",froot,tt)

  message(sprintf("Working on file %s",fn))

  if(file.exists(fn)){

    #Get the data...
    data <- rconf_rdata(fn,summarize = F, verbose = F)#T)

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

    #TODO: PLOT THE UNBOUND MOLECULES


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
}


#froot="D:/sjh/stringmol/out3smsp/"
#froot="D:/sjh/stringmol/smsp/1705smsp/out5/"
froot="~/Desktop/paulien/smsp/1705smsp/out3/"
froot="~/Desktop/paulien/smsp/1705smsp/out2/"
froot="~/Desktop/paulien/smsp/1705smsp/out1/"


#froot="~/Desktop/paulien/smsp/1705smspr/out3/"
froot="~/Desktop/paulien/smsp/1705smspr/out2/"

#TODO: Need to redo the figure for 'letterbox' arena shapes
#froot="~/Desktop/paulien/smsp/1705sm250r/out3/"

#pdf(file = "~/smsp5types.pdf",width = 12,height = 10)

tstep <- 20000 #sometimes its 10000
tt <-   20000
sps <- list()
nsps <- 1
gotdata<-T
imno = 1
while(gotdata){

<<<<<<< HEAD
  #This is how we'll do it in a function:
  png(file = sprintf("%sFrtypes%07d.png",froot,imno),width = 500,height=500)
  fn <- sprintf("%sout1_%d.conf",froot,tt)
  if(file.exists(fn)){
    rtype_plot(fn)
  }
  dev.off()

  png(file = sprintf("%srtypes%07d.png",froot,imno),width = 500,height=500)
=======
  png(file = sprintf("%srtypes%07d.png",froot,imno),width = 900,height=400)
>>>>>>> 5c6d970f44149d50348c3b84edc4952b07f83af7

  par(mfrow=c(1,2),  bg = "grey40", fg = "white", col.lab = "white", col.axis = "white", col.main = "white")

  fn <- sprintf("%sout1_%d.conf",froot,tt)


  if(file.exists(fn)){

    message(sprintf("Working on file %s",fn))

    #Get the data...
    data <- rconf_rdata(fn, verbose = F,summarize = F)#T)

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

    #TODO: PLOT THE UNBOUND MOLECULES


    # Now do the plot...
    ps <- rconf_params(fn)

<<<<<<< HEAD
    plot(NA,xlim=c(0,ps$gridx+50),ylim=c(0,ps$gridy),asp=1,main=sprintf("T= %d",tt))
=======
    plot(NA,xlim=c(0,125),ylim=c(0,100),asp=1,main=sprintf("T= %d",tt))
>>>>>>> 5c6d970f44149d50348c3b84edc4952b07f83af7
    rect(xleft = 0,xright=ps$gridx,ybottom = 0,ytop=ps$gridy,col="black",lty = "blank")

    for(mm in 1:nrow(data)){
      xx <- data$actx[mm]
      yy <- data$acty[mm]
      rc <- rcols[which(rtypes == data$type[mm])]
      rect(xleft = xx,xright = xx+1,ybottom = yy,ytop=yy+1,col=rc,lty = "blank")
      xx <- data$passx[mm]
      yy <- data$passy[mm]
      rect(xleft = xx,xright = xx+1,ybottom = yy,ytop=yy+1,col=rc,lty = "blank")
<<<<<<< HEAD
      legend("topright",pch=15,col = rcols,legend=rtypes)

    }
=======
    }
    #legend("topright",pch=15,col = rcols,legend=rtypes)


    # NOW THE BARPLOT:
    par(mar=c(10,6,1,0))
    ttt <- table(data$type)
    ttt <- ttt[rtypes]
    barplot(ttt,las=2,col=rcols,ylim=c(0,2000))
>>>>>>> 5c6d970f44149d50348c3b84edc4952b07f83af7
  }
  else{
    gotdata<-F
  }
  dev.off()
  tt <- tt+tstep
  imno <- imno+1
}

#dev.off()
