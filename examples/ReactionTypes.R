
#This reloads bacollite...
devtools::test()
require(png)


#froot="D:/sjh/stringmol/out3smsp"
#froot="D:/sjh/stringmol/smsp/1705smsp/out5"
froot="~/Desktop/paulien/smsp/1705smsp/out5"
#froot="~/Desktop/paulien/smsp/1705smsp/out4"
#froot="~/Desktop/paulien/smsp/1705smsp/out3"
#froot="~/Desktop/paulien/smsp/1705smsp/out2"
#froot="~/Desktop/paulien/smsp/1705smsp/out1"


#froot="~/Desktop/paulien/smsp/1705smspr/out3"
#froot="~/Desktop/paulien/smsp/1705smspr/out2"



#TODO: Need to redo the figure for 'letterbox' arena shapes
#froot="~/Desktop/paulien/smsp/1705sm250r/out3"

#pdf(file = "~/smsp5types.pdf",width = 12,height = 10)

tstep <- 20000 #sometimes its 10000
tt <-   20000
sps <- list()
nsps <- 1
gotdata<-T
imno = 1
while(gotdata){

  png(file = sprintf("%s/rtypes%07d.png",froot,imno),width = 900,height=800)

  par(mfrow=c(2,2),  bg = "grey40", fg = "white", col.lab = "white", col.axis = "white", col.main = "white")

  fn <- sprintf("%s/out1_%d.conf",froot,tt)


  if(file.exists(fn)){

    message(sprintf("Working on file %s",fn))

    #Get the data...
    data <- rconf_rdata(fn, verbose = F,summarize = F)#T)

    actset <- unique(data$actseq)
    if(length(actset)>0){
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
    }
    #TODO: PLOT THE UNBOUND MOLECULES


    # Now do the plot...
    ps <- rconf_params(fn)

    plot(NA,xlim=c(0,125),ylim=c(0,100),asp=1,main=sprintf("Reactions at T= %d",tt))
    rect(xleft = 0,xright=ps$gridx,ybottom = 0,ytop=ps$gridy,col="black",lty = "blank")

    for(mm in 1:nrow(data)){
      xx <- data$actx[mm]
      yy <- data$acty[mm]
      rc <- rcols[which(rtypes == data$type[mm])]
      rect(xleft = xx,xright = xx+1,ybottom = yy,ytop=yy+1,col=rc,lty = "blank")
      xx <- data$passx[mm]
      yy <- data$passy[mm]
      rect(xleft = xx,xright = xx+1,ybottom = yy,ytop=yy+1,col=rc,lty = "blank")
    }
    #legend("topright",pch=15,col = rcols,legend=rtypes)


    # NOW THE BARPLOT:
    oldpar <- par(mar=c(10,6,1,0))
    ttt <- table(data$type)
    ttt <- ttt[rtypes]
    barplot(ttt,las=2,col=rcols,ylim=c(0,2000))


    # Empty plot for Length image
    par(oldpar) 
    plot(NA,xlim=c(0,125),ylim=c(0,100),asp=1,main=sprintf("Length at T= %d",tt))
    rect(xleft = 0,xright=ps$gridx,ybottom = 0,ytop=ps$gridy,col="black",lty = "blank")
    img<- readPNG(sprintf("%spng/lenframe%07d.png",froot,tt))
    rasterImage(img,0,ps$gridy,ps$gridx,0,interpolate=F)


    # Empty plott for spno image
    par(oldpar) 
    plot(NA,xlim=c(0,125),ylim=c(0,100),asp=1,main=sprintf("Species at T= %d",tt))
    rect(xleft = 0,xright=ps$gridx,ybottom = 0,ytop=ps$gridy,col="black",lty = "blank")
    img<- readPNG(sprintf("%spng/sppframe%07d.png",froot,tt))
    rasterImage(img,0,ps$gridy,ps$gridx,0,interpolate=F)

  }
  else{
    gotdata<-F
  }
  dev.off()
  tt <- tt+tstep
  imno <- imno+1
}

#dev.off()
