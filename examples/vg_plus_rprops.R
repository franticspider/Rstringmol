

library(Rstringmol)
library(stringr)
source("~/git/Rstringmol/Rstringmol/R/reaction_type.R")


sm_barplot <- function(pdata,xlim=NA,ylim=NA,title="",xlab="",ylab=""){

  dt <- as.data.frame(table(pdata))
  dt$pdata <- as.integer(as.character(dt$pdata))

  #if(is.na(xlim))
  #  xlim = range(dt$pdata)

  if(is.na(ylim[1]))
    ylim = range(dt$Freq)

  plot(NA,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,bty="l")
  title(title,line = -2)

  rect(xleft = dt$pdata-0.5,ybottom = 0,xright = dt$pdata+0.5,ytop=dt$Freq)

  return(max(dt$Freq))

}


draw_reaction <- function(data,acty,pasy,pcols){

  ja <- data$actlen #jitter(data$actlen,amount = 0)
  jp <- data$paslen #jitter(data$paslen)
  mindot <- 0.25

  points(x=ja,y=rep(acty,length(data$actlen)),col=pcols,pch=19,cex=(mindot+log(data$nobs)))
  points(x=jp,y=rep(pasy,length(data$actlen)),col=pcols,pch=19,cex=(mindot+log(data$nobs)))
  segments(x0=ja,x1=jp,y0=rep(acty,length(data$actlen)),y1=rep(pasy,length(data$actlen)),col=pcols,lwd = (max(1,(mindot+log(data$nobs))/2)))


}

sm_lap_plot <- function(data,xlim=range(c(data$actlen,data$paslen)),title=""){

 acty <- 20
 pasy <- 10

 data <- data[order(data$actlen),]


 pcols <- adjustcolor(rainbow(n = nrow(data)), alpha.f = 0.4)


  plot(NA,xlim=xlim,ylim=c(5,25),yaxt="n",xlab="length",ylab ="",axes = F)
  axis(1)

  # have to do a for loop otherwise the segments are drawn out of sequence to the points:
  for(dd in 1:nrow(data))  draw_reaction(data[dd,],acty,pasy,pcols[dd])

  text("Active",x=(xlim[1]+xlim[2]/2),y=acty+3)
  text("Passive",x=(xlim[1]+xlim[2]/2),y=pasy-3)
  title(title)
}




parasite_plot_old <- function(infroot,outfn,title="",xlim = c(0,100),from = 20000, to = 2000000, by = 20000,type = "Parasite"){

	xmax <- 100
        maxfreq <- 0

	pdf(file=outfn,width = 9, height = 9)
	#for(tt in seq(from = 20000, to = 2000000, by = 20000)){
	for(tt in seq(from = from, to = to, by = by)){
	  message(sprintf("%d",tt))
	  #fdata <- rconf_rdata(sprintf("~/Desktop/paulien/smsp/1705smsp/out1/out1_%d.conf",tt))

	  fn <- sprintf("%s%d.conf",infroot,tt)
	  if(file.exists(fn)){

		  fdata <- rconf_rdata(fn)

		  ur <- unique(data.frame(actseq = fdata$actseq, passeq = fdata$passeq))
		  ur$actseq <- as.character(ur$actseq)
		  ur$passeq <- as.character(ur$passeq)
		  ur$type <- "unknown"
		  ur$product <- ""

		  ur$actlen <- str_length(ur$actseq)
		  ur$paslen <- str_length(ur$passeq)

		  par(mfrow=c(7,1))
		  par(mar=c(3,5,1,1))
		  afreq <- sm_barplot(ur$actlen,xlim=xlim,xlab="length",ylab="frequency",
                                     title=sprintf("%s,T=%d\nAll Active Spp",title,tt),
                                     ylim=c(0,500))
		  pfreq <- sm_barplot(ur$paslen,xlim=xlim,xlab="length",ylab="frequency",title="All Passive Spp",ylim=c(0,500))
		  message(sprintf("maxfreq was %d, afreq =  %d, pfreq = %d",maxfreq,afreq,pfreq))
		  maxfreq <- max(maxfreq,afreq,pfreq)

		  for (rr in 1:nrow(ur)){
		    rt <- reaction_type(ur$actseq[rr],ur$passeq[rr])
		    #message(sprintf("%d: %s %s - Reaction Type = %s",rr,ur$actseq[rr],ur$passeq[rr],rt$type))
		    ur$type[rr] <- rt$type
		    ur$product[rr] <- rt$product
		    ur$dexec[rr] <- rt$deterministicExec
		    ur$dbind[rr] <- rt$deterministicBind
		    ur$nobs[rr] <- nrow(fdata[fdata$actseq == ur$actseq[rr] & fdata$passeq == ur$passeq[rr],])
		  }


		  #ur$actlen <- as.integer(str_length(pur$actseq))
		  #ur$paslen <- as.integer(str_length(pur$passeq))

		  pur <- rbind(ur[ur$type=="SelfSelfReplicator",] , ur[ur$type=="NonSelfReplicator" ,])
		  sm_lap_plot(pur,xlim=xlim,title=sprintf("Replicator reactions"))
		  smlur <- pur[str_length(pur$passeq)<20,]
                  #if(nrow(smlur>0)){
                  #  message("SHORT REPLICATORS!")
                  #  for(rr in 1:nrow(smlur)){
                  #    print_mols("SHORT",smlur$actseq[rr],smlur$passeq[rr],"","")
                  #  }
                  #}



		  pur <- ur[ur$type=="Parasite",]
		  sm_lap_plot(pur,xlim=xlim,title=sprintf("Parasites"))

		  #nur <- ur[ur$type!="Parasite",]
		  #sm_lap_plot(nur,xlim=xlim)
		  #title(sprintf("Non-parasitic reactions"))

		  pur <- rbind(ur[ur$type=="SelfSelfNoProduct",] ,  ur[ur$type=="NonSelfNoProduct" ,])
		  sm_lap_plot(pur,xlim=xlim,title=sprintf("Reactions with no product"))


		  pur <- ur[ur$type=="Macromutation",]
		  sm_lap_plot(pur,xlim=xlim,title=sprintf("NonCatalytic reactions"))

		  pur <- rbind(ur[ur$type=="SelfSelfDifferentProduct",] , ur[ur$type=="NonSelfDifferentProduct" ,])
		  sm_lap_plot(pur,xlim=xlim,title=sprintf("Reactions with a different product"))



		  rm(ur)
		  rm(fdata)
		  rm(pur)
		  #rm(nur)
	  }
	  else{
	    message(sprintf("file %s not found, exiting...",fn))
            break
	  }

	}
	dev.off()
}







parasite_plot <- function(infroot,outfn,title="",xlim = c(0,100),from = 20000, to = 2000000, by = 20000,type = "Parasite"){

	xmax <- 100
    maxfreq <- 0

	pdf(file=outfn,width = 9, height = 9)
	for(tt in seq(from = from, to = to, by = by)){
	  message(sprintf("%d",tt))
	  
	  fn <- sprintf("%s%d.conf",infroot,tt)
	  if(file.exists(fn)){

		  fdata <- rconf_rdata(fn)

          #create a data frame of unique reactions for this file:
		  ur <- unique(data.frame(actseq = fdata$actseq, passeq = fdata$passeq))
		  ur$actseq <- as.character(ur$actseq)
		  ur$passeq <- as.character(ur$passeq)
		  ur$type <- "unknown"
		  ur$product <- ""
		  
		  
		  ur$pp_obligate <- F
		  
		  
		  
		  
		  

		  ur$actlen <- str_length(ur$actseq)
		  ur$paslen <- str_length(ur$passeq)

		  par(mfrow=c(7,1))
		  par(mar=c(3,5,1,1))
		  afreq <- sm_barplot(ur$actlen,xlim=xlim,xlab="length",ylab="frequency",
                                     title=sprintf("%s,T=%d\nAll Active Spp",title,tt),
                                     ylim=c(0,500))
		  pfreq <- sm_barplot(ur$paslen,xlim=xlim,xlab="length",ylab="frequency",title="All Passive Spp",ylim=c(0,500))
		  message(sprintf("maxfreq was %d, afreq =  %d, pfreq = %d",maxfreq,afreq,pfreq))
		  maxfreq <- max(maxfreq,afreq,pfreq)

		  for (rr in 1:nrow(ur)){
		    rt <- reaction_type(ur$actseq[rr],ur$passeq[rr])
		    #message(sprintf("%d: %s %s - Reaction Type = %s",rr,ur$actseq[rr],ur$passeq[rr],rt$type))
		    ur$type[rr] <- rt$type
		    ur$product[rr] <- rt$product
		    ur$dexec[rr] <- rt$deterministicExec
		    ur$dbind[rr] <- rt$deterministicBind
		    ur$nobs[rr] <- nrow(fdata[fdata$actseq == ur$actseq[rr] & fdata$passeq == ur$passeq[rr],])
		  }


		  #ur$actlen <- as.integer(str_length(pur$actseq))
		  #ur$paslen <- as.integer(str_length(pur$passeq))

		  pur <- rbind(ur[ur$type=="SelfSelfReplicator",] , ur[ur$type=="NonSelfReplicator" ,])
		  sm_lap_plot(pur,xlim=xlim,title=sprintf("Replicator reactions"))
		  smlur <- pur[str_length(pur$passeq)<20,]
                  #if(nrow(smlur>0)){
                  #  message("SHORT REPLICATORS!")
                  #  for(rr in 1:nrow(smlur)){
                  #    print_mols("SHORT",smlur$actseq[rr],smlur$passeq[rr],"","")
                  #  }
                  #}



		  pur <- ur[ur$type=="Parasite",]
		  sm_lap_plot(pur,xlim=xlim,title=sprintf("Parasites"))

		  #nur <- ur[ur$type!="Parasite",]
		  #sm_lap_plot(nur,xlim=xlim)
		  #title(sprintf("Non-parasitic reactions"))

		  pur <- rbind(ur[ur$type=="SelfSelfNoProduct",] ,  ur[ur$type=="NonSelfNoProduct" ,])
		  sm_lap_plot(pur,xlim=xlim,title=sprintf("Reactions with no product"))

		  pur <- ur[ur$type=="Macromutation",]
		  sm_lap_plot(pur,xlim=xlim,title=sprintf("NonCatalytic reactions"))

		  pur <- rbind(ur[ur$type=="SelfSelfDifferentProduct",] , ur[ur$type=="NonSelfDifferentProduct" ,])
		  sm_lap_plot(pur,xlim=xlim,title=sprintf("Reactions with a different product"))

		  rm(ur)
		  rm(fdata)
		  rm(pur)
		  #rm(nur)
	  }
	  else{
	    message(sprintf("file %s not found, exiting...",fn))
            break
	  }

	}
	dev.off()
}






testing <- T
#parasite_plot("~/Desktop/paulien/smsp/1705smsp/out1/out1_","test_smsp_out1.pdf",title="box-box1",to=100000)
parasite_plot("~/Desktop/paulien/smsp/1705smspr/out5/out1_","smspr_out5.pdf",title="box-rand5")

if(!testing){
	#box-box
	parasite_plot("~/Desktop/paulien/smsp/1705smsp/out1/out1_","smsp_out1.pdf",title="box-box1")
	parasite_plot("~/Desktop/paulien/smsp/1705smsp/out2/out1_","smsp_out2.pdf",title="box-box2")
	parasite_plot("~/Desktop/paulien/smsp/1705smsp/out3/out1_","smsp_out3.pdf",title="box-box3")
	parasite_plot("~/Desktop/paulien/smsp/1705smsp/out5/out1_","smsp_out5.pdf",title="box-box5")

	#box-rand
	parasite_plot("~/Desktop/paulien/smsp/1705smspr/out2/out1_","smspr_out2.pdf",title="box-rand2")
	parasite_plot("~/Desktop/paulien/smsp/1705smspr/out3/out1_","smspr_out3.pdf",title="box-rand3")

	#slot-box
	parasite_plot("~/Desktop/paulien/smsp/1705sm250/out4/out1_","sm250_out4.pdf",title="slot-box4")
	parasite_plot("~/Desktop/paulien/smsp/1705sm250/out5/out1_","sm250_out5.pdf",title="slot-box5")

	#slot-rand
	parasite_plot("~/Desktop/paulien/smsp/1705sm250r/out3/out1_","sm250r_out3.pdf",title="slot-box4",to=340000)

	#evoevo
	parasite_plot("~/Desktop/paulien/smsp/evoevo/out1/out1_","evoevo_out1.pdf",title="evoevo1",from=1000000)
	parasite_plot("~/Desktop/paulien/smsp/evoevo/out2/out1_","evoevo_out2.pdf",title="evoevo2",from=1000000)
	parasite_plot("~/Desktop/paulien/smsp/evoevo/out3/out1_","evoevo_out3.pdf",title="evoevo3",from=1000000)
	parasite_plot("~/Desktop/paulien/smsp/evoevo/out4/out1_","evoevo_out4.pdf",title="evoevo4",from=1000000)
	parasite_plot("~/Desktop/paulien/smsp/evoevo/out5/out1_","evoevo_out5.pdf",title="evoevo5",from=1000000)
}

