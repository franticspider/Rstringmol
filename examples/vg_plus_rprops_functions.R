

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

 pcols <- rainbow(n = nrow(data), alpha = 0.4)
 #pcols <- adjustcolor(pcols, alpha.f = 0.4)


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

	pdf(file=outfn,width = 9, height = 12)
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

		  par(mfrow=c(9,1))
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





# Find a reaction in a list from the active and passive sequences.
# If single, this will return NA or a single row.
lookup.reaction <- function(rlist,actseq,passeq,single = T){

  lookup<-rlist[rlist$actseq == actseq,]
  if(nrow(lookup)>0){
    lookup<-lookup[lookup$passeq == passeq,]
    if(nrow(lookup)>0){
      if(nrow(lookup)==1){
        return(lookup)
      }else{
        if(single){
          message(sprintf("WARNING: %d: >1 row for %s and %s",rr,actseq,passeq))
          return(NULL)
        }
      }
    }
  }
  #message("no rows in lookup")
  return(NULL)
}




network.properties <- function(pwp,mustexist=F){

  pwp$np_MutualRepl <- F
  pwp$np_Parasite <- F
  pwp$np_Hypercycle <- F

  for(rr in 1:nrow(pwp)){

    # Make sure sequences are not the same:
    if(pwp$actseq[rr]!=pwp$passeq[rr]){

      # Declare it as a parasite *unless* it does mutual replication
      if(pwp$pp_Replicator[rr]){

        pwp$np_Parasite[rr]<- T
        #pwp$np_MutualRepl[rr] <- F # Not needed - already set above

        if(mustexist){
          flip <- lookup.reaction(pwp, actseq = pwp$passeq[rr], passeq = pwp$actseq[rr])
        }
        else{
          #TODO: tidy this up!
          flip <- runReactionFP(c(pwp$passeq[rr],pwp$actseq[rr]))
          flur <- setup.pwprops(pwp$passeq[rr],pwp$actseq[rr])
          flip <- calc.pwprops(flur,1,flip)
        }



        if(!is.null(flip)){
          if(flip$pp_Replicator[1]){
            pwp$np_Parasite[rr]<- F
            pwp$np_MutualRepl[rr] <- T

            # Gather the info to detect a hypercycle
            actact <- lookup.reaction(pwp, actseq = pwp$actseq[rr], passeq = pwp$actseq[rr],single=T)
            paspas <- lookup.reaction(pwp, actseq = pwp$passeq[rr], passeq = pwp$passeq[rr],single=T)

            #booleans to hold whether self-copy happens:
            asc <- F
            psc <- F


            if(!is.null(actact)){
              #message(sprintf("%d: actact has %d rows",rr,nrow(actact)))
              if(actact$pp_Replicator[1])
                asc <- T
            }

            if(!is.null(paspas)){
              #message(sprintf("%d: paspas has %d rows",rr,nrow(paspas)))
              if(paspas$pp_Replicator[1])
                psc <- T
            }

            if(!asc && !psc)
              pwp$np_Hypercycle[rr] <- T
          }
        }
      }
    }
  }
  return(pwp)
}


print.props<- function(ur){

  if(ur$pp_obligate)
    message("Reaction is Obligate")
  if(ur$pp_NoProduct)
    message("Reaction has No Product")
  if(ur$pp_NewProduct)
    message("Reaction generates New Product")
  if(ur$pp_Replicator)
    message("Reaction is a Replicator")
  if(ur$pp_SelfReplicator)
    message("Reaction is a Self Replicator")
  if(ur$pp_AutoReplicator)
    message("Reaction Auto Replicates")
  if(ur$pp_Complex)
    message("Reaction is Complex")
  if(ur$pp_Jumper)
    message("Reaction is a Jumper")
  if(ur$np_Hypercycle)
    message("Reaction is a Hypercycle")
  if(ur$np_MutualRepl)
    message("Reaction is a mutual replicator")
  if(ur$np_Parasite)
    message("Reaction is parasitic")

}




calc.pwprops <- function(ur,rr,reaction){
  ur$bprob[rr]  <- reaction$bprob
  ur$nsteps[rr] <- reaction$count

  if(reaction$product == "empty")
    ur$pp_NoProduct[rr] <- T
  else{
    if(reaction$product != ur$actseq[rr] & reaction$product != ur$passeq[rr])
      ur$pp_NewProduct[rr] <- T
  }

  #Replicator types
  if(ur$actseq[rr] == ur$passeq[rr]){
    if(reaction$product == reaction$mPassive){
      ur$pp_SelfReplicator[rr] <- T
    }
  }else{
    if(reaction$product == reaction$mPassive)
      ur$pp_Replicator[rr] <- T
    if(reaction$product == reaction$mActive)
      ur$pp_AutoReplicator[rr] <- T
  }

  #Complex
  if(ur$actseq[rr] != reaction$mActive | ur$passeq[rr] != reaction$mPassive)
    ur$pp_Complex[rr] <- T

  #Jumper
  if(reaction$mActive == "" | reaction$mPassive == "")
    ur$pp_Jumper[rr] <- T


  #################################################
  #ur$type[rr] <- "notyp
  ur$product[rr] <- reaction$product
  ur$dexec[rr]   <- reaction$deterministicExec
  ur$dbind[rr]   <- reaction$deterministicBind

  #TODO: make sure nobs is set!
  #ur$nobs[rr] <- nrow(fdata[fdata$actseq == ur$actseq[rr]
  #                          & fdata$passeq == ur$passeq[rr],])

  #message("next")

  return(ur)
}




setup.pwprops <- function(actseq,passeq){

  ur <- unique(data.frame(actseq = actseq, passeq = passeq, stringsAsFactors = F))
  #ur$actseq <- as.character(ur$actseq)
  #ur$passeq <- as.character(ur$passeq)

  ur$type <- "unknown"
  ur$product <- ""

  #These are the properties agreed with Susan
  ur$pp_obligate <- F

  ur$pp_NoProduct <-F
  ur$pp_NewProduct <-F

  ur$pp_Replicator <-F
  ur$pp_SelfReplicator <-F
  ur$pp_AutoReplicator <-F

  #ur$pp_Hypercycle <-F
  #ur$pp_Parasitic <-F
  ur$pp_Complex <-F
  ur$pp_Jumper <-F

  ur$bprob <-0
  ur$nsteps <- 0

  ur$actlen <- str_length(ur$actseq)
  ur$paslen <- str_length(ur$passeq)

  return(ur)
}


pairwise.properties <- function(fn){

  fdata <- rconf_rdata(fn)
  fdata$actseq <- as.character(fdata$actseq)
  fdata$passeq <- as.character(fdata$passeq)


  #TODO: figure out how to do this for both a file of reactions and a calculated reaction
  #create a data frame of unique reactions for this file:
  ur <- unique(data.frame(actseq = fdata$actseq, passeq = fdata$passeq, stringsAsFactors = F))
  #ur$actseq <- as.character(ur$actseq)
  #ur$passeq <- as.character(ur$passeq)


  ur$type <- "unknown"
  ur$product <- ""

  #These are the properties agreed with Susan
  ur$pp_obligate <- F

  ur$pp_NoProduct <-F
  ur$pp_NewProduct <-F

  ur$pp_Replicator <-F
  ur$pp_SelfReplicator <-F
  ur$pp_AutoReplicator <-F

  #ur$pp_Hypercycle <-F
  #ur$pp_Parasitic <-F
  ur$pp_Complex <-F
  ur$pp_Jumper <-F

  ur$bprob <-0
  ur$nsteps <- 0

  ur$actlen <- str_length(ur$actseq)
  ur$paslen <- str_length(ur$passeq)

  for (rr in 1:nrow(ur)){
    #################################################
    # Reaction properties

    #rt <- reaction_type(ur$actseq[rr],ur$passeq[rr])
    # Run the reaction(s):
    result <- runReactionFP(c(ur$actseq[rr],ur$passeq[rr]))
    # conv <- result
    # if(ur$actseq[rr] != ur$passeq[rr]){
    #   if(result$deterministicBind){
    #     ur$pp_Obligate <- T
    #     conv <- NULL
    #   }
    #   else{
    #     conv <- runReactionFP(c(ur$passeq[rr],ur$actseq[rr]))
    #   }
    # }

    ur$bprob[rr] <- result$bprob
    ur$nsteps[rr] <- result$count

    if(result$product == "empty")
      ur$pp_NoProduct[rr] <- T
    else{
      if(result$product != ur$actseq[rr] & result$product != ur$passeq[rr])
        ur$pp_NewProduct[rr] <- T
    }

    #Replicator types
    if(ur$actseq[rr] == ur$passeq[rr]){
      if(result$product == result$mPassive){
        ur$pp_SelfReplicator[rr] <- T
      }
    }else{
      if(result$product == result$mPassive)
        ur$pp_Replicator[rr] <- T
      if(result$product == result$mActive)
        ur$pp_AutoReplicator[rr] <- T
    }



    # if(!is.null(conv)){
    #   #Replicator / Hypercycle
    #   if((result$product == ur$passeq[rr] & conv$product == ur$actseq[rr])
    #      |
    #      (result$product == ur$actseq[rr] & conv$product == ur$passeq[rr])){
    #     ur$pp_Replicator[rr] <- T
    #     if(ur$actseq[rr] != ur$passeq[rr])
    #       ur$pp_Hypercycle[rr] <- T
    #   }
    #   #Parasite
    #   if(result$product == ur$passeq[rr] & conv$product != ur$actseq[rr])
    #     ur$pp_Parasitic[rr] <- T
    # }
    # else{
    #   if((result$product == ur$passeq[rr])
    #      #     |
    #      #   (result$product == ur$actseq[rr])
    #   ){
    #     ur$pp_Parasitic[rr] <- T
    #   }
    # }

    #message("Complex?")
    #Complex
    if(result$mActive != ur$actseq[rr] | result$mPassive != ur$passeq[rr])
      ur$pp_Complex[rr] <- T

    #Jumper
    if(result$mActive == "" | result$mPassive == "")
      ur$pp_Jumper[rr] <- T


    #message("product")

    #################################################
    #ur$type[rr] <- "notyp
    ur$product[rr] <- result$product
    ur$dexec[rr] <- result$deterministicExec
    ur$dbind[rr] <- result$deterministicBind
    ur$nobs[rr] <- nrow(fdata[fdata$actseq == ur$actseq[rr]
                              & fdata$passeq == ur$passeq[rr],])
    #message("next")
  }


  return(ur)

}




parasite_plot <- function(infroot,outfn,title="",xlim = c(0,100),from = 20000, to = 2000000, by = 20000,type = "Parasite"){

	xmax <- 100
    maxfreq <- 0

	pdf(file=outfn,width = 9, height = 11)
	for(tt in seq(from = from, to = to, by = by)){
	  message(sprintf("%d",tt))

	  fn <- sprintf("%s%d.conf",infroot,tt)
	  if(file.exists(fn)){


	    ur<-pairwise.properties(fn)


		  par(mfrow=c(9,1))
		  par(mar=c(3,5,1,1))
		  afreq <- sm_barplot(ur$actlen,xlim=xlim,xlab="length",ylab="frequency",
                                     title=sprintf("%s,T=%d\nAll Active Spp",title,tt),
                                     ylim=c(0,500))
		  pfreq <- sm_barplot(ur$paslen,xlim=xlim,xlab="length",ylab="frequency",title="All Passive Spp",ylim=c(0,500))
		  message(sprintf("maxfreq was %d, afreq =  %d, pfreq = %d",maxfreq,afreq,pfreq))
		  maxfreq <- max(maxfreq,afreq,pfreq)

      #xlim <- 300

      pur <- rbind(ur[ur$pp_NoProduct,])
      sm_lap_plot(pur,xlim=xlim,title=sprintf("Reactions with no product (n=%d)",sum(pur$nobs)))

		  pur <- rbind(ur[ur$pp_Replicator,])
		  sm_lap_plot(pur,xlim=xlim,title=sprintf("Replicator reactions (n=%d)",sum(pur$nobs)))
		  #smlur <- pur[str_length(pur$passeq)<20,]

		  pur <- ur[ur$pp_SelfReplicator,]
		  sm_lap_plot(pur,xlim=xlim,title=sprintf("Self-Replicator reactions (n=%d)",sum(pur$nobs)))


		  pur <- ur[ur$pp_AutoReplicator,]
		  sm_lap_plot(pur,xlim=xlim,title=sprintf("AutoReplicator reactions (n=%d)",sum(pur$nobs)))

		  #nur <- ur[ur$type!="Parasite",]
		  #sm_lap_plot(nur,xlim=xlim)
		  #title(sprintf("Non-parasitic reactions"))


		  pur <- rbind(ur[ur$pp_NewProduct,])
		  sm_lap_plot(pur,xlim=xlim,title=sprintf("New product (n=%d)",sum(pur$nobs)))

		  pur <- ur[ur$pp_Complex,]
		  sm_lap_plot(pur,xlim=xlim,title=sprintf("Complex reactions (n=%d)",sum(pur$nobs)))



		  pur <- rbind(ur[ur$pp_Jumper,])
		  sm_lap_plot(pur,xlim=xlim,title=sprintf("Jumper (n=%d)",sum(pur$nobs)))


		  rm(ur)
		  #rm(fdata)
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



######################################

require(igraph)
require(Rstringmol)
require(stringr)


radius.from.freq <- function(freq){

  minr <- 3

  return(minr+sqrt((freq)/pi))
}

plot.rnet<- function(fn,onlyreplicators=T){

  pur <- pairwise.properties(fn)
  seqs <- as.data.frame(table(c(pur$actseq,pur$passeq)),stringsAsFactors = F)
  seqs <- seqs[order(seqs$Freq,decreasing = T),]
  seqs <- data.frame(seq = seqs$Var1, freq = seqs$Freq, len = str_length(seqs$Var1),stringsAsFactors = F)
  seqs$name <- sprintf("%03d",1:nrow(seqs))

  edges <- NULL
  for(ss in 1:nrow(seqs)){
    ae <- pur[pur$actseq == seqs$seq[ss],]
    if(nrow(ae)>0){
      ae$from <- seqs$name[ss]
      ae$to <- ""
      for(pp in 1:nrow(ae))
        ae$to[pp] <- seqs$name[seqs$seq == ae$passeq[pp]]

      if(is.null(edges))
        edges<-ae
      else
        edges <- rbind(edges,ae)
    }
  }

  #rearrange the cols:
  edges <- edges[,c(18,19,1:17)]
  seqs <- seqs[,c(4,1:3)]

  #We only want the replicator reactions
  if(onlyreplicators)
    edges<-edges[(edges$pp_Replicator | edges$pp_SelfReplicator | edges$pp_AutoReplicator),]

  # Give the seqs type numbers for colouring
  seqs$type <- 1
  sr <- pur$actseq[pur$pp_SelfReplicator]
  seqs$type[seqs$seq %in% sr] <- 2

  # Give the edges type numbers for colouring
  edges$type <- 1
  edges$type[edges$pp_Replicator] <- 2
  edges$type[edges$pp_SelfReplicator] <- 3
  edges$type[edges$pp_AutoReplicator] <- 4
  edges$type[edges$pp_Jumper] <- 5
  ecols <- c("grey","blue","green","red","pink")

  #Create the net
  net <- graph_from_data_frame(d=edges, vertices=seqs[seqs$name %in% c(edges$from, edges$to),], directed=T)
  l <- layout_with_fr(net,niter = 5000) # Good but nodes are close together - occasional L-shaped layout
  #l <- layout_with_drl(net) # majority of layouts are l-shaped
  #l <- layout_with_lgl(net,root=1)

  #Format the edges
  E(net)$color <- ecols[E(net)$type]
  E(net)$width <- E(net)$nobs

  #Format the vertices
  cols <- c("tomato","blue")
  V(net)$color <- cols[V(net)$type]
  V(net)$size <- 4 + radius.from.freq(V(net)$freq)#V(net)$freq

  #Plot the net with formatting
  plot(net, layout = l, vertex.label.cex = radius.from.freq(V(net)$freq)#(3+V(net)$freq) / 4
      )

}




runproplist <- function(froot,from=20000,to=2000000,step=20000,outfn=NULL){

  rundata <- list()
  dp <- 1
  for(rr in seq(from,to,step)){



    #testinfn <- sprintf("~/Desktop/paulien/smsp/1705smsp/out3/out1_%d.conf",rr)
    testinfn <- sprintf("%s%d.conf",froot,rr)
    message(sprintf("%d: file is %s",rr,testinfn))

    ug <- pairwise.properties(testinfn)

    ng <- network.properties(ug)

    rundata[[dp]] <- ng

    dp <- dp + 1

  }

  if(!is.null(outfn)){
    save(rundata,file=outfn)
  }

  return(rundata)

}
