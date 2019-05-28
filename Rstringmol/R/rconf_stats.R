

####Function 'arguments'
#fn <- sprintf("%sout1_200000.conf",froot)
#verbose <- T

rconf_stats <- function(fn,spfn,verbose=F){

  spldata <-  splist_stats(spfn)


  if(!is.data.frame(spldata)){
    return(NA)
  }


  if(verbose)message("Loaded spldata using splist_stats")


  ####Function 'body'
  reactions <- rconf_rdata(fn,verbose)

  #Get the product molecule from the splist:
  for(rr in 1:nrow(reactions)){

    if(verbose)message(sprintf("Reaction %d seen %d times",rr,reactions$count[rr]))


    #NA reactions - this needs debugging....missing data in config file, possibly after a cleave...
    if(is.na(reactions$actno[rr]) || is.na(reactions$pasno[rr])){

      if(verbose){
        message(sprintf("Reaction %d has NAs  (NA)",rr))
        message(sprintf("   Active  is %d: %s",reactions$actno[rr],reactions$actseq[rr]))
        message(sprintf("   Passive is %d: %s",reactions$pasno[rr],reactions$passeq[rr]))
      }
      reactions$type[rr]= "NA"
      next
    }

    #TODO: we need to fix the spldata structure we get from a restart file..
    if(verbose)message("\n1\n")
    dd <- spldata[(spldata$act == reactions$actno[rr]) & (spldata$pass == reactions$pasno[rr]), ]
    if(nrow(dd)==0){
      if(verbose){message(sprintf("Reaction %d has no product (NP)",rr))
        message(sprintf("   Active  is %d: %s",reactions$actno[rr],reactions$actseq[rr]))
        message(sprintf("   Passive is %d: %s",reactions$pasno[rr],reactions$passeq[rr]))
      }
      reactions$type[rr]= "NP"
      next
    }

    if(verbose)message("\n2\n")
    #if identical species in the reaction
    if(reactions$actno[rr] == reactions$pasno[rr]){
      if(nrow(dd)==1){
        if(dd$spp == reactions$actno[rr]){
          if(verbose){message(sprintf("Reaction %d is an exact self:self replicator (ER)*",rr))
              message(sprintf("   Active  is %d: %s",reactions$actno[rr],reactions$actseq[rr]))
              message(sprintf("   Passive is %d: %s",reactions$pasno[rr],reactions$passeq[rr]))
            }
          reactions$type[rr]= "ER"
          next
        }
        else{
          if(verbose){message(sprintf("Self:self reaction with different product (DP)*",rr))
            message(sprintf("   Active  is %d: %s",reactions$actno[rr],reactions$actseq[rr]))
            message(sprintf("   Passive is %d: %s",reactions$pasno[rr],reactions$passeq[rr]))
          }
          reactions$type[rr]= "DP"
          next
        }
      }
      if(nrow(dd)>1){
        if(verbose){message(sprintf("Reaction %d is a pathological replicator (PR)",rr))
          message(sprintf("   Active  is %d: %s",reactions$actno[rr],reactions$actseq[rr]))
          message(sprintf("   Passive is %d: %s",reactions$pasno[rr],reactions$passeq[rr]))
        }
        reactions$type[rr]= "PR"
        next
      }
    }
    else{
      if(nrow(dd)==1){
        if(dd$spp == reactions$pasno[rr]){
          if(verbose){message(sprintf("Reaction %d is parsitic (PA)*",rr))
            message(sprintf("   Active  is %d: %s",reactions$actno[rr],reactions$actseq[rr]))
            message(sprintf("   Passive is %d: %s",reactions$pasno[rr],reactions$passeq[rr]))
          }
          reactions$type[rr]= "PA"
          next
        }
        else{
          if(verbose){message(sprintf("none-self reaction with different product (DN)*",rr))
            message(sprintf("   Active  is %d: %s",reactions$actno[rr],reactions$actseq[rr]))
            message(sprintf("   Passive is %d: %s",reactions$pasno[rr],reactions$passeq[rr]))
          }
          reactions$type[rr]= "DN"
          next
        }
      }
      if(nrow(dd)>1){
        if(verbose){message(sprintf("Reaction %d is a pathological non-self replicator (PN)",rr))
          message(sprintf("   Active  is %d: %s",reactions$actno[rr],reactions$actseq[rr]))
          message(sprintf("   Passive is %d: %s",reactions$pasno[rr],reactions$passeq[rr]))
        }
        reactions$type[rr]= "PN"
        next
      }
    }
    if(verbose){message(sprintf("OOPS: Reaction %d can't be classified",rr))
      message(sprintf("   Active  is %d: %s",reactions$actno[rr],reactions$actseq[rr]))
      message(sprintf("   Passive is %d: %s",reactions$pasno[rr],reactions$passeq[rr]))
    }
  }

  return (reactions)

}


#We can't do this because we can't see the *product* of a given reaction - we need some history.
# Conclusion: we need to consolidate all the splist files into a 'master' from the end of the run
rconf_stats_nosplist <- function(fn,verbose=F){

  ####Function 'body'
  reactions <- rconf_rdata(fn,verbose)

  #Get the product molecule from the splist:
  for(rr in 1:nrow(reactions)){

    if(verbose)message(sprintf("Reaction %d seen %d times",rr,reactions$count[rr]))


    #NA reactions - this needs debugging....missing data in config file, possibly after a cleave...
    if(is.na(reactions$actno[rr]) || is.na(reactions$pasno[rr])){

      if(verbose){
        message(sprintf("Reaction %d has NAs  (NA)",rr))
        message(sprintf("   Active  is %d: %s",reactions$actno[rr],reactions$actseq[rr]))
        message(sprintf("   Passive is %d: %s",reactions$pasno[rr],reactions$passeq[rr]))
      }
      reactions$type[rr]= "NA"
      next
    }

    #TODO: we need to fix the spldata structure we get from a restart file..
    if(verbose)message("\n1\n")
    #dd <- spldata[(spldata$act == reactions$actno[rr]) & (spldata$pass == reactions$pasno[rr]), ]
    #dd <- reactions[reactions$actno == rea
    if(nrow(dd)==0){
      if(verbose){message(sprintf("Reaction %d has no product (NP)",rr))
        message(sprintf("   Active  is %d: %s",reactions$actno[rr],reactions$actseq[rr]))
        message(sprintf("   Passive is %d: %s",reactions$pasno[rr],reactions$passeq[rr]))
      }
      reactions$type[rr]= "NP"
      next
    }

    if(verbose)message("\n2\n")
    if(reactions$actno[rr] == reactions$pasno[rr]){
      if(nrow(dd)==1){
        if(dd$spp == reactions$actno[rr]){
          if(verbose){message(sprintf("Reaction %d is an exact self:self replicator (ER)*",rr))
            message(sprintf("   Active  is %d: %s",reactions$actno[rr],reactions$actseq[rr]))
            message(sprintf("   Passive is %d: %s",reactions$pasno[rr],reactions$passeq[rr]))
          }
          reactions$type[rr]= "ER"
          next
        }
        else{
          if(verbose){message(sprintf("Self:self reaction with different product (DP)*",rr))
            message(sprintf("   Active  is %d: %s",reactions$actno[rr],reactions$actseq[rr]))
            message(sprintf("   Passive is %d: %s",reactions$pasno[rr],reactions$passeq[rr]))
          }
          reactions$type[rr]= "DP"
          next
        }
      }
      if(nrow(dd)>1){
        if(verbose){message(sprintf("Reaction %d is a pathological replicator (PR)",rr))
          message(sprintf("   Active  is %d: %s",reactions$actno[rr],reactions$actseq[rr]))
          message(sprintf("   Passive is %d: %s",reactions$pasno[rr],reactions$passeq[rr]))
        }
        reactions$type[rr]= "PR"
        next
      }
    }
    else{
      if(nrow(dd)==1){
        if(dd$spp == reactions$pasno[rr]){
          if(verbose){message(sprintf("Reaction %d is parsitic (PA)*",rr))
            message(sprintf("   Active  is %d: %s",reactions$actno[rr],reactions$actseq[rr]))
            message(sprintf("   Passive is %d: %s",reactions$pasno[rr],reactions$passeq[rr]))
          }
          reactions$type[rr]= "PA"
          next
        }
        else{
          if(verbose){message(sprintf("none-self reaction with different product (DN)*",rr))
            message(sprintf("   Active  is %d: %s",reactions$actno[rr],reactions$actseq[rr]))
            message(sprintf("   Passive is %d: %s",reactions$pasno[rr],reactions$passeq[rr]))
          }
          reactions$type[rr]= "DN"
          next
        }
      }
      if(nrow(dd)>1){
        if(verbose){message(sprintf("Reaction %d is a pathological non-self replicator (PN)",rr))
          message(sprintf("   Active  is %d: %s",reactions$actno[rr],reactions$actseq[rr]))
          message(sprintf("   Passive is %d: %s",reactions$pasno[rr],reactions$passeq[rr]))
        }
        reactions$type[rr]= "PN"
        next
      }
    }
    if(verbose){message(sprintf("OOPS: Reaction %d can't be classified",rr))
      message(sprintf("   Active  is %d: %s",reactions$actno[rr],reactions$actseq[rr]))
      message(sprintf("   Passive is %d: %s",reactions$pasno[rr],reactions$passeq[rr]))
    }
  }

  return (reactions)

}
