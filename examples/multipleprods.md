---
title: "R Notebook"
output: html_notebook
---

```{r}
require(Rstringmol)
source("../examples/vg_plus_rprops_functions.R")


pr <- function(rundata){
  
  for(rr in 1:length(rundata)){
    d<-rundata[[rr]]
    message(sprintf("\nMax nprod for run %d is %d",rr*20000,max(rundata[[rr]]$nprod)))
    md <- d[d$nprod>1,]
    if(nrow(md>0)){
      for(kk in 1:nrow(md)){
        message(sprintf("%d\nact=%s\npas=%s\n",kk,md$actseq[kk],md$passeq[kk]))    
      }
    }
  }
  
}

load("rp1705smsp1.RData");pr(rundata)

#rundata <- runproplist("~/Desktop/paulien/smsp/1705smsp/out2/out1_",outfn = "rp1705smsp2.RData")
#rundata <- runproplist("~/Desktop/paulien/smsp/1705smsp/out3/out1_",outfn = "rp1705smsp3.RData")
#rundata <- runproplist("~/Desktop/paulien/smsp/1705smsp/out5/out1_",outfn = "rp1705smsp5.RData")

#rundata <- runproplist("~/Desktop/paulien/smsp/1705smspr/out2/out1_",outfn = "rp1705smspr2.RData")
#rundata <- runproplist("~/Desktop/paulien/smsp/1705smspr/out3/out1_",outfn = "rp1705smspr3.RData")

#rundata <- runproplist("~/Desktop/paulien/smsp/1705sm250/out4/out1_",outfn = "rp1705sm2504.RData")
#rundata <- runproplist("~/Desktop/paulien/smsp/1705sm250/out5/out1_",outfn = "rp1705sm2505.RData")
```


