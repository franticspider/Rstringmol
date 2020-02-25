

library(Rstringmol)
library(stringr)
source("~/git/Rstringmol/Rstringmol/R/reaction_type.R")








#TODO: put the final functions into Rstringmol properly!
source("../examples/vg_plus_rprops_functions.R")



#parasite_plot("~/Desktop/paulien/smsp/1705smsp/out1/out1_","test_smsp_out1.pdf",title="box-box1",to=100000)
parasite_plot("~/Desktop/paulien/smsp/1705smspr/out5/out1_","smspr_out5c.pdf",title="box-rand5")


#if(!testing){
#box-box
parasite_plot("~/Desktop/paulien/smsp/1705smsp/out1/out1_","csmsp_out1n.pdf",title="box-box1")
parasite_plot("~/Desktop/paulien/smsp/1705smsp/out2/out1_","csmsp_out2n.pdf",title="box-box2")

parasite_plot("~/Desktop/paulien/smsp/1705smsp/out3/out1_","csmsp_out3n.pdf",title="box-box3")

parasite_plot("~/Desktop/paulien/smsp/1705smsp/out5/out1_","csmsp_out5n.pdf",title="box-box5")

#box-rand
parasite_plot("~/Desktop/paulien/smsp/1705smspr/out2/out1_","csmspr_out2n.pdf",title="box-rand2")
parasite_plot("~/Desktop/paulien/smsp/1705smspr/out3/out1_","csmspr_out3n.pdf",title="box-rand3")

#slot-box
parasite_plot("~/Desktop/paulien/smsp/1705sm250/out4/out1_","csm250_out4n.pdf",title="slot-box4")
parasite_plot("~/Desktop/paulien/smsp/1705sm250/out5/out1_","csm250_out5n.pdf",title="slot-box5")

#slot-rand
parasite_plot("~/Desktop/paulien/smsp/1705sm250r/out3/out1_","csm250r_out3n.pdf",title="slot-box4",to=340000)

#evoevo
parasite_plot("~/Desktop/paulien/smsp/evoevo/out1/out1_","cevoevo_out1n.pdf",title="evoevo1",from=1000000)
parasite_plot("~/Desktop/paulien/smsp/evoevo/out2/out1_","cevoevo_out2n.pdf",title="evoevo2",from=1000000)
parasite_plot("~/Desktop/paulien/smsp/evoevo/out3/out1_","cevoevo_out3n.pdf",title="evoevo3",from=1000000)
parasite_plot("~/Desktop/paulien/smsp/evoevo/out4/out1_","cevoevo_out4n.pdf",title="evoevo4",from=1000000)
parasite_plot("~/Desktop/paulien/smsp/evoevo/out5/out1_","cevoevo_out5n.pdf",title="evoevo5",from=1000000)



do_networks <- T

if(do_networks){

  #Do a basic test:
  testinfn <- "~/Desktop/paulien/smsp/1705smspr/out5/out1_40000.conf"
  if(file.exists(testinfn)){
    pdf(file = "testnet.pdf",width = 38,height = 38)
     plot.rnet(testinfn)
     title("\ntest title",cex.main=15)
    dev.off()
  }
  else{
    message(sprintf("File %s doesn't exist!", testinfn))
  }

  #Now let's try a proper run:
  #parasite_plot("~/Desktop/paulien/smsp/1705smsp/out1/out1_","csmsp_out1n.pdf",title="box-box1")
  froot <- "~/Desktop/paulien/smsp/1705smsp/out3/out1_"
  pdf(file="csmsp_out3net.pdf",width=38,height=38)
    for(tt in seq(from = 20000, to = 2000000, by = 20000)){
      message(sprintf("%d",tt))

      fn <- sprintf("%s%d.conf",froot,tt)
      if(file.exists(fn)){
        plot.rnet(fn,onlyreplicators = T)
        title(sprintf("\nt=%d",tt),cex.main=15)
      }
    }
  dev.off()


  #Biofabric alternative
  pdf(file = "testnetbf.pdf",width = 38,height = 38)
  library(RBioFabric)
  bioFabric(net)
  dev.off()

}

plotnetprops <- function(froot,outfn,from = 20000, to = 2000000, by = 20000){

  #Now let's try a proper run:
  #parasite_plot("~/Desktop/paulien/smsp/1705smsp/out1/out1_","csmsp_out1n.pdf",title="box-box1")
  pdf(file=outfn,width=38,height=38)
  #for(tt in seq(from = 20000, to = 2000000, by = 20000)){
  for(tt in seq(from = from, to = to, by = by)){
    message(sprintf("%d",tt))

    fn <- sprintf("%s%d.conf",froot,tt)
    if(file.exists(fn)){
      plot.rnet(fn,onlyreplicators = T)
      title(sprintf("\nt=%d",tt),cex.main=15)
    }
  }
  dev.off()
}


plotnetprops("~/Desktop/paulien/smsp/1705smsp/out5/out1_","csmsp_out5net.pdf")


plotnetprops("~/Desktop/paulien/smsp/1705smspr/out2/out1_","csmspr_out2net.pdf")


# TODO: put a consolidated visual format into here:
#network.plot(testinfn,"smspr_out5_network.pdf")




