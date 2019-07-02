





froot="~/Desktop/paulien/smsp/1705smsp/out3/"


fn <- sprintf("%sout1_1800000.conf",froot)



data <- rconf_rdata(fn,summarize = F,verbose = T)

cf <- read.table(fn,as.is = T, fill = T, sep="\n",comment.char="",stringsAsFactors = F,blank.lines.skip = F)

