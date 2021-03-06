#Rgraphviz is at bioconductor
#source("http://bioconductor.org/biocLite.R")
#biocLite("Rgraphviz")

#require("Rgraphviz")

#Following this at: http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r/rgraphviz/
#> library(Rgraphviz)
#> test.matrix<-matrix(rep(c(0,1,0,0), 9), ncol=6, nrow=6)
#> rownames(test.matrix)<-c("a", "b", "c", "d", "e", "f")
#> colnames(test.matrix)<-c("a", "b", "c", "d", "e", "f")
#> test.matrix
#  a b c d e f
#a 0 0 0 0 0 0
#b 1 0 1 0 1 0
#c 0 0 0 0 0 0
#d 0 1 0 1 0 1
#e 0 0 0 0 0 0
#f 1 0 1 0 1 0
#> am.graph<-new("graphAM", adjMat=test.matrix, edgemode="directed")
#> am.graph
#A graphAM graph with directed edges
#Number of Nodes = 6
#Number of Edges = 9
#> plot(am.graph, attrs = list(node = list(fillcolor = "lightblue"),
#                                edge = list(arrowsize=0.5)))






#This version was an attempt to make the tree directly in R
make_tree <- function(fn){

	#x<- read.table(fn)
	x <- read.table(fn,header=T,sep=",")

	spnames <- unique(x$spp)
	#spnames <- unique(x$seq)


	x$depth <- as.integer(x$depth)
	x$time <- as.integer(x$time)
	x$spp <- as.integer(x$spp)
	x$act<- as.integer(x$act)
	x$pass<-as.integer(x$pass)


	#get a list of unique species
	spp <- unique(x$spp)


	nspp<- length(spp)

	#Create an adjacency table....

	adj.matrix<-matrix(rep(0,nspp*nspp), ncol = nspp, nrow = nspp)

	rownames(adj.matrix)<-spnames
	colnames(adj.matrix)<-spnames

	for(i in 1:nrow(x)){
		message(sprintf("Linking %d to active %d and passive %d",x$spp[i],x$act[i],x$pass[i]))


		rowno = match(x$spp[i],spp)
		colac = match(x$act[i],spp)
		colpa = match(x$pass[i],spp)

		adj.matrix[colac,rowno] = 1
		adj.matrix[colpa,rowno] = 1
	}

	gr <- new("graphAM", adjMat=adj.matrix, edgemode="directed")
	plot(gr, attrs = list(node = list(fillcolor = "lightblue"),edge = list(arrowsize=0.5)))

}


sp_community<- function(){

  fn<-"~/Desktop/paulien/confs/out3/splist100000.dat"

  ofn <- "out.dot"

  x <- read.table(fn,header=T,sep=",")
  x$ACT <- as.integer(x$ACT)
  x$PASS <- as.integer(x$PASS)
  x$COUNT <- as.integer(x$COUNT)


  #DOTFILE HEADER
  ofp <- file(ofn)
  writeLines(c("digraph simple_hierarchy {",
  #"size = \"50,50!\"; ratio=fill;",
  "nodesep=1; sep =0.5; len = 0.1;",
  "overlap = false"), ofp)
  close(ofp)


  min_num_binds = 4
  x <- x[x$COUNT>min_num_binds,]

  #WRITE THE NODES
  ##Create the node list
  sp1 <- array(c(as.integer(x$ACT),as.character(x$SEQA)),c(length(x$ACT),2))
  sp2 <- array(c(as.integer(x$PASS),as.character(x$SEQB)),c(length(x$PASS),2))
  splist <- rbind(sp1,sp2)
  spu <- unique(splist)
  #spu[,1] <- as.integer(spu[,1])
  write("",file = ofn, append=T)
  write("/* SPECIES */",file = ofn, append=T)
  for(i in 1:nrow(spu)){
  	message(sprintf("Writing spp %d with sequence %s",as.integer(spu[i,1]),spu[i,2]))
  	line <- sprintf("{spp%d [label=\"%d\n%s\",fontsize = 6]}",as.integer(spu[i,1]),as.integer(spu[i,1]),spu[i,2])
  	write(line, file = ofn, append=T)
  }
  write("",file = ofn, append=T)
  message("finished writing nodes")





  #WRITE THE GRAPH EDGES
  write("",file = ofn, append=T)
  write("/* DESCENT */",file = ofn, append=T)
  for(i in 1:nrow(x)){
  		cc <- "black"
  		if(x$COUNT[i]>30)
  			cc<- "red"

  		line <- sprintf("\tspp%d -> spp%d [penwidth = %d, color = %s, len = %0.3f]",x$ACT[i],x$PASS[i],
  		#min(15,x$COUNT[i]),
  		floor((min(15,x$COUNT[i])-min_num_binds+1)*1.5),
  		cc,100/x$COUNT[i])
  	if(x$COUNT[i] > min_num_binds){
  		message(sprintf("Writing %s",line))
  		write(line, file = ofn, append=T)
  	}
  }
  write("",file = ofn, append=T)
  message("finished writing edges")



  #DOTFILE FOOTER
  write(";}",file = ofn, append=T)

  message("Writing dot to pdf")
  #Use a system call to make the pdf:

  syscall <- sprintf("neato -Gepsilon=.000000000001 -Gstart=5 -Tpdf out.dot -o community_thr%d.pdf",min_num_binds)
  system(syscall)
}








