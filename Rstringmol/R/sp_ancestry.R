

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


vizanc <- function(){
  fn<-"ancestry121.txt"
  ofn <- "out.dot"

  x <- read.table(fn,header=T,sep=",")
  x$depth <- as.integer(x$depth)
  x$time <- as.integer(x$time)
  x$spp <- as.integer(x$spp)
  x$act<- as.integer(x$act)
  x$pass<-as.integer(x$pass)
  x <- x[order(x$time),]
  times <- unique(x$time)


  #DOTFILE HEADER
  ofp <- file(ofn)
  writeLines(c("digraph simple_hierarchy {",
  #"size=\".83,1.17!\" ratio=fill;",
  "nodesep=0.1; sep =0.1; len = 0.1;",
  "overlap = scale"), ofp)
  close(ofp)


  #WRITE THE TIME LINE
  write("",file = ofn, append=T)
  write("/* TIME LINE */",file = ofn, append=T)
  write("{node [shape=plaintext, fontsize=16];",file = ofn, append=T)
  for(i in 1:length(times)){
  	message(sprintf("Writing %d to timeline",times[i]))
  	line <- ""
  	if(i<length(times))
  		line <- sprintf("%d ->",times[i])
  	else
  		line <- sprintf("%d;",times[i])
  	write(line, file = ofn,append=T)
  }
  write("}",file = ofn, append=T)
  write("",file = ofn, append=T)



  #WRITE THE GRAPH NODES
  spnodes <- x[!duplicated(x$spp),]
  write("",file = ofn, append=T)
  write("/* SPECIES */",file = ofn, append=T)
  for(i in 1:nrow(spnodes)){
  	message(sprintf("Writing spp %d at time %d",spnodes$spp[i],spnodes$time[i]))
  	line <- sprintf("{rank=same %d; anc%d [label=\"%d\n%s\",fontsize = 6]}",spnodes$time[i],spnodes$spp[i],spnodes$spp[i],spnodes$seq[i])
  	write(line, file = ofn, append=T)
  }
  write("",file = ofn, append=T)
  message("finished writing nodes")



  #WRITE THE GRAPH EDGES
  spnodes <- as.data.frame(unique(cbind(x$spp,x$act,x$pass)))
  colnames(spnodes) <- c("spp","act","pass")
  write("",file = ofn, append=T)
  write("/* DESCENT */",file = ofn, append=T)
  for(i in 1:nrow(spnodes)){
  	if(spnodes$act[i] > -1){
  		message(sprintf("Link %d: %d + %d -> %d",i,spnodes$act[i],spnodes$pass[i],spnodes$spp[i]))
  		line <- sprintf("\t anc%d -> anc%d [color = red]",spnodes$act[i],spnodes$spp[i])
  		write(line, file = ofn, append=T)
  		line <- sprintf("\t anc%d -> anc%d [color = blue]",spnodes$pass[i],spnodes$spp[i])
  		write(line, file = ofn, append=T)
  	}
  }
  write("",file = ofn, append=T)
  message("finished writing edges")



  #DOTFILE FOOTER
  write(";}",file = ofn, append=T)

  #Use a system call to make the pdf:
  system("dot -Tpdf out.dot -o ancestry1.pdf")

  #Use a system call to make the pdf:
  system("neato -Tpdf out.dot -o ancestry1_neato.pdf")

  #Use a system call to make the pdf:
  system("dot -Tpdf out.dot -o ancestry1.pdf")
}









