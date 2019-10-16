#!/bin/bash

OUTFILE="out.txt"

if [ "$#" -ne 1 ]; then
	echo "No arguments, writing to \"${OUTFILE}\""
else
	echo "Writing to $1"
	OUTFILE=$1
fi

{
	echo "TESTS" 
	
	echo "" 
	echo "Replicator 1"
	echo "./stst \$=?\>\$EBUB^B\$=?\>\$AVBO%}HOB \$=?\>\$EBUB^B\$=?\>\$AVBO%}HOB"
	./stst \$=?\>\$EBUB^B\$=?\>\$AVBO%}HOB \$=?\>\$EBUB^B\$=?\>\$AVBO%}HOB 
	
	echo ""
	echo "Parasite 1"
	echo "./stst \$=?\>\$EBUB^B\$=?\>\$AVBO%}HOB \$=?\>\$B"
	./stst \$=?\>\$EBUB^B\$=?\>\$AVBO%}HOB \$=?\>\$B 
	
	echo ""
	echo "Replicator 2"
	echo "./stst \$=?\>\$EBUB^A\$=?\>\$AVB\$=?\>\$EBUB^B\$=?\>\$AVBO%HOB  \$=?\>\$EBUB^A\$=?\>\$AVB\$=?\>\$EBUB^B\$=?\>\$AVBO%HOB"
	./stst \$=?\>\$EBUB^A\$=?\>\$AVB\$=?\>\$EBUB^B\$=?\>\$AVBO%HOB  \$=?\>\$EBUB^A\$=?\>\$AVB\$=?\>\$EBUB^B\$=?\>\$AVBO%HOB
	
	echo ""
	echo "Replicator 3"
	echo "./stst OYHOB\$=?\>\$EBUB^A\$=?\>\$AVB\$=?\>\$EBUB^B\$=?\>\$AVBO%HOB  BLUBO\$=?\>\$EBUB^A\$=?\>\$AVB\$=?\>\$EBUB^B\$=?\>\$AVBO%HOB"
	./stst OYHOB\$=?\>\$EBUB^A\$=?\>\$AVB\$=?\>\$EBUB^B\$=?\>\$AVBO%HOB  BLUBO\$=?\>\$EBUB^A\$=?\>\$AVB\$=?\>\$EBUB^B\$=?\>\$AVBO%HOB
	
	echo ""
	echo "Replicator 4"
	echo "./stst ./stst \$=?\>\$UO^A\$T\>D\$VOO^BC=?\>\$BLUO%}=YH \$=?\>\$UO^A\$T\>D\$VOO^BC=?\>\$BLUO%}=YH"
	./stst \$=?\>\$UO^A\$T\>D\$VOO^BC=?\>\$BLUO%}=YH \$=?\>\$UO^A\$T\>D\$VOO^BC=?\>\$BLUO%}=YH
} > $OUTFILE 
