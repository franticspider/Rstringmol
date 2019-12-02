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
	echo "1: Replicator 1"
	echo "./stst \$=?\>\$EBUB^B\$=?\>\$AVBO%}HOB \$=?\>\$EBUB^B\$=?\>\$AVBO%}HOB"
	./stst \$=?\>\$EBUB^B\$=?\>\$AVBO%}HOB \$=?\>\$EBUB^B\$=?\>\$AVBO%}HOB 
	
	echo ""
	echo "2: Parasite 1"
	echo "./stst \$=?\>\$EBUB^B\$=?\>\$AVBO%}HOB \$=?\>\$B"
	./stst \$=?\>\$EBUB^B\$=?\>\$AVBO%}HOB \$=?\>\$B 
	
	echo ""
	echo "3: Replicator 2"
	echo "./stst \$=?\>\$EBUB^A\$=?\>\$AVB\$=?\>\$EBUB^B\$=?\>\$AVBO%HOB  \$=?\>\$EBUB^A\$=?\>\$AVB\$=?\>\$EBUB^B\$=?\>\$AVBO%HOB"
	./stst \$=?\>\$EBUB^A\$=?\>\$AVB\$=?\>\$EBUB^B\$=?\>\$AVBO%HOB  \$=?\>\$EBUB^A\$=?\>\$AVB\$=?\>\$EBUB^B\$=?\>\$AVBO%HOB
	
	echo ""
	echo "4: Replicator 3"
	echo "./stst OYHOB\$=?\>\$EBUB^A\$=?\>\$AVB\$=?\>\$EBUB^B\$=?\>\$AVBO%HOB  BLUBO\$=?\>\$EBUB^A\$=?\>\$AVB\$=?\>\$EBUB^B\$=?\>\$AVBO%HOB"
	./stst OYHOB\$=?\>\$EBUB^A\$=?\>\$AVB\$=?\>\$EBUB^B\$=?\>\$AVBO%HOB  BLUBO\$=?\>\$EBUB^A\$=?\>\$AVB\$=?\>\$EBUB^B\$=?\>\$AVBO%HOB
	
	echo ""
	echo "5: Replicator 4"
	echo "./stst ./stst \$=?\>\$UO^A\$T\>D\$VOO^BC=?\>\$BLUO%}=YH \$=?\>\$UO^A\$T\>D\$VOO^BC=?\>\$BLUO%}=YH"
	./stst \$=?\>\$UO^A\$T\>D\$VOO^BC=?\>\$BLUO%}=YH \$=?\>\$UO^A\$T\>D\$VOO^BC=?\>\$BLUO%}=YH
	
	
	echo ""
	echo "6: Long reaction with no product"
	./stst \$=?\>$EBUB^A\$=?\>\$AVB\$=?\>\$EBUB^B\$=?\>\$AVB}%}HOB\$=?\>\$EBUB^A\$=?\>\$AVB\$=?\>\$EBUB^B\$=?\>\$AVB \$=?\>\$EBUB^A\$=?\>\$AVB\$=?\>\$EBUB^B\$=?\>\$AVBO%}HOBP 
} > $OUTFILE 
