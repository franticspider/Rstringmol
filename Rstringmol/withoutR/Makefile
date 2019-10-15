

CC = gcc
WARN = -Wall


all: stst


stst: main.o alignment.o randutil.o mt19937-2.o \
	rsmData.o memoryutil.o instructions.o stringmanip.o
	${CC} ${WARN} -o stst main.o alignment.o randutil.o mt19937-2.o \
	rsmData.o memoryutil.o instructions.o stringmanip.o \
	-lstdc++ -lm
  
main.o: main.cpp
	${CC} ${WARN} -c main.cpp
	
alignment.o: alignment.cpp alignment.h
	${CC} ${WARN} -c alignment.cpp
	
randutil.o: randutil.cpp randutil.h
	${CC} ${WARN} -c randutil.cpp
	
mt19937-2.o: mt19937-2.cpp mt19937-2.h
	${CC} ${WARN} -c mt19937-2.cpp
	
rsmData.o: rsmData.cpp rsmData.h
	${CC} ${WARN} -c rsmData.cpp
	
memoryutil.o: memoryutil.cpp memoryutil.h
	${CC} ${WARN} -c memoryutil.cpp
	
instructions.o: instructions.cpp instructions.h
	${CC} ${WARN} -c instructions.cpp
	
stringmanip.o: stringmanip.cpp stringmanip.h
	${CC} ${WARN} -c stringmanip.cpp
	
clean:
	rm *.o
	rm stst
