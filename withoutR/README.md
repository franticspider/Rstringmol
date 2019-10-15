

# Introduction

The stringmol "reaction" code is now in a self-contained program to allow experiments on different configurations and hardware to be carried out. 

# Tasks

In order to get the program to run on the FPGA, the *main* task is to get the program size down. I think we can do this by: 

- removing/simplifying `alignment.cpp` to use a different sequence alignment algorithm
- removing `memoryutil.cpp` and using standard malloc operations instead.
- removing `randutil.cpp` and `mt19937-2.cpp` and using standard random number generators instead
- replacing the class in `rsmData.cpp` with a couple of simple C functions
- considering reducing the instruction set to 16 instructions: 7 opcodes and 9 n-ops. 

# Compiling and Running the program

The programs are in `.cpp` files to make them easier to integrate with R, but they are mainly written in standard ANSI C. If g++ is installed, then the program should compile easily using `make`. 

The compiled program takes two strings as arguments. A simple reaction can be carried out like this: 

		./stst AAA NNN

This should yield the following output: 

-`in initmyrand, seed is 1 (1)`: states the RNG seed
-`bprob,1.000000`: The two strings are an exact match, so they bind with probability 1
-`deterministicBind,FALSE`: Since the strings align in the same position, it's a coin toss as to which string exectues (this is called the active string, and the non-executing string is called the passive string).
- `m0status,3`: Status flag for debugging
- `m1status,2`: ditto
- `mActive,AAA`: The sequence of the active string when the program terminates
- `mPassive,NNN`: The sequence of the passive string
- `product,empty`: The sequence of any new molecule. It is possible that more than one molecule will be produced. This is currently not recorded. 
- `count,3`: The number of opcodes executed
- `m0status,1`: Status flag for debugging
- `m1status,1`: Status flag for debugging
- `deterministicExec,TRUE`: Whether any stochasticity was detected in the program execution

## Special symbols in command line arguments

Unfortunately, the `$` and `>` symbols have special meaning on the unix command line, so they need to be escaped with a `\` if they appear in the strings! For example the strings `$=?>G^AQC$=?>G^BQC$=?>E$BLGUO%}` and `$=?>G^AQC$=?>G^BQC$=?>E$BLUO%}` have to be input as `\$=?\>G^AQC\$=?\>G^BQC\$=?\>E\$BLGUO\%}` and `\$=?\>G^AQC\$=?\>G^BQC\$=?\>E\$BLUO\%}` like this: 

		./stst \$=?\>G^AQC\$=?\>G^BQC\$=?\>E\$BLGUO\%} \$=?\>G^AQC\$=?\>G^BQC\$=?\>E\$BLUO\%}



# Program size

A quick `ls -l` on the .o and executable files shows: 

```
ls -l *.o
-rw-rw-r-- 1 sjh sjh 25408 Oct 15 08:58 alignment.o
-rw-rw-r-- 1 sjh sjh  4488 Oct 15 08:59 instructions.o
-rw-rw-r-- 1 sjh sjh 18520 Oct 15 08:58 main.o
-rw-rw-r-- 1 sjh sjh  1992 Oct 15 08:59 memoryutil.o
-rw-rw-r-- 1 sjh sjh  3040 Oct 15 08:58 mt19937-2.o
-rw-rw-r-- 1 sjh sjh  2616 Oct 15 08:58 randutil.o
-rw-rw-r-- 1 sjh sjh  2096 Oct 15 09:00 rsmData.o
-rw-rw-r-- 1 sjh sjh  1888 Oct 15 08:59 stringmanip.o
ls -l stst
-rwxrwxr-x 1 sjh sjh 49136 Oct 15 09:00 stst
```

Note that there is currently no optimisation for size in these figures. 

The file `alignment.cpp` is the biggest memory hog. Note that the alignment process builds an NxM array to do a sequence alignment, where N and M are the length of the two strings being compared. The algorithm also references a large table of substitution scores, which is of size Ax(A+1), where A is the alphabet size (currently 33).


## Unused functions

I used a program called `cppcheck` in linux to scan for unused functions. Before doing this the program size was 68K(!) I have now reduced this to 49K. 

Note that I have commented out the unused functions rather than delete them because they may be useful in future when we get to run a full grid -- they may be used in the 'master' program that is running on a conventional CPU for example.

# Symbolic link to library files. 

The initial build of this program uses symbolic links to the source files contained in the R sub-project at `Rstringmol/src`. This allowed the functions to be tested via R. I suggest that 


# Random numbers

stringmol uses the Mersenne Twister algorithm to generate random numbers. This is to make the code portable, since it doesn't depend on system-specific RNGs. 

# Inexact alignment

By far the largest footprint in memory is the sequence alignment function, which uses the smith waterman algorithm 


