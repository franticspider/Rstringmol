

# Unused functions

I used a program called `cppcheck` in linux to scan for unused functions. Before doing this the program size was 68K(!)



# Random numbers

stringmol uses the Mersenne Twister algorithm to generate random numbers. This is to make the code portable, since it doesn't depend on system-specific RNGs. 

# Inexact alignment

By far the largest footprint in memory is the sequence alignment function, which uses the smith waterman algorithm 


# Special symbols in command line
