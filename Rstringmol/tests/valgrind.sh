#!/bin/bash


# Use this script to run the tests in command-line R using valgrind 

# EXPLANATION:
# This calls R with the debug flag -d, with the string in parentheses after it describing the debug environment
# valgrind is the main debug envifonment, but there are other options within the quotes:
#     --tool=memcheck           - use the memcheck tool
#     --leak-chcek=full         - do a full check for leaks
#     --show-reachable=yes      - report when memory isn't freed that could have been
#     --log-file=~/valgrind.txt - write the output to a file
# -f testthat.R means that R will start and execute the file testthat.R. an alternative would be -e "source('testthat.R')"

R -d "valgrind --tool=memcheck --leak-check=full --show-reachable=yes --log-file=~/valgrind.txt" -f testthat.R 
