#include <Rcpp.h>
using namespace Rcpp;


#include "SMspp.h" //Contains the definition of s_ag - the basic species unit

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwelve(NumericVector x) {
  return x * 12;
}



/****************************************
Procedure: doReaction
*****************************************/
//' React 2 stringmols together - determine whether the run is deterministic or not
//'
//' @param seqVector the sequence of the two strings, active first, then passive.
// [[Rcpp::export]]
int doReaction(Rcpp::StringVector seqVector) {

  s_ag *m1,*m2;


  return 0;
}



