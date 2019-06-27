#include <Rcpp.h>
using namespace Rcpp;

#include "rsmData.h"

rsmData::rsmData(){




  char status[20];
  float bprob;
  int count;
  int m0status;
  int m1status;
  char mActive[MAXL0];
  char mPassive[MAXL0];
  int errcode;
  bool deterministicBind;
  bool deterministicExec;
  char product[MAXL0];
}
