#include <Rcpp.h>
using namespace Rcpp;

#include <stdio.h>
#include "rsmData.h"

rsmData::rsmData(){


  sprintf(status,"none");
  bprob = 0;
  count = 0;
  m0status = -1;
  m1status = -1;
  sprintf(mActive,"unset");
  sprintf(mPassive,"unset");
  errcode = 0;
  deterministicBind = true; //we have to prove it's false!
  deterministicExec = true; //we have to prove it's false!
  sprintf(product,"unset");
}
