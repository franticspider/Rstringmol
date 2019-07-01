#include <Rcpp.h>
using namespace Rcpp;

#include <stdio.h>
#include "rsmData.h"

rsmData::rsmData(){

  sprintf(string0,"");
  sprintf(string1,"");

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

void rsmData::setString0(char *s){
  sprintf(string0,"%s",s);
}

void rsmData::setString1(char *s){
  sprintf(string1,"%s",s);
}

void rsmData::toFile(char *fn){
  FILE *fp;
  fp = fopen(fn,"w");
  if(fp == NULL) return;



}


rsmData::~rsmData(){



}
