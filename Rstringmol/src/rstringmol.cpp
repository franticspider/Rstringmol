#include <Rcpp.h>
using namespace Rcpp;

//Error codes
#define SM_ERR_BADNSTRINGS  1

//These header files were for .c files, renamed to .cpp in Rstringmol!
#include "memoryutil.h"
#include "mt19937-2.h"
#include "randutil.h"

#include "stringmanip.h"
#include "alignment.h"

#include "SMspp.h" //Contains the definition of s_ag - the basic species unit

//set in stringPM.cpp
unsigned int maxl0 = 2001;



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


/* TODO: Move this outside stringPM..
 * Create an 'agent', which is a string 'molecule'
 * NB: The string is not allocated here - done outside the function
 *
 * parameters:
 *
 * alab: a single-character label for the string (e.g. 'C')
 *
 * original version:
 *
 *  s_ag * stringPM::make_ag(int alab)
 *
 *  extra program arguments have been added to remove the need to use the stringPM class
 */
s_ag * make_ag(int alab, int agct = 0){

  s_ag *ag;

  //printf("Spatial make_ag called\n");fflush(stdout);

  if((ag = (s_ag *) mymalloc(1,sizeof(s_ag)))!=NULL){
    ag->label=alab;
    ag->next = NULL;
    ag->prev = NULL;
    ag->exec = NULL;
    ag->pass = NULL;
    ag->S = NULL;
    ag->spp = NULL;
    ag->status = B_UNBOUND;
    ag->idx = agct++;
    ag->nbind=0;
    ag->ect=0;
    ag->biomass=0;
    ag->x=-1;
    ag->y=-1;
    return ag;
  }
  else{
    printf("mymalloc error\n");fflush(stdout);
    getchar();
    return NULL;
  }
}

//This is a common way to create molecules, e.g. in stringPM::load_agents()
s_ag * make_mol(std::string seq){
  s_ag * pag;
  pag = make_ag('X');
  pag->S =(char *) malloc(maxl0*sizeof(char));
  memset(pag->S,0,maxl0*sizeof(char));
  strncpy(pag->S,seq.c_str(),strlen(seq.c_str()));
  pag->len = strlen(pag->S);
  return pag;
}



//float stringPM::get_sw(s_ag *a1, s_ag *a2, align *sw){
float get_sw(s_ag *a1, s_ag *a2, align *sw,
             s_sw *swlist){

  float bprob;
  char *comp;
  s_sw *swa;
  swt	*blosum;

  blosum =  default_table();



  //SUGGEST: pass in pointer to the species - not its index
  swa = read_sw(swlist,a1->spp->spp,a2->spp->spp);

  if(swa==NULL){

    //get_string_comp(a1);
    comp = string_comp(a1->S);

    bprob = SmithWatermanV2(comp,a2->S,sw,blosum,0);
    //bprob = SmithWaterman(comp,a2->S,sw,blosum,0);

    free(comp);

    align sw2;

    bprob = SmithWatermanV2(a1->S,a2->S,&sw2,blosum,0);

    //TODO: SUGGEST: pass in pointer to the species - not its index
    store_sw(&swlist,sw,a1->spp->spp,a2->spp->spp);
  }
  else{
    load_sw(swa,sw);
  }

  bprob = get_bprob(sw);

  //if(verbose_bind){
  //  printf("Alignment:\nm1: %d to %d\nm2: %d to %d\nscore = %f\nProb = %f = %E\n",sw->s1,sw->e1,sw->s2,sw->e2,sw->score,bprob,bprob);
  //}

  return bprob;
}





/****************************************
Procedure: doReaction
*****************************************/
//' React 2 stringmols together - determine whether the run is deterministic or not
//'
//' @param seqVector the sequence of the two strings, active first, then passive.
//' @export
// [[Rcpp::export]]
Rcpp::List doReaction(Rcpp::StringVector seqVector) {

  s_ag *m0,*m1;
  align sw;

  Rcpp::List Lresult = Rcpp::List::create(Rcpp::Named("product") = "empty", _["status"] = "none", _["errcode"] = 0);
  if(seqVector.length() != 2){
    Rprintf("ERROR: 2 stringmols required, %d given\n",seqVector.length());
    Lresult["status"] = "bad number of input strings";
    Lresult["errcode"] = SM_ERR_BADNSTRINGS;
    return Lresult;
  }

  std::string string0 = Rcpp::as<std::string>(seqVector[0]);
  std::string string1 = Rcpp::as<std::string>(seqVector[1]);


  //create the agents from the strings
  /* To Do this, we need to call 'make_ag' which is a function of the StringPM class
  * However, this function should be moved!
  * So for now, we'll copy the function into the code and just give a warning
  */
  Rprintf("REMINDER: make_ag function is a member of stringPM but shouldn't be!\n");
  m0 = make_mol(string0.c_str());
  m1 = make_mol(string1.c_str());

  Rprintf("Mol 0 has seq %s and length %d\n",m0->S,m0->len);
  Rprintf("Mol 1 has seq %s and length %d\n",m1->S,m1->len);

 //run get_sw() to get bind prob - see stringPM::testbind()
 //Carries out the Smith-Waterman alignment and gets the binding probability
 float bprob = get_sw(m0,m1,&sw);


  //use set_exec to determine the active and passive strings



  Lresult["status"] = "finished";


  return Lresult;
}



