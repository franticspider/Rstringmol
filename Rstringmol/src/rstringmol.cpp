#include <Rcpp.h>
#include <string.h>

using namespace Rcpp;

//Error codes
#define SM_ERR_BADNSTRINGS  1

//Threshold on classifying an exec_step as deterministic
//0.001 means if exec is the same in 999/1000 times, then we can say it's deterministic
#define EXEC_DET_THR (0.001)


//These header files were for .c files, renamed to .cpp in Rstringmol!
#include "memoryutil.h"
#include "mt19937-2.h"
#include "randutil.h"

#include "stringmanip.h"
#include "alignment.h"

#include "SMspp.h" //Contains the definition of s_ag - the basic species unit

#include "rsmData.h"
#include "instructions.h"

#include "sw_strip.h"

//set in stringPM.cpp
//const unsigned int maxl = 2000;
//const unsigned int maxl0 = maxl+1;
const int granular_1 = 0; //Whether we are doing granular stringmol or not

//function headers for things we've borrowed from stringPM

int     append_ag(s_ag **list, s_ag *ag);
int     rewind_bad_ptrs(s_ag* act);
int     check_ptrs(s_ag* act);
//int     h_pos(s_ag *pag, char head);
//int     hcopy(s_ag *act);
int     cleave(s_ag *act,s_ag **nexthead);
s_ag *  make_ag(int alab, int agct = 0);
s_ag *  make_mol(std::string seq);
float   get_bprob(align *sw);
float   get_sw(s_ag *a1, s_ag *a2, align *sw, swt *blosum, bool usestrip);
bool    set_exec(s_ag *A, s_ag *B, align *sw);
int     unbind_ag(s_ag * pag, char sptype, int update, l_spp *pa, l_spp *pp);
bool    exec_step(s_ag *act, s_ag *pass, swt *blosum, s_ag **nexthead, icount *ic);
void    print_ptr_offset(FILE *fp, char *S, char *p,int F, char c);
void    print_exec(FILE *fp, s_ag *act);
int     free_ag(s_ag *pag);
void    free_swt(swt *table);



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


//Let's see what happens if we do this...
// [[Rcpp::export]]
double mtrand(){
  float rng = rand0to1();
  return (double) rng;
}

////////////////////////////////////////////////////////////////////////

//int stringPM::append_ag(s_ag **list, s_ag *ag){
int             append_ag(s_ag **list, s_ag *ag){
  s_ag *pag;

  //printf("appending: list = %p, ag = %p\n",*list,ag);
  if(*list==NULL){
    *list=ag;
    //printf("appended: list = %p, ag = %p\n",*list,ag);
  }
  else{
    pag = *list;
    while(pag->next != NULL){
      pag = pag->next;
    }
    pag->next = ag;
    ag->prev = pag;
  }

  return 0;
}



//int stringPM::rewind_bad_ptrs(s_ag* act){
int             rewind_bad_ptrs(s_ag* act){

  int plen,alen,pdist;
  char *ps;

  //PUT DANGLING POINTERS AT THE *END* OF THE STRINGS:
  //DO THE PASSIVE POINTERS FIRST:
  plen = strlen(act->pass->S);
  if(plen){
    ps = act->pass->S;

    pdist = act->i[0]-ps;
    if(pdist>plen || pdist<0)
      act->i[0]=ps+plen;

    pdist = act->r[0]-ps;
    if(pdist>plen || pdist<0)
      act->r[0]=ps+plen;

    pdist = act->w[0]-ps;
    if(pdist>plen || pdist<0)
      act->w[0]=ps+plen;

    pdist = act->f[0]-ps;
    if(pdist>plen || pdist<0)
      act->f[0]=ps+plen;
  }
  else{//Toggle everything off this string...
    act->i[0]=act->pass->S;
    act->r[0]=act->pass->S;
    act->w[0]=act->pass->S;
    act->f[0]=act->pass->S;

    act->it=1;
    act->rt=1;
    act->wt=1;
    act->ft=1;
  }


  //DO THE ACTIVE POINTERS NOW
  alen = strlen(act->S);
  ps = act->S;

  if(alen){
    pdist = act->i[1]-ps;
    if(pdist>alen || pdist<0)
      act->i[1]=ps+alen;

    pdist = act->r[1]-ps;
    if(pdist>alen || pdist<0)
      act->r[1]=ps+alen;

    pdist = act->w[1]-ps;
    if(pdist>alen || pdist<0)
      act->w[1]=ps+alen;

    pdist = act->f[1]-ps;
    if(pdist>alen || pdist<0)
      act->f[1]=ps+alen;
  }
  else{//Toggle everything off this string...
    act->i[1]=act->S;
    act->r[1]=act->S;
    act->w[1]=act->S;
    act->f[1]=act->S;

    if(plen){
      act->it=0;
      act->rt=0;
      act->wt=0;
      act->ft=0;
    }
  }

  //TODO: error checking on this!
  return 0;
}



//int stringPM::check_ptrs(s_ag* act){
int             check_ptrs(s_ag* act){

#ifdef VERBOSE
  print_exec(stdout,act,act->pass);
#endif
  //int len,pdist;
  //char *ps;

  //Sort pointers out first - even if there's going to be an error!
  rewind_bad_ptrs(act);

  //Step 1: make sure act and pass *have* strings...
  if(!strlen(act->S)){
#ifdef VERBOSE
    printf("Zero length active string - dissoc\n");
#endif
    return 1;
  }
  if(!strlen(act->pass->S)){
#ifdef VERBOSE
    printf("Zero length passive string - dissoc\n");
#endif
    return 2;
  }
#ifdef VERBOSE
  print_exec(stdout,act,act->pass);
#endif

  return 0;

}













//int stringPM::cleave(s_ag *act){
int             cleave(s_ag *act,s_ag **nexthead){

  int dac = 0,cpy;
  s_ag *c,*pass,*csite;

  pass = act->pass;

  //pick the mol containing the cleave site:
  csite = act->ft?act:pass;

  if(act->f[act->ft]-csite->S < csite->len){

    //1: MAKE THE NEW MOLECULE FROM THE CLEAVE POINT

    //Can't really say what the label is easily - for ECAL, it's always pass
    c = make_ag(pass->label);//,1);

    //Copy the cleaved string to the agent
    char *cs;

    c->S =(char *) malloc(MAXL0*sizeof(char));
    memset(c->S,0,MAXL0*sizeof(char));

    cs = csite->S;
    cpy = strlen(cs);
    //Check that we aren't creating a zero-length molecule:
    if(!cpy){
      printf("WARNING: Zero length molecule being created!\n");
    }

    //Make the parent structure: ALL DONE NOW IN update_lineage
    //c->pp = splist->make_parents(act->spp,pass->spp);

    cpy -= act->f[act->ft]-cs;

    strncpy(c->S,act->f[act->ft],cpy);
    c->len = strlen(c->S);
#ifdef VERBOSE
    printf("String %d created:\n%s\n",c->idx,c->S);
#endif

    //ALL DONE NOW IN update_lineage
    //Fill in the birth certificate:
    //Parents now set in update_lineage
    //c->paspp = act->spp;
    //c->ppspp = pass->spp;

    //Check the lineage
    //update_lineage(c,'C',1,act->spp,pass->spp,act->biomass);
    act->biomass=0; //reset this; we might continue to make stuff!

    //append the agent to nexthead
    append_ag(nexthead,c);

    //2: HEAL THE PARENT

    memset(act->f[act->ft],0,cpy*sizeof(char));

    csite->len = strlen(csite->S);

    if((dac = check_ptrs(act))){
      switch(dac){
      case 1://Destroy active - only append passive
        unbind_ag(pass,'P',1,act->spp,pass->spp);
        //append_ag(&nexthead,pass);
        //free_ag(act);
        break;
      case 2://Destroy passive - only append active
        unbind_ag(act,'A',1,act->spp,pass->spp);
        //append_ag(&nexthead,act);
        //free_ag(pass);
        break;
      case 3://Destroy both
        Rprintf("destroying both molecules - This should never happen\n");
        unbind_ag(act,'A',1,act->spp,pass->spp);
        unbind_ag(pass,'P',1,act->spp,pass->spp);
        //free_ag(act);
        //free_ag(pass);
        break;
      default://This can't be right can it?
        if(act->ft == act->it){
          act->i[act->it]--;
        }
        break;
      }
    }
  }
  if(!dac){
    act->i[act->it]++;
    //dac=-1;
  }

  return dac;
}








//int stringPM::free_ag(s_ag *pag){
int free_ag(s_ag *pag){

  if(pag->S != NULL){
    //printf("destroying agent %d, code = %s\n",pag->idx,pag->S);
    free(pag->S);
  }

  free(pag);
  pag = NULL;

  return 0;
}






/* TODO: Declare 's_ag' without using malloc - need to modify the struct to do this!
 *
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
s_ag * make_ag(int alab, int agct){

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

//This is a common way to create molecules, e.g. in stringPM::load_agents(), although there isn't a function with this name
s_ag * make_mol(std::string seq){
  s_ag * pag;
  pag = make_ag('X');
  pag->S =(char *) malloc(MAXL0*sizeof(char));
  memset(pag->S,0,MAXL0*sizeof(char));
  strncpy(pag->S,seq.c_str(),strlen(seq.c_str()));
  pag->len = strlen(pag->S);
  return pag;
}







//float stringPM::get_bprob(align *sw){
float get_bprob(align *sw){
  float bprob = 0.;
  //This is the old bind prob, with a modifier for short strings:

  //l is the shorter alignment
  int l = sw->e1-sw->s1 < sw->e2-sw->s2 ? sw->e1-sw->s1 : sw->e2-sw->s2;
  if(l<=2)
    bprob=0;
  else{
    //bprob = pow(sw->score,l)/pow(l,l);
    //BRUTAL HACK:
    float s = sw->score<l-1.124? sw->score : l-1.124;
    bprob = s/(l-1.124);
  }

  return bprob;
}







//float stringPM::get_sw(s_ag *a1, s_ag *a2, align *sw){
//This is to be replaced with get_sw_stip()!
float             get_sw(s_ag *a1, s_ag *a2, align *sw, swt *blosum, bool useStrip = false, bool verbose = false){  //}, s_sw *swlist){

  float bprob;
  char *comp;

  comp = string_comp(a1->S);

  if(useStrip)
    bprob = SmithWatermanStrip(comp,a2->S,sw,blosum,verbose);
  else
    bprob = SmithWatermanV2(comp,a2->S,sw,blosum,verbose);

  free(comp);

  bprob = get_bprob(sw);

  return bprob;
}


//void stringPM::set_exec(s_ag *A, s_ag *B, align *sw){
bool set_exec(s_ag *A, s_ag *B, align *sw){

  s_ag *active,*passive;
  int active_idx,passive_idx;


  bool deterministic = true;
  if(sw->s1 == sw->s2)
    deterministic = false;

  /*
   * The following if statement shows that assignment of
   * active and passive roles is deterministic within this function
   *
   * The only way asssignment can be 'random' is if sw->s1 == sw->s2
   * in this case, whichever molecule is 'A' will be the active one,
   * so if A is picked randomly from the current list, that's where
   * the randomness in the assignment comes from
   */
  if(sw->s1>=sw->s2){
    active = A;
    passive = B;
    active_idx = sw->s1;
    passive_idx = sw->s2;
  }
  else{
    active = B;
    passive = A;
    active_idx = sw->s2;
    passive_idx = sw->s1;
  }

  active->status = B_ACTIVE;
  active->exec = NULL;
  active->pass = passive;

  active->biomass = 0;

  passive->status = B_PASSIVE;
  passive->exec = active;
  passive->pass = NULL;



  active->f[0] = active->i[0] = active->r[0] = active->w[0] = &(passive->S[passive_idx]);//&(passive->S[sw->s1]);
  active->f[1] = active->i[1] = active->r[1] = active->w[1] = &(active->S[active_idx]);//&(active->S[sw->s2]);
  active->ft   = active->it   = active->rt   = active->wt = 1;


  passive->f[0] = passive->i[0] = passive->r[0] = passive->w[0] = 0;
  passive->f[1] = passive->i[1] = passive->r[1] = passive->w[1] = 0;
  passive->ft   = passive->it   = passive->rt   = passive->wt = 0;

#ifdef V_VERBOSE
  printf("Bind finished - looks like:\n");
  print_exec(stdout,active,passive);
#endif

  return deterministic;

}




//int stringPM::unbind_ag(s_ag * pag, char sptype, int update, l_spp *pa, l_spp *pp){
int             unbind_ag(s_ag * pag, char sptype, int update, l_spp *pa, l_spp *pp){

  int found=0;
  //int mass=0;

  if(pag->status==B_ACTIVE){
    //mass = pag->biomass;
    pag->biomass = 0;
  }


  pag->status = B_UNBOUND;
  pag->pass = NULL;
  pag->exec = NULL;

  pag->ect=0;

  pag->f[0] = pag->i[0] = pag->r[0] = pag->w[0] = 0;
  pag->f[1] = pag->i[1] = pag->r[1] = pag->w[1] = 0;
  pag->ft   = pag->it   = pag->rt   = pag->wt = 0;

  //Not needed in Rstringmol!:
  //found = update_lineage(pag,sptype,update,pa,pp,mass);
  return found;
}


//int stringPM::exec_step(s_ag *act, s_ag *pass){
bool             exec_step(s_ag *act, s_ag *pass, swt *blosum, s_ag **nexthead, icount *ic){

  char  *tmp;
  int   dac=0;
  //int   safe_append=1;
  float pbprob=1.0;

  switch(*(act->i[act->it])){//*iptr[it]){


  /*************
   *   H_SEARCH  *
   *************/
  case '$'://h-search
    //act->ft = act->it;
    char *cs;
    if(act->ft)
      cs = act->S;
    else
      cs = act->pass->S;
    tmp = HSearch(act->i[act->it],cs,blosum,&(act->it),&(act->ft),&pbprob);
    act->f[act->ft] = tmp;
    act->i[act->it]++;
    break;

    /*
    //put the active flow head on the string where the active ip is:
    act->ft = act->it;
    char *cs;
    if(act->it)
    cs = act->S;
    else
    cs = act->pass->S;
    act->f[act->ft] = HSearch(act->i[act->it],cs,blosum);
    act->i[act->it]++;
    break;
    */

    /*************
    *   P_MOVE  *
    *************/
  case '>':
    tmp=act->i[act->it];
    tmp++;
    ic->c_move++;
    switch(*tmp){
    case 'A':
      act->it = act->ft;
      act->i[act->it] = act->f[act->ft];
      act->i[act->it]++;
      break;
    case 'B':
      act->rt = act->ft;
      act->r[act->rt] = act->f[act->ft];
      act->i[act->it]++;
      break;
    case 'C':
      act->wt = act->ft;
      act->w[act->wt] = act->f[act->ft];
      act->i[act->it]++;
      break;
    default:
      act->it = act->ft;
    act->i[act->it] = act->f[act->ft];
    act->i[act->it]++;
    break;
    }
    break;


    /************
    *   HCOPY  *
    ************/
  case '='://h-copy
    //Rprintf("executing h-copy\n");
    if(hcopy(act,ic)<0){
      unbind_ag(act,'A',1,act->spp,pass->spp);
      unbind_ag(pass,'P',1,act->spp,pass->spp);
    }
    break;


    /************
    *   INC_R  * GRANULAR STRINGMOL ONLY!
    ************/
  case '+'://h-copy
    if(granular_1==1){
      //printf("Incrementing read \n");
      /* Select the modifier */
      tmp=act->i[act->it];
      tmp++;
      switch(*tmp){
      case 'A':
        act->i[act->it]++;
        break;
      case 'B':
        act->r[act->rt]++;
        break;
      case 'C':
        act->w[act->wt]++;
        break;
      default:
        act->f[act->ft]++;
      break;
      }
    }
    act->i[act->it]++;
    break;



    /************
    *  TOGGLE  *
    ************/
  case '^'://p-toggle: toggle active pointer
    tmp=act->i[act->it];
    tmp++;
    ic->c_togg++;
    switch(*tmp){
    case 'A':
      act->it = 1-act->it;
      break;
    case 'B':
      act->rt = 1-act->rt;
      break;
    case 'C':
      act->wt = 1-act->wt;
      break;
    default:
      act->ft = 1-act->ft;
    break;
    }
    act->i[act->it]++;
    break;

    /************
    *  IFLABEL *
    ************/
  case '?'://If-label
    act->i[act->it]=IfLabel(act->i[act->it],act->r[act->rt],act->S,blosum,&pbprob);
    break;


    /************
    *  CLEAVE  *
    ************/
  case '%':
    if((dac = cleave(act,nexthead))){
      //safe_append=0;	//extract_ag(&nowhead,p);
    }
    break;

    /**************
     *  TERMINATE *
     **************/
  case 0:
  case '}'://ex-end - finish execution

//#ifdef V_VERBOSE
    //printf("Unbinding...\n");
    //fflush(stdout);
//#endif
    unbind_ag(act,'A',1,act->spp,pass->spp);
    unbind_ag(pass,'P',1,act->spp,pass->spp);

    break;

  default://Just increment the i-pointer
    act->i[act->it]++;
  break;
  }
#ifdef V_VERBOSE
  printf("Exec step - looks like:\n");
  print_exec(stdout,act,pass);
#endif

/*Not needed in Rstringmol (but vital in stringPM):
  if(safe_append){
    act->ect++;
    append_ag(&nexthead,act);
    append_ag(&nexthead,pass);
  }
  energy--;
*/

  if(pbprob > EXEC_DET_THR && pbprob < (1.0-EXEC_DET_THR))
    return true; //step was nondeterministic
  else
    return false;//step was deterministic
}





//void stringPM::print_ptr_offset(FILE *fp, char *S, char *p,int F, char c){
void             print_ptr_offset(FILE *fp, char *S, char *p,int F, char c){
  int i,n=p-S;
  if(n<0){
    Rprintf("Problem calculating pointer location\n");
    fflush(stdout);
  }
  for(i=0;i<n;i++)
    Rprintf(" ");
  Rprintf("%c\n",F?c-32:c);
}






//void stringPM::print_exec(FILE *fp, s_ag *act, s_ag *pas){
void             print_exec(FILE *fp, s_ag *act){

  s_ag *pas = act->pass;

  //Diagnostics to screen for active:
  if(act->status == B_ACTIVE){
    //Diagnostics to screen for passive:
    if(pas != NULL){
      if(!strlen(pas->S))
        Rprintf("Zero length passive string\n");
      Rprintf("%6d:\n%s\n",pas->idx,pas->S);
      print_ptr_offset(fp,pas->S,act->i[0],1-act->it,'i');
      print_ptr_offset(fp,pas->S,act->f[0],1-act->ft,'f');
      print_ptr_offset(fp,pas->S,act->r[0],1-act->rt,'r');
      print_ptr_offset(fp,pas->S,act->w[0],1-act->wt,'w');
    }

    if(!strlen(act->S))
      Rprintf("Zero length active string\n");
    Rprintf("%6d:\n%s\n",act->idx,act->S);
    print_ptr_offset(fp,act->S,act->i[1],act->it,'i');
    print_ptr_offset(fp,act->S,act->f[1],act->ft,'f');
    print_ptr_offset(fp,act->S,act->r[1],act->rt,'r');
    print_ptr_offset(fp,act->S,act->w[1],act->wt,'w');
  }
  else{
    Rprintf("\"Active\" molecule isn't Active!\n");
  }

/*
  act->len = strlen(act->S);
  pas->len = strlen(pas->S);


  if(pas->len<=(int) maxl){
    //printf("Passive string length = %d\n",pas->len);
  }
  else
    Rprintf("Passive string length = %d - TOO LONG\n",pas->len);

  if(act->len<=(int) maxl){
    //printf("Active  string length = %d\n",act->len);
  }
  else
    Rprintf("Active  string length = %d - TOO LONG\n",act->len);
*/

}





void free_swt(swt *table){
  //commented bits are reverse order allocations from default table..
  int i;

  if(table != NULL){

    for(i=0;i<table->N;i++){
      free(table->adj[i]);
    }
    free(table->adj);

    for(i=0;i<table->N+1;i++){
      free(table->T[i]);
    }
    free(table->T);

    free(table->key);

    free(table);
  }
}







/****************************************
 Procedure: doComplement
 *****************************************/
//' Carry out a Smith-Waterman alignment
//'
//' @param input the input string (Assume error checking in R...)
//' @param verbose whether to print output or not
//' @export
// [[Rcpp::export]]
Rcpp::String doComplement(Rcpp::String input, bool verbose = false) {

  Rprintf("Inside doCmplement\n");

  std::string string0 = input.get_cstring();

  char *instr,*comp;

  instr = (char *) malloc(MAXL0*sizeof(char));
  memset(instr,0,MAXL0*sizeof(char));

  strncpy(instr,string0.c_str(),strlen(string0.c_str()));

  comp = string_comp(instr);

  free(instr);

  Rcpp::String output(comp);
  return(output);

}






/****************************************
 Procedure: doSWAlign
 *****************************************/
//' Carry out a Smith-Waterman alignment
//'
//' @param seqVector the sequence of the two strings, active first, then passive.
//' @export
// [[Rcpp::export]]
List doSWAlign(Rcpp::StringVector seqVector, bool strip = false, bool verbose = false) {

  if(verbose)
    Rprintf("Inside doSWAlign\n");

  s_ag *m0,*m1;

  align sw;
  swt	*blosum;
  blosum = NULL;

  /*
   int match;		// the number of matching characters.
   float score; 	// the score of the match
   float prob;		// the probability of the match - used for determining events based on the score/match
   int s1;			// start of the match in string 1
   int e1;			// end of the match in string 1
   int s2;			// start of the match in string 2
   int e2;			// end of the match in string 2
   */

  List Lresult = List::create(_["status"] = ((String) "none"),
                              _["match"] = 0,
                              _["score"] = 0.0,
                              _["prob"] =  0.0,
                              _["s1"] = -1,
                              _["e1"] = -1,
                              _["s2"] = -1,
                              _["e2"] = -1,
                              _["bprob"] = 0.0
  );
  if(seqVector.length() != 2){
    Rprintf("ERROR: 2 stringmols required, %d given\n",seqVector.length());
    Lresult["status"] = ((String) "bad number of input strings");
    Lresult["errcode"] = SM_ERR_BADNSTRINGS;

    return Lresult;
  }

  blosum =  default_table();

  std::string string0 = Rcpp::as<std::string>(seqVector[0]);
  std::string string1 = Rcpp::as<std::string>(seqVector[1]);


  //create the agents from the strings
  /* To Do this, we need to call 'make_ag' which is a function of the StringPM class
   * However, this function should be moved!
   * So for now, we'll copy the function into the code and just give a warning
   */
  m0 = make_mol(string0.c_str());
  m1 = make_mol(string1.c_str());

  if(verbose){
    Rprintf("Mol 0 has seq %s and length %d\n",m0->S,m0->len);
    Rprintf("Mol 1 has seq %s and length %d\n",m1->S,m1->len);
  }

  //Do the equivalent of stringPM::testbind():
  //run get_sw() to get bind prob - see stringPM::testbind()
  //Carries out the Smith-Waterman alignment and gets the binding probability
  float bprob = get_sw(m0,m1,&sw,blosum,strip,verbose);
  if(verbose){
    Rprintf("Bind probability for these molecules is %f\n",bprob);
    Rprintf("Match runs from (i,j) %d,%d to %d,%d\n",sw.s1,sw.s2,sw.e1,sw.e2);
    Rprintf("Match, Score and Prob values are %d, %0.3f and %0.3f respectively\n",sw.match,sw.score,sw.prob);
  }



  Lresult["bprob"] = (bprob);

  Lresult["status"] = ((String) "Aligned");
  Lresult["match"] = sw.match;
  Lresult["score"] = sw.score;
  Lresult["prob"] =  sw.prob;
  Lresult["s1"] = sw.s1;
  Lresult["e1"] = sw.e1;
  Lresult["s2"] = sw.s2;
  Lresult["e2"] = sw.e2;


  /* Clean up and return */

  //Tidy up the memory:
  free_ag(m0);
  free_ag(m1);

  //align sw; (no pointers in this, so assume it doesn't need deallocating)
  //swt	*blosum;
  free_swt(blosum);

  return Lresult;
}




/****************************************
 Procedure: count_spp
 *****************************************/

int nagents(s_ag *head, int state){
  s_ag *pag;
  int count=0;
  pag = head;
  while(pag!=NULL){
    switch(state){
    case -1:
      count++;
      break;
    default:
      if(pag->status == state)
        count++;
      /* no break */
    }
    pag=pag->next;
  }
  return count;
}


/****************************************
Procedure: doReaction
 WARNING: this function causes segfaults in R if called too many times - must be a memory leak
 in Lresult list - valgrind doesn't seem to spot it, so it'll take some work to track it down
*****************************************/
//' React 2 stringmols together - determine whether the run is deterministic or not
//'
//' @param seqVector the sequence of the two strings, active first, then passive.
//' @export
// [[Rcpp::export]]
List doReaction(Rcpp::StringVector seqVector, bool verbose = false, const int climit = 1000) {

  s_ag *m0,*m1;

  s_ag *product; //This is used to hold any new molecules that are produced.. it's called 'nexthead' in the functions because that's what it's called in stringPM
  product = NULL;

  align sw;
  swt	*blosum;
  blosum = NULL;

  icount ic;
  init_counter(&ic);


  List Lresult = List::create(_["product"] = "empty",
            _["status"] = ((String) "none"),
            _["bprob"] = 0.0,
            _["count"] = 0,
            _["m0status"] = -1,
            _["m1status"] = -1,
            _["mActive"] = ((String)"unset"),
            _["mPassive"] = ((String)"unset"),
            _["errcode"] = 0,
					  _["deterministicBind"] = false,
					  _["deterministicExec"] = true,
					  _["product"] = ((String) ""),
					  _["nprod"] = 0,
					  //
					  _["ccopy"] = 0,
            _["cmove"] = 0,
            _["cover"] = 0,
            _["ctogg"] = 0
				);


  if(seqVector.length() != 2){
    Rprintf("ERROR: 2 stringmols required, %d given\n",seqVector.length());
    Lresult["status"] = "bad number of input strings";
    Lresult["errcode"] = SM_ERR_BADNSTRINGS;

    return Lresult;
  }

  blosum =  default_table();

  std::string string0 = Rcpp::as<std::string>(seqVector[0]);
  std::string string1 = Rcpp::as<std::string>(seqVector[1]);


  //create the agents from the strings
  /* To Do this, we need to call 'make_ag' which is a function of the StringPM class
  * However, this function should be moved!
  * So for now, we'll copy the function into the code and just give a warning
  */
  m0 = make_mol(string0.c_str());
  m1 = make_mol(string1.c_str());

  if(verbose){
    Rprintf("Mol 0 has seq %s and length %d\n",m0->S,m0->len);
    Rprintf("Mol 1 has seq %s and length %d\n",m1->S,m1->len);
  }

  //Do the equivalent of stringPM::testbind():
  //run get_sw() to get bind prob - see stringPM::testbind()
  //Carries out the Smith-Waterman alignment and gets the binding probability
  float bprob = get_sw(m0,m1,&sw,blosum);
  if(verbose)
    Rprintf("Bind probability for these molecules is %f\n",bprob);

  Lresult["bprob"] = (bprob);

  //

  //use set_exec to determine the active and passive strings
  Lresult["deterministicBind"] = set_exec(m0,m1,&sw);

  Lresult["m0status"] = ((int) m0->status);
  Lresult["m1status"] = ((int) m1->status);

  int count = 0;
  //const int climit = 1000;
  //run exec_setp until the reactants dissassociate

  bool NDStep = false;
  Lresult["deterministicExec"] = true;
  while( count <= climit){
    if(verbose)Rprintf("\n========== STEP %d ==========\n",count);
    if(m0->status == B_ACTIVE){
      NDStep = exec_step(m0,m0->pass,blosum,&product,&ic);
      if(verbose)print_exec(stdout,m0);
      Lresult["mActive"]  = ((String) m0->S);
      Lresult["mPassive"] = ((String) m1->S);
    }

    if(m1->status == B_ACTIVE){
      NDStep = exec_step(m1,m1->pass,blosum,&product,&ic);
      if(verbose)print_exec(stdout,m1);
      Lresult["mActive"]  = ((String) m1->S);
      Lresult["mPassive"] = ((String) m0->S);
    }

    if(m0->status == B_UNBOUND || m1->status == B_UNBOUND ){
      if(verbose)Rprintf("Reaction has ended\n");
      break;
    }

    if(NDStep)
      Lresult["deterministicExec"] = false;

    count++;
  }

  if(product != NULL){
    Lresult["product"] = ((String) product->S);
  }
  Lresult["count"] = count;
  Lresult["m0status"] = ((int) m0->status);
  Lresult["m1status"] = ((int) m1->status);

  Lresult["status"] = "finished";
  Lresult["ccopy"] = ic.c_copy;
  Lresult["cmove"] = ic.c_move;
  Lresult["cover"] = ic.c_over;
  Lresult["ctogg"] = ic.c_togg;

  Lresult["nprod"] = nagents(product,-1);

  //Tidy up the memory:
  free_ag(m0);
  free_ag(m1);
  if(product != NULL)free_ag(product);

  //align sw; (no pointers in this, so assume it doesn't need deallocating
  //swt	*blosum;
  free_swt(blosum);


  return Lresult;
}



//////////////////////
// ALTERNATIVE APPROACH: WRITE OUTPUT TO FILE - LOAD IT IN R
//////////////////////

/****************************************
 Procedure: doReactionFP
 *****************************************/
//' React 2 stringmols together - determine whether the run is deterministic or not - write to a file
//'
//' @param seqVector the sequence of the two strings, active first, then passive.
//' @export
// [[Rcpp::export]]
void doReactionFP(Rcpp::StringVector seqVector,  Rcpp::StringVector fnVector, bool verbose = false) {

  s_ag *m0,*m1;

  s_ag *product; //This is used to hold any new molecules that are produced.. it's called 'nexthead' in the functions because that's what it's called in stringPM
  product = NULL;

  align sw;
  swt	*blosum;
  blosum = NULL;

  rsmData result;

  icount ic;
  init_counter(&ic);
  //Rprintf("Testing fnVector stuff\n");


  //Rprintf("fnVector size is %d:\n",fnVector.size());

  for(int ff=0;ff<fnVector.size();ff++){
    std::string fileName = Rcpp::as<std::string>(fnVector[ff]);
    //Rprintf("fnVector[%d] is %s\n",ff,fileName.c_str());
  }

  if(fnVector.size()!=1){
    Rprintf("ERROR: Too many tmpfiles, exiting doReactionFP\n");
    return;
  }


  std::string fileName = Rcpp::as<std::string>(fnVector[0]);

  FILE * fp;
  if((fp=fopen(fileName.c_str(),"w")) == NULL){

    Rprintf("ERROR: unable to open file %s, exiting doReactionFP\n", fileName.c_str());
    return;

  }

  //Rprintf("Writing reaction data to file %s \n",fileName.c_str());


  /*
  Rcpp::List Lresult = Rcpp::List::create(_["product"] = "empty",
                                          _["status"] = "none",
                                          _["bprob"] = 0.0,
                                          _["count"] = 0,
                                          _["m0status"] = -1,
                                          _["m1status"] = -1,
                                          _["mActive"] = "unset",
                                          _["mPassive"] = "unset",
                                          _["errcode"] = 0);
  if(seqVector.length() != 2){
    Rprintf("ERROR: 2 stringmols required, %d given\n",seqVector.length());
    Lresult["status"] = "bad number of input strings";
    Lresult["errcode"] = SM_ERR_BADNSTRINGS;

    return Lresult;
  }
  */

  blosum =  default_table();

  std::string string0 = Rcpp::as<std::string>(seqVector[0]);
  std::string string1 = Rcpp::as<std::string>(seqVector[1]);



  //create the agents from the strings
  /* To Do this, we need to call 'make_ag' which is a function of the StringPM class
  * However, this function should be moved!
  * So for now, we'll copy the function into the code and just give a warning
  */
  m0 = make_mol(string0.c_str());
  m1 = make_mol(string1.c_str());

  result.setString0(m0->S);
  result.setString1(m1->S);
  //fprintf(fp,"string0,%s\n",string0.c_str());
  //fprintf(fp,"string1,%s\n",string1.c_str());

  if(verbose){
    Rprintf("Mol 0 has seq %s and length %d\n",m0->S,m0->len);
    Rprintf("Mol 1 has seq %s and length %d\n",m1->S,m1->len);
  }

  //Do the equivalent of stringPM::testbind():
  //run get_sw() to get bind prob - see stringPM::testbind()
  //Carries out the Smith-Waterman alignment and gets the binding probability
  float bprob = get_sw(m0,m1,&sw,blosum);
  if(verbose)
    Rprintf("Bind probability for these molecules is %f\n",bprob);


  fprintf(fp,"bprob,%0.6f\n",bprob);

  bool dtb = set_exec(m0,m1,&sw);
  bool dtexec = true;

  if(dtb)
    fprintf(fp,"deterministicBind,TRUE\n");
  else
    fprintf(fp,"deterministicBind,FALSE\n");

  fprintf(fp,"m0status,%d\n",m0->status);
  fprintf(fp,"m1status,%d\n",m1->status);


  int count = 0;
  const int climit = 1000;
  //run exec_setp until the reactants dissassociate

  char mActive[MAXL0];
  char mPassive[MAXL0];

  bool NDStep = false;
  while( count <= climit){
    if(verbose)Rprintf("\n========== STEP %d ==========\n",count);
    if(m0->status == B_ACTIVE){
      NDStep = exec_step(m0,m0->pass,blosum,&product,&ic);
      if(verbose)print_exec(stdout,m0);
      //Lresult["mActive"]  = ((String) m0->S);
      //Lresult["mPassive"] = ((String) m1->S);
      sprintf(mActive,"%s",m0->S);
      sprintf(mPassive,"%s",m1->S);
    }

    if(m1->status == B_ACTIVE){
      NDStep = exec_step(m1,m1->pass,blosum,&product,&ic);
      if(verbose)print_exec(stdout,m1);
      //Lresult["mActive"]  = ((String) m1->S);
      //Lresult["mPassive"] = ((String) m0->S);
      sprintf(mActive,"%s",m1->S);
      sprintf(mPassive,"%s",m0->S);
    }

    if(m0->status == B_UNBOUND || m1->status == B_UNBOUND ){
      if(verbose)Rprintf("Reaction has ended\n");
      break;
    }

    if(NDStep){
      //  Lresult["deterministicExec"] = false;
      dtexec = false;
    }

    count++;
  }


  fprintf(fp,"mActive,%s\n",mActive);
  fprintf(fp,"mPassive,%s\n",mPassive);
  if(product != NULL)
    fprintf(fp,"product,%s\n",product->S);
  else
    fprintf(fp,"product,empty\n");

  fprintf(fp,"count,%d\n",count);
  fprintf(fp,"m0status,%d\n",m0->status);
  fprintf(fp,"m1status,%d\n",m1->status);
  if(dtexec)
    fprintf(fp,"deterministicExec,TRUE\n");
  else
   fprintf(fp,"deterministicExec,FALSE\n");

  fprintf(fp,"ccopy,%d\n",ic.c_copy);
  fprintf(fp,"cmove,%d\n",ic.c_move);
  fprintf(fp,"cover,%d\n",ic.c_over);
  fprintf(fp,"ctogg,%d\n",ic.c_togg);

  fprintf(fp,"nprod,%d",nagents(product,-1));

  //Tidy up the memory:
  free_ag(m0);
  free_ag(m1);
  if(product != NULL)free_ag(product);

  //align sw; (no pointers in this, so assume it doesn't need deallocating
  //swt	*blosum;
  free_swt(blosum);


  //return Lresult;
  fclose(fp);


}





//////////////////////
// EXPERIMENTAL MODULE
//////////////////////
/*

class Foo {
public:
  Foo(double x_, double y_, double z_ ):
  x(x_), y(y_), z(z_) {}
  double x;
  double y;
  double get_z() { return z; }
  void set_z( double z_ ) { z = z_; }
private:
  double z;
};

RCPP_MODULE(mod_foo) {
  class_<Foo>( "Foo" )
  .constructor<double,double,double>()
  .field( "x", &Foo::x )
  .field_readonly( "y", &Foo::y )
  .property( "z", &Foo::get_z, &Foo::set_z )
  ;
}
*/



