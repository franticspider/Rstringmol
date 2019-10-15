
///////////////////////////////////////////////////////////////////////////////
// Standard libs:
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


///////////////////////////////////////////////////////////////////////////////
// Hash defines: 
//Error codes
#define SM_ERR_BADNSTRINGS  1

//Threshold on classifying an exec_step as deterministic
//0.001 means if exec is the same in 999/1000 times, then we can say it's deterministic
#define EXEC_DET_THR (0.001)


///////////////////////////////////////////////////////////////////////////////
// project libs:
//These header files were for .c files, renamed to .cpp in Rstringmol!
#include "memoryutil.h"
#include "randutil.h"
#include "SMspp.h"
#include "alignment.h"
#include "rsmData.h"
#include "instructions.h"
#include "stringmanip.h"

///////////////////////////////////////////////////////////////////////////////
// CONSTANTS
//set in stringPM.cpp
const unsigned int maxl = 2000;
const unsigned int maxl0 = maxl+1;
const int granular_1 = 0; //Whether we are doing granular stringmol or not


///////////////////////////////////////////////////////////////////////////////
//function headers for things we've borrowed from the stringPM class

int     append_ag(s_ag **list, s_ag *ag);
int     rewind_bad_ptrs(s_ag* act);
int     check_ptrs(s_ag* act);
int     h_pos(s_ag *pag, char head);
int     hcopy(s_ag *act);
int     cleave(s_ag *act,s_ag **nexthead);
s_ag *  make_ag(int alab, int agct = 0);
s_ag *  make_mol(char * seq);
float   get_bprob(align *sw);
float   get_sw(s_ag *a1, s_ag *a2, align *sw, swt *blosum);
bool    set_exec(s_ag *A, s_ag *B, align *sw);
int     unbind_ag(s_ag * pag, char sptype, int update, l_spp *pa, l_spp *pp);
bool    exec_step(s_ag *act, s_ag *pass, swt *blosum, s_ag **nexthead);
void    print_ptr_offset(FILE *fp, char *S, char *p,int F, char c);
void    print_exec(FILE *fp, s_ag *act);
int     free_ag(s_ag *pag);
void    free_swt(swt *table);


///////////////////////////////////////////////////////////////////////////////

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






//int stringPM::h_pos(s_ag *pag, char head){
int h_pos(s_ag *pag, char head){

  char *ph;
  char *ps;

  if(pag->status != B_ACTIVE)
    printf("ERROR: attempting head position for inactive string");

  switch(head){
  case 'w':
    ph = pag->w[pag->wt];
    if(pag->wt)
      ps = pag->S;
    else
      ps = pag->pass->S;
    break;
  case 'f':
    ph = pag->f[pag->ft];
    if(pag->ft)
      ps = pag->S;
    else
      ps = pag->pass->S;
    break;
  case 'i':
    ph = pag->i[pag->it];
    if(pag->it)
      ps = pag->S;
    else
      ps = pag->pass->S;
    break;
  case 'r':
    ph = pag->r[pag->rt];
    if(pag->rt)
      ps = pag->S;
    else
      ps = pag->pass->S;
    break;
  }

  return ph-ps;

}




//int stringPM::hcopy(s_ag *act){
int hcopy(s_ag *act){

  //s_ag *pass;
  //pass = act->pass;
  //int cidx;
  //float rno;
  int safe = 1;// this gets set to zero if any of the tests fail..

  act->len = strlen(act->S);
  act->pass->len = strlen(act->pass->S);

  int p;
  if( (p = h_pos(act,'w'))>=(int) maxl){

    //#ifdef DODEBUG
    printf("Write head out of bounds: %d\n",p);
    //#endif

    //just to make sure no damage is done:
    if(act->wt)
      act->S[maxl]='\0';
    else
      act->pass->S[maxl]='\0';

    act->i[act->it]++;

    safe = 0;
    return -1;
  }
  if(h_pos(act,'r')>=(int) maxl){
    printf("Read head out of bounds\n");

    act->i[act->it]++;

    safe = 0;
    return -2;
  }



  if(*(act->r[act->rt]) == 0){
    //possibly return a negative value and initiate a b
    safe = 0;
    //return -3;
  }

  //if(h_pos(act,'w')>=maxl){
  //	//We are at the end of a copy...
  //	//...so just increment *R
  //	act->r[act->rt]++;
  //	safe = 0;
  //}

  if(safe){
    //MUTATION: (NOT NEEDED IN RSTRINGMOL ATM)
    /*
    rno=rand0to1();
    if(rno<indelrate){//INDEL

      //should follow the blosum table for this....
      rno=rand0to1();
      if(rno<0.5){//insert
        //first do a straight copy..
        *(act->w[act->wt])=*(act->r[act->rt]);

        //no need to test for granular here since we are inserting...
        act->w[act->wt]++;

        //Then pick a random instruction:
        cidx = (float) rand0to1() * blosum->N;

        //insert the random instruction
        *(act->w[act->wt])=blosum->key[cidx];
        if(granular_1==0){
          act->w[act->wt]++;
        }
      }
      else{//delete
        act->i[act->it]++;
      }

      if(granular_1==0){
        act->r[act->rt]++;
      }

    }
    else{
      if(rno<subrate+indelrate){//INCREMENTAL MUTATION
        cidx = sym_from_adj(*(act->r[act->rt]),blosum);
        *(act->w[act->wt])=cidx;
      }
      else{//NO MUTATION
        *(act->w[act->wt])=*(act->r[act->rt]);
      }
      if(granular_1==0){
        act->w[act->wt]++;
        act->r[act->rt]++;
      }
    }
     */


    //If no mutation, it's easy:
    *(act->w[act->wt])=*(act->r[act->rt]);
    act->w[act->wt]++;
    act->r[act->rt]++;



  }
  //update lengths
  act->len = strlen(act->S);
  act->pass->len = strlen(act->pass->S);




  //Increment the instruction pointer
  act->i[act->it]++;


#ifdef VERBOSE
  if(mut)
    printf("Mutant event %d. new string is:\n%s\n\n",mut,act->wt?act->S:act->pass->S);
#endif
  act->biomass++;
  ///biomass++;
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

    c->S =(char *) malloc(maxl0*sizeof(char));
    memset(c->S,0,maxl0*sizeof(char));

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
        printf("destroying both molecules - This should never happen\n");
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
s_ag * make_mol(char * seq){
  s_ag * pag;
  pag = make_ag('X');
  pag->S =(char *) malloc(maxl0*sizeof(char));
  memset(pag->S,0,maxl0*sizeof(char));
  strncpy(pag->S,seq,strlen(seq));
  pag->len = strlen(pag->S);
  return pag;
}







//float stringPM::get_bprob(align *sw){
float get_bprob(align *sw){
  float bprob = 0.;
  //This is the old bind prob, with a modifier for short strings:

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
float             get_sw(s_ag *a1, s_ag *a2, align *sw, swt *blosum){//}, s_sw *swlist){

  float bprob;
  char *comp;
  s_sw *swa;

  //SUGGEST: pass in pointer to the species - not its index
  //this attempts to load a cached sw alignment - which won't happen here!
  //swa = read_sw(swlist,a1->spp->spp,a2->spp->spp);
  swa = NULL;

  if(swa==NULL){

    //get_string_comp(a1);
    comp = string_comp(a1->S);

    bprob = SmithWatermanV2(comp,a2->S,sw,blosum,0);
    //bprob = SmithWaterman(comp,a2->S,sw,blosum,0);

    free(comp);

    align sw2;

    bprob = SmithWatermanV2(a1->S,a2->S,&sw2,blosum,0);

    //TODO: SUGGEST: pass in pointer to the species - not its index
    //store_sw(&swlist,sw,a1->spp->spp,a2->spp->spp);
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
bool             exec_step(s_ag *act, s_ag *pass, swt *blosum, s_ag **nexthead){

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
    tmp = HSearch(act->i[act->it],cs,blosum,&(act->it),&(act->ft),maxl,&pbprob);
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
    //printf("executing h-copy\n");
    if(hcopy(act)<0){
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
    act->i[act->it]=IfLabel(act->i[act->it],act->r[act->rt],act->S,blosum,maxl,&pbprob);
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
    printf("Problem calculating pointer location\n");
    fflush(stdout);
  }
  for(i=0;i<n;i++)
    printf(" ");
  printf("%c\n",F?c-32:c);
}






//void stringPM::print_exec(FILE *fp, s_ag *act, s_ag *pas){
void             print_exec(FILE *fp, s_ag *act){

  s_ag *pas = act->pass;

  //Diagnostics to screen for active:
  if(act->status == B_ACTIVE){
    //Diagnostics to screen for passive:
    if(pas != NULL){
      if(!strlen(pas->S))
        printf("Zero length passive string\n");
      printf("%6d:\n%s\n",pas->idx,pas->S);
      print_ptr_offset(fp,pas->S,act->i[0],1-act->it,'i');
      print_ptr_offset(fp,pas->S,act->f[0],1-act->ft,'f');
      print_ptr_offset(fp,pas->S,act->r[0],1-act->rt,'r');
      print_ptr_offset(fp,pas->S,act->w[0],1-act->wt,'w');
    }

    if(!strlen(act->S))
      printf("Zero length active string\n");
    printf("%6d:\n%s\n",act->idx,act->S);
    print_ptr_offset(fp,act->S,act->i[1],act->it,'i');
    print_ptr_offset(fp,act->S,act->f[1],act->ft,'f');
    print_ptr_offset(fp,act->S,act->r[1],act->rt,'r');
    print_ptr_offset(fp,act->S,act->w[1],act->wt,'w');
  }
  else{
    printf("\"Active\" molecule isn't Active!\n");
  }

/*
  act->len = strlen(act->S);
  pas->len = strlen(pas->S);


  if(pas->len<=(int) maxl){
    //printf("Passive string length = %d\n",pas->len);
  }
  else
    printf("Passive string length = %d - TOO LONG\n",pas->len);

  if(act->len<=(int) maxl){
    //printf("Active  string length = %d\n",act->len);
  }
  else
    printf("Active  string length = %d - TOO LONG\n",act->len);
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

















int main(int argc, char *argv[]){



  //Additional structures to make this main work like rstringmol.cpp:
  bool verbose = false; 
  FILE *fp;
  fp = stdout;

  if(verbose)printf("Hello!\n");fflush(stdout);

  initmyrand(1);
  
  s_ag *m0,*m1;
  s_ag *product; //This is used to hold any new molecules that are produced.. it's called 'nexthead' in the functions because that's what it's called in stringPM
  product = NULL;

  align sw;
  swt	*blosum;
  blosum = NULL;

  rsmData result;


  if(verbose)printf("Making blosum table\n");fflush(stdout);

  blosum =  default_table();
  
  //create the agents from the strings
  /* To Do this, we need to call 'make_ag' which is a function of the StringPM class
  * However, this function should be moved!
  * So for now, we'll copy the function into the code and just give a warning
  */
  if(verbose)printf("Making molecules\n");fflush(stdout);
  m0 = make_mol(argv[1]);
  m1 = make_mol(argv[2]);

  result.setString0(m0->S);
  result.setString1(m1->S);

  if(verbose){
    printf("Mol 0 has seq %s and length %d\n",m0->S,m0->len);
    printf("Mol 1 has seq %s and length %d\n",m1->S,m1->len);
  }

  //Do the equivalent of stringPM::testbind():
  //run get_sw() to get bind prob - see stringPM::testbind()
  //Carries out the Smith-Waterman alignment and gets the binding probability
  float bprob = get_sw(m0,m1,&sw,blosum);
  if(verbose)
    printf("Bind probability for these molecules is %f\n",bprob);


  fprintf(fp,"bprob,%0.6f\n",bprob);

  bool dtb = set_exec(m0,m1,&sw);
  bool dtexec = true;

  if(dtb)
    fprintf(fp,"deterministicBind,TRUE\n");
  else
    fprintf(fp,"deterministicBind,FALSE\n");

  fprintf(fp,"m0status,%d\n",m0->status);
  fprintf(fp,"m1status,%d\n",m1->status);


  /*
  Lresult["bprob"] = (bprob);


  //use set_exec to determine the active and passive strings
  Lresult["deterministicBind"] = set_exec(m0,m1,&sw);

  Lresult["m0status"] = ((int) m0->status);
  Lresult["m1status"] = ((int) m1->status);

  Lresult["deterministicExec"] = true;
  */





  int count = 0;
  const int climit = 1000;
  //run exec_setp until the reactants dissassociate

  char mActive[MAXL0];
  char mPassive[MAXL0];

  bool NDStep = false;
  while( count <= climit){
    if(verbose)printf("\n========== STEP %d ==========\n",count);
    if(m0->status == B_ACTIVE){
      NDStep = exec_step(m0,m0->pass,blosum,&product);
      if(verbose)print_exec(stdout,m0);
      //Lresult["mActive"]  = ((String) m0->S);
      //Lresult["mPassive"] = ((String) m1->S);
      sprintf(mActive,"%s",m0->S);
      sprintf(mPassive,"%s",m1->S);
    }

    if(m1->status == B_ACTIVE){
      NDStep = exec_step(m1,m1->pass,blosum,&product);
      if(verbose)print_exec(stdout,m1);
      //Lresult["mActive"]  = ((String) m1->S);
      //Lresult["mPassive"] = ((String) m0->S);
      sprintf(mActive,"%s",m1->S);
      sprintf(mPassive,"%s",m0->S);
    }

    if(m0->status == B_UNBOUND || m1->status == B_UNBOUND ){
      if(verbose)printf("Reaction has ended\n");
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
