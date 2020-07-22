/* Copyright (C) 2009-2012 Simon Hickinbotham                           */
/* When you use this, send an email to: sjh436@gmail.com                */
/* with an appropriate reference to your work.                          */

/* This file is part of STRINGMOL										*/

/* STRINGMOL is free software: you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */

/* This program is distributed in the hope that it will be useful,      */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */

/* You should have received a copy of the GNU General Public License    */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>.*/



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <Rcpp.h>
using namespace Rcpp;

#include "memoryutil.h"
#include "randutil.h"

//string stuff
#include "stringmanip.h"
#include "alignment.h"

#include "SMspp.h"
#include "rsmData.h"

#include "instructions.h"

//extern const int  maxl;

//Use this to control whether h-search is stochastic or not
#define SOFT_SEARCH



//Forward function declarations
int     h_pos(s_ag *pag, char head);


////////counting functions

void init_counter(icount *ct){
  ct->c_copy=0;
  ct->c_move=0;
  ct->c_over=0;
  ct->c_togg=0;
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
int hcopy(s_ag *act, icount * ct){

  //s_ag *pass;
  //pass = act->pass;
  //int cidx;
  //float rno;
  int safe = 1;// this gets set to zero if any of the tests fail..

  act->len = strlen(act->S);
  act->pass->len = strlen(act->pass->S);

  int p;

  //If wptr is too far, just increment iptr
  if( (p = h_pos(act,'w'))>=(int) MAXL){

    //#ifdef DODEBUG
    Rprintf("Write head out of bounds: %d\n",p);
    //#endif

    //just to make sure no damage is done:
    if(act->wt)
      act->S[MAXL]='\0';
    else
      act->pass->S[MAXL]='\0';

    act->i[act->it]++;

    safe = 0;
    return -1;
  }

    //If rptr is too far, just increment iptr
  if(h_pos(act,'r')>=(int) MAXL){
    Rprintf("Read head out of bounds\n");

    act->i[act->it]++;

    safe = 0;
    return -2;
  }

  //If rptr is not pointing at an opcode
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

    if(ct !=NULL){
      ct->c_copy++;
      if(*(act->w[act->wt])!=0){
        ct->c_over++;
      }
    }

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





int LabLength(char *ip, const int maxl){

	int len=0;
	ip = ip+1;
	while(*ip > 64 && *ip <91){
		len++;
		ip++;
	}
	if(len>maxl){
		printf("Label = %d, longer than maxl (= %d!!\n",len,maxl);
	}
	return len;
}


////////////////////
// $: H-Search    //
////////////////////
char * HSearch(char *iptr, char *sp, swt *T, int *itog, int *ftog, float *pbprob){

	char *ip,*tp,tmp[MAXL];
	ip = iptr;
	int i,len=0;
	align A;

	//create space for the match
	memset(tmp,0,MAXL*sizeof(char));

	/*NOTE: We are currently searching from the start of the string with the active flow pointer.
	Perhaps we should start AT the flow pointer, and "loop around" to the beginning of the string if
	there is no match in the first part.

	So the string:

	ABCDEFGHIJKLMNOPQRSTUVWXYZ
	          f

	would be serached as if it was:

	JKLMNOPQRSTUVWXYZABCDEFGHI

	..the best match along this line would be returned. Note it is possible that we could position <f> at the
	start of the line using this technique with a little modification. But probably better to implement a "decrement"
	operator */

	/*
	//first get the length of the string
	while(*ip > 64 && *ip <91){
		len++;
		ip++;
	}
	*/

	//get the length of the label
	len = LabLength(ip, MAXL);

	//set the default positions if there's no match:
	tp = iptr+len;
	ip=iptr+1;

  //special case if len is zero - different behaviour to if len is 1 or 2, which is also different to if len is >2
	if(!len){
		//Ensure that the toggles are set:
		*ftog = *itog;
		return iptr;
	}

	memset(tmp,0,128*sizeof(char));
	strncpy(tmp,iptr+1,len);
	//generate the complement:
	for(i=0;i<len;i++)
		tmp[i] = AlphaComp(tmp[i]);

	//SmithWaterman(tmp,sp,&A,T,0);
	float bprob = SmithWatermanV2(tmp,sp,&A,T,0);

#ifndef SOFT_SEARCH
	//TODO: this will always match if any symbols match. There is no stochastic element..
	//match will only be zero if NO symbols match
	if(A.match){
#else
	int l = A.e1-A.s1 < A.e2-A.s2 ? A.e1-A.s1 : A.e2-A.s2;

	//if l is too short, bprob is zero which means the soft search will never match
	if(l<=2){
	  //printf("l<=2, so bprob = 0\n");
		bprob=0;
	}
	else{
	  //printf("A.score = %f; l = %d\n",A.score,l);
		bprob = pow(A.score,l)/pow(l,l);
	}

	float s = A.score<l-1.124? A.score : l-1.124;
	bprob = s/(l-1.124);
	*pbprob = bprob;
	//printf("Hsearch: s = %f;  bprob = %f\n",s,bprob);

	float rno = rand0to1();
	if(rno<bprob){//search success!
#endif
		return tp + A.e2 - (tp-sp);
	}

	// if 'l' is 1 or 2, we *always* get to here
	//Ensure that the toggles are set - if no match found, we currenty move *F to *I - might leave it on the opposite string if it was there in the 1st place...:
	*ftog = *itog;
	return tp;
}



////////////////////
// ?: If-Label    //
////////////////////
char * IfLabel(char *ip, char *rp, char *sp, swt *T, float *iprob){

	char tmp[MAXL],tmp2[MAXL];
	int i,len = LabLength(ip, MAXL);
	ip++;
	align A;

        align_init(&A);

	switch(len){

	case 0:
		if(!*rp)
			return ip+len+1;
		return ip+len;

		break;


	case 1: //Possibly switch to look at different heads here....not implemented yet....
		if(!*rp)
			return ip+len+1;
		return ip+len;
		break;

	default: //use an alignment to move the pointers

		memset(tmp ,0,MAXL*sizeof(char));
		memset(tmp2,0,MAXL*sizeof(char));
		strncpy(tmp,ip,len);
		//generate the complement:
		for(i=0;i<len;i++)
			tmp[i] = AlphaComp(tmp[i]);

		strncpy(tmp2,rp,len);

		//SmithWaterman(tmp,tmp2,&A,T,0);
		SmithWatermanV2(tmp,tmp2,&A,T,0);

		bool ae = align_event(&A,len);
                *iprob = A.prob;
		if(ae)
		  return ip+len+1;
		else
		  return ip+len;
		//original formulation:
		//if(align_event(&A,len))
		//	return ip+len+1;
		//return ip+len;
		break;
	}

	//TODO: We should never reach this point, but check:
	return 0;
}

