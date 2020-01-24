#include <Rcpp.h>
#include <string.h>


#include "alignment.h"

typedef struct {
  int pos_ii;
  int pos_jj;
  float maxscore;
} swPos;



bool procH(char l1, char *s2, int l2, float *H, swt *swT, swPos *swp, bool forward = true){

  int jj,sj,oo = 0;
  float ms,is,ds,pmax=0.;

  float Hjminus1 = 0;
  bool updated = false;

  if(forward){
    for(jj=1;jj<=l2;jj++){
      sj=jj-1;

      //Calculate the three scores:

      //Match/Mismatch score
      ms = Hjminus1 + wT(l1,s2[sj],swT);

      //Deletion score
      ds = H[jj] + wT(l1,0,swT);

      //Insertion score
      is = H[jj-1] + wT(0,s2[sj],swT);

      //record the proxy for H[ii-1][jj] - it'll be needed next time round as H[ii-1][jj-1]!
      //DEBUG: Rprintf("s2[jj] = %c, H[jj] = %0.2f, hjminus1 = %0.2f\n",s2[sj],H[jj],Hjminus1);
      Hjminus1 = H[jj];

      //update H for new row at jj
      H[jj] = 0;
      if(ms > H[jj])H[jj] = ms;
      if(ds > H[jj])H[jj] = ds;
      if(is > H[jj])H[jj] = is;


      //record the posn of the highest score:
      if(H[jj]>swp->maxscore){
        swp->maxscore = H[jj];
        swp->pos_jj=jj;
        updated = true;
      }
    }
  }
  else{
    for(jj=l2;jj>=1;jj--){
      sj=jj+1;

      //Calculate the three scores:

      //Match/Mismatch score
      ms = Hjminus1 + wT(l1,s2[jj-1],swT);

      //Deletion score
      ds = H[jj] + wT(l1,0,swT);

      //Insertion score
      is = H[sj] + wT(0,s2[jj-1],swT);

      //record the proxy for H[ii-1][jj] - it'll be needed next time round as H[ii-1][jj-1]!
      //DEBUG: Rprintf("s2[jj] = %c, H[jj] = %0.2f, hjminus1 = %0.2f\n",s2[sj],H[jj],Hjminus1);
      Hjminus1 = H[jj];

      //update H for new row at jj
      H[jj] = 0;
      if(ms > H[jj]){
        H[jj] = ms;
      }
      if(ds > H[jj])H[jj] = ds;
      if(is > H[jj])H[jj] = is;


      //record the posn of the highest score:
      if(H[jj]>swp->maxscore){
        swp->maxscore = H[jj];
        swp->pos_jj=jj;
        updated = true;
      }
    }
  }

  return updated;
}


swPos procStrip(char *s1, char *s2, align *A, swt *swT, int verbose){

  int i,si,j,sj;   int l1 = strlen(s1);
  int l2 = strlen(s2);
  float *H;

  swPos pos,ros;
  pos.pos_ii = 0;
  pos.pos_jj = 0;
  pos.maxscore = 0.;

  ros.pos_ii = 0;
  ros.pos_jj = 0;
  ros.maxscore = 0.;

  //Set up the array - TODO: declare globally if H is likely to be reused.
  H =(float *) malloc((l2+1)*sizeof(float));
  memset(H,0,(l2+1)*sizeof(float));


  /*R does it like this:
   for(ii in 1:str_length(seq1)){
     oo <- sw_strip(str_sub(seq1,ii,ii),seq2,oo)
     maxl <- max(oo)
     if(maxl > maxsc){
       maxi <- ii
       maxsc <- max(oo)
       maxj <- min(which(oo == maxsc))
     }
   }
   */
  if(verbose){
    /* print the "header row" of the alignment table (=s2) */
    print_header(s2,l2);
  }

  //Process the "forward" matrix to get the end point
  for(i=1;i<=l1;i++){
    si=i-1;

    if(procH(s1[i-1],s2,l2,H,swT,&pos)){
      pos.pos_ii = i;
    }

    if(verbose)print_Hrow(s1[i-1],&(H[1]),l2);
  }

  if(verbose){
    Rprintf("After forward, i,j and score are %d, %d and %0.2f\n",pos.pos_ii,pos.pos_jj,pos.maxscore);

    /* print the "header row" of the alignment table (=s2) */
    print_header(s2,l2);
  }

  //Now process in "reverse" from the end point to get the start point..
  memset(H,0,(l2+1)*sizeof(float));
  for(i=pos.pos_ii;i>0;i--){
    if(procH(s1[i-1],s2,pos.pos_jj,H,swT,&ros,false)){
      ros.pos_ii = i;
    }
    if(verbose)print_Hrow(s1[i-1],&(H[1]),l2);
  }

  if(verbose)Rprintf("After reverse, i,j and score are %d, %d and %0.2f\n",ros.pos_ii,ros.pos_jj,ros.maxscore);

  A->s1 = ros.pos_ii-1;
  A->s2 = ros.pos_jj-1;
  A->e1 = pos.pos_ii;
  A->e2 = pos.pos_jj;
  A->score = pos.maxscore;
  //Note that this is an upper bound on the match, but we don't use it anyway!
  A->match = (A->e1 - A->s1) > (A->e2 - A->s2) ? (A->e1 - A->s1) : (A->e2 - A->s2);
  A->prob = 0.0; //TODO: should really use get_bprob() here - would be clearer!



  free(H);
  return pos;
}



/*
 * THIS VERSION USES A PROPER TRACE BACK MATRIX
 */
int SmithWatermanStrip(char *s1, char *s2, align *A, swt *swT, int verbose){

  //Forward processing to find the end
  procStrip(s1,s2,A,swT,verbose);

  //Find the end, then reverse the strings


  //Reverse processing to find the start


  return 0;
}
