/*
    This part of dtwKnn is derived from the UCR Suite.  The original license is below.
 */
    
/***********************************************************************/
/************************* DISCLAIMER **********************************/
/***********************************************************************/
/** This UCR Suite software is copyright protected © 2012 by          **/
/** Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,            **/
/** Gustavo Batista and Eamonn Keogh.                                 **/
/**                                                                   **/
/** Unless stated otherwise, all software is provided free of charge. **/
/** As well, all software is provided on an "as is" basis without     **/
/** warranty of any kind, express or implied. Under no circumstances  **/
/** and under no legal theory, whether in tort, contract,or otherwise,**/
/** shall Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,      **/
/** Gustavo Batista, or Eamonn Keogh be liable to you or to any other **/
/** person for any indirect, special, incidental, or consequential    **/
/** damages of any character including, without limitation, damages   **/
/** for loss of goodwill, work stoppage, computer failure or          **/
/** malfunction, or for any and all other damages or losses.          **/
/**                                                                   **/
/** If you do not agree with these terms, then you you are advised to **/
/** not use this software.                                            **/
/***********************************************************************/
/***********************************************************************/

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

#include "deque.h"

/*Initial the queue at the begining step of envelop calculation*/
void deq_init(struct deque *d, int capacity)
{
  d->capacity = capacity;
  d->size = 0;
  d->dq = Calloc(d->capacity,int);
  d->f = 0;
  d->r = d->capacity-1;
}

/* Destroy the queue*/
void deq_destroy(struct deque *d)
{
  Free(d->dq);
}

 /*Insert to the queue at the back*/
void deq_push_back(struct deque *d, int v)
{
  d->dq[d->r] = v;
  d->r--;
  if (d->r < 0)
    d->r = d->capacity-1;
  d->size++;
}

 /*Delete the current (front) element from queue*/
void deq_pop_front(struct deque *d)
{
  d->f--;
  if (d->f < 0)
    d->f = d->capacity-1;
  d->size--;
}

 /*Delete the last element from queue*/
void deq_pop_back(struct deque *d)
{
  d->r = (d->r+1)%d->capacity;
  d->size--;
}

 /*Get the value at the current position of the circular queue*/
int deq_front(struct deque *d)
{
  int aux = d->f - 1;

  if (aux < 0)
    aux = d->capacity-1;
  return d->dq[aux];
}

 /*Get the value at the last position of the circular queueint back(struct deque *d)*/
int deq_back(struct deque *d)
{
  int aux = (d->r+1)%d->capacity;
  return d->dq[aux];
}

/* Check whether or not the queue is empty*/
int deq_empty(struct deque *d)
{
  return d->size == 0;
}
 
/*
 Finding the envelop of min and max value for LB_Keogh
 Implementation idea is intoruduced by Danial Lemire in his paper
 "Faster Retrieval with a Two-Pass Dynamic-Time-Warping Lower Bound", Pattern Recognition 42(9), 2009.
 */
void lower_upper_lemire(double *t, int len, int r, double *l, double *u)
{
  struct deque du, dl;

  deq_init(&du, 2*r+2);
  deq_init(&dl, 2*r+2);

  deq_push_back(&du, 0);
  deq_push_back(&dl, 0);

  for (int i = 1; i < len; i++)
  {
    if (i > r)
    {
      u[i-r-1] = t[deq_front(&du)];
      l[i-r-1] = t[deq_front(&dl)];
    }
    if (t[i] > t[i-1])
    {
      deq_pop_back(&du);
      while (!deq_empty(&du) && t[i] > t[deq_back(&du)])
        deq_pop_back(&du);
    }
    else
    {
      deq_pop_back(&dl);
      while (!deq_empty(&dl) && t[i] < t[deq_back(&dl)])
        deq_pop_back(&dl);
    }
    deq_push_back(&du, i);
    deq_push_back(&dl, i);
    if (i == 2 * r + 1 + deq_front(&du))
      deq_pop_front(&du);
    else if (i == 2 * r + 1 + deq_front(&dl))
      deq_pop_front(&dl);
  }
  for (int i = len; i < len+r+1; i++)
  {
    u[i-r-1] = t[deq_front(&du)];
    l[i-r-1] = t[deq_front(&dl)];
    if (i-deq_front(&du) >= 2 * r + 1)
      deq_pop_front(&du);
    if (i-deq_front(&dl) >= 2 * r + 1)
      deq_pop_front(&dl);
  }
  deq_destroy(&du);
  deq_destroy(&dl);
}

SEXP lower_upper_env_C(SEXP x, SEXP win){
  SEXP l, u;
  int Win = asInteger(win);
  int i, P=0;
  R_len_t nr  = nrows(x);
  R_len_t nc  = ncols(x);
  R_len_t len = nr*nc;
  if(nr<Win*2+1) error(" win must be < half the length of x ");


  double * _x;
  if(!isReal(x))  { 
   x =  PROTECT( coerceVector(x,REALSXP)); P++;
  } 
  _x = REAL(x);

  PROTECT(l     = allocVector(REALSXP,len));       P++;
  PROTECT(u     = allocVector(REALSXP,len));       P++;

  double  *_l   = REAL(l);
  double  *_u   = REAL(u);

  if(nc==1){ 
    lower_upper_lemire(_x, len, Win, _l, _u);
  } else {
    for(int i = 0; i<nc; i++)
      lower_upper_lemire(_x+nr*i, nr, Win, _l+nr*i, _u+nr*i);
    DUPLICATE_ATTRIB(l,x);
    DUPLICATE_ATTRIB(u,x);
  }
  

  SEXP ans, names;
  names = PROTECT(allocVector(STRSXP,2));       P++;
  SET_STRING_ELT(names,0,mkChar("upper")); 
  SET_STRING_ELT(names,1,mkChar("lower"));
  ans = PROTECT(allocVector(VECSXP, 2));        P++;
  SET_VECTOR_ELT(ans, 0, u); 
  SET_VECTOR_ELT(ans, 1, l);
  setAttrib(ans, R_NamesSymbol, names);
  UNPROTECT(P);
  return(ans);
}
 
