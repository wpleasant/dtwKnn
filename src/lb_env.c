
/*
 *   Author: William Pleasant 
 *   
 *   Based on "ascending minima" algorithm. original code by Richarda Harter na svetu.
 *   http:richardhartersworld.com/cri/2001/slidingmin.html 
 *   
 *   lb_env doubles the speed and replaces Danial Lemire's veresion "lower_upper_lemire" 
 *   finding the upper and lower envelope for the LB_Keogh. 
 *
 */


#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

#include "lb_env.h"


void lb_env(double *in, double *lo, double *hi, int n, int k)
{   
  int i,j,ii;
  int kk = k*2+1;
  if(n<1) error("n must be > 0");

  /* structs  */
  struct pairs_double * ring_lo;
  struct pairs_double * minpair;
  struct pairs_double * end_lo;
  struct pairs_double * last_lo;

  struct pairs_double * ring_hi;
  struct pairs_double * maxpair;
  struct pairs_double * end_hi;
  struct pairs_double * last_hi;


  /* init lower env */
  ring_lo = malloc(kk * sizeof *ring_lo);
  if (!ring_lo) error("malloc error");
  end_lo  = ring_lo + kk;
  last_lo = ring_lo;
  minpair = ring_lo;
  minpair->value = in[0];
  minpair->death = kk;


  /* init upper env */
  ring_hi = malloc(kk * sizeof *ring_hi);
  if (!ring_hi) error("malloc error");
  end_hi  = ring_hi + kk;
  last_hi = ring_hi;
  maxpair = ring_hi;
  maxpair->value = in[0];
  maxpair->death = kk;

  /* start and main window  */
  ii = 0;
  for (i=1;i<=n+k;i++) {
    if(ii<n-1) ii++;
    if(i>k){
      lo[i-k-1] = minpair->value;
      hi[i-k-1] = maxpair->value;
    }

    /* lower */
    if (minpair->death == i) {
      minpair++;
      if (minpair >= end_lo) minpair = ring_lo;
    }
    if (in[ii] <= minpair->value) {
      minpair->value = in[ii];
      minpair->death = i+kk;
      last_lo = minpair;
    } else {
      while (last_lo->value >= in[ii]) {
        if (last_lo == ring_lo) last_lo = end_lo;
        --last_lo;
      }
      ++last_lo;
      if (last_lo == end_lo) last_lo = ring_lo;
      last_lo->value = in[ii];
      last_lo->death = i+kk;
    }

    /* upper */
    if (maxpair->death == i) {
      maxpair++;
      if (maxpair >= end_hi) maxpair = ring_hi;
    }
    if (in[ii] >= maxpair->value) {
      maxpair->value = in[ii];
      maxpair->death = i+kk;
      last_hi = maxpair;
    } else {
      while (last_hi->value <= in[ii]) {
        if (last_hi == ring_hi) last_hi = end_hi;
        --last_hi;
      }
      ++last_hi;
      if (last_hi == end_hi) last_hi = ring_hi;
      last_hi->value = in[ii];
      last_hi->death = i+kk;
    }
  }
  free(ring_lo);
  free(ring_hi);
}

SEXP lb_env_C(SEXP x, SEXP win, SEXP dropAttr){
  SEXP l,u;
  int Win = asInteger(win);
  int i, P=0;
  R_len_t nr  = nrows(x);
  R_len_t nc  = ncols(x);
  R_len_t len = nr*nc;
  int drop = asLogical(dropAttr);

  double * _x, *_l, *_u;
  if(!isReal(x))  { 
    x =  PROTECT( coerceVector(x,REALSXP)); P++;
  } 
  _x = REAL(x);

  if(nc>1 && drop) {
    PROTECT(l     = allocMatrix(REALSXP,nr,nc));     P++;
    PROTECT(u     = allocMatrix(REALSXP,nr,nc));     P++;
  } else { 
    PROTECT(l     = allocVector(REALSXP,len));       P++;
    PROTECT(u     = allocVector(REALSXP,len));       P++;
  }

  _l   = REAL(l);
  _u   = REAL(u);

  /*for(int i = 0; i<nc*nr; i++) { // set NA ???
    _l[i] = NA_REAL;
    _u[i] = NA_REAL;
    }*/


  /* Work here.  apply to each column of matrix if any */
  for(int i = 0; i<nc; i++){
    lb_env(_x+nr*i, _l+nr*i, _u+nr*i, nr, Win);
  }

  if(!drop){
    DUPLICATE_ATTRIB(l,x);
    DUPLICATE_ATTRIB(u,x);
  }
  /* returns a named list */
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



