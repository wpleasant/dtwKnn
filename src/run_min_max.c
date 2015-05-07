

/* The "ascending minima" algorithm. original code by Richarda Harter na svetu.
 * http://richardhartersworld.com/cri/2001/slidingmin.html 
 *
 * Adapted by William Pleasant  2015-05-05 
 */


#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

#include "run_min_max.h"

void minwindow(double *in, double *out, int n, int k)
{
  int i;
  struct pairs_d * ring;
  struct pairs_d * minpair;
  struct pairs_d * end;
  struct pairs_d * last;

  ring = malloc(k * sizeof *ring);
  if (!ring) error("malloc error");
  end  = ring + k;
  last = ring;
  minpair = ring;
  minpair->value = in[0];
  minpair->death = k;
  out[0] = in[0];

  for (i=1;i<n;i++) {
    if (minpair->death == i) {
      minpair++;
      if (minpair >= end) minpair = ring;
    }
    if (in[i] <= minpair->value) {
      minpair->value = in[i];
      minpair->death = i+k;
      last = minpair;
    } else {
      while (last->value >= in[i]) {
        if (last == ring) last = end;
        --last;
      }
      ++last;
      if (last == end) last = ring;
      last->value = in[i];
      last->death = i+k;
    }
    out[i] = minpair->value;
  }
  free(ring);
}

void maxwindow(double *in, double *out, int n, int k)
{
  int i;
  struct pairs_d * ring;
  struct pairs_d * maxpair;
  struct pairs_d * end;
  struct pairs_d * last;

  ring = malloc(k * sizeof *ring);
  if (!ring) error("malloc error");
  end  = ring + k;
  last = ring;
  maxpair = ring;
  maxpair->value = in[0];
  maxpair->death = k;
  out[0] = in[0];

  for (i=1;i<n;i++) {
    if (maxpair->death == i) {
      maxpair++;
      if (maxpair >= end) maxpair = ring;
    }
    if (in[i] >= maxpair->value) {
      maxpair->value = in[i];
      maxpair->death = i+k;
      last = maxpair;
    } else {
      while (last->value <= in[i]) {
        if (last == ring) last = end;
        --last;
      }
      ++last;
      if (last == end) last = ring;
      last->value = in[i];
      last->death = i+k;
    }
    out[i] = maxpair->value;
  }
  free(ring);
}

SEXP run_min_max_C(SEXP x, SEXP win, SEXP wflag, SEXP dropAttr){
  SEXP l;
  int Win = asInteger(win);
  int i, P=0;
  R_len_t nr  = nrows(x);
  R_len_t nc  = ncols(x);
  R_len_t len = nr*nc;
  int wf = asInteger(wflag);
  int drop = asLogical(dropAttr);

  double * _x, *_l, *_u;
  if(!isReal(x))  { 
    x =  PROTECT( coerceVector(x,REALSXP)); P++;
  } 
  _x = REAL(x);

  if(nc>1 && drop) {
    PROTECT(l     = allocMatrix(REALSXP,nr,nc));     P++;
  } else { 
    PROTECT(l     = allocVector(REALSXP,len));       P++;
  }
  _l   = REAL(l);

  if(wf==1) {
    for(int i = 0; i<nc; i++)
      minwindow(_x+nr*i, _l+nr*i, nr, Win);
    if(!drop) DUPLICATE_ATTRIB(l,x);
    UNPROTECT(P);
    return(l);
  } else 
    if(wf==2) {
      for(int i = 0; i<nc; i++)
        maxwindow(_x+nr*i, _l+nr*i, nr, Win);
      if(!drop) DUPLICATE_ATTRIB(l,x);
      UNPROTECT(P);
      return(l);
    } else 
      if(wf==3) {
        SEXP u;
        if(nc>1 && drop){
          PROTECT(u     = allocMatrix(REALSXP,nr,nc));     P++;
        } else {
          PROTECT(u     = allocVector(REALSXP,len));       P++;
        }
        _u   = REAL(u);

        for(int i = 0; i<nc; i++){
          minwindow(_x+nr*i, _l+nr*i, nr, Win);
          maxwindow(_x+nr*i, _u+nr*i, nr, Win);
        }
        if(!drop){
          DUPLICATE_ATTRIB(l,x);
          DUPLICATE_ATTRIB(u,x);
        }
        SEXP ans, names;
        names = PROTECT(allocVector(STRSXP,2));       P++;
        SET_STRING_ELT(names,0,mkChar("max")); 
        SET_STRING_ELT(names,1,mkChar("min"));
        ans = PROTECT(allocVector(VECSXP, 2));        P++;
        SET_VECTOR_ELT(ans, 0, u); 
        SET_VECTOR_ELT(ans, 1, l);
        setAttrib(ans, R_NamesSymbol, names);
        UNPROTECT(P);
        return(ans);
      } else {
        error("wflag must be 1 for min, 2 for max, or 3 for both min and max.");
      } 
}



