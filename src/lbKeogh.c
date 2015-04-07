#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <limits.h>
#include "deque.h"


#define dist(x,y) ((x-y)*(x-y))

static double lb_keogh_no_sort(double *t, double *uo, double *lo,
    int len,  double best_so_far)
{
  double lb = 0;
  double x, d;
  for (int i = 0; i < len && lb < best_so_far; i++)
  {
    x = t[i];
    d = 0;
    if (x > uo[i])
      d = dist(x,uo[i]);
    else if(x < lo[i])
      d = dist(x,lo[i]);
    lb += d;
  }
  return lb;
}

static double lb_keogh_no_sort_cum(double *t, double *uo, double *lo,
    double *cb, int len, double best_so_far, int csumcb)
{
  double lb = 0.0;
  double x, d;
  for (int i = 0; i < len && lb < best_so_far; i++)
  {
    x = t[i];
    d = 0;
    if (x > uo[i])
      d = dist(x,uo[i]);
    else if(x < lo[i])
      d = dist(x,lo[i]);
    lb += d;
    cb[i] = d;
  }
  if(csumcb) {
    long double sum = 0.0;
    for(int i = len-1; i >= 0; i--) 
    {
      sum   += cb[i];
      cb[i] =  sum;
    }
  }
  return lb;
}


SEXP lb_keogh_C(SEXP x, SEXP y, SEXP win, SEXP returnsum, SEXP csum, SEXP dropAttr){
  SEXP  ans;
  if(isNewList(x) || isNewList(y)) error("x and y can't be lists or data.frames");
  int P      = 0;
  int nr     = nrows(x);
  int nc     = ncols(x);
  int ncY    = ncols(y);
  if(nrows(y) != nr) error("x and y should have the same number of rows");
  if(nc > 1 && !isArray(x)) error("x must be numeric vector or an array ");
  if(ncY!=1 && ncY != nc ) error(" y must have num columns == 1 or num colums x");
  int Win       = asInteger(win);
  int csumcb    = asLogical(csum);
  int retSum    = asLogical(returnsum);
  int drop      = asLogical(dropAttr);
  double  lsum  = 0.0;
  double *_ans;
  if(retSum && csumcb) csumcb=0;

  if(!isReal(x))  { 
    x =  PROTECT( coerceVector(x,REALSXP)); P++;
  } 
  double * _x = REAL(x);
  if(!isReal(y))  { 
    y =  PROTECT( coerceVector(y,REALSXP)); P++;
  } 
  double * _y = REAL(y);

  /* allocate for lower and upper bounds */
  double *_l, *_u; 
  _l  = (double*) R_alloc(nr,sizeof(double));
  _u  = (double*) R_alloc(nr,sizeof(double));
  memset(_l,0,(size_t)(nr)*sizeof(double));
  memset(_u,0,(size_t)(nr)*sizeof(double));

  /*
     If ncY > 1 then we return the  multivariate version inspired by 
     Toni M. Rath and R. Manmatha " Lower-Bounding of Dynamic Time Warping  
     for Multivariate Time Series ". Otherwise we check y to each column in x
  */

  if(ncY==1) lower_upper_lemire(_y, nr, Win, _l, _u);
  if(!retSum){
    PROTECT(ans   = allocMatrix(REALSXP,nr,nc));       P++;
    _ans = REAL(ans);
    memset(_ans,0,(size_t)(nr*nc)*sizeof(double));
    for(int j = 0; j<nc; j++){
      if(ncY > 1) lower_upper_lemire(_y+(nr*j), nr, Win, _l, _u);
      lsum = lb_keogh_no_sort_cum(_x+(nr*j), _u, _l, _ans+(nr*j), nr, INFINITY, csumcb);
    }
  } else {
    PROTECT(ans  = allocVector(REALSXP,nc));           P++;
    _ans = REAL(ans);
    memset(REAL(ans),0,(size_t)(nc)*sizeof(double));
    for(int j = 0; j<nc; j++){
      if(ncY > 1) lower_upper_lemire(_y+(nr*j), nr, Win, _l, _u);
      _ans[j] = lb_keogh_no_sort( _x+(nr*j), _u, _l, nr, INFINITY);
    }
  }
  /* duplicate all attributes if not returning sum and checking for one to many */
  if(nc>1 && !retSum && !drop) 
    DUPLICATE_ATTRIB(ans, x);
  else
    if(nc>1 && retSum && !drop){
      /* duplicate colnames if returning sum from each column for one to many */
      SEXP colnames = GetColNames(getAttrib(x, R_DimNamesSymbol));
      if(!isNull(colnames))
        namesgets(ans,colnames);
    }
  UNPROTECT(P);
  return(ans);
}

