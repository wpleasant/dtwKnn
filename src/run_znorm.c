/*
    This part of dtwKnn is derived from the UCR Suite.  
*/


#include <R.h>
#include <Rinternals.h>

SEXP run_znorm_C(SEXP x, SEXP win, SEXP epoch, SEXP backfill, SEXP dropAttr)
{
  SEXP ans;
  R_len_t  ii=0, i=0, j=0, ep=0, sp=0, ed=0;
  int P=0;
  int W     = asInteger(win);
  int WW    = W-1;
  int EPOCH = asInteger(epoch);
  R_len_t nr    = nrows(x);
  int nc    = ncols(x);
  int nj    = 0;
  int bf    = asLogical(backfill);
  int DA    = asLogical(dropAttr);
  if(nr<2) return(x);

  R_len_t len  = (R_len_t)nr*nc;


  double ex=0.0, ex2=0.0, mean=0.0, std=0.0, d;
 
  double * _x;
  if(!isReal(x))  { 
   x =  PROTECT( coerceVector(x,REALSXP)); P++;
  } 
  _x = REAL(x);

  ans = PROTECT(allocVector(REALSXP,len));       P++;
  double  *_ans   = REAL(ans);

  /* if win > nr reset W to nr*/
  if( W > nr )  W = nr;
  if( W > EPOCH) { 
    EPOCH = W;
    warning("EPOCH %d was increased to %d\n", EPOCH, W);
  }
  if(EPOCH>nr) EPOCH=nr;
  
  double w =  (double)W; 
  
  //memset(_ans,0,(size_t)len*sizeof(double));
  
  for(j=0;j<nc;j++){
    nj = nr*j;
    ex=0; ex2=0;
    d = _x[nj];
    ex += d;
    ex2 += d*d;
    ep = EPOCH;
    for(i=1; i<nr; i++)
    {
      /* If EPOCH data have been computed refresh ex and ex2 */
      if(i == ep)
      { 
        ex  = 0; 
        ex2 = 0;
        for(ii = i-WW; ii < i; ii++)
        { 
          d = _x[ii+nj];
          ex += d;
          ex2 += d*d;
        }
        ep += EPOCH;
      }
      d = _x[i+nj];
      ex += d;
      ex2 += d*d;
      if( i >= WW )
      {
        mean =  ex/w;
        std  =  ex2/w;
        std  =  sqrt(std-mean*mean);
        _ans[i+nj]  =  (d-mean)/std;
        sp   = (i-WW)+nj;
        ex  -=  _x[sp];
        ex2 -=  _x[sp]*_x[sp];
      }
    }
    if(bf){
      for(ii=WW+nj-1; ii>=nj; ii--)
        _ans[ii] = _ans[W+nj];
    } else {
      for(ii=WW+nj-1; ii>=nj; ii--)
        _ans[ii] = NA_REAL;
    }
  }
  if(!DA) DUPLICATE_ATTRIB(ans,x);
  UNPROTECT(P);
  return(ans);
}


SEXP run_MeanSD_C(SEXP x, SEXP win, SEXP epoch, SEXP backfill, SEXP dropAttr){
  SEXP mu, se;
  R_len_t  ii=0, i=0, j=0, ep=0, sp=0, ed=0;
  //R_len_t  len = LENGTH(x);
  int P=0;
  int W     = asInteger(win);
  int WW    = W-1;
  int EPOCH = asInteger(epoch);
  int nr    = nrows(x);
  int nc    = ncols(x);
  int nj    = 0;
  int bf    = asLogical(backfill);
  int DA    = asLogical(dropAttr);
  if(nr<2) return(x);
  R_len_t len  = (R_len_t)nr*nc;

  double ex=0.0, ex2=0.0, mean=0.0, std=0.0, d;
 
  double * _x;
  if(!isReal(x))  { 
   x =  PROTECT( coerceVector(x,REALSXP)); P++;
  } 
  _x = REAL(x);
  
  
  mu = PROTECT(allocVector(REALSXP,len));       P++;
  se = PROTECT(allocVector(REALSXP,len));       P++;
  double  *_mu   = REAL(mu);
  double  *_se   = REAL(se);
  

  /* if win > nr reset W to nr*/
  if( W > nr )  W = nr;
  if( W > EPOCH) { 
    EPOCH = W;
    warning("EPOCH %d was increased to %d\n", EPOCH, W);
  }
  if(EPOCH>nr) EPOCH = nr;

  
  double w =  (double)W; 
 
  for(j=0;j<nc;j++){
    nj = nr*j;
    ex=0; ex2=0;
    d = _x[nj];
    ex += d;
    ex2 += d*d;
    ep = EPOCH;
    for(i=1; i<nr; i++)
    {
      /* If EPOCH data have been computed refresh ex and ex2 */
      if(i == ep)
      {
        ex  = 0; 
        ex2 = 0;
        for(ii = i-WW; ii < i; ii++)
        { 
          d   = _x[ii+nj];
          ex  += d;
          ex2 += d*d;
        }
        ep += EPOCH;
      }
      d = _x[i+nj];
      ex += d;
      ex2 += d*d;
      if( i >= WW )
      {
        mean =  ex/w;
        std  =  ex2/w;
        std  =  sqrt(std-mean*mean);
        _mu[i+nj] = mean;
        _se[i+nj] = std;
        sp   = (i-WW)+nj;
        ex  -=  _x[sp];
        ex2 -=  _x[sp]*_x[sp];
      }
    }
    if(bf){
      for(ii=WW+nj-1; ii>=nj; ii--){
        _mu[ii] = _mu[WW+nj];
        _se[ii] = _se[WW+nj];
      }
    } else {
      for(ii=WW+nj-1; ii>=nj; ii--){
        _mu[ii] = NA_REAL;
        _se[ii] = NA_REAL;
      }
    }
  }

  if(!DA){
    DUPLICATE_ATTRIB(mu,x);
    DUPLICATE_ATTRIB(se,x);
  }
  /*mkNamed needs null termination*/
  const char *nms[] = {"mu","sd",""};
  SEXP ans = mkNamed(VECSXP, nms);
  SET_VECTOR_ELT(ans, 0, mu);
  SET_VECTOR_ELT(ans, 1, se);

  UNPROTECT(P);
  return(ans);
}


