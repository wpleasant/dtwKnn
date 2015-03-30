/*
    This part of dtwKnn is derived from the UCR Suite.  
*/


#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP run_znorm_C(SEXP x, SEXP win, SEXP epoch, SEXP backfill){
  SEXP ans;
  R_len_t  ii=0, i=0, j=0, ep=0, sp=0, ed=0;
  R_len_t  len = LENGTH(x);
  int P=0;
  int W     = asInteger(win);
  int WW    = W-1;
  int EPOCH = asInteger(epoch);
  int nr    = nrows(x);
  int nc    = ncols(x);
  int nj    = 0;
  int bf    = asLogical(backfill);


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
  
  double w =  (double)W; 
  
  memset(_ans,0,(size_t)len*sizeof(double));
   
  for(j=0;j<nc;j++){
    nj = nr*j;
    ii = 0;
    /* run through x restarting every epoch */
    while(ii<nr){
      ed = ((ii+EPOCH)<nr) ? (ii+EPOCH) : (nr);
      ex  = 0.0;
      ex2 = 0.0;
      if(ii==0){
       /* init and backfill first w rows */
        for(i=0; i<W; i++) {
          d       = _x[i+nj];
          ex      += d;
          ex2     += d*d;
        }
        mean = ex/w;
        std  = ex2/w;
        std  = sqrt(std-mean*mean);
        if(bf){
        for(i = 0 ; i < W; i++) 
          _ans[i+nj] = (_x[i+nj] - mean)/std;
        } else {
          for(i = 0 ; i < W; i++) 
          _ans[i+nj] = NA_REAL;
        }
        sp       = nj;
        ex      -= _x[sp];
        ex2     -= _x[sp]*_x[sp];
        ii = W;
      } else 
      if(ii < nr ){
        for(i=ii-WW; i<ii; i++) {
          d = _x[i+nj];
          ex+= d;
          ex2+= d*d;
        }
      }
      while(ii < ed) {
        d    =  _x[ii+nj];
        ex   += d;
        ex2  += d*d;
        mean =  ex/w;
        std  =  ex2/w;
        std  =  sqrt(std-mean*mean);
        _ans[ii+nj]  =  (d-mean)/std;
        sp   = (ii-WW)+nj;
        ex-=  _x[sp];
        ex2-=  _x[sp]*_x[sp];
        ii++;
      }
    }
  }
  DUPLICATE_ATTRIB(ans,x);
  UNPROTECT(P);
  return(ans);
}

SEXP run_MeanSD_C(SEXP x, SEXP win, SEXP epoch, SEXP backfill){
  SEXP mu, se;
  R_len_t  ii=0, i=0, j=0, ep=0, sp=0, ed=0;
  R_len_t  len = LENGTH(x);
  int P=0;
  int W     = asInteger(win);
  int WW    = W-1;
  int EPOCH = asInteger(epoch);
  int nr    = nrows(x);
  int nc    = ncols(x);
  int nj    = 0;
  int bf    = asLogical(backfill);

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
  
  memset(_mu,0,(size_t)len*sizeof(double));
  memset(_se,0,(size_t)len*sizeof(double));

  /* if win > nr reset W to nr*/
  if( W > nr )  W = nr;
  if( W > EPOCH) { 
    EPOCH = W;
    warning("EPOCH %d was increased to %d\n", EPOCH, W);
  }
  
  double w =  (double)W; 
  
  
  for(j=0;j<nc;j++){
    ii = 0;
    nj = nr*j;
    /* run through x restarting every epoch */
    while(ii<nr){
      ed = ((ii+EPOCH)<nr) ? (ii+EPOCH) : (nr);
      ex  = 0.0;
      ex2 = 0.0;
      if(ii==0){
       /* init and backfill first w rows */
        for(i=0; i<W; i++) {
          d       = _x[i+nj];
          ex      += d;
          ex2     += d*d;
        }
        mean = ex/w;
        std  = ex2/w;
        std  = sqrt(std-mean*mean);
        if(bf){
          for(i = 0 ; i < W; i++){
            _mu[i+nj] = mean;
            _se[i+nj] = std;
          }
        } else {
          for(i = 0 ; i < W; i++){ 
            _mu[i+nj] = NA_REAL;
            _se[i+nj] = NA_REAL;
          }
        }
        sp       = nj;
        ex      -= _x[sp];
        ex2     -= _x[sp]*_x[sp];
        ii = W;
      } else 
      if(ii < nr ){
        for(i=ii-WW; i<ii; i++) {
          d = _x[i];
          ex+= d;
          ex2+= d*d;
        }
      }
      while(ii < ed) {
        d    =  _x[ii+nj];
        ex   += d;
        ex2  += d*d;
        mean =  ex/w;
        std  =  ex2/w;
        _se[ii+nj]  =  sqrt(std-mean*mean);
        _mu[ii+nj]  =  mean;
        sp   = (ii-WW)+nj;
        ex   -=  _x[sp];
        ex2  -=  _x[sp]*_x[sp];
        ii++;
      }
     }
  }
  DUPLICATE_ATTRIB(mu,x);
  DUPLICATE_ATTRIB(se,x);
  /*mkNamed needs null termination*/
  const char *nms[] = {"mu","sd",""};
  SEXP ans = mkNamed(VECSXP, nms);
  SET_VECTOR_ELT(ans, 0, mu);
  SET_VECTOR_ELT(ans, 1, se);

  UNPROTECT(P);
  return(ans);
}


