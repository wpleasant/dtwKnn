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
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include "deque.h"

#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))
#define dist(x,y) ((x-y)*(x-y))

#define INF 1e20       //Pseudo Infitinte number for this code

/* Data structure for sorting the query*/
typedef struct ucr_Index
{   
  double value;
  int    index;
}  ucr_Index;

/* Sorting function for the query, sort by abs(z_norm(q[i])) from high to low */
static int comp(const void *a, const void* b)
{   
  ucr_Index* x = (ucr_Index*)a;
  ucr_Index* y = (ucr_Index*)b;
  return fabs(y->value) - fabs(x->value);   // high to low
}

/* Calculate quick lower bound
   Usually, LB_Kim take time O(m) for finding top,bottom,fist and last.
   However, because of z-normalization the top and bottom cannot give siginifant benefits.
   And using the first and last points can be computed in constant time.
   The prunning power of LB_Kim is non-trivial, especially when the query is not long, say in length 128.
   */
static double lb_kim_hierarchy(double *t, double *q, int j, int len, 
    double mean, double std, double bsf)
{
  /// 1 point at front and back
  double d, lb;
  double x0 = (t[j] - mean) / std;
  double y0 = (t[(len-1+j)] - mean) / std;
  lb = dist(x0,q[0]) + dist(y0,q[len-1]);
  if (lb >= bsf)   return lb;

  /// 2 points at front
  double x1 = (t[(j+1)] - mean) / std;
  d = min(dist(x1,q[0]), dist(x0,q[1]));
  d = min(d, dist(x1,q[1]));
  lb += d;
  if (lb >= bsf)   return lb;

  /// 2 points at back
  double y1 = (t[(len-2+j)] - mean) / std;
  d = min(dist(y1,q[len-1]), dist(y0, q[len-2]) );
  d = min(d, dist(y1,q[len-2]));
  lb += d;
  if (lb >= bsf)   return lb;

  /// 3 points at front
  double x2 = (t[(j+2)] - mean) / std;
  d = min(dist(x0,q[2]), dist(x1, q[2]));
  d = min(d, dist(x2,q[2]));
  d = min(d, dist(x2,q[1]));
  d = min(d, dist(x2,q[0]));
  lb += d;
  if (lb >= bsf)   return lb;

  /// 3 points at back
  double y2 = (t[(len-3+j)] - mean) / std;
  d = min(dist(y0,q[len-3]), dist(y1, q[len-3]));
  d = min(d, dist(y2,q[len-3]));
  d = min(d, dist(y2,q[len-2]));
  d = min(d, dist(y2,q[len-1]));
  lb += d;

  return lb;
}

/*
   LB_Keogh 1: Create Envelop for the query
   Note that because the query is known, envelop can be created once at the begenining.

   Variable Explanation,
order : sorted indices for the query.
uo, lo: upper and lower envelops for the query, which already sorted.
t     : a circular array keeping the current data.
j     : index of the starting location in t
cb    : (output) current bound at each position. It will be used later for early abandoning in DTW.
*/
static double lb_keogh_cumulative(int* order, double *t, double *uo, double *lo,
    double *cb, int j, int len, double mean, double std, double best_so_far)
{
  double lb = 0;
  double x, d;

  for (int i = 0; i < len && lb < best_so_far; i++)
  {
    x = (t[(order[i]+j)] - mean) / std;
    d = 0;
    if (x > uo[i])
      d = dist(x,uo[i]);
    else if(x < lo[i])
      d = dist(x,lo[i]);
    lb += d;
    cb[order[i]] = d;
  }
  return lb;
}

/* LB_Keogh 2: Create Envelop for the data
   Note that the envelops have been created (in main function) when each data point has been read.

   Variable Explanation,
tz: Z-normalized data
qo: sorted query
cb: (output) current bound at each position. Used later for early abandoning in DTW.
l,u: lower and upper envelop of the current data
*/
static double lb_keogh_data_cumulative(int* order, double *tz, double *qo, double *cb,
    double *l, double *u, int len, double mean, double std, double best_so_far)
{
  double lb = 0;
  double uu,ll,d;

  for (int i = 0; i < len && lb < best_so_far; i++)
  {
    uu = (u[order[i]]-mean)/std;
    ll = (l[order[i]]-mean)/std;
    d = 0;
    if (qo[i] > uu)
      d = dist(qo[i], uu);
    else
    {   
      if(qo[i] < ll)
        d = dist(qo[i], ll);
    }
    lb += d;
    cb[order[i]] = d;
  }
  return lb;
}

/// Calculate Dynamic Time Wrapping distance
/// A,B: data and query, respectively
/// cb : cummulative bound used for early abandoning
/// r  : size of Sakoe-Chiba warpping band

static double dtw(double* A, double* B, double *cb, int m, int r, double bsf)
{

  double *cost;
  double *cost_prev;
  double *cost_tmp;
  int i,j,k;
  double x,y,z,min_cost;

  /// Instead of using matrix of size O(m^2) or O(mr), we will reuse two array of size O(r).

  cost       = R_Calloc((2*r+1), double);
  cost_prev  = R_Calloc((2*r+1), double);
  for(k=0; k<2*r+1; k++) {
    cost[k]       = INF;
    cost_prev[k]  = INF;
  }

  for (i=0; i<m; i++)
  {
    k = max(0,r-i);
    min_cost = INF;

    for(j=max(0,i-r); j<=min(m-1,i+r); j++, k++)
    {
      /// Initialize all row and column
      if ((i==0)&&(j==0))
      {
        cost[k]=dist(A[0],B[0]);
        min_cost = cost[k];
        continue;
      }

      if ((j-1<0)||(k-1<0))     y = INF;
      else                      y = cost[k-1];
      if ((i-1<0)||(k+1>2*r))   x = INF;
      else                      x = cost_prev[k+1];
      if ((i-1<0)||(j-1<0))     z = INF;
      else                      z = cost_prev[k];

      /// Classic DTW calculation
      cost[k] = min( min( x, y) , z) + dist(A[i],B[j]);

      /// Find minimum cost in row for early abandoning (possibly to use column instead of row).
      if (cost[k] < min_cost)
      {   
        min_cost = cost[k];
      }
    }

    /// We can abandon early if the current cummulative distace with lower bound together are larger than bsf
    if (i+r < m-1 && min_cost + cb[i+r+1] >= bsf)
    {   
      Free(cost);
      Free(cost_prev);
      return min_cost + cb[i+r+1];
    }

    /// Move current array to previous array.
    cost_tmp = cost;
    cost = cost_prev;
    cost_prev = cost_tmp;
  }
  k--;

  /// the DTW distance is in the last cell in the matrix of size O(m^2) or at the middle of our array.
  double final_dtw = cost_prev[k];
  Free(cost);
  Free(cost_prev);
  return final_dtw;
}



static  double ucr_set_knn(double dist, double *_kvec, int *_lvec, int *wk, int *_wk,
    int *last_loc, int loc, int K, int JS){
  double bsf = 0.0;
  if(loc - *last_loc > JS){
    _kvec[*wk]  = dist;
    _lvec[*wk]  = loc;
    *_wk = *wk;
    *last_loc   =  loc;
  } else {
    if(dist < _kvec[*_wk]){
      _kvec[*_wk]  = dist;
      _lvec[*_wk]  = loc;
      *last_loc    = loc;
    }
  }
  /* Now for K-NN redefine bsf to the "worst of the" "best so far" and save it's index at wk */
  for(int k=0; k<K; k++) { 
    if(_kvec[k]>bsf) {
      bsf = _kvec[k];
      *wk  = k;
    }
  }
  return bsf;
}

SEXP ucr_dtw_knn_C( SEXP query, SEXP reference, SEXP window, SEXP sizeK
    ,SEXP epoch, SEXP sortd, SEXP sizeJ, SEXP Depth,SEXP verbose)
{
  double bsf;          /// best-so-far
  double *t, *q;       /// data array and query array
  int *order;          ///new order of the query
  double *u, *l, *qo, *uo, *lo,*tz,*cb, *cb1, *cb2,*u_d, *l_d;


  double d, m_d;
  R_len_t i , j, ii,k,P=0;
  double ex , ex2 , mean, std;
  int m=-1, r=-1;
  R_len_t loc = 0, lenR, last_loc=0;
  int kim = 0,keogh = 0, keogh2 = 0;
  double dist=0, lb_kim=0, lb_k=0, lb_k2=0;
  double  *u_buff, *l_buff, * ref;
  ucr_Index *Q_tmp;

  /* read size of the query */
  m       = LENGTH(query);
  m_d     = (double)m;
  lenR    = LENGTH(reference);

  if(m<4) error("query length must be > 3");
  if(lenR<m) error("reference length must be >= query length\n");

  /* For every EPOCH points, all cummulative values
      , such as ex (sum), ex2 (sum square), will be restarted for reducing the floating point error.*/
  int EPOCH = asInteger(epoch);
  if(EPOCH > lenR)  EPOCH = lenR;
  
  /* K as in the max number of neighbors */
  int K     = asInteger(sizeK);
  if(K<1) K = 1;

  int wk    = 0;  /* index position of worst of the KNN */
  int _wk   = 0;  /* the index position of the previous worst of the KNN. Used if JS > 1 */
  
  SEXP kvec;
  PROTECT(kvec   = allocVector(REALSXP,K));       P++;
  double *_kvec  = REAL(kvec);

  SEXP lvec;
  PROTECT(lvec   = allocVector(INTSXP,K));       P++;
  int *_lvec  = INTEGER(lvec);
  
  for(k = 0; k<K; k++){
    _kvec[k] = INF;
    _lvec[k] = 0;
  }


  /*  get max jump size */
  int JS  = asInteger(sizeJ);
  if(K==1)
    JS = 1;
  else
  if(JS > lenR) 
    JS = lenR-1;

  int depth  = asInteger(Depth);
  if(depth < 1) depth =1;


  if(!isReal(query))  { 
    query = PROTECT(coerceVector(query,REALSXP)); P++;
  }
  
  SEXP qu;
  PROTECT( qu = duplicate(query)) ;                       P++;
  q   = REAL(qu);

  
  if(!isReal(reference))  { 
    reference = PROTECT(coerceVector(reference,REALSXP)); P++;
  }
  ref   = REAL(reference);
  


  /* read warping windows */
  double R = asReal(window);
  if (R<=1)
    r = floor(R*m_d);
  else
    r = floor(R);
  if(r>m_d) r = m_d;

  /// malloc everything here
  qo      = (double*) R_alloc(m,sizeof(double));
  uo      = (double*) R_alloc(m,sizeof(double));
  lo      = (double*) R_alloc(m,sizeof(double));
  order   = (int*)    R_alloc(m,sizeof(int));
  Q_tmp   = (ucr_Index *) R_alloc(m,sizeof(ucr_Index));
  u       = (double*) R_alloc(m,sizeof(double));
  l       = (double*) R_alloc(m,sizeof(double));
  cb      = (double*) R_alloc(m,sizeof(double));
  cb1     = (double*) R_alloc(m,sizeof(double));
  cb2     = (double*) R_alloc(m,sizeof(double));
  u_d     = (double*) R_alloc(m,sizeof(double));
  l_d     = (double*) R_alloc(m,sizeof(double));
  t       = (double*) R_alloc(m*2,sizeof(double));
  tz      = (double*) R_alloc(m,sizeof(double));
  u_buff  = (double*) R_alloc(lenR,sizeof(double));
  l_buff  = (double*) R_alloc(lenR,sizeof(double));


  /// Read query file
  bsf = INF;
  i = 0;
  j = 0;
  ex = ex2 = 0;
  
  for(i = 0; i < m; i++)
  {
    d   = q[i];
    ex += d;
    ex2 += d*d;
    q[i] = d;
  }

  /// Do z-normalize the query, keep in same array, q
  mean = ex/m_d;
  std = ex2/m_d;
  std = sqrt(std-mean*mean);
  for( i = 0 ; i < m ; i++ )
    q[i] = (q[i] - mean)/std;

  /// Create envelop of the query: lower envelop, l, and upper envelop, u
  lower_upper_lemire(q, m, r, l, u);

  /// Sort the query one time by abs(z-norm(q[i]))
  for( i = 0; i<m; i++)
  {
    Q_tmp[i].value = q[i];
    Q_tmp[i].index = i;
  }
  qsort(Q_tmp, m, sizeof(ucr_Index),comp);

  /// also create another arrays for keeping sorted envelop
  for( i=0; i<m; i++)
  {   
    int o = Q_tmp[i].index;
    order[i] = o;
    qo[i] = q[o];
    uo[i] = u[o];
    lo[i] = l[o];
  }

  /// Initial the cummulative lower bound
  for( i=0; i<m; i++)
  {   
    cb[i]=0;
    cb1[i]=0;
    cb2[i]=0;
  }
  
  i = 0;          // current index of the data in current chunk of size EPOCH
  j = 0;          // the starting index of the data in the circular array, t
  ii = 0;
  k  = 0;
  ex = ex2 = 0;
  R_len_t  I;    // the starting index of the data in current chunk of size EPOCH
  R_len_t  ep  = EPOCH;

  
  
  lower_upper_lemire(ref, lenR, r, l_buff, u_buff);

  /* Do main task here.. */
  ex=0; ex2=0;
  /* Init index zero and start loop at one */ 
  d = ref[i];
  ex += d;
  ex2 += d*d;
  t[i%m] = d;
  t[(i%m)+m] = d;


  for(i=1; i<lenR; i++)
  {
    /* If EPOCH data have been computed refresh ex and ex2 */
    if(i == ep)
    {
      if(i==lenR-1) break;
      ex = 0; 
      ex2 = 0;
      for(ii = i-m+1; ii < i; ii++)
      { 
        d = ref[ii];
        ex += d;
        ex2 += d*d;
      }
      ep += EPOCH;
    }

    /* A bunch of data has been read and pick one of them at a time to use */
    d = ref[i];

    /* Calcualte sum and sum square */
    ex += d;
    ex2 += d*d;


    /* t is a circular array for keeping current data */
    t[i%m] = d;

    /* Double the size for avoiding using modulo "%" operator */
    t[(i%m)+m] = d;

    /* Start the task when there are more than m-1 points in the current chunk*/
    if( i >= m-1 )
    {
      mean = ex/m_d;
      std  = ex2/m_d;
      std  = sqrt(std-mean*mean);

      /* compute the start location of the data in the current circular array, t*/
      j = (i+1)%m;
      /* the start location of the data in the current chunk */
      I = i-(m-1);

      /* Use a constant lower bound to prune the obvious subsequence */
      lb_kim = lb_kim_hierarchy(t, q, j, m, mean, std, bsf);

      if (lb_kim < bsf)
      {
        if(depth==1){
          dist = lb_kim;
          goto checkDist;

        }
        /* Use a linear time lower bound to prune; z_normalization of t will be computed on the fly.
            uo, lo are envelop of the query.*/
        lb_k = lb_keogh_cumulative(order, t, uo, lo, cb1, j, m, mean, std, bsf);
        if (lb_k < bsf)
        {
          if(depth==2){
            dist = lb_k;
            goto checkDist;

          }
          /* Take another linear time to compute z_normalization of t.
             Note that for better optimization, this can merge to the previous function.*/
          for(k=0;k<m;k++)  {   tz[k] = (t[(k+j)] - mean)/std; }

          /* Use another lb_keogh to prune
             qo is the sorted query. tz is unsorted z_normalized data.
             l_buff, u_buff are big envelop for all data in this chunk */
          lb_k2 = lb_keogh_data_cumulative(order, tz, qo, cb2, l_buff+I, u_buff+I, m, mean, std, bsf);
          if (lb_k2 < bsf)
          {
            if(depth==3){
              dist = min(lb_k,lb_k2);
              goto checkDist;
            }

            /* Choose better lower bound between lb_keogh and lb_keogh2 to be used in early abandoning DTW
                 Note that cb and cb2 will be cumulative summed here.*/
            if (lb_k > lb_k2)
            {
              cb[m-1]=cb1[m-1];
              for(k=m-2; k>=0; k--)
                cb[k] = cb[k+1]+cb1[k];
            }
            else
            {
              cb[m-1]=cb2[m-1];
              for(k=m-2; k>=0; k--)
                cb[k] = cb[k+1]+cb2[k];
            }
            /* Compute DTW and early abandoning if possible */
            dist = dtw(tz, q, cb, m, r, bsf);
 checkDist: 
            if( dist < bsf )
            { 
              loc       = I+1; /* R's starting index of reference data */
              if(K==1) /* orginal case UCR Suite K=1 */
              {
                _kvec[0] = dist;
                _lvec[0] = loc;
                bsf      = dist;
              } 
              else /* Update K nearest neighbor subject to jump size */
              {
                bsf = ucr_set_knn(dist, _kvec, _lvec, &wk, &_wk, &last_loc, loc, K, JS);
              }
            }
          } else
            keogh2++;
        } else
          keogh++;
      } else
        kim++;

      /* Reduce obsolute points from sum and sum square */
      ex  -= t[j];
      ex2 -= t[j]*t[j];
    }
  }


  /* Use R's sort with index to sort in decreasing order */
  if(asLogical(sortd)&& K>1) rsort_with_index(_kvec,  _lvec, K);
  
  SEXP ans, names, pruned;
  /* fill the ans list */
  PROTECT(ans  = allocVector(VECSXP,  2));              P++;
  SET_VECTOR_ELT(ans, 0, kvec);
  SET_VECTOR_ELT(ans, 1, lvec);
  /* name the list */
  PROTECT(names = allocVector(STRSXP, 2));              P++;
  SET_STRING_ELT(names, 0, mkChar("distance"));
  SET_STRING_ELT(names, 1, mkChar("location"));
  setAttrib(ans, R_NamesSymbol, names);
  /* print if verbose */
  if(asLogical(verbose)){
    Rprintf("\n");
    Rprintf("Pruned by LB_Kim    : %6.2f%%\n", ((double) kim / i)*100);
    if(depth>1) Rprintf("Pruned by LB_Keogh  : %6.2f%%\n", ((double) keogh / i)*100);
    if(depth>2) Rprintf("Pruned by LB_Keogh2 : %6.2f%%\n", ((double) keogh2 / i)*100);
    if(depth>3) Rprintf("DTW Calculation     : %6.2f%%\n", 100-(((double)kim+keogh+keogh2)/i*100));
  }
  UNPROTECT(P);
  return ans;
}

