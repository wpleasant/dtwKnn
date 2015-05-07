/*
 * 
 *  Author: William Pleasant 
 *   
 *   Based on "ascending minima" algorithm. original code by Richarda Harter na svetu.
 *   http://richardhartersworld.com/cri/2001/slidingmin.html 
 *   
 *   lb_env doubles the speed and replaces Danial Lemire's veresion "lower_upper_lemire" 
 *   finding the upper and lower envelope for the LB_Keogh. 
 *
 */


struct pairs_double {
  double value;
  int death;
};


void lb_env(double *in, double *lo, double *hi, int n, int k);

