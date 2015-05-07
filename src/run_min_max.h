

/* The "ascending minima" algorithm. original code by Richarda Harter na svetu.
 * http://richardhartersworld.com/cri/2001/slidingmin.html */

struct pairs_d {
  double value;
  int death;
};

void minwindow(double *in, double *out, int n, int k);
void maxwindow(double *in, double *out, int n, int k);





