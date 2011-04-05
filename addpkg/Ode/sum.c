// K M Briggs 98 Dec 29
// was cplusplus/sum.cc KMB 98 Dec 28
// Plain C version 99 Jan 04
/*
* The code realises the Priest-algorithm found in Nicholas J. Higham,
* Accuracy and Stability of Numerical Algorithms, SIAM, 1996 on page 97.
* (The computed sum is accurate virtually to full precision under
* assumptions satisfied by the IEEE arithmetic.)
*/

#include "sum.h"

#define TRYSUM
#define DBGSUM
#undef TRYSUM
#undef DBGSUM

#ifdef TRYSUM
int main() {
  double eps=1e-12;
  printf("%22.18g\n",1.0e-9+eps-sum0(3,1.0e-9,-1.0-3*eps,1.0+4*eps));
  printf("%22.18g\n",1.0e-9+eps-sum2(3,1.0e-9,-1.0-3*eps,1.0+4*eps));
  return 0;
}
#endif

int comparedown(const void *x, const void *y) {
  /* for non-increasing |x[i]|... */
  double fx = fabs(*(double *) x), fy = fabs(*(double *) y);
  return fx>fy?-1:(fx==fy?0:1);
}

double sum0(const unsigned n, ...) {
  int i; double s=0;
  va_list xx;
  va_start(xx,n);
  for (i = 0; i < n; i++) s+=va_arg(xx,double);
  va_end(xx);
  return s;
}

double sum2(const unsigned n, ...) {
  double y, u, t, v, z, s1, s0, c0 = 0;
  int i;
  double *x;
  va_list xx;
  #ifdef X86
  x86_FIX;
  #endif
  if (n < 1) {
    fprintf(stderr, "n=%d in sum2!\n", n);
    exit(1);
  }
  x=(double*)malloc(n*sizeof(double));
  va_start(xx,n);
  for (i = 0; i < n; i++) x[i]=va_arg(xx,double);
  va_end(xx);
  qsort((void *) x, n, sizeof(double), comparedown);
  #ifdef DBGSUM
    for (i = 0; i < n; i++) printf("%22.18g ",x[i]);
    printf("\n");
  #endif
  s0=x[0];
  for (i = 1; i < n; i++) {
    if (0.0 == x[i]) break; /* return s0; */ /* no more nonzero values */
    y = c0 + x[i];
    u = x[i] - (y - c0);
    t = y + s0;
    v = y - (t - s0);
    z = u + v;
    s1 = t + z;
    c0 = z - (s1 - t);
    s0 = s1;
  }
  #ifdef X86
  END_x86_FIX;
  #endif
  free(x);
  return s0;
}

/*
double sum1(double *x, const unsigned n) {
  double y, u, t, v, z, s1, s0, c0 = 0;
  int i;
  if (n < 1) {
    fprintf(stderr, "n= %d in sum1!\n", n);
    exit(1);
  }
  qsort((void *) x, n, sizeof(double), comparedown);
  s0=x[0];
  for (i = 1; i < n; i++) {
    if (0.0 == x[i]) return s0; // no more nonzero values
    y = c0 + x[i];
    u = x[i] - (double) (y - c0);
    t = y + s0;
    v = y - (double) (t - s0);
    z = u + v;
    s1 = t + z;
    c0 = z - (double) (s1 - t);
    s0 = s1;
  }
  return s0;
}
*/
