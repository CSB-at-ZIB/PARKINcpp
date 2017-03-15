/* Ode.c 
   K M Briggs 
   98 Oct 09, 13 first version
   99 Jan 04,1 8 using dop853b
   99 Aug 12 new argument nout
   00 April 13 renamed Ode
   01 Aug 08 ode.py python preprocessor
   01 Aug 13,17 -r and -a options
   01 Sep 13 Improved nout control: the number of output lines
             should always be nout+1, the first line having x value x0
	     and the last line having x value x1
 
   ode solver: "Ode [-rrtoler] [-aatoler] rhs1 rhs2 ... iv1 iv2 ... x0 x1 nout"
   independent var: x
   dependent vars:  a,b,c,d,...
   integration range x0...x1 (x1 may be <= x0)
   rtoler=relative tolerance
   atoler=absolute tolerance
   nout+1 lines are output (including IV).

   Examples:
   Exponential growth:
     ./Ode a 1 0 1 100 | p
     (NB: "Ode -a 1 0 1 100" will have -a interpreted as a option.  Use "0-a")
   Sine: 
     ./Ode b -a 0 1 0 6.28 100 | cols -c1,2 | p
   Damped sine: 
     ./Ode b-a/5 -a 0 1 0 62.8 100 | cols -c1,2 | p
   Logistic:
     ./Ode "a*(1-a)" 0.1 0 5 100 | p
   Richards: nu=2 beta=3 kappa=4 y0=0.1
     ./Ode "1.5*a*(1-a^2/16)" 0.1 0 4 100 | p
   Lotka-Volterra: 
     ./Ode "a-1.3*a*b" "1.5*a*b-2*b" 1 2 0 20 100 | cols -c1,2 | p
   Gause:
     ./Ode "a*(0.8-0.00016*a-0.00048*b)" "b*(0.6-0.0003*b-0.00015*a)" 20 20 0 20 100 | colex 1 3 | p
   Lorenz: 
     ./Ode "3*(b-a)" "-a*c+26*a-b" "a*b-c" 0 1 0 0 30 100 | cols -c1,3 | p
     ./Ode "16*(b-a)" "a*(45.92-c)-b" "a*b-4*c" 0 1 0 0 30 100 | cols -c1,2 | p # dlia
   Adam's model:
     ./Ode -r1e-4 -a1e-4 "-0.05*sin(2*pi()*x/18.2)-0.06*sin(2*pi()*x/19)-0.05*sin(2*pi()*x/23)+0.6*a-4*a^3" 0.4 40 56 1000 | p
*/

#define EIGHT /* if defined, uses eighth-order dop853, else dopri5 */

#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>
#include "formulc.h"
#ifdef EIGHT
  #include "dop853b.h"
  #define contd5 contd8
  #define dopri5 dop853
  void fcn(unsigned n,double x, double *y, double *dy, double *cd);
#else
  #include "dopri5.h"
  void fcn (unsigned n, double x, double *y, double *f);
#endif
void solout(long nr,double xold, double x, double* y, unsigned n, int* irtrn);

formu* f;
char* vars;
int neq,nout,nlines=0;
double dx;

/* extra formulc functions... */
double sqr(double x) { return x*x; }
double cub(double x) { return x*x*x; }
double sin2pi(double x) { return sin(6.283185307179586477*x); }
double cos2pi(double x) { return cos(6.283185307179586477*x); }
double smoothiside(double x, double g) {  // smoothiside
  double a,b;
  if (x<0.0) return 0.0;
  a=g*x;
  if (a>=1.0) return 1.0;
  b=a/(1.0-a);
  return -expm1(-b*b);
}
double pulse(double x, double g) { /* ___|-|___ */
  double a,b;
  if (x<=0.0) return 0.0;
  if (x>=1.0) return 0.0;
  a=g*x;
  if (a>1.0) { /* 1/g < x < 1 */
    if (a<g-1) return 1.0; /* 1/g < x < 1-1/g */
    b=g*(1-x)/(1.0-g*(1-x)); return -expm1(-b*b); /* 1-1/g < x < 1  */
  }
  b=a/(1.0-a); return -expm1(-b*b); /* 0 < x < 1/g  */
}

int main(int argc, char *argv[]) {
  int i,j,inform,len,err;
  double x,*y;
  double xend;           /* final x-value (xend-x may be positive or negative) */
#ifdef EIGHT
  double rtoler=1.0e-12;  /* relative error tolerance */
  double atoler=1.0e-12;  /* absolute error tolerance */
  double *cd;            /* not used */
#else
  double rtoler=1.0e-9;  /* relative error tolerance */
  double atoler=1.0e-9;  /* absolute error tolerance */
#endif
  int itoler=0;          /* switch for rtoler and atoler */
  int iout=2;            /* switch for calling solout */
  double uround=0.0;     /* rounding unit */
  double safe=0.0;       /* safety factor */
  double fac1=0.0;       /* parameters for step size selection */
  double fac2=0.0;
  double beta=0.0;       /* for stabilized step size control */
  double hmax;           /* maximal step size */
  double h;              /* initial step size */
  long nmax=10000;       /* maximal number of allowed steps */
  int meth=0;            /* switch for the choice of the coefficients */
  long nstiff=0;         /* test for stiffness */
  unsigned nrdens;       /* number of components for which dense output is required */
  unsigned licont=0;     /* declared length of icont */
  unsigned* icont=NULL;  /* indexes for which dense output is required  NULL=all */
  FILE* fileout=stderr;  /* messages stream */
  const char *OdeVersion=" "__DATE__" "__TIME__;
  if (argv[1][0]=='-' && argv[1][1]=='r') { /* relative error */
    rtoler=atof(argv[1]+2);
    if (rtoler<1.0e-15 || rtoler>0.1) rtoler=1e-9;
    argc--; argv++;
  }
  if (argv[1][0]=='-' && argv[1][1]=='a') { /* absolute error */
    atoler=atof(argv[1]+2);
    if (atoler<1.0e-15 || atoler>0.1) atoler=1e-9;
    argc--; argv++;
  }
  nrdens=neq=(argc-4)/2;
  if (2*neq!=argc-4) { 
    fprintf(fileout,"| Ode: error, wrong number of arguments: argc=%d\n",argc); 
    fprintf(fileout,"| Usage: Ode rhs1 rhs2 ... iv1 iv2 ... xstart xend nout\n"); 
    fprintf(fileout,"| Example: compute e: Ode a 1 0 1 1\n"); 
    return 1; 
  }
  /* extra formulc functions... */
  fnew("smoothiside",(Func)smoothiside,2,0);
  fnew("step",(Func)smoothiside,2,0);
  fnew("pulse",(Func)pulse,2,0);
  fnew("expm1",(Func)expm1,1,0);
  fnew("log1p",(Func)log1p,1,0);
  fnew("sin2pi",(Func)sin2pi,1,0);
  fnew("cos2pi",(Func)cos2pi,1,0);
  fnew("sqr",(Func)sqr,1,0);
  fnew("cub",(Func)cub,1,0);
  y= (double*)malloc(neq*sizeof(double));
  f=  (formu*)malloc(neq*sizeof(formu));
  vars=(char*)malloc((neq+2)*sizeof(char));
  if (neq<1 || neq>20) { fprintf(fileout,"| ode: error, neq=%d\n",neq); return 1; }
  fprintf(fileout,"|----- Ode 1.1 by Keith Briggs --- compiled%s -------\n",OdeVersion);
  if (neq==1) fprintf(fileout,"|   one equation,");
         else fprintf(fileout,"|   %d equations,",neq);
  vars[0]='x'; vars[neq+1]=0; /* terminate string */
  fprintf(fileout," variables: independent x, dependent ");
  for(i=0; i<neq; i++) { vars[i+1]='a'+i; fprintf(fileout," %c ",vars[i+1]); }
  fprintf(fileout,"\n");
  for(i=0; i<neq; i++) {
    fprintf(fileout,"|     d%c/dx=%s\n",vars[i+1],argv[i+1]);
    f[i]=translate(argv[i+1],vars,&len,&err);
    if (!fnot_empty(f[i]) || err>0) { 
      fprintf(fileout,"|");  
      for(i=0; i<err+12; i++) fprintf(fileout," "); fprintf(fileout,"^\n");
      fprintf(fileout,"|   syntax error near marked character, aborting.\n");  
      return 1;
    }
    y[i]=atof(argv[i+neq+1]);
  }
  fprintf(fileout,"|   initial values ");
  for(i=0; i<neq; i++) fprintf(fileout,"%c=%g ",vars[i+1],y[i]);
  fprintf(fileout,"\n");
  x=atof(argv[argc-3]); 
  xend=atof(argv[argc-2]); 
  nout=atoi(argv[argc-1]);
  if (nout<=0) {
    fprintf(fileout,"|   nout = %d is invalid\n",nout);
    return 1;
  }
  dx=(hmax=(xend-x))/nout; h=0.1*hmax;
  fprintf(fileout,"|   x range = %g ... %g\n",x,xend);
  fprintf(fileout,"|   relative tolerance=%g\n",rtoler); 
  fprintf(fileout,"|   absolute tolerance=%g\n",atoler); 
  fprintf(fileout,"|-----------------------------------------------------------------------\n");
  printf("%e ",x); 
  for (j=0; j<neq; j++) printf("%e ",y[j]); 
  printf("\n");
  if (dx==0.0) return 0;
  inform=dopri5(neq,fcn,x,y,xend,&rtoler,&atoler,itoler,solout,iout,fileout,
    uround,safe,fac1,fac2,beta,hmax,h,nmax,meth,nstiff,nrdens,icont,licont
#ifdef EIGHT
    ,cd
#endif
    );
  if (nlines<nout && (dx>0?x<xend:x>xend)) {
    printf("%e ",xend); 
    for (j=0; j<neq; j++) printf("%e ",y[j]); 
    printf("\n");
  }
  if (inform!=1) return inform+1; /* inform=1 => ok */
  return 0;
} 

void fcn(unsigned n, double x, double *y, double *dy
#ifdef EIGHT
  ,double *cd
#endif
) { 
  int i;
  make_var('x',x); for (i=0; i<n; i++) make_var('a'+i,y[i]);
  for (i=0; i<n; i++) dy[i]=fval_at(f[i]);
  return;
}

void solout(long nr, double xold, double x, double* y, unsigned n, int* irtrn) {
  int i,i1,i2,j; double xx;
  i1=xold/dx; i2=x/dx;
  if (dx>0) for (i=i1; i<=i2; i++) {
    xx=dx*i;
    if ((xx>xold) && (xx<x)) {
      printf("%e ",xx); 
      for (j=0; j<n; j++) printf("%e ",contd5(j,xx)); 
      printf("\n");
      nlines++;
    }
  }
  else /* dx<0 */ for (i=i1; i<=i2; i++) {
    xx=dx*i;
    if ((xx<xold) && (xx>=x)) {
      printf("%e ",xx); 
      for (j=0; j<n; j++) printf("%e ",contd5(j,xx)); 
      printf("\n");
      nlines++;
    }
  }
}
