/* KMB 98 Oct 14 added Func4,Func5,etc */
/* KMB 99 Jan 05 allowed upper case variables */

/*FORMULC.C 2.22           as of 2/19/98        */
/*A byte-compiler for mathematical functions */
/*Copyright (c) 1993-98 by Harald Helfgott        */

/* 	This program is provided "as is", without any explicit or */
/* implicit warranty. */
/* This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.*/
/* Programmer's Address:
	    Harald Helfgott
	    MB 1807, Brandeis University
	    P.O. Box 9110
	    Waltham, MA 02254-9110
	    U.S.A.
	    hhelf@cs.brandeis.edu */

/* Many thanks to Dr. Ralf Grosse-Kunstleve (ralf@kristall.erdw.ethz.ch),
and other contributors */

/* Needed for HP-UX 10.2 */
#ifdef HPUX 
#include <sys/_inttypes.h>
#endif

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include <ctype.h>
#ifndef HPUX
#include <time.h>
#endif
#include "formulc.h"

static double pi(void);
static double value(formu function);
static const char *i_error; /*pointer to the character in source[]
                        that causes an error */
#define Max_ctable 255
   /*maximum number of items in a table of constants */
   /* Max_ctable must be less than 256 */
static int i_pctable; /* number of items in a table of constants -
                         used only by the translating functions */
static double *i_ctable; /*current table of constants -
                           used only by the translating functions */
static UCHAR *i_trans(UCHAR *function, char *begin, char *end);
static char  *my_strtok(char *s);
static UCHAR *comp_time(UCHAR *function, UCHAR *fend, int npars);

static char *errmes = NULL;
static void fset_error(char *);

//KMB static double param['z'-'a'+1];
static double param['z'-'A'+1];

typedef struct {
  char *name;
  Func f;    /* pointer to function*/
  int n_pars; /* number of parameters (0, 1, 2 or 3) */
  int varying; /* Does the result of the function vary
                  even when the parameters stay the same?
                  varying=1 for e.g. random-number generators. */
} formu_item;
#define TABLESIZE 256
#define STD_LIB_NUM 13
static formu_item ftable[TABLESIZE]=
{
  {"exp", exp,1,0},
  {"ln",  log,1,0},
  {"sin", sin,1,0},
  {"cos", cos,1,0},
  {"tan", tan,1,0},
  {"asin", asin,1,0},
  {"acos", acos,1,0},
  {"atan", atan,1,0},
  {"atan2",(Func) atan2,2,0},
  {"abs",  fabs,1,0},
  {"sqrt",  sqrt,1,0},
  {"pi", (Func) pi,0,0},
  {"rnd", (Func) rnd, 0, 1}, /*returns a random number from 0 to 1 */
  {NULL,NULL,0}};
/* please, don't use this array directly;
   methods are provided for its manipulation */

/*********************************************************/
/* The following routines manipulate the table of functions */

int read_table(int i, char *name, int *n_pars, int *varying)
/* returns 1 if succesful */
/* returns 0 otherwise */
{
 if(!ftable[i].f) {
  fset_error("index out of bounds");
  return 0;
 }
 else {
  strcpy(name,ftable[i].name);
  *n_pars = ftable[i].n_pars;
  *varying = ftable[i].varying;
  fset_error(NULL);
  return 1;
 }
}

int where_table(char *name)
/* If the function exists, where_table() returns the index of its name
    in the table. Otherwise, it returns -1. */
{
 formu_item *table_p;

 for(table_p=ftable; table_p->f != NULL &&
        strcmp(name,table_p->name); table_p++)
   ;
 if(table_p->f == NULL) /*The end of the table has been reached,
                 but name[] is not there. */
  {
    fset_error("function not found");
    return -1;
  }
 else {
   fset_error(NULL);
   return table_p - ftable;
 }
}

int fdel(char *name)
/* If the function exists, it is deleted and a non-negative value
    is returned. */
/* Otherwise, -1 is returned. */
/* Original library functions may not be deleted. */
{
 int place;
 formu_item *scan;

 if((place=where_table(name)) == -1)
  return -1; /* there is an error message already */
 if(place<STD_LIB_NUM) {
  fset_error("original functions may not be deleted");
  return -1;
 }
 free(ftable[place].name);
 for(scan = &ftable[place]; scan->f!=NULL; scan++) {
  scan->name  =  (scan+1)->name;
  scan->f     =  (scan+1) -> f;
  scan->n_pars = (scan+1) -> n_pars;
 }
 fset_error(NULL);
 return scan-ftable;
} /*end of fdel */

int fnew(char *name, Func f, int n_pars, int varying)
/* 0 is rendered if there is an error */
/* 1 is rendered otherwise */
{
 formu_item *where;

 if(n_pars<0 || n_pars>9/*KMB*/) {
  fset_error("invalid number of parameters");
  return 0;
 }
 for(where=ftable; where->f != NULL && strcmp(name,where->name); where++);
 if(where->f != NULL) {
  where->f=f;
  where->varying = varying;
  where->n_pars = n_pars;   /*old function is superseded */
  fset_error(NULL);
  return 1;
 } else if((where-ftable) >= TABLESIZE-1) {
  fset_error("function table full");
  return 0;
 }
 else {
  where->name = (char *) calloc(strlen(name)+1, sizeof(char));
  if(where->name==NULL) {
    fset_error("no memory");
    return 0;
  }
  strcpy(where->name,name);
  where->f=f;
  where->varying = varying;
  where->n_pars = n_pars;
  fset_error(NULL);
  return 1;
 }
}  /* end of fnew */

/***********************************************************/
/* Error functions                                         */

static void fset_error(char *s)
/* fset_error(NULL) and fset_error("") erase
   any previous error message */
{
 if (s == NULL || *s == '\0') errmes = NULL; /* an empty error message means
                                                   that there is no error */
 else errmes = s;
}

const char *fget_error(void)
{
 return errmes;
}

/**********************************************************/
/* Evaluating functions                                   */

double fval_at(formu function)
{
  fset_error(NULL);
  return value(function);
}

void make_var(char var, double value)
/*for use with fval_at */
/* make_var('x',3); makes x=3 */
{
  param[var-'a']=value;
}

double f_x_val(formu function, double x)
{
 fset_error(NULL);
 param['x'-'a']=x;
 return value(function);
}

double fval(formu function, char *args, ...)
{
 va_list ap;
 double result;

 fset_error(NULL);
 va_start(ap, args);
 while(*args)
  param[(*args++)-'a'] = va_arg(ap, double);
 va_end(ap);
 result=value(function);
 return result;
}


#define BUFSIZE 500
/* bufsize is the size of the stack for double value(formu func) */
static double value(formu func)
{
 double buffer[BUFSIZE];
 register double *bufp = buffer;
          /* points to the first free space in the buffer */
 double x,y,z,w1,w2,w3,w4,w5,w6; /*KMB*/
 register double result;
 register UCHAR *function=func.code;
 register double *ctable=func.ctable;

 if(!function) {
   fset_error("empty coded function");
   return 0; /* non-existent function; result of
                            an unsuccesful call to translate */
 }
 for(;;) {
   switch(*function++) {
    case '\0':goto finish; /* there is a reason for this "goto":
                              this function must be as fast as possible */
    case 'D': *bufp++ = ctable[*function++];
              break;
    case 'V': *bufp++ = param[(*function++)-'a'];
              break;
    case 'M':result = -(*--bufp);
             *bufp++ = result;
             break;
    case '+':y = *(--bufp);
             result = y + *(--bufp);
             *bufp++ = result;
          break;
    case '-':y = *--bufp;
             result= *(--bufp) - y;
             *bufp++ = result;
             break;
    case '*':y = *(--bufp);
             result = *(--bufp) * y;
             *bufp++ = result;
             break;
    case '/':y = *--bufp;
             result = *(--bufp) / y;
             *bufp++ = result;
             break;
    case '^':y = *--bufp;
             result = pow(*(--bufp),y);
             *bufp++ = result;
             break;
    case 'F':switch(ftable[*function].n_pars) {
               case 0:*bufp++ = ((Func0)ftable[*function++].f)();
                      break;
               case 1:x = *--bufp;
                      *bufp++ = ftable[*function++].f(x);
                      break;
               case 2:y = *--bufp;
                      x = *--bufp;
                      *bufp++ = ((Func2)ftable[*function++].f)(x,y);
                      break;
               case 3:z = *--bufp;
                      y = *--bufp;
                      x = *--bufp;
                      *bufp++ = ((Func3)ftable[*function++].f)(x,y,z);
                      break;
               case 4: /*KMB*/
	              w1 = *--bufp;
	              z = *--bufp;
                      y = *--bufp;
                      x = *--bufp;
                      *bufp++ = ((Func4)ftable[*function++].f)(x,y,z,w1);
                      break;
               case 5: /*KMB*/
	              w2 = *--bufp;
	              w1 = *--bufp;
	              z = *--bufp;
                      y = *--bufp;
                      x = *--bufp;
                      *bufp++ = ((Func5)ftable[*function++].f)(x,y,z,w1,w2);
                      break;
               case 6: /*KMB*/
	              w3 = *--bufp;
	              w2 = *--bufp;
	              w1 = *--bufp;
	              z = *--bufp;
                      y = *--bufp;
                      x = *--bufp;
                      *bufp++ = ((Func6)ftable[*function++].f)(x,y,z,w1,w2,w3);
                      break;
               case 7: /*KMB*/
	              w4 = *--bufp;
	              w3 = *--bufp;
	              w2 = *--bufp;
	              w1 = *--bufp;
	              z = *--bufp;
                      y = *--bufp;
                      x = *--bufp;
                      *bufp++ = ((Func7)ftable[*function++].f)(x,y,z,w1,w2,w3,w4);
                      break;
               case 8: /*KMB*/
	              w5 = *--bufp;
	              w4 = *--bufp;
	              w3 = *--bufp;
	              w2 = *--bufp;
	              w1 = *--bufp;
	              z = *--bufp;
                      y = *--bufp;
                      x = *--bufp;
                      *bufp++ = ((Func8)ftable[*function++].f)(x,y,z,w1,w2,w3,w4,w5);
                      break;
               case 9: /*KMB*/
	              w6 = *--bufp;
	              w5 = *--bufp;
	              w4 = *--bufp;
	              w3 = *--bufp;
	              w2 = *--bufp;
	              w1 = *--bufp;
	              z = *--bufp;
                      y = *--bufp;
                      x = *--bufp;
                      *bufp++ = ((Func9)ftable[*function++].f)(x,y,z,w1,w2,w3,w4,w5,w6);
                      break;
               default:fset_error("I2: too many parameters\n");
                       return 0;
              }
             break;
    default:fset_error("I1: unrecognizable operator");
         return 0;
   }
 }
 finish: if((bufp-buffer)!=1)
   fset_error("I3: corrupted buffer");
 return buffer[0];
} /* end of value */

/**********************************************************/
/* Manipulation of data of type formu                     */

void destrf(formu old)
{
 fset_error(NULL);
 free(old.code);
 free(old.ctable);
}

void make_empty(formu f)
{
 fset_error(NULL);
 f.code=NULL;
 f.ctable=NULL;
}

int fnot_empty(formu f)
{
 fset_error(NULL);
 return(f.code!=NULL);
}

/*********************************************************/
/* Interpreting functions                                */

static int isoper(char c)
{
 return ((c == '+') || (c == '-') || (c == '*') || (c == '/')
                    || (c == '^'));
}

static int is_code_oper(UCHAR c)
{
 return ((c == '+') || (c == '-') || (c == '*') || (c == '/')
                    || (c == '^') || (c == 'M'));
}
static int isin_real(char c)
/* + and - are not included */
{
 return (isdigit(c) || c=='.' || c=='E');
}

size_t max_size(const char *source)
/* gives an upper estimate of the size required for
   the coded form of source (including the final '\0') */
/* Take care when modifying: the upper estimate
   returned by max_size must not also accomodate
   *proper* output, but also *improper* output
   which takes place before the translator detects an error. */
{
 int numbers=0;
 int functions=0;
 int operators=0;
 int variables=0;

/* const size_t func_size=2*sizeof(UCHAR); */ /* not needed */
 const size_t var_size=2*sizeof(UCHAR);
 const size_t num_size=sizeof(UCHAR)+sizeof(double);
 const size_t op_size=sizeof(UCHAR);
 const size_t end_size=sizeof('\0');

 const char *scan;

 for(scan=source; *scan; scan++)
  if(isalpha(*scan) && (*scan != 'E'))
  {
    if(isalpha(*(scan+1))) ; /* it is a function name,
                                it will be counted later on */
    else
     if(*(scan+1) == '(')  functions++;
     else variables++;
  }

 if(isoper(*source)) operators++;
 if(*source != '\0')
  for(scan = source+1; *scan; scan++)
   if(isoper(*scan) && *(scan-1) != 'E') operators++;

 /* counting numbers.. */
 scan=source;
 while(*scan)
  if(isin_real(*scan) || ((*scan == '+' || *scan == '-') &&
                           scan>source && *(scan-1)=='E'))
   {numbers++;
    scan++;
    while(isin_real(*scan) || ((*scan == '+' || *scan == '-') &&
                                scan>source && *(scan-1)=='E'))
     scan++;
   }
  else scan++;

 return(numbers*num_size + operators*op_size + functions*num_size
                         + variables*var_size + end_size);
 /*Do you wonder why "function" is multiplied with "num_size"
   and not with func_size? This function calculates an upper-bound
   (i.e. pessimistic) estimate. It supposes that all functions are
   converted into doubles by comp_time. For example, pi() actually
   becomes a double. */
}

/***********************************************************/
/* Interface for interpreting functions                     */

formu translate(const char *sourc, const char *args, int *leng,
                                                        int *error)
{
 UCHAR *result;
 char *source;
 const char *scan, *scarg;
 UCHAR *function;
 UCHAR *nfunc; /* used to free unused heap space */
 size_t size_estim; /* upper bound for the size of the
                                        coded function */
 double *ctable;
 formu returned; /*the value to be returned by the function
                   is stored here */


 i_error=NULL;

 source = (char *) malloc(strlen(sourc) + 1);
 if(source==NULL) {
   fset_error("no memory");
   *leng = 0;
   *error = 0; /* first character */
   returned.code = NULL;
   returned.ctable = NULL;
   return(returned);
 }
 strcpy(source,sourc);
 /* FORMULC's routines must have their own copy of sourc
    because the copy could be temporarily modified.
    Modifying a string constant can make some computers crash. */

 /* search for undeclared parameters */
 for(scan=source; *scan != '\0'; scan++) {
  //KMB if(islower(*scan) && !isalpha(*(scan+1)) &&
  if((islower(*scan)||isupper(*scan)) && !isalpha(*(scan+1)) &&
      (scan==source || !isalpha(*(scan-1))) ) {
   for(scarg=args; *scarg != '\0' && *scarg != *scan; scarg++)
     ;
   if(*scarg == '\0') /*parameter not found */
    {
     i_error = scan;

     fset_error("undeclared parameter");
     *leng = 0;
     *error = i_error - source;
     returned.code=NULL;
     returned.ctable=NULL;
     free(source);
     return(returned);
    }
  }
 }  /* end of search for undeclared... */

 size_estim=max_size(source); /* upper estimate of the size
                                 of the coded function,
                                 which doesn't exist yet */

 if(!(function = (UCHAR *) malloc(size_estim))) {
  /* out of memory */
  fset_error("no memory");
  *leng = 0;
  *error = -1;
  returned.code=NULL;
  returned.ctable=NULL;
  free(source);
  return (returned);
 }

 /*table of memory is cleaned: */
 i_pctable=0;
 if(!(i_ctable = (double *) malloc(Max_ctable * sizeof(double)) )) {
  /* out of memory */
  fset_error("no memory");
  free(function);
  *leng = 0;
  *error = -1;
  returned.code=NULL;
  returned.ctable=NULL;
  free(source);
  return (returned);
 }
 ctable = i_ctable;

 fset_error(NULL);
 /* THIS IS THE CORE STATEMENT */
 result=i_trans(function,(char *) source,(char *) source+strlen(source));

 if(!result || fget_error()) {
  free(function);
  free(i_ctable);
  *leng = 0;
  if(i_error)
   *error = i_error-source;
  else *error = -1; /* internal error or out of memory */
  returned.code=NULL;
  returned.ctable=NULL;
  free(source);
  return (returned);
 }
 else { /* OK */
  *result = '\0';
  *error = -1;
  *leng = result-function;

  /* free unused heap space.. */
  if(((*leng)+1) * sizeof(UCHAR) > size_estim)
   /* one must use (*leng)+1 instead of *leng because '\0'
      has not been counted */
   {
    fset_error("I4: size estimate too small");
    returned.code=NULL;
    returned.ctable=NULL;
    free(source);
    return (returned);
   }
  else if(((*leng)+1) * sizeof(UCHAR) < size_estim) {
    nfunc = (UCHAR *) malloc(((*leng)+1) * sizeof(UCHAR));
      if(nfunc) {
        memcpy( nfunc, function, ((*leng)+1) * sizeof(UCHAR) );
        free(function);
        function=nfunc;
      }
  } /* end of if-else stairs */

  /* free heap space hoarded by i_ctable.. */
  if(i_pctable<Max_ctable) {
    ctable = (double *) malloc(i_pctable * sizeof(double));
    if(ctable) {
      memcpy(ctable, i_ctable, i_pctable * sizeof(double));
      free(i_ctable);
    } else ctable = i_ctable;
  } else ctable = i_ctable;

  returned.code=function;
  returned.ctable=ctable;
  fset_error(NULL);
  free(source);
  return(returned);
 } /* end of OK */
}  /* end of translate */

static UCHAR *comp_time(UCHAR *function, UCHAR *fend, int npars)
  /* calculates at "compile time" */
  /* Postconditions: If the coded expression in *function..*(fend-1)
      can be calculated, its value is stored in *function..*(fend-1) */
  /* comp_time returns a pointer to the first character after the
     end of the coded function; if this function cannot be evaluated
     at compile time, comp_time returns fend, of course.  */
  /* Only memory positions from *function to *comp_time are touched. */
{
  UCHAR *scan;
  UCHAR temp;
  double tempd;
  int i;
  formu trans;

  scan=function;
  for(i=0; i<npars; i++) {
   if(*scan++ != 'D') return fend;
   scan++;
  }

  if(!( ( scan == fend - (sizeof((UCHAR) 'F')+sizeof(UCHAR))
           && *(fend-2) == 'F' && ftable[*(fend-1)].varying == 0) ||
         ( scan == fend - sizeof(UCHAR)
           && is_code_oper(*(fend-1)) ) )
    )
    /* compile-time evaluation is done only if
       1) everything but the ending function consists of doubles
       AND
       2) the function does not vary when its parameters remain the same
          (i.e. random-number generators are not evaluated at compile time)
          */
   return fend;

  temp = *fend;
  *fend = '\0';

  trans.code=function;
  trans.ctable=i_ctable;
  tempd = value(trans);
  *fend = temp;
  *function++ = 'D';
  i_pctable -= npars;
  *function++ = (UCHAR) i_pctable;
  i_ctable[i_pctable++] = tempd;

  return function;
} /* end of comp_time */

static char *my_strtok(char *s)
/* a version of strtok that respects parentheses */
/* token delimiter = comma */
{
 int pars;
 static char *token=NULL;
 char *next_token;

 if(s!=NULL) token=s;
 else if(token!=NULL) s=token;
 else return NULL;

 for(pars=0; *s != '\0' && (*s != ',' || pars!=0); s++) {
   if(*s == '(') ++pars;
   if(*s == ')') --pars;
 }
 if(*s=='\0') {
  next_token=NULL;
  s=token;

  token=next_token;
  return s;
 } else {
  *s = '\0';
  next_token=s+1;
  s=token;

  token=next_token;
  return s;
 }
} /* end of my_strtok */


/************************************************************/
/* Here begins the core of interpretation                   */

#define TWO_OP {                                 \
    if((tempu=i_trans(function,begin,scan)) &&      \
       (temp3=i_trans(tempu,scan+1,end)) ) {       \
    *temp3++ = *scan; /* copies operator */                 \
    temp3 = comp_time(function,temp3,2); /*tries to simplify expression*/ \
   if(fget_error()) return NULL; /* internal error in comp_time */  \
   else return temp3; /* expression has been translated */ \
  } else return NULL; /* something is wrong with the operands */ \
 }

#define ERROR_MEM {    \
   fset_error("no memory"); \
   i_error=NULL;  \
   return NULL;   \
  }
static UCHAR *i_trans(UCHAR *function, char *begin, char *end)
 /* the source is *begin .. *(end-1) */
 /* returns NULL if a normal error or an internal error occured;
    otherwise, returns a pointer
    to the first character after the end of function[] */
 /* i_trans() does not write a '\0' at the end of function[], */
 /* but it MAY touch its end (i.e. *i_trans) without changing it.*/
{
 int pars;     /* parentheses */
 char *scan;
 UCHAR *tempu, *temp3;
 char *temps;
 char tempch;
 double tempd;
 char *endf;     /* points to the opening
                    parenthesis of a function (e.g. of sin(x) ) */
 int n_function;
 int space;
 int i;

 char *paramstr[MAXPAR];
 char *par_buf;

 if(begin>=end) {
  fset_error("missing operand");
  i_error = begin;
  return NULL;
 }

 /* test paired parentheses */
 for(pars=0, scan=begin; scan<end && pars>=0; scan++) {
  if(*scan == '(') pars++;
  else if(*scan == ')') pars--;
 }
 if(pars<0 || pars>0) {
  fset_error("unmatched parentheses");
  i_error = scan-1;
  return NULL;
 }

 /* plus and binary minus */
 for(pars=0, scan=end-1; scan>=begin; scan--) {
  if(*scan == '(') pars++;
  else if(*scan == ')') pars--;
  else if(!pars && (*scan == '+' || ((*scan == '-') && scan!=begin))
                                          /* recognizes unary
                                             minuses */
             && (scan==begin || *(scan-1)!='E') )
          /* be wary of misunderstanding exponential notation */
   break;
 }

 if(scan >= begin) TWO_OP

 /* multiply and divide */
 for(pars=0, scan=end-1; scan>=begin; scan--) {
  if(*scan == '(') pars++;
  else if(*scan == ')') pars--;
  else if(!pars && (*scan == '*' || *scan == '/' ))
   break;
 }

 if(scan >= begin) TWO_OP

 /* unary minus */
 if(*begin == '-') {
   tempu=i_trans(function,begin+1,end);
   if(tempu) {
     *tempu++ = 'M';
     tempu=comp_time(function,tempu,1); /*tries to simplify
                                          expression*/
     if(fget_error()) return NULL; /* internal error in comp_time */
     else return tempu;
   } else return NULL;
 }

 /* power */
 for(pars=0, scan=end-1; scan>=begin; scan--) {
  if(*scan == '(') pars++;
  else if(*scan == ')') pars--;
  else if(!pars && (*scan == '^'))
   break;
 }

 if(scan >= begin) TWO_OP

 /* erase white space */
 while(isspace(*begin))
  begin++;
 while(isspace(*(end-1)))
  end--;

 /* parentheses around the expression */
 if(*begin == '(' && *(end-1) == ')')
  return i_trans(function,begin+1,end-1);

 /* variable */
 //KMB if(end == begin+1 && islower(*begin)) {
 if (end == begin+1 && (islower(*begin)||isupper(*begin))) {
  *function++ = 'V';
  *function++ = *begin;
  return function;
 }

 /* number */
 tempch = *end;
 *end = '\0';
 tempd=strtod(begin,(char**) &tempu);
 *end = tempch;
 if((char*) tempu == end) {
  *function++ = 'D';
  if (i_pctable < Max_ctable)
    {
      i_ctable[i_pctable] = tempd;
      *function++ = (UCHAR) i_pctable++;
    }
  else
    {
      fset_error("too many constants");
      i_error=begin;
      return NULL;
    }
  return function;
 }

 /*function*/
 if(!isalpha(*begin) && *begin != '_')
                        /* underscores are allowed */
 {
  fset_error("syntax error");
  i_error=begin;
  return NULL;
 }
 for(endf = begin+1; endf<end && (isalnum(*endf) || *endf=='_');
                                                           endf++);
 tempch = *endf;
 *endf = '\0';
 if((n_function=where_table(begin)) == -1) {
  *endf = tempch;
  i_error=begin;
  /* error message has already been created */
  return NULL;
 }
 *endf = tempch;
 if(*endf != '(' || *(end-1) != ')') {
  fset_error("improper function syntax");
  i_error=endf;
  return NULL;
 }
 if(ftable[n_function].n_pars==0) {
  /*function without parameters (e.g. pi() ) */
   space=1;
   for(scan=endf+1; scan<(end-1); scan++)
    if(!isspace(*scan)) space=0;
   if(space) {
    *function++ = 'F';
    *function++ = n_function;
    function = comp_time(function-2,function,0);
    if(fget_error()) return NULL; /* internal error in comp_time */
    else return function;
   } else {
    i_error=endf+1;
    fset_error("too many parameters");
    return NULL;
   }
 } else {    /*function with parameters*/
    tempch = *(end-1);
    *(end-1) = '\0';
    par_buf = (char *) malloc(strlen(endf+1)+1);
    if(!par_buf)
     ERROR_MEM;
    strcpy(par_buf, endf+1);
    *(end-1) = tempch;
    /* look at the first parameter */
    for(i=0; i<ftable[n_function].n_pars; i++) {
     if( ( temps=my_strtok((i==0) ? par_buf : NULL) ) == NULL )
      break; /* too few parameters */
     paramstr[i]=temps;
    }
    if(temps==NULL) {
     /* too few parameters */
     free(par_buf);
     i_error=end-2;
     fset_error("too few parameters");
     return NULL;
    }
    if((temps=my_strtok(NULL))!=NULL) {
     /* too many parameters */
     free(par_buf);
     i_error=(temps-par_buf)+(endf+1); /* points to the first character
                                          of the first superfluous
                                          parameter */
     fset_error("too many parameters");
     return NULL;
    }

    tempu=function;
    for(i=0; i<ftable[n_function].n_pars; i++)
     if(!(tempu=i_trans( tempu, paramstr[i],
                                 paramstr[i]+strlen(paramstr[i]) ) ) )
     {
      i_error=(i_error-par_buf)+(endf+1); /* moves i_error to
                                           the permanent copy of the
                                           parameter */
      free(par_buf);
      return NULL; /* error in one of the parameters */
     }
    /* OK */
    free(par_buf);
    *tempu++ = 'F';
    *tempu++ = n_function;
    tempu = comp_time(function,tempu,ftable[n_function].n_pars);
    if(fget_error()) return NULL; /* internal error in comp_time */
    else return tempu;
 }
}

/* Here is the definition of some functions in the FORMULC standard
   library */
static double pi(void)
{
 return 3.14159265358979323846264;
}

#ifndef MY_RND
void r250_init(int seed);
double dr250(void);

double rnd(void)
{
  return dr250();
}

void rnd_init(void)
{
  r250_init(time(NULL));
}
#endif


