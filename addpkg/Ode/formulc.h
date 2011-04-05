/* KMB 98 Oct 14 added Func4,5 */
/* FORMULC.H (c) 1993-97 Harald A. Helfgott
 */
/* This program must be distributed with its corresponding FORMULC.DOC */
/* The full copyright and availability notice is in FORMULC.DOC	      */
/* 	This program is provided "as is", without any explicit or */
/* implicit warranty. */
/* This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.*/

#define UCHAR unsigned char
#define MAXPAR 9 /*KMB*/
                    /* maximum number of parameters */
typedef struct {
 UCHAR *code;
 double *ctable;
} formu;

typedef double (*Func)(double);
typedef double (*Func2)(double,double);
typedef double (*Func3)(double,double,double);
typedef double (*Func4)(double,double,double,double);
typedef double (*Func5)(double,double,double,double,double);
typedef double (*Func6)(double,double,double,double,double,double);
typedef double (*Func7)(double,double,double,double,double,double,double);
typedef double (*Func8)(double,double,double,double,double,double,double,double);
typedef double (*Func9)(double,double,double,double,double,double,double,double,double);
typedef double (*Func0)(void);

formu translate(const char *source, const char *args, int *length,
                                                        int *error);
void destrf(formu);
void make_empty(formu);
int fnot_empty(formu);
const char *fget_error(void);

double fval_at(formu function);
void make_var(char var, double value);
double fval(formu function, char *args, ...);
double f_x_val(formu function, double x);

int fnew(char *name, Func f, int n_of_pars, int varying);
int read_table(int i, char *name, int *n_of_pars, int *varying);
int where_table(char *name);
int fdel(char *name);

double rnd(void);
void rnd_init(void);
/* If MY_RND is defined, rnd() and rnd_init() must be defined by the user.*/
/* Otherwise, formulc.c uses the random-number generator r250  */
/* (written by W. L. Maier, S. Kirkpatrick and E. Stoll) */

/* rnd_init is used by formulc.c only if STAND_ALONE is defined.   */
/* If FORMULC is compiled without STAND_ALONE, it is the user's    */
/* responsibility to initialize her random-number generator.       */
