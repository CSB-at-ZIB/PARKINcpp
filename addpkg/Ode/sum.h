#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <math.h>
int comparedown(const void *x, const void *y);
double sum0(const unsigned n, ...);
double sum2(const unsigned n, ...);
#define x86_FIX \
  unsigned short __old_cw, __new_cw; \
  asm volatile ("fnstcw %0":"=m" (__old_cw)); \
  __new_cw = (__old_cw & ~0x300) | 0x200; \
  asm volatile ("fldcw %0": :"m" (__new_cw));
#define END_x86_FIX  asm volatile ("fldcw %0": :"m" (__old_cw));
