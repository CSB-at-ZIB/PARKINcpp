cmake_minimum_required (VERSION 2.6)

# This include path is necessary on Mac due to the non standard 
# include<malloc.h>
# TODO: replace malloc.h by stdlib.h in the source code and remove this line.
INCLUDE_DIRECTORIES(/usr/include/sys/)

add_library(Ode
            dop853b.c 
            sum.c)