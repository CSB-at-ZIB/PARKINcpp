// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-03-30 td
// last changed:
//
#include <cmath>
#include <cstdlib> // for rand(), RAND_MAX, ...

#include "UserFunc.h"

using namespace PARKIN;

//---------------------------------------------------------------------------

double
PARKIN::randu()
// return a random number ( uniformly distributed in [0,1[ )
{
    return ( std::rand() / (1.0 + RAND_MAX) );
}

//---------------------------------------------------------------------------

double
PARKIN::randn()
{
    double u = 0.0, v = 0.0;
    while ( u == 0.0 ) u = randu();
    while ( v == 0.0 ) v = randu();
    return std::sqrt(-2*std::log(u)) * std::cos(2*M_PI*v);
    // or: std::sqrt(-2*std::log(u)) * std::sin(2*M_PI*v);
}

//---------------------------------------------------------------------------
