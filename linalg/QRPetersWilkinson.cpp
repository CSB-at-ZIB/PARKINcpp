// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-01-19 td
// last changed:
//

#include "QRPetersWilkinson.h"

using namespace PARKIN;

//---------------------------------------------------------------------

bool
QRPetersWilkinson::prepare(Matrix const& qrA, Vector const& diag, long rank)
{
    long n = qrA.nc();

    // _rank = 0;

    if ( (0 < rank) && (rank < n) )
    {
        // _rank = rank;
        return false;
    }

    return false;
}

//---------------------------------------------------------------------

bool
QRPetersWilkinson::solvePInv(long rank, Vector& v)
{
    return false;
}

//---------------------------------------------------------------------

bool
QRPetersWilkinson::solveR(Matrix const& qrA, Vector const& diag,
                          long rank, Vector const& b,
                          Vector& v)
{
    return false;
}

//---------------------------------------------------------------------

Real
QRPetersWilkinson::projectIntoSubspace(long rank, Vector& w)
{
    Real del = 0.0;

    return del;
}

//---------------------------------------------------------------------
