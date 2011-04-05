// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-01-19 td
// last changed:
//

#include "QRMoorePenrose.h"

using namespace PARKIN;

//---------------------------------------------------------------------

bool
QRMoorePenrose::prepare(Matrix const& qrA, Vector const& diag, long rank)
{
    long n = qrA.nc();

    // _rank = rank;

//std::cerr << "entering:  ### QRMoorePenrose::prepare() ###\n";

    if ( (0 < rank) && (rank < n) )
    {
	    long    ker = n - rank;
	    long    kr1 = ker + 1;
	    long    i2 = 0, rk1 = rank + 1;
        Real    sh, s = 0.0;
        Vector  v;

        _qrAH.zeros(n,n);
        v.zeros(n);

        // compute I + V V^T = L L^T
        for (long j = 1; j <= ker; ++j)
        {
            // compute V(j,1..p) = R^(-1) S(1..p,j) (the j-th row of V)
            for (long ii = 1; ii <= rank; ++ii)
            {
                long i = rk1 - ii;
                // index i running from  _rank  to  1.
                // i = rk1 - ii;
                s = qrA(i,j+rank);          // note the shift s(i,j) := _qrA(i,j+rank)
                if (ii != 1)
                {
                    sh = 0.0;
                    for (long l = i2; l <= rank; ++l)
                    {
                        sh += qrA(i,l) * v(l);
                    }
                    s = s - sh;
                }
                i2 = i;
                v(i) = s / diag(i);
                _qrAH(j,n-ii+1) = v(i);     // note the swap j,i in the index order !!!
            }
        }

        v.zero();

        for (long j = kr1; j <= n; ++j)   // now work with W := V^T
        {
            // v.zero();
            for (long l = 1; l <= ker; ++l)
            {
                v(l) = _qrAH(l,j);
            }

            // compute L^T(1..j,j), i.e. the j-th column of L^T
            // such that   I + W^T W = I + V V^T = L L^T   (see above)
            for (long i = kr1; i <= j; ++i)
            {
                s = 0.0;
                for (long l = 1; l < i; ++l)
                {
                    s += _qrAH(l,i) * v(l);
                }
                if (i != j)
                {
                    v(i) = - s / _qrAH(i,i);
                    _qrAH(i,j) = -v(i);
                }
            }
            if (s > -1.0)
            {
                _qrAH(j,j) = std::sqrt(s + 1.0);
            }
            else
            {
                return false;
            }
        }

        // _rank = rank;
        return true;
    }

    return false;
}

//---------------------------------------------------------------------

bool
QRMoorePenrose::solvePInv(long rank, Vector& v)
{
    long n = _qrAH.nc();
    long ker = n - rank;
    // long kr1 = ker + 1;
    // long rk1 = rank + 1;
    Real s = 0.0;

//std::cerr << "entering:  ### QRMoorePenrose::solvePInv() ###\n";

    for (long j = 1; j <= rank; ++j)
    {
        s = v(j);
        long jj = ker + j;
        for (long l = 1; l < j; ++l)
            s += _qrAH(ker+l,jj) * v(l);
        v(j) = -s / _qrAH(jj,jj);
    }

    for (long j = rank; j > 0; --j)
    {
        s = v(j);
        long jj = ker + j;
        for (long l = j+1; l <= rank; ++l)
            s += _qrAH(jj,ker+l) * v(l);
        v(j) = -s / _qrAH(jj,jj);
    }

    for (long j = 1; j <= ker; ++j)
    {
        s = 0.0;
        for (long l = 1; l <= rank; ++l)
            s += _qrAH(j,ker+l) * v(l);
        v(rank+j) = s;
    }

    return true;
}

//---------------------------------------------------------------------

bool
QRMoorePenrose::solveR(Matrix const& qrA, Vector const& diag,
                       long rank, Vector const& b,
                       Vector& v)
{
    long rk1 = rank+1, j1 = rk1;
    Real sh, s = 0.0;

//std::cerr << "entering:  ### QRMoorePenrose::solveR() ###\n";

    for (long jj = 1; jj <= rank; ++jj)
    {
        long j = rk1 - jj;
        s = b(j);
        if (jj != 1)
        {
            sh = 0.0;
            for (long l = j1; l <= rank; ++l)
                sh += qrA(j,l) * v(l);
            s -= sh;
        }
        j1 = j;
        v(j) = s / diag(j);
    }

    return true;
}

//---------------------------------------------------------------------

Real
QRMoorePenrose::projectIntoSubspace(long rank, Vector& w)
{
    long    n = _qrAH.nc();
    long    ker = n - rank;
    long    kr1 = ker + 1;
    Real    s2, s = 0.0;
    Real    del = 0.0;

    for (long j = kr1; j <= n; ++j)
    {
        s2 = 0.0;
        for (long l = 1; l < j; ++l) s2 += _qrAH(l,j)*w(l);

        s    = ( w(j) - s2 ) / _qrAH(j,j);
        del += s*s;
        w(j) = s;
    }

    return del;
}

//---------------------------------------------------------------------
