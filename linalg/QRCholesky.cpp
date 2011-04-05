// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-01-19 td
// last changed:
//

#include "QRCholesky.h"

using namespace PARKIN;

//---------------------------------------------------------------------

bool
QRCholesky::prepare(Matrix const& qrA, Vector const& diag, long rank)
{
    long    n = qrA.nc();

    // _rank = rank;

    if ( (0 < rank) && (rank < n) )
    {
        long    i2 = 0, rk1 = rank + 1;
        Real    sh, s = 0.0;
        Vector  v; // , d = diag;

        _qrAH.zeros(n,n);
        v.zeros(n);

        // compute I + V^T V = L L^T
        for (long j = rk1; j <= n; ++j)
        {
            // compute V(1..p,j) = R^(-1) S(1..p,j) (the j-th column of V)
            for (long ii = 1; ii <= rank; ++ii)
            {
                long i = rk1 - ii;
                // index i running from  _rank  to  1.
                // i = rk1 - ii;
                s = qrA(i,j);
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
                _qrAH(i,j) = v(i);
            }
            // compute L^T(1..j,j), i.e. the j-th column of L^T
            // such that   L L^T = I + V^T V   (see above)
            for (long i = rk1; i <= j; ++i)
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
QRCholesky::solvePInv(long rank, Vector& v)
{
    long n = _qrAH.nr();
    long rk1 = rank+1;
    Real sh, s = 0.0;

    for (long j = rk1; j <= n; ++j)
    {
        s = 0.0;
        for (long l = 1; l < j; ++l)
            s += _qrAH(l,j) * v(l);
        v(j) = -s / _qrAH(j,j);
    }
    long j1 = n+1;
    for (long jj = 1; jj <= n; ++jj)
    {
        long j = n - jj + 1;
        s = 0.0;
        if ( jj != 1 )
        {
            sh = 0.0;
            for (long l = j1; l <= n; ++l)
                sh += _qrAH(j,l) * v(l);
            s = sh;
        }
        if ( (jj != 1) && (j <= rank) )
        {
            v(j) -= s;
        }
        else
        {
            j1 = j;
            v(j) = -(s + v(j)) / _qrAH(j,j);
        }
    }

    return true;
}

//---------------------------------------------------------------------

bool
QRCholesky::solveR(Matrix const& qrA, Vector const& diag,
                   long rank, Vector const& b,
                   Vector& v)
{
    long rk1 = rank+1, j1 = rk1;
    Real sh, s = 0.0;

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
QRCholesky::projectIntoSubspace(long rank, Vector& w)
{
    long    n = _qrAH.nc();
    long    rk1 = rank + 1;
    Real    s2, s = 0.0;
    Real    del = 0.0;

    for (long j = rk1; j <= n; ++j)
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
