// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-08-24 td
// last changed:
//

#include "LIMEXTrajectory.h"

using namespace PARKIN;

//----------------------------------------------------------------------------
LIMEXTrajectory::LIMEXTrajectory(int n) :
    ODETrajectory(),
    _iLow(0), _iHigh(0), _hermite()
{
    _n = n;
}
//----------------------------------------------------------------------------
LIMEXTrajectory::~LIMEXTrajectory()
{
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void
LIMEXTrajectory::clear()
{
    _trajectory.clear();
    _hermite.clear();
}
//----------------------------------------------------------------------------
void
LIMEXTrajectory::insert(double t, int n, double* y)
{
    if ( n != _n ) return;

    TrajData data;

    data.t = t;
    data.y.clear();
    for (int j = 0; j < n; ++j)
    {
        data.y.push_back( *y++ );
    }

    if ( _trajectory.empty() )
    {
        _iLow = _iHigh = 0;
        _trajectory.push_back( data );
        _hermite.push_back( HermiteData() );

        return;
    }

    if ( (_trajectory[_iLow].t <= t) && (t <= _trajectory[_iHigh].t) )
    {
        if ( (_trajectory[_iLow].t == t) || (_trajectory[_iHigh].t == t) )
        {
            return;
        }
    }
    else
    {
        _iLow = 0;
        _iHigh = _trajectory.size()-1;

        while ( (_iHigh - _iLow) > 1 )
        {
            unsigned idx = (_iHigh + _iLow) / 2;

            if ( t < _trajectory[idx].t)
            {
                _iHigh = idx;
            }
            else
            {
                _iLow = idx;
            }
        }
    }

    std::vector<TrajData>::iterator    itTra = _trajectory.begin();
    std::vector<HermiteData>::iterator itHer = _hermite.begin();

    if ( 0 < _iHigh )
    {
        itTra += _iHigh;
        itHer += _iHigh;
    }

    _trajectory.insert(itTra , data);
    _hermite.insert(itHer , HermiteData());
}
//----------------------------------------------------------------------------
void
LIMEXTrajectory::append(double t, int n, double* y,
                        int k, int N, double coeff[], double t1, double t2)
{
    if ( n != _n ) return;

    TrajData    dataTra;
    HermiteData dataHer;

    dataTra.t = t;
    dataTra.y.clear();
    for (int i = 0; i < n; ++i)
    {
        dataTra.y.push_back( *y++ );
    }

    dataHer.t1 = t1;
    dataHer.t2 = t2;
    dataHer.kOrder = k;
    dataHer.nDim = n;
    dataHer.coeff.clear();
    for (int j = 0; j < k+2; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            dataHer.coeff.push_back( coeff[N*j + i] );
        }
    }

    _trajectory.push_back( dataTra );
    _hermite.push_back( dataHer );
}
//----------------------------------------------------------------------------
std::vector<Real>
LIMEXTrajectory::eval(double t)
{
    if ( _trajectory.empty() )
    {
        return std::vector<Real>();
    }

    if ( (_trajectory[_iLow].t <= t) && (t <= _trajectory[_iHigh].t) )
    {
        if ( _trajectory[_iLow].t == t )
        {
            return _trajectory[_iLow].y;
        }

        if ( _trajectory[_iHigh].t == t )
        {
            return _trajectory[_iHigh].y;
        }
    }
    else
    {
        _iLow = 0;
        _iHigh = _trajectory.size()-1;

        while ( (_iHigh - _iLow) > 1 )
        {
            unsigned idx = (_iHigh + _iLow) / 2;

            if ( t < _trajectory[idx].t )
            {
                _iHigh = idx;
            }
            else
            {
                _iLow = idx;
            }
        }
    }

    return evalHermite(t, _hermite[_iHigh]);
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
std::vector<Real>
LIMEXTrajectory::evalHermite(double t, HermiteData const& herm)
{
    unsigned          n = herm.nDim;
    unsigned          k = herm.kOrder;
    double            tFrac = (t - herm.t1) / (herm.t2 - herm.t1);
    double            tmp = tFrac - 1.0;
    std::vector<Real> yEval(n);

    // -----------------------------------------------------------------------
    //
    //    Evaluate  the Hermite  interpolation polynomials  by the  Horner
    //    scheme.
    //
    // -----------------------------------------------------------------------

    for (unsigned i = 0; i < n; ++i)
    {
        yEval[i] = herm.coeff[n*(k+1) + i];
    }

    for (unsigned j = 0; j < k; ++j)
    {
        for (unsigned i = 0; i < n; ++i)
        {
            yEval[i] *= tmp;
            yEval[i] += herm.coeff[n*(k-j) + i];
        }
    }

    for (unsigned i = 0; i < n; ++i)
    {
        yEval[i] *= tFrac;
        yEval[i] += herm.coeff[i];
    }

    return yEval;
}
//----------------------------------------------------------------------------

/*
    c
    c-----------------------------------------------------------------------
    c
    c     Evaluate  the Hermite  interpolation polynomials  by the  Horner
    c     scheme.
    c
    c-----------------------------------------------------------------------
    c
          tmp = tFac - one
    c
          call dcopy ( n, Dense(1,k+2), 1, y_Interp, 1 )
    c
          do j = 1, k
             do i = 1, n
                y_Interp(i) = Dense(i,k+2-j) + tmp * y_Interp(i)
             end do
          end do
    c
          do i = 1, n
             y_Interp(i) = Dense(i,1) + tFac * y_Interp(i)
          end do
*/
