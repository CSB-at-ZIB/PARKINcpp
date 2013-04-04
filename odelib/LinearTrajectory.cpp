// Copyright (C) 2010 - 2013
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2013-04-04 td
// last changed:
//

#include "LinearTrajectory.h"

using namespace PARKIN;

//----------------------------------------------------------------------------
LinearTrajectory::LinearTrajectory(int n) :
    ODETrajectory(),
    _iLow(0), _iHigh(0)
{
    _n = n;
}
//----------------------------------------------------------------------------
LinearTrajectory::~LinearTrajectory()
{
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void
LinearTrajectory::clear()
{
    _trajectory.clear();
}
//----------------------------------------------------------------------------
void
LinearTrajectory::insert(double t, int n, double* y)
{
  insertLin(t,n,y);
}
//----------------------------------------------------------------------------
void
LinearTrajectory::insertLin(double t, int n, double* y)
{
    if ( n != _n ) return;

    TrajData data;

    data.t = t;
    data.y.clear();
    data.dy.clear();
    for (int j = 0; j < n; ++j)
    {
        data.y.push_back( *y++ );
    }

    if ( _trajectory.empty() )
    {
        _iLow = _iHigh = 0;
        _trajectory.push_back( data );

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

    if ( 0 < _iHigh )
    {
        itTra += _iHigh;
    }

    _trajectory.insert(itTra , data);
}
//----------------------------------------------------------------------------
void
LinearTrajectory::appendLin(double t, int n, double* y)
{
    if ( n != _n ) return;

    TrajData  dataTra;

    dataTra.t = t;
    dataTra.y.clear();
    dataTra.dy.clear();
    for (int i = 0; i < n; ++i)
    {
        dataTra.y.push_back( *y++ );
    }

    _trajectory.push_back( dataTra );
}
//----------------------------------------------------------------------------
std::vector<Real>
LinearTrajectory::eval(double t)
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

    if ( _iHigh <= 0 )
    {
        return std::vector<Real>(_n);
    }

    return evalLinear(t, _trajectory[_iHigh-1],
                          _trajectory[_iHigh]);
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
std::vector<Real>
LinearTrajectory::evalLinear(double t,
                              TrajData const& trajLeft,
                              TrajData const& trajRight)
{
    unsigned          n        = _n;

    double            tFrac    = (t - trajLeft.t) /
                                   (trajRight.t - trajLeft.t);
    double            tFrac_1  = 1.0 - tFrac;

    std::vector<Real> yEval(n);

    for (unsigned j = 0; j < n; ++j)
    {
        double tmp = tFrac_1*trajLeft.y[j] + tFrac*trajRight.y[j];

        yEval[j] = tmp;
    }

    return yEval;
}
//----------------------------------------------------------------------------

