// Copyright (C) 2010 - 2013
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2013-02-18 td
// last changed:
//

#include "CubicHermiteTrajectory.h"

using namespace PARKIN;

//----------------------------------------------------------------------------
CubicHermiteTrajectory::CubicHermiteTrajectory(int n) :
    ODETrajectory(),
    _iLow(0), _iHigh(0)
{
    _n = n;
}
//----------------------------------------------------------------------------
CubicHermiteTrajectory::~CubicHermiteTrajectory()
{
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void
CubicHermiteTrajectory::clear()
{
    _trajectory.clear();
}
//----------------------------------------------------------------------------
void
CubicHermiteTrajectory::insert(double t, int n, double* y)
{
  double dy[n];

  insertHerm(t,n,y,dy);
}
//----------------------------------------------------------------------------
void
CubicHermiteTrajectory::insertHerm(double t, int n, double* y, double* dy)
{
    if ( n != _n ) return;

    TrajData data;

    data.t = t;
    data.y.clear();
    data.dy.clear();
    for (int j = 0; j < n; ++j)
    {
        data.y.push_back( *y++ );
        data.dy.push_back( *dy++ );
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
CubicHermiteTrajectory::appendHerm(double t, int n, double* y, double* dy)
{
    if ( n != _n ) return;

    TrajData  dataTra;

    dataTra.t = t;
    dataTra.y.clear();
    dataTra.dy.clear();
    for (int i = 0; i < n; ++i)
    {
        dataTra.y.push_back( *y++ );
        dataTra.dy.push_back( *dy++ );
    }

    _trajectory.push_back( dataTra );
}
//----------------------------------------------------------------------------
std::vector<Real>
CubicHermiteTrajectory::eval(double t)
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

    return evalHermite(t, _trajectory[_iHigh-1],
                           _trajectory[_iHigh]);
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
std::vector<Real>
CubicHermiteTrajectory::evalHermite(double t,
                                    TrajData const& trajLeft,
                                    TrajData const& trajRight)
{
    unsigned          n        = _n;

    double            tFrac    = (t - trajLeft.t) /
                                   (trajRight.t - trajLeft.t);
    double            tFrac2   = tFrac * tFrac;
    double            tFrac_1  = tFrac - 1.0;
    double            tFrac_12 = tFrac_1 * tFrac_1;

    double            Herm0 = (2*tFrac + 1.0) * tFrac_12;
    double            Herm1 = tFrac * tFrac_12;
    double            Herm2 = - tFrac2 * (2*tFrac - 3.0);
    double            Herm3 = tFrac2 * tFrac_1;

    std::vector<Real> yEval(n);

    for (unsigned j = 0; j < n; ++j)
    {
        double tmp;

        tmp  =  trajLeft.y[j]*Herm0 +  trajLeft.dy[j]*Herm1;
        tmp += trajRight.y[j]*Herm2 + trajRight.dy[j]*Herm3;

        yEval[j] = tmp;
    }

    return yEval;
}
//----------------------------------------------------------------------------

