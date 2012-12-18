// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-12-13 td
// last changed:
//
#include <algorithm>
#include "DOP853.h"

using namespace PARKIN;

DOP853* DOP853Wrapper::_obj = 0;

//----------------------------------------------------------------------------
DOP853::DOP853() :
    ODESolver(),
    _n(0), _fcn(0), _x0(0.0), _y0(0), _y(0), _xend(0.0),
    _itol(0),

    _solout(0), _iout(0), _fileout(0),

    _uround(0.0), _safe(0.9), _fac1(0.333), _fac2(6.0), _beta(0.0),
    _hmax(0.0), _h(0.0), _nmax(100000), _meth(1), _nstiff(1000),
    _nrdens(0), _icont(0), _licont(0), _cd(0)
{ }
//----------------------------------------------------------------------------
DOP853::~DOP853()
{
    delete[] _y; _y = 0;
    delete[] _y0; _y0 = 0;
}
//----------------------------------------------------------------------------
int
DOP853::integrate()
{
    _nrdens = _n;
    for (unsigned j = 0; j < _n; ++j) _y[j] = _y0[j];

    _solPoints.clear();
    _solution.clear();
    _data.clear();

    DOP853Wrapper::setObj(*this);

    return dop853(
                    _n, _fcn, _x0, _y, _xend,
                    &_rtol, &_atol, _itol,
                    _solout, _iout, _fileout,
                    _uround, _safe, _fac1, _fac2, _beta,
                    _hmax, _h,
                    _nmax, _meth, _nstiff,
                    _nrdens, _icont, _licont,
                    _cd
                 );
}
//----------------------------------------------------------------------------
int
DOP853::integrate(unsigned n, double* yIni,
                  double xLeft, double xRight)
{
    if (_n != n) return -99;

    _nrdens = _n;

    _solPoints.clear();
    _solution.clear();
    // _data.clear();

    DOP853Wrapper::setObj(*this);

    return dop853(
                    _n, _fcn, xLeft, yIni, xRight,
                    &_rtol, &_atol, _itol,
                    _solout, _iout, _fileout,
                    _uround, _safe, _fac1, _fac2, _beta,
                    _hmax, _h,
                    _nmax, _meth, _nstiff,
                    _nrdens, _icont, _licont,
                    _cd
                 );
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
ODESolver::Grid&
DOP853::getAdaptiveGridPoints()
{
    return _solPoints;
}
//----------------------------------------------------------------------------
ODESolver::Trajectory&
DOP853::getAdaptiveSolution()
{
    return _solution;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void
DOP853::setODESystem(
                     FcnEqDiff      fcn,
                     double         x0,
                     Grid const&    y0,
                     Grid const&    refGrid,
                     double         xend
                    )
{
    long j = 0;
    if ( _n != y0.size() )
    {
        _n = y0.size();
        delete[] _y0; _y0 = new double[_n];
        delete[] _y;  _y  = new double[_n];
    }
    _fcn  = fcn;
    _x0   = x0;
    _xend = xend;

    for (GridIterConst it = y0.begin(); it != y0.end(); ++it) _y0[j++] = *it;

    _datPoints = refGrid;
    std::sort( _datPoints.begin(), _datPoints.end() );
    _data.clear();
    _solPoints.clear();
    _solution.clear();

    DOP853Wrapper::setObj(*this);

    _solout    = DOP853Wrapper::solout;
    _iout      = 2;
    _nrdens    = _n;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
extern "C"
void
DOP853Wrapper::solout(
                      long      nr,
                      double    xold,
                      double    x,
                      double*   y,
                      unsigned  n,
                      int*      irtrn
                     )
{
    // std::stringstream       s;
    double*                 yy = y;
    ODESolver::Grid         dpt = _obj->getDataGridPoints();
    ODESolver::Trajectory&  dat = _obj->getDataTrajectory();
    ODESolver::Grid&        xpt = _obj->getSolutionGridPoints();
    ODESolver::Trajectory&  sol = _obj->getSolutionTrajectory();
    GridIterConst           beg = std::lower_bound( dpt.begin(), dpt.end(), xold);
    GridIterConst           end = std::upper_bound( dpt.begin(), dpt.end(), x);

    for (GridIterConst it = beg; it != end; ++it)
    {
        double xx = (double) *it;

        if ( (xold < xx) && (xx < x) )
        {
            for (unsigned j = 0; j < n; ++j)
            {
                // s << j;
                dat[j].push_back( contd8(j,xx) );
            }
        }
        else if ( xx == x )
        {
            for (unsigned j = 0; j < n; ++j)
            {
                // s << j;
                dat[j].push_back( *(yy++) );
            }
        }
    }

    if ( xpt.empty() || (xpt.back() == xold) )
    {
        xpt.push_back(x);
        for (unsigned j = 0; j < n; ++j) sol[j].push_back( *(y++) );
    }

    *irtrn = 0;
}
//----------------------------------------------------------------------------
