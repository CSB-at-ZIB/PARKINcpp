// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-12-13 td
// last changed:
//
#include <algorithm>
#include "DOP853.h"

using namespace PARKIN;

//----------------------------------------------------------------------------
FirstOrderODESystem*    DOP853Wrapper::_ode = 0;
DOP853*                 DOP853Wrapper::_obj = 0;
//----------------------------------------------------------------------------
extern "C"
void
DOP853Wrapper::xfcn(unsigned n, double t, double* y, double* f,
                    double* cd)
{
    int info;

    if ( _ode == 0 ) { return; }

    _ode -> computeDerivatives( t, y, f, &info );
}
//----------------------------------------------------------------------------
extern "C"
void
DOP853Wrapper::solout(
                      long      nr,
                      double    tOld,
                      double    t,
                      double*   y,
                      unsigned  n,
                      int*      irtrn
                     )
{
    // std::stringstream       s;
    double*                 yy = y;
    ODESolver::Grid         dpt = _obj -> getDataGridPoints();
    ODESolver::Trajectory&  dat = _obj -> getDataTrajectory();
    ODESolver::Grid&        tpt = _obj -> getSolutionGridPoints();
    ODESolver::Trajectory&  sol = _obj -> getSolutionTrajectory();
    GridIterConst           beg = std::lower_bound( dpt.begin(), dpt.end(), tOld);
    GridIterConst           end = std::upper_bound( dpt.begin(), dpt.end(), t);

    for (GridIterConst it = beg; it != end; ++it)
    {
        double xx = (double) *it;

        if ( (tOld < xx) && (xx < t) )
        {
            for (unsigned j = 0; j < n; ++j)
            {
                // s << j;
                dat[j].push_back( contd8(j,xx) );
            }
        }
        else if ( xx == t )
        {
            for (unsigned j = 0; j < n; ++j)
            {
                // s << j;
                dat[j].push_back( *(yy++) );
            }
        }
    }

    if ( tpt.empty() || (tpt.back() == tOld) )
    {
        tpt.push_back(t);
        for (unsigned j = 0; j < n; ++j) sol[j].push_back( *(y++) );
    }

    *irtrn = 0;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
DOP853::DOP853() :
    ODESolver(ODE_SOLVER_DOP853),
    _n(0), _fcn(DOP853Wrapper::xfcn),
    _t0(0.0), _y0(0), _y(0), _tEnd(0.0), _itol(0),

    _solout(DOP853Wrapper::solout), _iout(0), _fileout(0),

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
    int rc = 0;

    _nrdens = _n;
    for (unsigned j = 0; j < _n; ++j) _y[j] = _y0[j];

    _solPoints.clear();
    _solution.clear();
    _data.clear();

    DOP853Wrapper::setObj(*this);

    rc = dop853(
                    _n, _fcn, _t0, _y, _tEnd,
                    &_rtol, &_atol, _itol,
                    _solout, _iout, _fileout,
                    _uround, _safe, _fac1, _fac2, _beta,
                    _hmax, _h,
                    _nmax, _meth, _nstiff,
                    _nrdens, _icont, _licont,
                    _cd
                 );

    return (rc == 1) ? 0 : rc;
}
//----------------------------------------------------------------------------
int
DOP853::integrate(unsigned n, double* yIni,
                  double tLeft, double tRight)
{
    int rc = 0;

    if (_n != n) return -99;

    _nrdens = _n;

    _solPoints.clear();
    _solution.clear();
    // _data.clear();

    DOP853Wrapper::setObj(*this);

    rc = dop853(
                    _n, _fcn, tLeft, yIni, tRight,
                    &_rtol, &_atol, _itol,
                    _solout, _iout, _fileout,
                    _uround, _safe, _fac1, _fac2, _beta,
                    _hmax, _h,
                    _nmax, _meth, _nstiff,
                    _nrdens, _icont, _licont,
                    _cd
                 );

    return (rc == 1) ? 0 : rc;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
int
DOP853::integrateSensitivitySystem(unsigned nDAE)
{
    int rc = 0;
/*
std::cerr << std::endl;
std::cerr << "*** DOP853::integrateSensitivitySystem() called ***";
std::cerr << std::endl;
*/
    _nrdens = _n;
    for (unsigned j = 0; j < _n; ++j) _y[j] = _y0[j];

    _solPoints.clear();
    _solution.clear();
    _data.clear();

    DOP853Wrapper::setObj(*this);

    rc = dop853(
                    _n, _fcn, _t0, _y, _tEnd,
                    &_rtol, &_atol, _itol,
                    _solout, _iout, _fileout,
                    _uround, _safe, _fac1, _fac2, _beta,
                    _hmax, _h,
                    _nmax, _meth, _nstiff,
                    _nrdens, _icont, _licont,
                    _cd
                 );

    return (rc == 1) ? 0 : rc;
}
//----------------------------------------------------------------------------
int
DOP853::integrateSensitivitySystem(unsigned nDAE,
                                    unsigned n, double* yIni,
                                    double tLeft, double tRight
                                  )
{
    int rc = 0;
/*
std::cerr << std::endl;
std::cerr << "*** DOP853::integrateSensitivitySystem( n, yIni, ... ) called ***";
std::cerr << std::endl;
*/
    if (_n != n) return -99;

    _nrdens = _n;

    _solPoints.clear();
    _solution.clear();
    // _data.clear();

    DOP853Wrapper::setObj(*this);

    rc = dop853(
                    _n, _fcn, tLeft, yIni, tRight,
                    &_rtol, &_atol, _itol,
                    _solout, _iout, _fileout,
                    _uround, _safe, _fac1, _fac2, _beta,
                    _hmax, _h,
                    _nmax, _meth, _nstiff,
                    _nrdens, _icont, _licont,
                    _cd
                 );

    return (rc == 1) ? 0 : rc;
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
                        FirstOrderODESystem&   ode,
                        double                 t0,
                        Grid const&            y0,
                        double                 tEnd,
                        int                    bandwidth
                    )
{
    Grid no_refGrid;

    setODESystem( ode, t0, y0, no_refGrid, tEnd, bandwidth );
}
//----------------------------------------------------------------------------
void
DOP853::setODESystem(
                        FirstOrderODESystem& ode,
                        double              t0,
                        Grid const&         y0,
                        Grid const&         refGrid,
                        double              tEnd,
                        int                 bandwidth
                    )
{
    setODESystem(
                    DOP853Wrapper::xfcn,
                    t0, y0, refGrid, tEnd,
                    bandwidth
                );

    DOP853Wrapper::setODE(ode);
}
//----------------------------------------------------------------------------
void
DOP853::setODESystem(
                     FcnEqDiff      fcn,
                     double         t0,
                     Grid const&    y0,
                     Grid const&    refGrid,
                     double         tEnd,
                     int             bandwidth
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
    _t0   = t0;
    _tEnd = tEnd;

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
