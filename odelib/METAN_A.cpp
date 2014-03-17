// Copyright (C) 2010 -2013
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2013-04-02 td
// last changed:
//
#include <algorithm>

#include "METAN_A.h"
#include "CubicHermiteTrajectory.h"
#include "LinearTrajectory.h"

using namespace PARKIN;

//----------------------------------------------------------------------------
FirstOrderODESystem*    METANWrapper::_ode = 0;
METAN_A*                METANWrapper::_obj = 0;
//----------------------------------------------------------------------------
extern "C"
void
METANWrapper::xfcn(int* n, double* t, double* y, double* dy, int* ifail)
{
    if ( _ode == 0 ) { *ifail = -999; return; }

    //

    _ode -> computeDerivatives( *t, y, dy, ifail );
}
//----------------------------------------------------------------------------
extern "C"
void
METANWrapper::xsout(
                      int*      n,
                      double*   tOld,
                      double*   t,
                      double*   y,
                      double*   dy
                   )
{
    // std::stringstream       s;
    // double*                 yy = y;
    ODESolver::Grid&        tpt = _obj -> getSolutionGridPoints();
    ODESolver::Trajectory&  sol = _obj -> getSolutionTrajectory();
    ODETrajectory*          tra = _obj -> getRawTrajectory();

    if ( tpt.empty() || (tpt.back() == *tOld) )
    {
        dynamic_cast<CubicHermiteTrajectory*>(tra) -> appendHerm(*t,*n,y,dy);
        /// dynamic_cast<LinearTrajectory*>(tra) -> appendLin(*t,*n,y);

        tpt.push_back(*t);
        for (int j = 0; j < *n; ++j)
        {
            sol[j].push_back( *(y++) );
        }
    }

}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
METAN_A::METAN_A() :
    ODESolver(ODE_SOLVER_METAN_A),
    _n(0), _fcn(METANWrapper::xfcn),
    _t0(0.0), _y0(0), _y(0), _tEnd(0.0),
    _hMax(0.0), _h(_inistep), _kFlag(4),
    _sout(METANWrapper::xsout)
{
    _trajectory = new CubicHermiteTrajectory(0);
    /// _trajectory = new LinearTrajectory(0);
}
//----------------------------------------------------------------------------
METAN_A::~METAN_A()
{
    delete[] _y0; _y0 = 0;
    delete[] _y; _y = 0;
    delete _trajectory;
}
//----------------------------------------------------------------------------
int
METAN_A::integrate()
{
    int             rc[1];
    GridIterConst   dBeg = _datPoints.begin();
    GridIterConst   dEnd = _datPoints.end();
    double          h[1];
    double          t[1];

    *rc = _kFlag;
    *h  = _h;
    *t  = _t0;

    for (int j = 0; j < _n; ++j) _y[j] = _y0[j];

    _solPoints.clear();
    _solution.clear();
    _data.clear();
    _trajectory->clear();

    METANWrapper::setObj(*this);

    //

    metan1_(
                &_n, _fcn, t, _y, &_tEnd,
                &_rtol, &_hMax, h, rc,
                _sout
           );

    //
    // Evaluate and save the _data at measurement points
    //
    for (GridIterConst it = dBeg; it != dEnd; ++it)
    {
        double tEval = (double) *it;

        std::vector<Real> yEval = _trajectory->eval( tEval );

        if ( yEval.size() == (unsigned) _n )
        {
            /// for (long j = 0; j < (_n - 1); ++j)
            /// {
            ///     _data[j].push_back( yEval[j+1] );
            /// }
            for (long j = 0; j < _n; ++j)
            {
                 _data[j].push_back( yEval[j] );
            }
        }
    }

    //

    return (*rc >= 0) ? 0 : rc[0];
}
//----------------------------------------------------------------------------
int
METAN_A::integrate(unsigned n, double* yIni,
                  double tLeft, double tRight)
{
    int             rc[1];
    GridIterConst   dBeg = _datPoints.begin();
    GridIterConst   dEnd = _datPoints.end();
    double          h[1];
    double          t[1];
    double          T = tRight;

    if ((unsigned)_n != n) return -99;

    *rc = _kFlag;
    *h  = _h;
    *t  = tLeft;

    _solPoints.clear();
    _solution.clear();
    // _data.clear();

    METANWrapper::setObj(*this);

    //

    metan1_(
                &_n, _fcn, t, yIni, &T,
                &_rtol, &_hMax, h, rc,
                _sout
           );

    //

    GridIterConst gBeg = std::lower_bound(dBeg, dEnd, tLeft);    // gBeg pointing to first element in [dBeg, dEnd[ that does *not* compare less than tLeft
    GridIterConst gEnd = std::upper_bound(dBeg, dEnd, tRight);   // gEnd pointing to first element in [dBeg, dEnd[ that compares greater than tRight

    //
    // Evaluate and save the _data at measuremet points
    //
//std::cerr << std::endl;
//std::cerr << "*** Data point eval ***" << std::endl;
    for (GridIterConst it = gBeg; it != gEnd; ++it)
    {
        double tEval = (double) *it;

        std::vector<Real> yEval = _trajectory->eval( tEval );

//std::cerr << "     t = " << tEval << std::endl;
//std::cerr << " yEval = " << std::endl;
//for (unsigned j = 0; j < yEval.size(); ++j)
//{
//std::cerr <<  "  " << yEval[j];
//}
//std::cerr << std::endl << std::endl;

        if ( yEval.size() == (unsigned) _n )
        {
            /// for (long j = 0; j < (_n - 1); ++j)
            /// {
            ///     _data[j][long(it-dBeg)] = yEval[j+1];
            /// }
            for (long j = 0; j < _n; ++j)
            {
                _data[j][long(it-dBeg)] = yEval[j];
            }
        }
    }
//std::cerr << "***" << std::endl;
    //

    return (*rc >= 0) ? 0 : rc[0];
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
int
METAN_A::integrateSensitivitySystem(unsigned nDAE)
{
    int             rc[1];
    GridIterConst   dBeg = _datPoints.begin();
    GridIterConst   dEnd = _datPoints.end();
    double          h[1];
    double          t[1];
/*
std::cerr << std::endl;
std::cerr << "*** METAN_A::integrateSensitivitySystem() called ***";
std::cerr << std::endl;
*/
    *rc = _kFlag;
    *h  = _h;
    *t  = _t0;

    for (int j = 0; j < _n; ++j) _y[j] = _y0[j];

    _solPoints.clear();
    _solution.clear();
    _data.clear();
    _trajectory->clear();

    METANWrapper::setObj(*this);

    //

    metan1_(
                &_n, _fcn, t, _y, &_tEnd,
                &_rtol, &_hMax, h, rc,
                _sout
           );

    //
    // Evaluate and save the _data at measurement points
    //
    for (GridIterConst it = dBeg; it != dEnd; ++it)
    {
        double tEval = (double) *it;

        std::vector<Real> yEval = _trajectory->eval( tEval );

        if ( yEval.size() == (unsigned) _n )
        {
            /// for (long j = 0; j < (_n - 1); ++j)
            /// {
            ///     _data[j].push_back( yEval[j+1] );
            /// }
            for (long j = 0; j < _n; ++j)
            {
                 _data[j].push_back( yEval[j] );
            }
        }
    }

    //

    return (*rc >= 0) ? 0 : rc[0];
}
//----------------------------------------------------------------------------
int
METAN_A::integrateSensitivitySystem(unsigned nDAE,
                                    unsigned n, double* yIni,
                                    double tLeft, double tRight
                                  )
{
    int             rc[1];
    GridIterConst   dBeg = _datPoints.begin();
    GridIterConst   dEnd = _datPoints.end();
    double          h[1];
    double          t[1];
    double          T = tRight;

/*
std::cerr << std::endl;
std::cerr << "*** METAN_A::integrateSensitivitySystem( n, yIni, ... ) called ***";
std::cerr << std::endl;
*/
    if ((unsigned)_n != n) return -99;

    *rc = _kFlag;
    *h  = _h;
    *t  = tLeft;

    _solPoints.clear();
    _solution.clear();
    // _data.clear();

    METANWrapper::setObj(*this);

    //

    metan1_(
                &_n, _fcn, t, yIni, &T,
                &_rtol, &_hMax, h, rc,
                _sout
           );

    //

    GridIterConst gBeg = std::lower_bound(dBeg, dEnd, tLeft);    // gBeg pointing to first element in [dBeg, dEnd[ that does *not* compare less than tLeft
    GridIterConst gEnd = std::upper_bound(dBeg, dEnd, tRight);   // gEnd pointing to first element in [dBeg, dEnd[ that compares greater than tRight

    //
    // Evaluate and save the _data at measurement points
    //
    for (GridIterConst it = gBeg; it != gEnd; ++it)
    {
        double tEval = (double) *it;

        std::vector<Real> yEval = _trajectory->eval( tEval );

        if ( yEval.size() == (unsigned) _n )
        {
            /// for (long j = 0; j < (_n - 1); ++j)
            /// {
            ///     _data[j][long(it-dBeg)] = yEval[j+1];
            /// }
            for (long j = 0; j < _n; ++j)
            {
                 _data[j][long(it-dBeg)] = yEval[j];
            }
        }
    }

    //

    return (*rc >= 0) ? 0 : rc[0];
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
ODESolver::Grid&
METAN_A::getAdaptiveGridPoints()
{
    return _solPoints;
}
//----------------------------------------------------------------------------
ODESolver::Trajectory&
METAN_A::getAdaptiveSolution()
{
    return _solution;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
std::string
METAN_A::getErrorMessage(int rc)
{
    std::string message;

    message.clear();

    switch(rc)
    {
        case 0:
            message = "Computation successful.";
            break;

        case -1:
            message = "Interval error:  TEND < T.";
            break;

        case -2:
            message = "More than NSTMAX basic integration steps per interval have been performed.";
            break;

        case -3:
            message = "More than JRMAX stepsize reductions occurred per basic integration step.";
            break;

        case -4:
            message = "Stepsize proposal for next basic integration too small.";
            break;

        default:
            message = "Unknown rc number.";
            break;
    }

    return message;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void
METAN_A::setODESystem(
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
METAN_A::setODESystem(
                        FirstOrderODESystem& ode,
                        double              t0,
                        Grid const&         y0,
                        Grid const&         refGrid,
                        double              tEnd,
                        int                 bandwidth
                    )
{
    setODESystem(
                    METANWrapper::xfcn,
                    t0, y0, refGrid, tEnd,
                    bandwidth
                );

    METANWrapper::setODE(ode);
}
//----------------------------------------------------------------------------
void
METAN_A::setODESystem(
                     FcnMetan        fcn,
                     double         t0,
                     Grid const&    y0,
                     Grid const&    refGrid,
                     double         tEnd,
                     int             bandwidth
                    )
{
    if ( (unsigned)_n != y0.size() )
    {
        _n = y0.size();
        delete[] _y0; _y0 = new double[_n];
        delete[] _y;  _y  = new double[_n];
    }
    _fcn  = fcn;
    _t0   = t0;
    _tEnd = tEnd;
    _hMax = (_maxstep <= 0.0) ? std::fabs( tEnd - t0 ) : _maxstep;
    _h    = _inistep;  // _rtol;

    _datPoints = refGrid;
    std::sort( _datPoints.begin(), _datPoints.end() );
    _data.clear();

    long j = 0;
    for (GridIterConst it = y0.begin(); it != y0.end(); ++it)
    {
        _y0[j] = *it;

        _data[j].resize( _datPoints.size() );

        ++j;
    }

    _solPoints.clear();
    _solution.clear();

    _trajectory->clear();
    _trajectory->setDim(_n);

    METANWrapper::setObj(*this);

    _sout = METANWrapper::xsout;
}
//----------------------------------------------------------------------------
