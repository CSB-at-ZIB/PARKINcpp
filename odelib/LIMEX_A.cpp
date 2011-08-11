// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-02-04 td
// last changed:
//
#include <algorithm>

#include "LIMEX_A.h"

using namespace PARKIN;

// int _iOptStandard[30] = { 0 };

//----------------------------------------------------------------------------
LIMEX_A::LIMEX_A() :
    ODESolver(),
    _n(0), _fcn(0), _jac(0),
    _t0(0.0), _T(0.0),
    _y0(0), _dy0(0),
    _h(0.0),
    _iOpt(), _rOpt(), _iPos(0),
    _ifail(),
    _kOrder(0), _Dense(), _t1(0.0), _t2(0.0)
{
    initOpt();
}
//----------------------------------------------------------------------------
LIMEX_A::~LIMEX_A()
{
    delete[] _y0;    _y0   = 0;
    delete[] _dy0;   _dy0  = 0;
    delete[] _iPos;  _iPos = 0;
}
//----------------------------------------------------------------------------
void
LIMEX_A::initOpt()
{
    ///

    // iOpt[0] - iOpt[17] must be set by caller on entry.
    // !!! iOpt[15] may be modified; all others are not modified !!!

    _iOpt[0]  =  0;     // Integration monitoring: 0 no output, 1 standard, 2 additional
    _iOpt[1]  =  0;     // Unit number for monitor ( == 6 if iOpt[0] > 0 )
    _iOpt[2]  =  0;     // Solution output: 0 no output, 1 initial&final vaules, 2 additional
    _iOpt[3]  =  0;     // Unit number for solution ( == 6 if iOpt[2] > 0 )
    _iOpt[4]  =  1;     // Singular or non-singualar matrix B: 0 sing, 1 non-sing
    _iOpt[5]  =  0;     // Determination of consistent initial values (CIV): 0 no, 1 determ
    _iOpt[6]  =  0;     // Numerical or analytical Jacobian: 0 num diff approx, 1 analytical

    _iOpt[7]  = -1;     // Lower bandwidth of Jacobian: 0 <= iOpt[7] <= n <= Max_Lower_Diags
    _iOpt[8]  = -1;     // Upper bandwidth of Jacobian: 0 <= iOpt[8] <= n <= Max_Upper_Diags

    _iOpt[9]  =  1;     // Re-use of Jacobian: 0 no re-use, 1 re-use of Jacobian in the following steps
    _iOpt[10] =  0;     // Switch for error tolerances: 0 rTol&aTol scalar, 1 rTol&aTol are vectors
    _iOpt[11] =  1;     // Switch for one step mode: 0 off, 1 return from each step, 2 return only from prescribed steps

    _iOpt[12] = -1;     // Dense output option: 0 off, 1 on equidist pts within interval, 2 on equidist pts within step (# in iOpt[13]), 3 on additional pts
    _iOpt[13] =  0;     // Number of equidistant points if iOpt[12] == 1 or 2
    _iOpt[14] =  0;     // Unit number for dense output (iOpt[14] == 0 suppresses dense output)

    _iOpt[15] = -1;     // Type of call, may be modified! 0 initial call, 1 successive call

    _iOpt[16] =  1;     // Behaviour at t_End: 0 stop exactly at t_End, 1 may compute/use values also for t > t_End
    _iOpt[17] =  0;     // PostScript plot of Jacobian: 0 no plot, j plot at j-th step, -1 plot for initial step (step 0)

    _iOpt[18] =
    _iOpt[19] =
    _iOpt[20] =
    _iOpt[21] =
    _iOpt[22] = -1;     // Not used in LIMEX_A (relevant in LIMEX_B, sparse Jacobians)

    _iOpt[23] =         // on return: Number of function evaluations
    _iOpt[24] =         // on return: Number of fcn evaltions for Jacobian computation
    _iOpt[25] =         // on return: Number of LU decompositions
    _iOpt[26] =         // on return: Number of back-substitions
    _iOpt[27] =         // on return: Number of integration steps
    _iOpt[28] =  0;     // on return: Number of Jacobian evaluations

    _iOpt[29] = -1;     // Not used in LIMEX_A (relevant in LIMEX_B, sparse Jacobians)

    _iOpt[30] =  0;     // !!! Only available in LIMD !!!
                        // Type of left-hand side B: 0 B=id, 1 B=const., 2 variable B
    _iOpt[31] =  1;     // !!! Only available in LIMDHERM !!!
                        // Interpolation mode: 0 no addition output, 1 give additional output (switched on)

    ///

    // rOpt[0] - rOpt[2] must be set by caller
    // all values are NOT modified

    _rOpt[0] =  0.0;    // Maximum allowed stepsize (default t_End - t_Begin if set to 0.0)
    _rOpt[1] = -1.0;    // Maximal distance between two dense output points (0 < rOpt[1] if iOpt[12] == 3)
    _rOpt[2] =  0.0;    // Upper limit for t (used only if iOpt[16] == 1; no upper limit will be considered if rOpt[2] < t_End)

    _rOpt[3] =
    _rOpt[4] = -1.0;    // Not used in LIMEX_A (relevant only in LIMEX_B, sparse Jacobian)

    ///

    return;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
int
LIMEX_A::integrate()
{
    int           maxEqns = (int) MAX_NO_EQNS;
    GridIterConst dBeg = _datPoints.begin();
    GridIterConst dEnd = _datPoints.end();
    double        t[1], z[_n], yEval[_n];
    double        T = 0.0;
    double*       ztmp = 0;

    *t = _t0;
    T = _T;

    for (long j = 0; j < _n; ++j)
    {
        z[j] = _y0[j];
        _dy0[j] = 0.0;
    }

    _solPoints.clear();
    _solution.clear();
    _data.clear();

    //

    _iOpt[11] = 1;      // single step mode ON
    _iOpt[12] = 0;      // no dense output
    _iOpt[13] = 0;
    _iOpt[14] = 0;

    _iOpt[15] = 0;      // type of call of limex_(): 0 initial call

    _iOpt[16] = 0;      // integration for t > T internally switched off
    _iOpt[31] = 1;      // switch on interpolation (for single step mode), only available in LIMDHERM !!!

    //

    _rOpt[0] = _rOpt[1] = _rOpt[2] = 0.0;

    //

    while ( (_ifail[0] == 0) && (*t < T) )
    {
        limdherm_(
                    &_n,
                    _fcn, _jac,
                    t, &T,
                    z, _dy0,
                    &_rtol, &_atol, &_h,
                    _iOpt, _rOpt, _iPos,
                    _ifail,
                    &_kOrder, _Dense, &_t1, &_t2
                 );

        // in single step mode: save the new, adaptive time point
        if (*t <= T)
        {
            _solPoints.push_back( *t );
            ztmp = z;
            for (long j = 0; j < _n; ++j)
            {
                _solution[j].push_back( *ztmp++ );
            }
        }

/*
        // NOTE: This disabled snippet is possible only with adjusted limdherm_()
        //       Basically, the '_Dense' array produced below is then already computed
        //       by a call to comp_herm_() inside limdherm_() !

        int    maxRowTab = MAX_ROW_TAB;
        int    nj[MAX_ROW_TAB] = { 1, 2, 3, 4, 5, 6, 7 };
        int    ipt[MAX_ROW_TAB];
        double tempWork[_n*(MAX_ROW_TAB+1)];

        ipt[0] = 3;
        for (int i = 1; i < maxRowTab; ++i) ipt[i] = ipt[i-1] + nj[i];

        comp_herm_(
                    &_n, &maxEqns, _Dense,
                    &_kOrder,
                    ipt, nj,
                    tempWork
                  );
*/

        // interpolate the measurement time points in ] t1, t2 ] subste of [dBeg, dEnd[
        // requirement: search range [ dBeg, dEnd [ has to be sorted!
        GridIterConst gBeg = std::lower_bound( dBeg, dEnd, _t1);    // gBeg pointing to first element in [dBeg, dEnd[ that does *not* compare less than _t1
        GridIterConst gEnd = std::upper_bound( dBeg, dEnd, _t2);    // gEnd pointing to first element in [dBeg, dEnd[ that compares greater than _t2

        for (GridIterConst it = gBeg; it != gEnd; ++it)
        {
            double tEval = (double) *it;

            if ( (_t1 < tEval) && (tEval < _t2) )
            {
                double tFrac = (tEval - _t1)/(_t2 - _t1);

            /*
             * Note: In the current version unfortunately NOT re-entrant
             *
                hermine_(
                            &_n, &_kOrder,
                            _Dense, &_t1, &_t2,
                            &tEval, yEval
                        );
            */
                eval_herm_(
                            &_n, &maxEqns, _Dense, &_kOrder,
                            &tFrac, yEval
                          );

                ztmp = yEval;
                for (long j = 0; j < _n; ++j)
                {
                    _data[j].push_back( *ztmp++ );
                }
            }
            else if ( tEval == _t2 )
            {
                ztmp = z;
                for (long j = 0; j < _n; ++j)
                {
                    _data[j].push_back( *ztmp++ );  // a simple append is sufficient here presumingly
                }
            }

        } // end for gBeg, gEnd

    } // end while limdherm_


    return _ifail[0];
}
//----------------------------------------------------------------------------
int
LIMEX_A::integrate( unsigned n, double* yIni,
                    double tLeft, double tRight)
{
    if ( n != (unsigned)_n ) return -99;

    int           maxEqns = (int) MAX_NO_EQNS;
    GridIterConst dBeg = _datPoints.begin();
    GridIterConst dEnd = _datPoints.end();
    double        t[1], z[_n], yEval[_n];
    double        T = 0.0;
    double*       ztmp = 0;

    if ( yIni == 0 )
    {
        for (long j = 0; j < _n; ++j)
        {
            z[j] = _y0[j];
            _dy0[j] = 0.0;
        }
    }
    else
    {
        for (long j = 0; j < _n; ++j)
        {
            z[j] = yIni[j];
        }
    }

    *t = tLeft;
    T  = tRight;

    _solPoints.clear();
    _solution.clear();

//std::cerr << "*** Data time points where measurements are taken:" << std::endl;
//    for (GridIterConst it = dBeg; it != dEnd; ++it)
//    {
//std::cerr << *it << ", ";
//    }
//std::cerr << std::endl;

    //

    _iOpt[11] = 1;      // single step mode ON
    _iOpt[12] = 0;      // no dense output (!), but interpolation will be used below
    _iOpt[13] = 0;
    _iOpt[14] = 0;

    if ( _iOpt[15] == -1 ) {
        _iOpt[15] = 0;  // type of call of limex_(): 0 initial call, 1 successive call

        _rOpt[0] = _rOpt[1] = _rOpt[2] = 0.0;
        _h = _rtol;
    }

    _iOpt[16] = 0;      // integration for t > T internally forbidden
    _iOpt[31] = 1;      // switch on interpolation (for single step mode), only available in LIMDHERM !!!

    //

    while ( (_ifail[0] == 0) && (*t < T) )
    {
        limdherm_(
                    &_n,
                    _fcn, _jac,
                    t, &T,
                    z, _dy0,
                    &_rtol, &_atol, &_h,
                    _iOpt, _rOpt, _iPos,
                    _ifail,
                    &_kOrder, _Dense, &_t1, &_t2
                 );

        // in single step mode: save the new, adaptive time point
        if (*t <= T)
        {
//std::cerr << "***" << std::endl;
//std::cerr << "*** Next adaptive time point: " << *t << std::endl;
            _solPoints.push_back( *t );
            ztmp = z;
            for (long j = 0; j < _n; ++j)
            {
                _solution[j].push_back( *ztmp++ );
            }
        }

/*
        // NOTE: This disabled snippet is possible only with adjusted limdherm_()
        //       Basically, the '_Dense' array produced below is then already computed
        //       by a call to comp_herm_() inside limdherm_() !

        int    maxRowTab = MAX_ROW_TAB;
        int    nj[MAX_ROW_TAB] = { 1, 2, 3, 4, 5, 6, 7 };
        int    ipt[MAX_ROW_TAB];
        double tempWork[_n*(MAX_ROW_TAB+1)];

        ipt[0] = 3;
        for (int i = 1; i < maxRowTab; ++i) ipt[i] = ipt[i-1] + nj[i];

        comp_herm_(
                    &_n, &maxEqns, _Dense,
                    &_kOrder,
                    ipt, nj,
                    tempWork
                  );
*/

        // interpolate the measurement time points in ] t1, t2 ] subste of [dBeg, dEnd[
        // requirement: search range [ dBeg, dEnd [ has to be sorted!
        GridIterConst gBeg = std::lower_bound( dBeg, dEnd, _t1);    // gBeg pointing to first element in [dBeg, dEnd[ that does *not* compare less than _t1
        GridIterConst gEnd = std::upper_bound( dBeg, dEnd, _t2);    // gEnd pointing to first element in [dBeg, dEnd[ that compares greater than _t2

//std::cerr << "*** Proceeding with subinterval ] " << _t1 << ", " << _t2 << " ]" << std::endl;
        for (GridIterConst it = gBeg; it != gEnd; ++it)
        {
            double tEval = (double) *it;

            if ( (_t1 < tEval) && (tEval < _t2) )
            {
                double tFrac = (tEval - _t1)/(_t2 - _t1);

//std::cerr << "*** Interpolation at time = " << tEval
//          << " ( tFraction: " << tFrac << "),"
//          << " ( #" << long(it-dBeg) << " ),"
//          << " Order: " << _kOrder << std::endl;
            /*
             * Note: a call to hermine_() destroys _Dense; thus not re-entrant...
             *
                hermine_(
                            &_n, &_kOrder,
                            _Dense, &_t1, &_t2,
                            &tEval, yEval
                        );
            */
                eval_herm_(
                            &_n, &maxEqns, _Dense, &_kOrder,
                            &tFrac, yEval
                          );

//std::cerr << "***    ";
                ztmp = yEval;
                for (long j = 0; j < _n; ++j)
                {
//std::cerr << *ztmp << ", ";
                    _data[j][long(it-dBeg)] = *ztmp++;  // see comment right below
                }
//std::cerr << std::endl;

            }
            else if ( tEval == _t2 )
            {
//std::cerr << "*** Right boundary time = " << _t2 << " ( #" << long(it-dBeg) << " )" << std::endl;
//std::cerr << "***    ";
                ztmp = z;
                for (long j = 0; j < _n; ++j)
                {
//std::cerr << *ztmp << ", ";
                    _data[j][long(it-dBeg)] = *ztmp++;  // allowing overwriting of data points!
                    // _data[j].push_back( z[j] );      // a simple append eventually produces double data points...
                }
//std::cerr << std::endl;

            }

        } // end for gBeg, gEnd

    } // end while limdherm_

    return _ifail[0];
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void
LIMEX_A::computeAndSaveLimexTrajectory(double* t, double T, double* y)
{
    double* ztmp;
    double* z;

    if ( y == 0 ) z = _y0; else z = y;

    while ( (_ifail[0] == 0) && (*t < T) )
    {
        // limex_(
        limd_(
                &_n,
                _fcn, _jac,
                t, &T,
                z, _dy0,
                &_rtol, &_atol, &_h,
                _iOpt, _rOpt, _iPos,
                _ifail
             );

        ztmp = z;
        for (long j = 0; j < _n; ++j)
        {
            if ( std::fabs( *ztmp ) < EPMACH ) *ztmp = 0.0;
            ++ztmp;
        }

        if (*t <= T)
        {
            _solPoints.push_back( *t );
            ztmp = z;
            for (long j = 0; j < _n; ++j)
            {
                _solution[j].push_back( *ztmp++ );
            }
        }
    }

    return;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
int
LIMEX_A::integrateWithoutInterpolation()
{
    GridIterConst gBeg = _datPoints.begin();
    GridIterConst gEnd = _datPoints.end();
    double        t = _t0;
    double        T = _T;
    double        z[_n];

    for (long j = 0; j < _n; ++j)
    {
        z[j] = _y0[j];
        _dy0[j] = 0.0;
    }

    _solPoints.clear();
    _solution.clear();
    _data.clear();

    //

    _iOpt[11] = 1;      // single step mode ON
    _iOpt[12] = 0;      // no dense output
    _iOpt[13] = 0;
    _iOpt[14] = 0;

    _iOpt[15] = 0;      // type of call of limex_(): 0 initial call

    _iOpt[16] = 1;      // integration for t > T internally allowed

    //

    _rOpt[0] = _rOpt[1] = _rOpt[2] = 0.0;

    //

    for (GridIterConst it = gBeg; it != gEnd; ++it)
    {
        T = *it;

        computeAndSaveLimexTrajectory( &t, T, z );

        for (long j = 0; j < _n; ++j)
        {
            _data[j].push_back( z[j] );
        }
    }

    if ( t != _T )
    {
        T = _T;

        computeAndSaveLimexTrajectory( &t, T, z );
    }

    return _ifail[0];
}
//----------------------------------------------------------------------------
int
LIMEX_A::integrateWithoutInterpolation( unsigned n, double* yIni,
                                        double tLeft, double tRight )
{
    // const Real reduce = 0.001;

    GridIterConst dBeg = _datPoints.begin();
    GridIterConst dEnd = _datPoints.end();
    GridIterConst gBeg = std::lower_bound( dBeg, dEnd, tLeft);
    GridIterConst gEnd = std::upper_bound( dBeg, dEnd, tRight);
    double        t = tLeft;
    double        T = tRight;

    if ( n != (unsigned)_n ) return -99;

//    for (long j = 0; j < _n; ++j)
//    {
//        _y0[j]  = yIni[j];
//        _dy0[j] = 0.0;
//    }

    _solPoints.clear();
    _solution.clear();

    //

    _iOpt[11] = 1;      // single step mode ON
    _iOpt[12] = 0;      // no dense output
    _iOpt[13] = 0;
    _iOpt[14] = 0;

    if ( _iOpt[15] == -1 ) {
        _iOpt[15] = 0;  // type of call of limex_(): 0 initial call, 1 successive call

        _rOpt[0] = _rOpt[1] = _rOpt[2] = 0.0;
        _h = _rtol;
    }

    _iOpt[16] = 0;      // integration for t > T internally forbidden

    //

    for (GridIterConst it = gBeg; it != gEnd; ++it)
    {
        T = *it;

        computeAndSaveLimexTrajectory( &t, T, yIni );

        // if ( t != tRight )
        {
            for (long j = 0; j < _n; ++j)
            {
                _data[j][long(it-dBeg)] = yIni[j];  // allowing overwriting of data points!
                // _data[j].push_back( yIni[j] );   // a simple append eventually produces double data points...
            }
        }
    }

    if ( t != tRight )
    {
        T = tRight;

        computeAndSaveLimexTrajectory( &t, T, yIni );
    }

    //

//    for (long j = 0; j < _n; ++j)
//    {
//        yIni[j] = _y0[j];
//    }

    return _ifail[0];
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
ODESolver::Trajectory&
LIMEX_A::getSimulatedData()
{
    return _data;
}
//----------------------------------------------------------------------------
void
LIMEX_A::setODESystem(
                     Fcn            fcn,
                     Jac            jac,
                     double         t0,
                     Grid const&    y0,
                     Grid const&    refGrid,
                     double         tEnd,
                     int            bandwidth
                    )
{
    if ( (unsigned)_n != y0.size() )
    {
        _n = (int)y0.size();
        delete[] _y0;    _y0   = new double[_n];
        delete[] _dy0;   _dy0  = new double[_n];
        delete[] _iPos;  _iPos = new int[_n];
        // delete[] _Dense; _Dense = new double[4000*30];
    }

    _fcn = fcn;
    _jac = jac;
    _t0  = t0;
    _T   = tEnd;
    _h   = _rtol; // 0.0;

    _datPoints = refGrid;
    std::sort( _datPoints.begin(), _datPoints.end() );
    _data.clear();

    long j = 0;
    for (GridIterConst it = y0.begin(); it != y0.end(); ++it)
    {
        _y0[j]   = *it;
        _dy0[j]  = 0.0;
        _iPos[j] = 0;

        _data[j].resize( _datPoints.size() );

        ++j;
    }

    _solPoints.clear();
    _solution.clear();

    _iOpt[6] = 0;       // NO analytic Jacobian (supplied by routine 'jac')

    int bw = (bandwidth == 0) ? _n : bandwidth;

    _iOpt[7] = bw;      // lower band of Jacobian set to max. in relation to 'fcn'
    _iOpt[8] = bw;      // upper band of Jacobian set to max. in relation to 'fcn'

//    _iOpt[17] = 10;     // 'Jacobian.ps' output in step #10
//    _iOpt[9]  = 0;      // ... but only if Jacobian re-use is switched off!

    _iOpt[15] = -1;      // type of limex call: 0 initial, 1 successive

    return;
}
//----------------------------------------------------------------------------
